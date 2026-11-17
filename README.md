# BigClipper

A tool to identify pile-ups of clipped reads whose clipped portions align elsewhere in the genome, useful for identifying long-range rearrangements that may be missed by other SV callers. Written for long-read genome sequencing data, e.g. from PacBio and Oxford Nanopore.

## Requirements
Bedtools

Python >= 3

**Python libraries:**
  * pysam
  * google-crc32c (only for verified streaming from Google Cloud Storage, see below)

The cohort script (`bigclipper_cohort.py`, see below) also needs GNU sort (`gsort` from coreutils on macOS) and htslib's `bgzip` and `tabix`.

## Installation
To create a conda environment with all the prerequisites for running BigClipper locally:

```
git clone https://github.com/yuliamostovoy/bigclipper.git
cd bigclipper
conda env create -f environment.yml
conda activate bigclipper
```

Alternatively, a docker image and WDL workflow are available for running BigClipper on the cloud.

## Usage
BigClipper runs in two phases. The first -- which is slower and only needs to be run once -- reads a BAM or CRAM file, detects soft-clipped reads with supplementary alignments, and outputs their information to an intermediate file. It also writes a second file listing clipped read ends whose clipped sequence is not aligned anywhere. The second phase parses the intermediate file with a set of parameters and outputs pile-ups of clipped reads that pass the filtering parameters, optionally flagging pile-ups explained by mobile element insertions. The second phase can be run multiple times with different parameters.

Phase 1: process the BAM/CRAM file
```
python bigclipper/scripts/bigclipper_processbam.py [-h] -o OUTPUT_PREFIX [-d OUTPUT_DIR] [-t THREADS] [-r REFERENCE] [-m MIN_CLIP] bamfile

positional arguments:
  bamfile               BAM/CRAM file to read: a local path, '-' for stdin, or a URL htslib can read (gs://, https://, s3://; not integrity-checked, see below)

optional arguments:
  -o OUTPUT_PREFIX, --output_prefix OUTPUT_PREFIX
                        Prefix for output file
  -d OUTPUT_DIR, --output_dir OUTPUT_DIR
                        Directory for output file [default=current directory]
  -t THREADS, --threads THREADS
                        Number of threads for BAM/CRAM decompression [default=4]
  -r REFERENCE, --reference REFERENCE
                        Reference FASTA (needed for CRAM input)
  -m MIN_CLIP, --min_clip MIN_CLIP
                        Minimum length (bp) of a clipped read end without a split-read partner to write to the unaligned-clips file [default=20]
```
Phase 1 outputs (both described [below](#intermediate-files)):
* `[OUTPUT_DIR]/[OUTPUT_PREFIX]_intermediate.bed`: split-read breakends, the input to phase 2
* `[OUTPUT_DIR]/[OUTPUT_PREFIX]_unaligned_clips.bed`: clipped read ends whose clipped sequence is not aligned anywhere

Most of phase 1's run time is spent decompressing the BAM, so it is I/O- and decompression-bound; 4 threads is usually enough. CRAM input is read without decoding read sequences and base qualities, which makes it faster still.

To process a BAM/CRAM in a Google bucket without copying it to disk first, stream it through `bigclipper_verifystream.py`. This script checks the streamed bytes against the size and CRC32C that the bucket stores for the object:
```
set -o pipefail
BAM=gs://bucket/sample.bam
read -r SIZE CRC < <(gcloud storage objects describe "$BAM" --format='value(size,crc32c_hash)')
gcloud storage cat "$BAM" \
    | python bigclipper/scripts/bigclipper_verifystream.py --size "$SIZE" --crc32c "$CRC" \
    | python bigclipper/scripts/bigclipper_processbam.py - -o sample
```
Don't skip the check or `pipefail`. Phase 1 reading from stdin can't tell a stream that stops between two BGZF blocks (or CRAM containers) from a complete file. It would exit successfully having read only part of the reads, and with BAM the only sign is an htslib warning that the EOF marker is absent. `gcloud storage cat` does not checksum what it streams. With the check in place, a short, long or corrupted stream makes the pipeline exit non-zero; discard the intermediate file if it does. Copying the file first with `gcloud storage cp` is the alternative: it validates checksums itself.

Passing a `gs://` (or `https://`, `s3://`) URL straight to phase 1 has htslib read it over the network. That path has no end-to-end checksum, and it has silently dropped records from GCS before ([htslib#1037](https://github.com/samtools/htslib/issues/1037), still open). htslib also reads `GCS_OAUTH_TOKEN` once and doesn't refresh it, and access tokens expire after about an hour ([htslib#803](https://github.com/samtools/htslib/issues/803)). Phase 1 prints a warning when given a URL. For large files, prefer the verified stream above.

Phase 2: process the intermediate file and output clusters
```
python bigclipper/scripts/bigclipper_getclusters.py [-h] [-d MIN_DIST] [-c MIN_CLUSTER_COUNT] [-s CLUSTER_DIST] [-u MAX_UNIQUE_BREAKENDS] [-o OUTPUT]
                                                   [--mei_vcf MEI_VCF] [--mei_bed MEI_BED] [--mei_info_key MEI_INFO_KEY] [--mei_classes MEI_CLASSES]
                                                   [--mei_exclude_svtypes MEI_EXCLUDE_SVTYPES]
                                                   [--mei_sample MEI_SAMPLE] [--mei_all_sites] [--mei_window MEI_WINDOW] [--mei_partner_frac MEI_PARTNER_FRAC]
                                                   [--pass_only]
                                                   intermediate_file

positional arguments:
  intermediate_file           Intermediate file produced by bigclipper_processbam.py

optional arguments:
  -h, --help            show help message and exit
  -d MIN_DIST, --min_dist MIN_DIST
                        Minimum reference distance between split alignment positions of a single read (increase to find longer-range split-reads) (bp) [default=1]
  -c MIN_CLUSTER_COUNT, --min_cluster_count MIN_CLUSTER_COUNT
                        Minimum number of reads in a cluster (counting only reads whose split partner passes --min_dist) [default=5]
  -s CLUSTER_DIST, --cluster_dist CLUSTER_DIST
                        Maximum distance between supplementary alignment positions to cluster them together as a single "supplementary alignment region" [default=50]
  -u MAX_UNIQUE_BREAKENDS, --max_unique_breakends MAX_UNIQUE_BREAKENDS
                        Maximum number of supplementary alignment regions within a cluster (reduces false positives by filtering high-copy-number repeats) [default=10]
  -o OUTPUT, --output OUTPUT
                        Output VCF path [default: <intermediate prefix>_d<d>_c<c>_s<s>_u<u>.vcf next to the intermediate file]

MEI filtering (see below):
  --mei_vcf MEI_VCF     SV VCF/BCF with MEI-annotated calls (can be given more than once)
  --mei_bed MEI_BED     BED file of MEI insertion positions; every interval is treated as an MEI (can be given more than once)
  --mei_info_key MEI_INFO_KEY
                        INFO key marking a VCF record as an MEI; the value is reported as the MEI class [default=MEI_type]
  --mei_classes MEI_CLASSES
                        Comma-separated values of --mei_info_key to treat as MEIs, case-insensitive (e.g. ALU,LINE,SVA) [default: any non-missing value]
  --mei_exclude_svtypes MEI_EXCLUDE_SVTYPES
                        Comma-separated SVTYPEs to ignore in the MEI VCF; NONE to use all.
                        Deletions are excluded by default because a deleted reference MEI is not a new insertion [default=DEL]
  --mei_sample MEI_SAMPLE
                        Sample name in the MEI VCF; only calls where this sample has a non-reference genotype are used.
                        Required if the VCF has more than one sample, unless --mei_all_sites is set.
  --mei_all_sites       Use every MEI site in the VCF regardless of genotype
  --mei_window MEI_WINDOW
                        Max distance (bp) between a cluster or partner breakend and an MEI call [default=100]
  --mei_partner_frac MEI_PARTNER_FRAC
                        Flag a cluster as MEI_PARTNER if at least this fraction of its supporting reads have their partner breakend at an MEI call [default=0.5]
  --pass_only           Write only PASS clusters (drop MEI-flagged ones)
```

## Output
VCF file containing filtered clusters of split reads, one record per cluster of clipped-read breakends:
* `POS`: leftmost breakend position in the cluster (`CIPOS` gives the range of breakend positions in the cluster)
* `QUAL`: highest mapping quality among the clipped reads
* `ORI`: orientation of the clipped reads. `+` if they are clipped on the right (the aligned sequence lies to the left of the breakend), `-` if they are clipped on the left
* `ALNS`: where the clipped portions align ("partner" breakends), as `chrom:start-end|ori|count`, most frequent first. The partner `ori` is the same as the cluster's `ORI` if the clipped portion aligns in the same orientation as the clipped read, and the opposite if it is inverted relative to it. So both ends of a deletion (or a tandem duplication) show matching orientations (`ORI=+` with partner `+`, `ORI=-` with partner `-`), while an inversion shows opposite ones
* `NREADS`: number of reads supporting the cluster
* `PARTNER_MAPQ`, `PARTNER_ALNLEN`: median MAPQ and aligned read length of the supplementary alignments in each `ALNS` region
* `FILTER`: `PASS`, or `MEI_SITE`/`MEI_PARTNER` if the cluster is explained by an MEI (see below)

The output is position-sorted, so it can be bgzipped and indexed with tabix.

## MEI filtering
Most false positives are caused by mobile element insertions (MEIs): the inserted copy of an Alu, L1 or SVA element in the sample aligns to a reference copy of that element elsewhere in the genome. Clipped reads then pile up both at the insertion site and at the reference copy. Given an SV callset for the same sample in which MEIs are annotated (e.g. with an `INFO/MEI_type` tag, as added by L1MEAID + IntactMEI), phase 2 flags:
* `MEI_SITE`: clusters within `--mei_window` bp of an MEI call (the insertion site itself)
* `MEI_PARTNER`: clusters where at least `--mei_partner_frac` of the supporting reads have their partner breakend within `--mei_window` bp of an MEI call (the reference copy that the inserted sequence aligned to)

Flagged clusters are kept in the output with a non-PASS `FILTER` (use `--pass_only` to drop them, or `bcftools view -f PASS`). `MEI_CLASS` and `MEI_PARTNER_FRAC` record what was matched. The log reports how many MEI calls were used, broken down by class, so check it once for unexpected values. If the annotation tag can also hold non-MEI values (e.g. a general repeat annotation with `TANDEM`), restrict it with `--mei_classes`. The callset must use the same reference and contig names as the BAM; BigClipper warns if no contig names are shared, but cannot detect two assemblies that use the same contig names. Since MEI filtering only involves phase 2, it can be applied to existing intermediate files without re-reading the BAMs.

For a large multi-sample callset, reading the whole VCF once per sample is slow (about 2.5 minutes for a 2.5 GB, 2,000-sample BCF). Extract the MEI records once and pass the smaller file instead, keeping genotypes for per-sample filtering:
```
bcftools view -i 'INFO/MEI_type!="."' -Ob -o cohort.MEIs.bcf cohort.bcf
```

## Intermediate files
Both phase 1 outputs are tab-separated, sorted by position, and use 0-based positions. Each breakend is at the start of the alignment for a left clip (`ori` `-`) and at its end for a right clip (`ori` `+`).

`_intermediate.bed` has one line per clipped alignment end with an adjacent split alignment:
`chrom, pos, pos+1, ., ., ori, mapq, aligned_length, read_name, partner breakend (chrom:pos<ori>), partner_mapq, partner_aligned_length, HP, PS`

`HP` and `PS` are the alignment's haplotype and phase set tags (e.g. from HiPhase or WhatsHap), or `.` if the BAM has none. HP values are only comparable between alignments with the same PS. Supplementary alignments are phased at their own locus, so they are tagged less often than primary alignments. Phase 2 doesn't use these columns; the cohort evidence table reports the HP composition of each sample's reads at each breakend site. Intermediate files from version 2.0 have only the first 12 columns, and files from before 2.0 have only the first 10 and partner positions that are 1 bp too high. Phase 2 reads both (and corrects the old partner positions), so they don't need to be regenerated.

`_unaligned_clips.bed` has one line per clipped read end of at least `--min_clip` bp whose clipped sequence is not aligned anywhere else: the soft clips of reads without supplementary alignments, and the outer ends of split reads. Each clipped end appears in only one of the two files.
`chrom, pos, pos+1, ., ., ori, mapq, aligned_length, read_name, clip_length, HP, PS`

These clips come from insertions too long to be aligned within the read, sequence missing from the reference, and adapters or low-quality read ends, so they are noisy one read at a time. They are meant for pile-up analysis across reads or samples. For a 30x HiFi genome the file is a few tens of MB.

## Running on Terra / Cromwell
`wdl/Clipper.wdl` runs both phases. By default it streams the BAM/CRAM from its bucket rather than localizing it (`stream_bam`; set `billing_project` for requester-pays buckets). The stream is verified as described above, and the task fails if the bytes don't match the object. `stream_bam = false` copies the file to disk with `gcloud storage cp` instead, which costs disk and time but uses gcloud's own checksum validation. Pass `intermediate_file` instead of `aligned_bam` to re-run phase 2 (e.g. with MEI filtering) on an existing intermediate file. The workflow outputs the VCF, the intermediate file and the unaligned-clips file (`min_clip` sets the minimum clip length for the latter).

To build the docker image (from the `docker` directory): `make build`, then `make push`.

## Cohort integration and filtering
`wdl/CohortClipper.wdl` (script `scripts/bigclipper_cohort.py`) combines the phase 1 outputs of a cohort into cohort-level tables and a VCF, and lists the rare breakends called in affected samples together with the read evidence for them in each parent. It runs on a Terra sample_set: `sample_ids`, `intermediate_files` and `unaligned_clips` are parallel arrays (e.g. `this.samples.sample_id`, `this.samples.intermediate_file`, `this.samples.unaligned_clips`), plus a PED file and, optionally, the cohort MEI-annotated SV callset (`mei_vcf`), BED files of regions to flag (`region_beds`) and the reference `.fai`. Sample IDs must be the IDs used in the PED file and the MEI VCF. The workflow first checks that the arrays have the same length, that no ID is repeated, that every sample is in the PED file (and in the MEI VCF, if given), and that each intermediate file and unaligned-clips file share a `<prefix>_intermediate.bed` / `<prefix>_unaligned_clips.bed` name, so that a shifted array can't attach one sample's clips to another sample.

### How samples are combined
1. **Clusters** (batches of `batch_size` samples): phase 2 on each sample with the same clustering as `bigclipper_getclusters.py` (default `min_dist` 5000, `min_cluster_count` 5) and per-sample MEI flags. Each cluster gives a breakend (chrom, ORI, position = median read breakend), and each of its partner groups with at least `junction_min_reads` (2) reads gives a junction.
2. **Breakend sites**: clusters on the same chromosome with the same ORI form a site if their positions are less than `site_window` (50) bp to the right of the site's first cluster. Sites are bounded and don't chain, so a rare breakend 60 bp from a common one stays separate. A site's range (`start`-`end`) spans the read breakends of its clusters. Cohort frequency and parental evidence are counted per site, independently of the partner: recurrent artifacts often have different partners in different samples, and a parent's reads may be clipped without an aligned partner.
3. **Junctions** (the VCF records): junctions from different samples are the same if ORI, partner chromosome and partner orientation are identical and both breakends are less than `site_window` bp from those of a representative junction (the one with the most reads). Matches are not chained.
4. **Evidence scan** (batches): each sample's split-read breakends (any partner, at any distance) and unaligned clips are counted at every site with the same orientation, within `site_window` - 1 bp of the site's range, whether or not the sample was called there. Reads are counted once per site.
5. **Tables and VCF**, then **filtering** (below).

### Filtering
* A sample is a **carrier** of a site if it was called there or has at least `af_min_reads` (2) reads there. **AF** = carriers / (2 x number of samples), counting each carrier as one allele (clipped reads don't give zygosity). Counting uncalled samples with 2 or more reads means recurrent noise just below the calling threshold still counts toward AF; single reads (about 0.3% of samples at a random position in HG003) don't. `AF_called` counts called samples only.
* A **candidate** is an affected sample (PED phenotype 2) called at a site with AF <= `max_af` (default 0.01, about 41 carriers in 2,097 samples). MEI-flagged candidates are kept and labeled unless `exclude_mei` is set.
* **Parental evidence**: a parent has evidence if it has any read (split or clip) or a call at the candidate's site, or at its mate site (the site of the other breakend of the candidate's own junction). `parental_evidence` is `paternal`, `maternal`, `both`, `no_parental_evidence`, or `parents_incomplete` (a parent isn't in the cohort and the other has no evidence). The de novo subset is the `no_parental_evidence` rows. A parent with low yield (`father_rel_yield`/`mother_rel_yield`, split-read lines relative to the cohort median) has less chance of showing reads at an inherited breakend.
* Both breakends of a junction become candidates if both were called, so most events appear twice (once per side).

To change the thresholds without redoing the integration, run the workflow again with `site_table`, `evidence_table` and `sample_stats` from the first run (plus `ped`, `max_af`, `af_min_reads`, `exclude_mei`): only the Filter task runs.

### Outputs
Positions in all cohort outputs are 1-based breakend bases as in VCF breakend notation: the last aligned base for an ORI `+` breakend (clipped on the right) and the first aligned base for ORI `-`. For `+` breakends this is 1 bp left of the `POS` written by `bigclipper_getclusters.py`.

* `cohort_site_table` (one row per breakend site, unfiltered): `pos` (median cluster position), `ori`, `start`/`end`; `n_called`, `n_called_families`, `n_called_affected`, `n_called_mei` (called samples, and how many had the cluster MEI-flagged); `n_called_nearby` (samples called at other sites of the same orientation within 2 x `site_window` bp); `n_carriers`, `AF`, `AF_called`; `n_samples_reads1`/`reads2` (samples with at least 1/2 reads, called or not), `n_samples_split1`, `n_samples_clip1`, `n_families_reads1`; `n_junctions`; the most common partner of the called samples (`top_partner_chrom`/`pos`/`ori`, `top_partner_samples`) and the site at its other breakend (`mate_site`); `mei_classes` (from the called samples' MEI flags), `mei_near` (MEI calls in `mei_vcf` near the site, any sample); `regions`/`partner_regions` (`region_beds` overlapping the site or its top partner).
* `cohort_evidence_table` (one row per site and sample with a call or any read): `called_reads` (0 if not called), `filter` (the sample's cluster FILTER), `SR` (reads with a split-read breakend at the site), `CL` (reads with an unaligned clip), `R` (distinct reads, either kind), `H1`/`H2`/`HU` (those reads on haplotype 1, 2 or unphased), and for called samples their own `partner` and `mate_site`.
* `cohort_sample_stats`: per sample, family, parents, affected status, cluster and junction counts, split-read and clip line counts, and their yield relative to the cohort median.
* `cohort_junctions_vcf` (bgzipped, tabix-indexed): one BND record per cohort junction, `ALT` in breakend notation (`N[chr:pos[` for ORI `+` with partner ori `+`, `N]chr:pos]` for `+`/`-`, `]chr:pos]N` for `-`/`-`, `[chr:pos[N` for `-`/`+`), `MATEID` linking the two sides when both were called. INFO: `NSAMPLES`, `NREADS`, `NMEI`, and the site-level `SITE`, `SITE_CALLED`, `SITE_CARRIERS`, `SITE_AF`, `REGIONS`. FILTER is `MEI` if more than half of the samples with the junction had it MEI-flagged. FORMAT `GT:NR:SR:CL`: `0/1` with the junction's read count `NR` for samples with the junction; `./.` with `SR`/`CL` for samples with reads at the site but no junction call; `.` otherwise. Clipped reads give no reference evidence, so there are no `0/0` genotypes.
* `cohort_junctions_table`: the same junctions as a TSV, with each sample's read count.
* `rare_candidates`, `denovo_candidates`: one row per candidate (site and affected sample) with the sample's evidence, partner and mate site, site AF and counts, `parental_evidence`, each parent's `R`/`SR`/`CL` at the site, `R` at the mate site (`father_mate_R`/`mother_mate_R`) and relative yield, other relatives with reads (`relatives_with_reads`, A/U = affected/unaffected), and the site's partner, MEI and region annotations.

### Resources

For 2,097 samples in batches of 50, with synthetic inputs of about 7,000 clusters and 11,000 junctions per sample (as in real samples):

| Step | Time | Peak memory | Output |
|---|---|---|---|
| `DefineSites` | 12 min | 4.1 GB | 6.9 M sites, 11 M junctions |
| `ScanEvidence` | 1.4–4.5 s per sample, plus 25 s to load the sites | 0.4 GB | |
| `CohortTables` | 18 min | 0.35 GB | 59 M evidence rows, 0.85 GB VCF |
| `Filter` | 7 min | 0.6 GB | |

The `DefineSites` (16 GB) and `CohortTables` (8 GB) memory defaults leave headroom for about twice these counts. Disk for every task is sized from its inputs. The synthetic evidence came from a noise model and may be half the real number of rows.

## Tests
```
python -m unittest discover -s tests -v
```
