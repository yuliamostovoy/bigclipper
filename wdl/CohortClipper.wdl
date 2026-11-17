version 1.0

######################################################################################
## BigClipper cohort integration and filtering
##
## From the per-sample outputs of Clipper.wdl (intermediate and unaligned-clips files) and a PED file:
##   - a cohort table of breakend sites (unfiltered), with the read evidence of every sample at each site
##   - a cohort VCF of junctions (BND records matched across samples)
##   - rare breakends called in affected samples (cohort AF <= max_af), with parental evidence
##   - the de novo subset of those (no reads at either breakend in either parent)
## Run on a sample_set: sample_ids, intermediate_files and unaligned_clips are parallel arrays (this.samples.*).
## To change the filter thresholds without redoing the integration, give site_table, evidence_table and
## sample_stats from a previous run: only the Filter task runs.
######################################################################################

workflow CohortClipper {
    input {
        Array[String] sample_ids = []
        Array[File] intermediate_files = []
        Array[File] unaligned_clips = []
        File ped
        File? mei_vcf
        Array[File] region_beds = []
        File? reference_fai
        String prefix = "cohort"

        Int min_dist = 5000
        Int min_cluster_count = 5
        Int max_cluster_dist = 50
        Int max_unique_breakends = 10
        Int site_window = 50
        Int junction_min_reads = 2
        Int batch_size = 50

        Float max_af = 0.01
        Int af_min_reads = 2
        Boolean exclude_mei = false

        String mei_info_key = "MEI_type"
        String? mei_classes
        String mei_exclude_svtypes = "DEL"
        Int mei_window = 100
        Float mei_partner_frac = 0.5

        File? site_table
        File? evidence_table
        File? sample_stats

        String docker = "quay.io/ymostovoy/clipper:2.2"
        RuntimeAttr? runtime_attr_call_clusters
        RuntimeAttr? runtime_attr_define_sites
        RuntimeAttr? runtime_attr_scan_evidence
        RuntimeAttr? runtime_attr_cohort_tables
        RuntimeAttr? runtime_attr_filter
    }

    parameter_meta {
        sample_ids:             "sample IDs (as in the PED file and the MEI VCF), e.g. this.samples.sample_id"
        intermediate_files:     "<prefix>_intermediate.bed from Clipper.wdl, in the same order as sample_ids"
        unaligned_clips:        "<prefix>_unaligned_clips.bed from Clipper.wdl, in the same order; each must have the same <prefix> as its intermediate file"
        ped:                    "PED file (family, sample, father, mother, sex, phenotype; phenotype 2 = affected) covering every sample"
        mei_vcf:                "cohort SV VCF/BCF with MEI-annotated calls (INFO/mei_info_key) and genotypes for every sample; clusters at a sample's own MEI calls are flagged MEI_SITE/MEI_PARTNER"
        region_beds:            "BED files of regions to flag in the tables (e.g. satellites, segdups); each is named by its file name"
        reference_fai:          "reference .fai, for contig lengths in the VCF header"
        min_dist:               "minimum distance (bp) between a read's clipped breakend and its split-read partner for the read to count toward a cluster (default=5000)"
        min_cluster_count:      "minimum number of reads in a cluster (default=5)"
        max_cluster_dist:       "max distance (bp) between breakends to cluster them together (default=50)"
        max_unique_breakends:   "maximum number of partner groups in a cluster (default=10)"
        site_window:            "clusters less than this many bp apart (same orientation) form one breakend site; junctions match if both breakends are less than this many bp apart (default=50)"
        junction_min_reads:     "partner groups with at least this many reads become junctions in the VCF (default=2)"
        batch_size:             "samples per CallClusters/ScanEvidence job (default=50)"
        max_af:                 "keep breakends with cohort AF at most this (default=0.01)"
        af_min_reads:           "a sample counts as a carrier if it was called at the site or has at least this many reads there (default=2)"
        exclude_mei:            "drop candidates whose cluster is MEI-flagged in the candidate's sample (default=false: keep and label them)"
        site_table:             "filter-only mode: site table from a previous run"
        evidence_table:         "filter-only mode: evidence table from a previous run"
        sample_stats:           "filter-only mode: sample stats from a previous run"
    }

    Boolean filter_only = defined(site_table) && defined(evidence_table) && defined(sample_stats)

    if (!filter_only) {
        call CheckInputs {
            input:
                sample_ids = sample_ids,
                intermediate_paths = intermediate_files,
                clips_paths = unaligned_clips,
                ped = ped,
                docker = docker
        }

        if (defined(mei_vcf)) {
            call ExtractMEIs {
                input:
                    mei_vcf = select_first([mei_vcf]),
                    sample_list = CheckInputs.sample_list,
                    mei_info_key = mei_info_key,
                    docker = docker
            }
        }

        Int n_samples = length(sample_ids)
        Int n_batches = (n_samples + batch_size - 1) / batch_size

        scatter (b in range(n_batches)) {
            scatter (k in range(batch_size)) {
                Int i = b * batch_size + k
                if (i < n_samples) {
                    String batch_sample = sample_ids[i]
                    File batch_intermediate = intermediate_files[i]
                    File batch_clip = unaligned_clips[i]
                }
            }
            Array[String] batch_samples = select_all(batch_sample)
            Array[File] batch_intermediates = select_all(batch_intermediate)
            Array[File] batch_clips = select_all(batch_clip)

            call CallClusters {
                input:
                    batch_name = "batch~{b}",
                    sample_list = CheckInputs.sample_list,
                    samples = batch_samples,
                    intermediate_files = batch_intermediates,
                    mei_bcf = ExtractMEIs.mei_bcf,
                    mei_bcf_index = ExtractMEIs.mei_bcf_index,
                    min_dist = min_dist,
                    min_cluster_count = min_cluster_count,
                    max_cluster_dist = max_cluster_dist,
                    max_unique_breakends = max_unique_breakends,
                    junction_min_reads = junction_min_reads,
                    mei_info_key = mei_info_key,
                    mei_classes = mei_classes,
                    mei_exclude_svtypes = mei_exclude_svtypes,
                    mei_window = mei_window,
                    mei_partner_frac = mei_partner_frac,
                    docker = docker,
                    runtime_attr_override = runtime_attr_call_clusters
            }
        }

        call DefineSites {
            input:
                clusters = CallClusters.clusters,
                junctions = CallClusters.junctions,
                site_window = site_window,
                prefix = prefix,
                docker = docker,
                runtime_attr_override = runtime_attr_define_sites
        }

        scatter (b in range(n_batches)) {
            call ScanEvidence {
                input:
                    batch_name = "batch~{b}",
                    sites = DefineSites.sites,
                    samples = batch_samples[b],
                    intermediate_files = batch_intermediates[b],
                    unaligned_clips = batch_clips[b],
                    site_window = site_window,
                    docker = docker,
                    runtime_attr_override = runtime_attr_scan_evidence
            }
        }

        call CohortTables {
            input:
                sample_list = CheckInputs.sample_list,
                ped = ped,
                sites = DefineSites.sites,
                calls = DefineSites.calls,
                junctions = DefineSites.junction_table,
                evidence = ScanEvidence.evidence,
                cluster_stats = CallClusters.cluster_stats,
                scan_stats = ScanEvidence.scan_stats,
                region_beds = region_beds,
                reference_fai = reference_fai,
                mei_bcf = ExtractMEIs.mei_bcf,
                mei_bcf_index = ExtractMEIs.mei_bcf_index,
                site_window = site_window,
                af_min_reads = af_min_reads,
                mei_info_key = mei_info_key,
                mei_classes = mei_classes,
                mei_exclude_svtypes = mei_exclude_svtypes,
                mei_window = mei_window,
                prefix = prefix,
                docker = docker,
                runtime_attr_override = runtime_attr_cohort_tables
        }
    }

    call Filter {
        input:
            site_table = select_first([site_table, CohortTables.site_table]),
            evidence_table = select_first([evidence_table, CohortTables.evidence_table]),
            sample_stats = select_first([sample_stats, CohortTables.sample_stats]),
            ped = ped,
            max_af = max_af,
            af_min_reads = af_min_reads,
            exclude_mei = exclude_mei,
            prefix = prefix,
            docker = docker,
            runtime_attr_override = runtime_attr_filter
    }

    output {
        File rare_candidates = Filter.rare_candidates
        File denovo_candidates = Filter.denovo_candidates
        File? cohort_site_table = CohortTables.site_table
        File? cohort_evidence_table = CohortTables.evidence_table
        File? cohort_sample_stats = CohortTables.sample_stats
        File? cohort_junctions_vcf = CohortTables.junctions_vcf
        File? cohort_junctions_vcf_index = CohortTables.junctions_vcf_index
        File? cohort_junctions_table = DefineSites.junction_table
    }
}

task CheckInputs {
    input {
        Array[String] sample_ids
        Array[String] intermediate_paths
        Array[String] clips_paths
        File ped
        String docker
    }

    command <<<
        set -euo pipefail
        python3 /opt/bigclipper/bigclipper_cohort.py check --samples ~{write_lines(sample_ids)} --ped ~{ped} \
            --intermediate_list ~{write_lines(intermediate_paths)} --clips_list ~{write_lines(clips_paths)}
        cp ~{write_lines(sample_ids)} samples.txt
    >>>

    output { File sample_list = "samples.txt" }

    runtime {
        cpu: 1
        memory: "2 GiB"
        disks: "local-disk 10 HDD"
        preemptible: 3
        maxRetries: 1
        docker: docker
    }
}

task ExtractMEIs {
    input {
        File mei_vcf
        File sample_list
        String mei_info_key
        String docker
    }

    Int disk_size = 2 * ceil(size(mei_vcf, "GB")) + 10

    command <<<
        set -euxo pipefail
        bcftools query -l ~{mei_vcf} | sort > vcf_samples.txt
        sort ~{sample_list} | comm -23 - vcf_samples.txt > missing.txt
        if [ -s missing.txt ]; then
            echo "Error: $(wc -l < missing.txt) samples are not in ~{basename(mei_vcf)}, e.g.:" >&2
            head missing.txt >&2
            exit 1
        fi
        bcftools view -i 'INFO/~{mei_info_key}!="."' -Ob -o meis.bcf ~{mei_vcf}
        bcftools index meis.bcf
        echo "$(bcftools index -n meis.bcf) MEI-annotated records" >&2
    >>>

    output {
        File mei_bcf = "meis.bcf"
        File mei_bcf_index = "meis.bcf.csi"
    }

    runtime {
        cpu: 1
        memory: "4 GiB"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 3
        maxRetries: 1
        docker: docker
    }
}

task CallClusters {
    input {
        String batch_name
        File sample_list
        Array[String] samples
        Array[File] intermediate_files
        File? mei_bcf
        File? mei_bcf_index
        Int min_dist
        Int min_cluster_count
        Int max_cluster_dist
        Int max_unique_breakends
        Int junction_min_reads
        String mei_info_key
        String? mei_classes
        String mei_exclude_svtypes
        Int mei_window
        Float mei_partner_frac
        String docker
        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = ceil(size(intermediate_files, "GB") + size(mei_bcf, "GB")) + 10

    command <<<
        set -euxo pipefail
        df -h .
        echo "~{batch_name}: ~{length(samples)} of $(wc -l < ~{sample_list}) samples" >&2
        paste ~{write_lines(samples)} ~{write_lines(intermediate_files)} > manifest.tsv
        /usr/bin/time -v python3 /opt/bigclipper/bigclipper_cohort.py clusters --manifest manifest.tsv --out_prefix ~{batch_name} \
            -d ~{min_dist} -c ~{min_cluster_count} -s ~{max_cluster_dist} -u ~{max_unique_breakends} \
            --junction_min_reads ~{junction_min_reads} \
            ~{"--mei_vcf " + mei_bcf} --mei_info_key ~{mei_info_key} ~{"--mei_classes " + mei_classes} \
            --mei_exclude_svtypes ~{mei_exclude_svtypes} --mei_window ~{mei_window} --mei_partner_frac ~{mei_partner_frac}
        df -h .
    >>>

    output {
        File clusters = "~{batch_name}.clusters.tsv.gz"
        File junctions = "~{batch_name}.junctions.tsv.gz"
        File cluster_stats = "~{batch_name}.cluster_stats.tsv"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             4,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  3,
        max_retries:        1,
        docker:             docker
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task DefineSites {
    input {
        Array[File] clusters
        Array[File] junctions
        Int site_window
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 6 * ceil(size(clusters, "GB") + size(junctions, "GB")) + 20

    command <<<
        set -euxo pipefail
        lscpu | grep -E '^CPU\(s\)|^Model name' || true
        THREADS=$(nproc)
        MEM_MB=$(awk '/MemTotal/ {print int($2 / 1024 / 2)}' /proc/meminfo)
        mkdir -p tmp
        /usr/bin/time -v python3 /opt/bigclipper/bigclipper_cohort.py sites \
            --clusters ~{sep=' ' clusters} --junctions ~{sep=' ' junctions} --window ~{site_window} \
            --out_prefix ~{prefix} --tmpdir tmp --sort_mem ${MEM_MB}M --threads ${THREADS}
        df -h .
    >>>

    output {
        File sites = "~{prefix}.sites.tsv.gz"
        File calls = "~{prefix}.calls.tsv.gz"
        File junction_table = "~{prefix}.junctions.tsv.gz"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             16,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        1,
        docker:             docker
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " SSD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task ScanEvidence {
    input {
        String batch_name
        File sites
        Array[String] samples
        Array[File] intermediate_files
        Array[File] unaligned_clips
        Int site_window
        String docker
        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = ceil(size(intermediate_files, "GB") + size(unaligned_clips, "GB") + 3 * size(sites, "GB")) + 10

    command <<<
        set -euxo pipefail
        df -h .
        paste ~{write_lines(samples)} ~{write_lines(intermediate_files)} ~{write_lines(unaligned_clips)} > manifest.tsv
        mkdir -p tmp
        /usr/bin/time -v python3 /opt/bigclipper/bigclipper_cohort.py scan --sites ~{sites} --manifest manifest.tsv \
            --window ~{site_window} --out_prefix ~{batch_name} --tmpdir tmp
        df -h .
    >>>

    output {
        File evidence = "~{batch_name}.evidence.tsv.gz"
        File scan_stats = "~{batch_name}.scan_stats.tsv"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             4,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  3,
        max_retries:        1,
        docker:             docker
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task CohortTables {
    input {
        File sample_list
        File ped
        File sites
        File calls
        File junctions
        Array[File] evidence
        Array[File] cluster_stats
        Array[File] scan_stats
        Array[File] region_beds
        File? reference_fai
        File? mei_bcf
        File? mei_bcf_index
        Int site_window
        Int af_min_reads
        String mei_info_key
        String? mei_classes
        String mei_exclude_svtypes
        Int mei_window
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 10 * ceil(size(sites, "GB") + size(calls, "GB") + size(junctions, "GB") + size(evidence, "GB")) + 20

    command <<<
        set -euxo pipefail
        lscpu | grep -E '^CPU\(s\)|^Model name' || true
        /usr/bin/time -v python3 /opt/bigclipper/bigclipper_cohort.py tables --samples ~{sample_list} --ped ~{ped} \
            --sites ~{sites} --calls ~{calls} --junctions ~{junctions} \
            --evidence ~{sep=' ' evidence} --cluster_stats ~{sep=' ' cluster_stats} --scan_stats ~{sep=' ' scan_stats} \
            ~{true="--region_bed " false="" length(region_beds) > 0}~{sep=' --region_bed ' region_beds} \
            ~{"--fai " + reference_fai} ~{"--mei_vcf " + mei_bcf} --mei_info_key ~{mei_info_key} ~{"--mei_classes " + mei_classes} \
            --mei_exclude_svtypes ~{mei_exclude_svtypes} --mei_window ~{mei_window} \
            --window ~{site_window} --af_min_reads ~{af_min_reads} --out_prefix ~{prefix} --threads $(nproc)
        df -h .
    >>>

    output {
        File site_table = "~{prefix}.site_table.tsv.gz"
        File evidence_table = "~{prefix}.evidence.tsv.gz"
        File sample_stats = "~{prefix}.sample_stats.tsv"
        File junctions_vcf = "~{prefix}.junctions.vcf.gz"
        File junctions_vcf_index = "~{prefix}.junctions.vcf.gz.tbi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             8,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        1,
        docker:             docker
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " SSD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task Filter {
    input {
        File site_table
        File evidence_table
        File sample_stats
        File ped
        Float max_af
        Int af_min_reads
        Boolean exclude_mei
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 3 * ceil(size(site_table, "GB") + size(evidence_table, "GB")) + 10

    command <<<
        set -euxo pipefail
        /usr/bin/time -v python3 /opt/bigclipper/bigclipper_cohort.py filter --site_table ~{site_table} \
            --evidence ~{evidence_table} --sample_stats ~{sample_stats} --ped ~{ped} \
            --max_af ~{max_af} --af_min_reads ~{af_min_reads} ~{if exclude_mei then "--exclude_mei" else ""} \
            --out_prefix ~{prefix}
    >>>

    output {
        File rare_candidates = "~{prefix}.rare_candidates.tsv"
        File denovo_candidates = "~{prefix}.denovo_candidates.tsv"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             4,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  3,
        max_retries:        1,
        docker:             docker
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
    String? docker
}
