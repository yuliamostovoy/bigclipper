import sys
import os
import time
import bisect
import argparse
import subprocess
from collections import Counter, defaultdict

VERSION = "2.2"
LEGACY_NCOLS = 10  # pre-2.0 intermediate files (partner positions 1 bp too high)
OLD_DEFAULT_SUFFIX = ".bed"


def main():
    parser = argparse.ArgumentParser(description="Filter intermediate file derived from bigclipper_processbam.py\n")
    parser.add_argument('intermediate_file', help='Intermediate file produced by bigclipper_processbam.py')
    parser.add_argument('-d', '--min_dist', dest='min_dist', help="Minimum reference distance between split alignment positions of a single read (increase to find longer-range split-reads) (bp) [default=1]", type=int, default=1)
    parser.add_argument('-c', '--min_cluster_count', dest='min_cluster_count', help="Minimum number of reads in a cluster (counting only reads whose split partner passes --min_dist) [default=5]", type=int, default=5)
    parser.add_argument('-s', '--cluster_dist', dest='cluster_dist', help="Maximum distance between supplementary alignment positions to cluster them together as a single \"supplementary alignment region\" [default=50]", type=int, default=50)
    parser.add_argument('-u', '--max_unique_breakends', dest='max_unique_breakends', help="Maximum number of supplementary alignment regions within a cluster (reduces false positives by filtering high-copy-number repeats) [default=10]", type=int, default=10)
    parser.add_argument('-o', '--output', dest='output', default=None, help="Output VCF path [default: <intermediate prefix>_d<d>_c<c>_s<s>_u<u>.vcf next to the intermediate file]")

    mei = parser.add_argument_group("MEI filtering", "Flag clusters explained by mobile element insertions (MEIs) in an SV callset for the same sample. "
                                    "A new MEI's inserted sequence aligns to a reference copy of the element elsewhere, creating split-read clusters both at the insertion site "
                                    "and at the reference copy. Flagged clusters keep their record but get a non-PASS FILTER.")
    mei.add_argument('--mei_vcf', dest='mei_vcf', action='append', default=[], help="SV VCF/BCF with MEI-annotated calls (can be given more than once)")
    mei.add_argument('--mei_bed', dest='mei_bed', action='append', default=[], help="BED file of MEI insertion positions; every interval is treated as an MEI (can be given more than once)")
    mei.add_argument('--mei_info_key', dest='mei_info_key', default='MEI_type', help="INFO key marking a VCF record as an MEI; the value is reported as the MEI class [default=MEI_type]")
    mei.add_argument('--mei_classes', dest='mei_classes', default=None, help="Comma-separated values of --mei_info_key to treat as MEIs, case-insensitive (e.g. ALU,LINE,SVA) [default: any non-missing value]")
    mei.add_argument('--mei_exclude_svtypes', dest='mei_exclude_svtypes', default='DEL', help="Comma-separated SVTYPEs to ignore in the MEI VCF; NONE to use all. "
                     "Deletions are excluded by default because a deleted reference MEI is not a new insertion [default=DEL]")
    mei.add_argument('--mei_sample', dest='mei_sample', default=None, help="Sample name in the MEI VCF; only calls where this sample has a non-reference genotype are used. "
                     "Required if the VCF has more than one sample, unless --mei_all_sites is set. Not needed for single-sample VCFs.")
    mei.add_argument('--mei_all_sites', dest='mei_all_sites', action='store_true', help="Use every MEI site in the VCF regardless of genotype (e.g. a cohort sites-only VCF)")
    mei.add_argument('--mei_window', dest='mei_window', type=int, default=100, help="Max distance (bp) between a cluster or partner breakend and an MEI call [default=100]")
    mei.add_argument('--mei_partner_frac', dest='mei_partner_frac', type=float, default=0.5, help="Flag a cluster as MEI_PARTNER if at least this fraction of its supporting reads have their partner breakend at an MEI call [default=0.5]")
    mei.add_argument('--pass_only', dest='pass_only', action='store_true', help="Write only PASS clusters (drop MEI-flagged ones)")

    args = parser.parse_args()

    if args.output:
        outpath = args.output
    else:
        prefix = args.intermediate_file[:-len(OLD_DEFAULT_SUFFIX)] if args.intermediate_file.endswith(OLD_DEFAULT_SUFFIX) else args.intermediate_file
        outpath = "%s_d%d_c%d_s%d_u%d.vcf" % (prefix, args.min_dist, args.min_cluster_count, args.cluster_dist, args.max_unique_breakends)

    mei_index = None
    if args.mei_vcf or args.mei_bed:
        mei_index = MEIIndex(args.mei_window)
        for path in args.mei_vcf:
            mei_index.load_vcf(path, args.mei_info_key, args.mei_classes, args.mei_exclude_svtypes, args.mei_sample, args.mei_all_sites)
        for path in args.mei_bed:
            mei_index.load_bed(path)

    records = []
    contigs = []  # in order of first appearance
    seen_contigs = set()
    n_legacy = 0
    for cluster in read_clusters(args.intermediate_file, args.cluster_dist):
        if cluster[0]['legacy']:
            n_legacy += 1
        for line in cluster:
            for contig in (line['ref'], line['partner_ref']):
                if contig not in seen_contigs:
                    seen_contigs.add(contig)
                    contigs.append(contig)
        if len(set(line['read'] for line in cluster)) < args.min_cluster_count:
            continue
        record = process_cluster(cluster, args.min_dist, args.min_cluster_count, args.max_unique_breakends, args.cluster_dist, mei_index, args.mei_partner_frac)
        if record is None:
            continue
        if args.pass_only and record['filters']:
            continue
        records.append(record)

    if n_legacy:
        sys.stderr.write("Note: intermediate file is in the pre-2.0 format; correcting its partner positions by -1 bp (SA-tag POS is 1-based)\n")
    if mei_index is not None:
        mei_index.check_contigs(seen_contigs)

    contig_order = {c: i for i, c in enumerate(contigs)}
    records.sort(key=lambda r: (contig_order[r['ref']], r['pos']))
    with open(outpath, 'w') as outfile:
        write_header(outfile, contigs, args, mei_index)
        for record in records:
            print_cluster(record, outfile)

    n_flagged = sum(1 for r in records if r['filters'])
    sys.stderr.write("Wrote %d clusters to %s%s\n" % (len(records), outpath,
                     " (%d flagged as MEI-related)" % n_flagged if mei_index is not None and not args.pass_only else ""))


def read_clusters(intermediate_file, cluster_dist):
    """Cluster clipped-read breakends with bedtools (by position, separately for each orientation) and
    yield each cluster as a list of parsed lines."""
    cmd = ['bedtools', 'cluster', '-d', str(cluster_dist), '-s', '-i', intermediate_file]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, universal_newlines=True)
    current_cluster = []
    current_cluster_ID = None
    for raw in proc.stdout:
        line = parse_line(raw)
        if line['cluster'] != current_cluster_ID:
            if current_cluster:
                yield current_cluster
            current_cluster = []
            current_cluster_ID = line['cluster']
        current_cluster.append(line)
    if current_cluster:
        yield current_cluster
    proc.stdout.close()
    if proc.wait() != 0:
        sys.exit("Error: '%s' failed with exit code %d" % (' '.join(cmd), proc.returncode))


def parse_line(raw):
    # columns: chrom, start, end, ., ., ori, mapq, aligned_length, read_name, partner, [partner_mapq, partner_aligned_length,
    # [HP, PS,]] cluster_ID  (pre-2.0 files have 10 columns, 2.0 files 12, later files 14; the cluster ID is always last)
    fields = raw.rstrip('\n').split('\t')
    legacy = len(fields) == LEGACY_NCOLS + 1
    partner = fields[9]
    partner_ref, partner_pos = partner[:-1].rsplit(':', 1)  # contig names may contain ':'
    partner_pos = int(partner_pos)
    if legacy:
        partner_pos -= 1
    return {'ref': fields[0], 'pos': int(fields[1]), 'ori': fields[5], 'mapq': int(fields[6]),
            'read': fields[8], 'partner_ref': partner_ref, 'partner_pos': partner_pos, 'partner_ori': partner[-1],
            'partner_mapq': None if legacy else int(fields[10]),
            'partner_aligned_length': None if legacy else int(fields[11]),
            'cluster': fields[-1], 'legacy': legacy}


def process_cluster(cluster, min_dist, min_cluster_count, max_unique_breakends, thresh, mei_index, mei_partner_frac):
    ref = cluster[0]['ref']
    pos = cluster[0]['pos']
    ID = int(cluster[0]['cluster'])

    # keep reads whose partner breakend is on another contig or at least min_dist away
    passing = [line for line in cluster
               if line['partner_ref'] != ref or abs(line['pos'] - line['partner_pos']) >= min_dist]
    if len(set(line['read'] for line in passing)) < min_cluster_count:
        return None

    # cluster the partner breakends (single linkage, separately per contig and orientation)
    by_partner = defaultdict(list)
    for line in passing:
        by_partner[(line['partner_ref'], line['partner_ori'])].append(line)
    cluster_info = []
    for (breakend_ref, ori), lines in by_partner.items():
        lines.sort(key=lambda x: x['partner_pos'])
        group = [lines[0]]
        for line in lines[1:]:
            if line['partner_pos'] - group[-1]['partner_pos'] > thresh:
                cluster_info.append(summarize_partner_group(breakend_ref, ori, group))
                group = []
            group.append(line)
        cluster_info.append(summarize_partner_group(breakend_ref, ori, group))

    if len(cluster_info) > max_unique_breakends:
        return None
    cluster_info.sort(key=lambda x: (-x['count'], x['ref'], x['pos1']))  # sort by frequency

    positions = [line['pos'] for line in cluster]
    record = {'ref': ref, 'pos': pos, 'ID': ID, 'ori': cluster[0]['ori'],
              'qual': max(line['mapq'] for line in cluster),
              'span': max(positions) - min(positions),
              'nreads': len(set(line['read'] for line in passing)),
              'cluster_info': cluster_info, 'filters': [], 'mei_classes': set(), 'mei_partner_frac': None}

    if mei_index is not None:
        site_classes = mei_index.query(ref, min(positions), max(positions))
        if site_classes:
            record['filters'].append('MEI_SITE')
            record['mei_classes'].update(site_classes)
        partner_hits = set()
        for line in passing:
            classes = mei_index.query(line['partner_ref'], line['partner_pos'], line['partner_pos'])
            if classes:
                partner_hits.add(line['read'])
                record['mei_classes'].update(classes)
        record['mei_partner_frac'] = len(partner_hits) / float(record['nreads'])
        if partner_hits and record['mei_partner_frac'] >= mei_partner_frac:
            record['filters'].append('MEI_PARTNER')
    return record


def summarize_partner_group(ref, ori, lines):
    positions = [line['partner_pos'] for line in lines]
    mapqs = [line['partner_mapq'] for line in lines if line['partner_mapq'] is not None]
    lengths = [line['partner_aligned_length'] for line in lines if line['partner_aligned_length'] is not None]
    # pos_med, local_pos: median partner and clipped-read breakends, for bigclipper_cohort.py (not in the VCF)
    return {'ref': ref, 'pos1': min(positions), 'pos2': max(positions), 'ori': ori, 'count': len(lines),
            'mapq': median(mapqs), 'aligned_length': median(lengths),
            'pos_med': median(positions), 'local_pos': median([line['pos'] for line in lines])}


def median(values):
    if not values:
        return None
    values = sorted(values)
    return values[len(values) // 2]


class MEIIndex:
    """Positions of MEI calls, per contig, for windowed lookups."""

    def __init__(self, window):
        self.window = window
        self.intervals = defaultdict(list)  # contig -> [(start, end, class)], 0-based breakend coordinates
        self.n_calls = 0
        self.sources = []

    def add(self, contig, start, end, mei_class):
        self.intervals[contig].append((start, end, mei_class))
        self.n_calls += 1

    def finalize(self):
        self.starts = {}
        self.max_len = {}
        for contig, ivs in self.intervals.items():
            ivs.sort()
            self.starts[contig] = [iv[0] for iv in ivs]
            self.max_len[contig] = max(iv[1] - iv[0] for iv in ivs)

    def query(self, contig, start, end):
        """Return the set of MEI classes of calls within self.window of [start, end]."""
        ivs = self.intervals.get(contig)
        if not ivs:
            return set()
        lo, hi = start - self.window, end + self.window
        starts = self.starts[contig]
        i = bisect.bisect_right(starts, hi)
        classes = set()
        min_start = lo - self.max_len[contig]
        while i > 0 and starts[i - 1] >= min_start:
            i -= 1
            iv_start, iv_end, mei_class = ivs[i]
            if iv_end >= lo:
                classes.add(mei_class)
        return classes

    def load_vcf(self, path, info_key, classes, exclude_svtypes, sample, all_sites):
        import pysam
        start_time = time.time()
        vcf = pysam.VariantFile(path)
        if info_key not in vcf.header.info:
            sys.exit("Error: INFO/%s is not defined in the header of %s (set --mei_info_key)" % (info_key, path))
        samples = list(vcf.header.samples)
        use_sample = None
        if all_sites:
            vcf.subset_samples([])
        elif sample is not None:
            if sample not in samples:
                sys.exit("Error: sample '%s' not found in %s (samples: %s)" % (sample, path, ', '.join(samples[:10]) + (' ...' if len(samples) > 10 else '')))
            vcf.subset_samples([sample])
            use_sample = sample
        elif len(samples) == 1:
            use_sample = samples[0]
        elif len(samples) == 0:
            pass  # sites-only VCF
        else:
            sys.exit("Error: %s has %d samples; set --mei_sample to the sample matching this BAM, or --mei_all_sites to use every MEI site" % (path, len(samples)))
        excluded_svtypes = set() if exclude_svtypes.upper() in ('NONE', '') else set(s.strip().upper() for s in exclude_svtypes.split(','))
        wanted_classes = None if not classes else set(c.strip().upper() for c in classes.split(','))

        n_before = self.n_calls
        n_tagged = 0
        skipped_class = Counter()
        skipped_svtype = Counter()
        n_skipped_gt = 0
        used_classes = Counter()
        for rec in vcf:
            mei_class = rec.info.get(info_key)
            values = mei_class if isinstance(mei_class, tuple) else (mei_class,)
            values = [str(v) for v in values if v is not None and str(v) not in MISSING_VALUES]
            if not values:
                continue
            n_tagged += 1
            mei_class = ','.join(values)
            if wanted_classes is not None and not any(v.upper() in wanted_classes for v in values):
                skipped_class[mei_class] += 1
                continue
            svtype = get_svtype(rec)
            if svtype in excluded_svtypes:
                skipped_svtype[svtype] += 1
                continue
            if use_sample is not None:
                gt = rec.samples[use_sample].get('GT')
                if not gt or not any(a is not None and a > 0 for a in gt):
                    n_skipped_gt += 1
                    continue
            if svtype in ('INS', 'BND', 'TRA') or rec.stop <= rec.pos:
                # for an insertion, POS is the base before the insertion, i.e. 0-based breakend position POS
                self.add(rec.chrom, rec.pos, rec.pos, mei_class)
            else:
                self.add(rec.chrom, rec.pos, rec.stop, mei_class)
            used_classes[mei_class] += 1
        vcf.close()
        self.finalize()
        n_used = self.n_calls - n_before
        criteria = ["SVTYPE not in %s" % ','.join(sorted(excluded_svtypes))] if excluded_svtypes else []
        if wanted_classes is not None:
            criteria.append("class in %s" % ','.join(sorted(wanted_classes)))
        criteria.append("non-ref in sample %s" % use_sample if use_sample else "all sites")
        desc = "%s: %d records with INFO/%s, %d used (%s)" % (path, n_tagged, info_key, n_used, '; '.join(criteria))
        self.sources.append(desc)
        sys.stderr.write("MEI calls from %s [%.1fs]\n" % (desc, time.time() - start_time))
        sys.stderr.write("  used, by class: %s\n" % format_tally(used_classes))
        if skipped_class:
            sys.stderr.write("  skipped, class not in --mei_classes: %s\n" % format_tally(skipped_class))
        if skipped_svtype:
            sys.stderr.write("  skipped, by SVTYPE: %s\n" % format_tally(skipped_svtype))
        if use_sample is not None:
            sys.stderr.write("  skipped, not carried by %s: %d\n" % (use_sample, n_skipped_gt))
        if n_used == 0:
            sys.stderr.write("WARNING: no MEI calls were loaded from %s; MEI filtering will have no effect\n" % path)

    def load_bed(self, path):
        import gzip
        opener = gzip.open if path.endswith('.gz') else open
        n_before = self.n_calls
        with opener(path, 'rt') as f:
            for line in f:
                if not line.strip() or line.startswith(('#', 'track', 'browser')):
                    continue
                fields = line.rstrip('\n').split('\t')
                mei_class = fields[3] if len(fields) > 3 and fields[3] else 'MEI'
                self.add(fields[0], int(fields[1]), int(fields[2]), mei_class)
        self.finalize()
        n_used = self.n_calls - n_before
        self.sources.append("%s: %d intervals" % (path, n_used))
        sys.stderr.write("MEI calls from %s: %d intervals\n" % (path, n_used))

    def check_contigs(self, intermediate_contigs):
        if self.n_calls and not (set(self.intervals) & set(intermediate_contigs)):
            sys.stderr.write("WARNING: none of the MEI call contigs (e.g. %s) match the intermediate file's contigs (e.g. %s); "
                             "check that both use the same reference and naming\n" % (
                                 ', '.join(sorted(self.intervals)[:3]), ', '.join(sorted(intermediate_contigs)[:3])))


MISSING_VALUES = ('', '.', 'NA', 'N/A', 'None', 'none', 'NONE')


def format_tally(counter, top=10):
    items = counter.most_common()
    text = ', '.join('%s=%d' % kv for kv in items[:top])
    return text + (', ... (%d more)' % (len(items) - top) if len(items) > top else '') if items else 'none'


def get_svtype(rec):
    svtype = rec.info.get('SVTYPE') if 'SVTYPE' in rec.header.info else None
    if isinstance(svtype, tuple):
        svtype = svtype[0] if svtype else None
    if svtype:
        return str(svtype).upper().split(':')[0]
    alt = rec.alts[0] if rec.alts else ''
    if alt.startswith('<'):
        return alt[1:-1].upper().split(':')[0]
    if len(alt) > len(rec.ref):
        return 'INS'
    if len(alt) < len(rec.ref):
        return 'DEL'
    return None


def write_header(outfile, contigs, args, mei_index):
    outfile.write("##fileformat=VCFv4.2\n")
    outfile.write("##source=BigClipper_v%s\n" % VERSION)
    outfile.write("##bigclipperCommand=%s\n" % ' '.join(sys.argv))
    outfile.write("##FILTER=<ID=PASS,Description=\"All filters passed\">\n")
    outfile.write("##FILTER=<ID=MEI_SITE,Description=\"Cluster is within %d bp of an MEI call (likely an MEI insertion site)\">\n" % args.mei_window)
    outfile.write("##FILTER=<ID=MEI_PARTNER,Description=\"At least %g of supporting reads have their partner breakend within %d bp of an MEI call (likely the reference copy of an inserted element)\">\n" % (args.mei_partner_frac, args.mei_window))
    outfile.write("##ALT=<ID=BND,Description=\"Breakend: pile-up of clipped reads whose clipped portions align elsewhere (see ALNS)\">\n")
    outfile.write("##INFO=<ID=ALNS,Number=.,Type=String,Description=\"Supplementary alignments for clipped reads, as chrom:start-end|ori|count. "
                  "ori is the same as the cluster's ORI if the clipped portion aligns in the same orientation as the clipped read, and the opposite if it is inverted relative to it\">\n")
    outfile.write("##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length of the SV\">\n")
    outfile.write("##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of the SV\">\n")
    outfile.write("##INFO=<ID=ORI,Number=1,Type=String,Description=\"Orientation of the clipped reads: + if clipped on the right (aligned sequence lies left of the breakend), - if clipped on the left\">\n")
    outfile.write("##INFO=<ID=NREADS,Number=1,Type=Integer,Description=\"Number of reads in the cluster whose partner breakend passes the min_dist filter\">\n")
    outfile.write("##INFO=<ID=CIPOS,Number=2,Type=Integer,Description=\"Range of clipped-read breakend positions in the cluster, relative to POS\">\n")
    outfile.write("##INFO=<ID=PARTNER_MAPQ,Number=.,Type=Integer,Description=\"Median MAPQ of the supplementary alignments in each ALNS region (same order as ALNS)\">\n")
    outfile.write("##INFO=<ID=PARTNER_ALNLEN,Number=.,Type=Integer,Description=\"Median aligned read length of the supplementary alignments in each ALNS region (same order as ALNS)\">\n")
    outfile.write("##INFO=<ID=MEI_CLASS,Number=.,Type=String,Description=\"Classes of MEI calls matched at the cluster or its partner breakends\">\n")
    outfile.write("##INFO=<ID=MEI_PARTNER_FRAC,Number=1,Type=Float,Description=\"Fraction of supporting reads whose partner breakend is within the MEI window of an MEI call\">\n")
    outfile.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
    if mei_index is not None:
        for source in mei_index.sources:
            outfile.write("##bigclipperMEICalls=%s\n" % source)
    for contig in contigs:
        outfile.write("##contig=<ID=%s>\n" % contig)
    outfile.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")


def print_cluster(record, outfile):
    # VCF format is 1-based!
    info = ["SVTYPE=BND", "ORI=%s" % record['ori'], "NREADS=%d" % record['nreads'], "CIPOS=0,%d" % record['span'],
            "ALNS=" + ','.join("%s:%d-%d|%s|%d" % (hit['ref'], hit['pos1'] + 1, hit['pos2'] + 1, hit['ori'], hit['count'])
                               for hit in record['cluster_info'])]
    if any(hit['mapq'] is not None for hit in record['cluster_info']):
        info.append("PARTNER_MAPQ=" + ','.join(fmt_opt(hit['mapq']) for hit in record['cluster_info']))
        info.append("PARTNER_ALNLEN=" + ','.join(fmt_opt(hit['aligned_length']) for hit in record['cluster_info']))
    if record['mei_classes']:
        info.append("MEI_CLASS=" + ','.join(sorted(c.replace(',', '|').replace(';', '|').replace(' ', '_') for c in record['mei_classes'])))
    if record['mei_partner_frac'] is not None:
        info.append("MEI_PARTNER_FRAC=%.3g" % record['mei_partner_frac'])
    filters = ';'.join(record['filters']) if record['filters'] else 'PASS'
    outfile.write("%s\t%d\t%d\tN\t<BND>\t%d\t%s\t%s\n" % (record['ref'], record['pos'] + 1, record['ID'], record['qual'], filters, ';'.join(info)))


def fmt_opt(value):
    return '.' if value is None else str(value)


if __name__ == "__main__": main()
