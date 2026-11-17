#!/usr/bin/env python3
"""
BigClipper cohort integration and filtering. Subcommands, in the order wdl/CohortClipper.wdl runs them:

  check     check sample IDs and file names against the PED file
  clusters  phase 2 clusters for a batch of samples, as cluster rows and junction rows (one per partner group)
  sites     cohort breakend sites and junctions
  scan      each sample's split-read and unaligned-clip reads at every site
  tables    site table, evidence table, sample stats and cohort VCF
  filter    rare breakends called in affected samples, with parental evidence

A breakend site (chrom, ORI, position) opens at a cluster and takes the same-ORI clusters less than --window bp
to its right, without chaining. Junctions from different samples match if ORI, partner chrom and partner ori are
equal and both positions are less than --window bp from a representative (the one with the most reads).

Positions are 1-based breakend bases as in VCF breakend notation: the last aligned base for '+' (clipped on the
right), the first for '-'.
"""
import sys
import os
import io
import gzip
import time
import bisect
import heapq
import shutil
import argparse
import itertools
import subprocess
import tempfile
from array import array
from collections import defaultdict, namedtuple, Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bigclipper_getclusters as gc  # noqa: E402

VERSION = "2.2"
FLIP = {'+': '-', '-': '+'}

SITES_COLS = ['site_id', 'chrom_key', 'chrom', 'pos', 'ori', 'start', 'end', 'n_clusters', 'n_called', 'n_called_mei',
              'mei_classes', 'top_partner_chrom', 'top_partner_pos', 'top_partner_ori', 'top_partner_samples']
CALLS_COLS = ['site_id', 'sample', 'called_reads', 'n_clusters', 'filter', 'mei_classes', 'pos', 'top_partner']
JUNCTION_COLS = ['junction_id', 'site_id', 'chrom_key', 'chrom', 'pos', 'ori', 'partner_chrom', 'partner_pos',
                 'partner_ori', 'n_samples', 'n_reads', 'n_mei', 'samples', 'mate_id']
SCAN_COLS = ['site_id', 'sample', 'SR', 'CL', 'R', 'H1', 'H2', 'HU']
EVIDENCE_COLS = ['site_id', 'sample', 'called_reads', 'filter', 'SR', 'CL', 'R', 'H1', 'H2', 'HU', 'partner', 'mate_site']
SITE_TABLE_COLS = ['site_id', 'chrom', 'pos', 'ori', 'start', 'end', 'n_clusters', 'n_called', 'n_called_families',
                   'n_called_affected', 'n_called_mei', 'n_called_nearby', 'n_carriers', 'AF', 'AF_called',
                   'n_samples_reads1', 'n_samples_reads2', 'n_samples_split1', 'n_samples_clip1', 'n_families_reads1',
                   'n_junctions', 'top_partner_chrom', 'top_partner_pos', 'top_partner_ori', 'top_partner_samples',
                   'mate_site', 'mei_classes', 'mei_near', 'regions', 'partner_regions']
SAMPLE_STATS_COLS = ['sample', 'family', 'father', 'mother', 'affected', 'n_clusters', 'n_clusters_mei', 'n_junctions',
                     'n_split_lines', 'n_split_reads', 'n_clip_lines', 'phased_frac', 'rel_split_lines', 'rel_clip_lines']
CANDIDATE_COLS = ['site_id', 'chrom', 'pos', 'ori', 'start', 'end', 'sample', 'family', 'called_reads', 'filter',
                  'SR', 'CL', 'R', 'H1', 'H2', 'HU', 'partner', 'mate_site', 'AF', 'n_carriers', 'AF_called', 'n_called',
                  'n_called_families', 'n_called_nearby', 'n_samples_reads1', 'parental_evidence',
                  'father', 'father_R', 'father_SR', 'father_CL', 'father_mate_R', 'father_rel_yield',
                  'mother', 'mother_R', 'mother_SR', 'mother_CL', 'mother_mate_R', 'mother_rel_yield', 'relatives_with_reads',
                  'top_partner_chrom', 'top_partner_pos', 'top_partner_ori', 'top_partner_samples',
                  'n_junctions', 'mei_classes', 'mei_near', 'regions', 'partner_regions']

Junction = namedtuple('Junction', 'chrom pos ori pchrom ppos pori sample cluster reads filter')


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0].strip(), formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--version', action='version', version='BigClipper cohort %s' % VERSION)
    sub = parser.add_subparsers(dest='command', required=True)

    p = sub.add_parser('check', help='check the sample list against the PED file')
    p.add_argument('--samples', required=True, help='file with one sample ID per line')
    p.add_argument('--ped', required=True)
    p.add_argument('--intermediate_list', help='intermediate file paths, one per line, in sample order')
    p.add_argument('--clips_list', help='unaligned-clips file paths, one per line, in sample order; each must have the same '
                   'prefix as the intermediate file on its line (<prefix>_intermediate.bed, <prefix>_unaligned_clips.bed)')
    p.add_argument('--n_intermediate', type=int, help='number of intermediate files given (must equal the number of samples)')
    p.add_argument('--n_clips', type=int, help='number of unaligned-clips files given (must equal the number of samples)')
    p.set_defaults(func=cmd_check)

    p = sub.add_parser('clusters', help='phase 2 clusters for a batch of samples')
    p.add_argument('--manifest', required=True, help='TSV: sample, intermediate file (one sample per line)')
    p.add_argument('--out_prefix', required=True)
    p.add_argument('-d', '--min_dist', type=int, default=5000, help='as in bigclipper_getclusters.py [default=5000]')
    p.add_argument('-c', '--min_cluster_count', type=int, default=5, help='[default=5]')
    p.add_argument('-s', '--cluster_dist', type=int, default=50, help='[default=50]')
    p.add_argument('-u', '--max_unique_breakends', type=int, default=10, help='[default=10]')
    p.add_argument('--junction_min_reads', type=int, default=2, help='write a junction row for partner groups with at least this many reads [default=2]')
    add_mei_args(p)
    p.add_argument('--mei_partner_frac', type=float, default=0.5, help='[default=0.5]')
    p.set_defaults(func=cmd_clusters)

    p = sub.add_parser('sites', help='cohort breakend sites and junctions')
    p.add_argument('--clusters', nargs='+', required=True, help='cluster rows from the clusters subcommand')
    p.add_argument('--junctions', nargs='+', required=True, help='junction rows from the clusters subcommand')
    p.add_argument('--window', type=int, default=50, help='site and junction matching distance (bp) [default=50]')
    p.add_argument('--out_prefix', required=True)
    add_sort_args(p)
    p.set_defaults(func=cmd_sites)

    p = sub.add_parser('scan', help="count each sample's reads at every site")
    p.add_argument('--sites', required=True, help='sites file from the sites subcommand')
    p.add_argument('--manifest', required=True, help='TSV: sample, intermediate file, unaligned-clips file')
    p.add_argument('--window', type=int, default=50, help='count reads less than this many bp outside the range of a site [default=50]')
    p.add_argument('--out_prefix', required=True)
    p.add_argument('--tmpdir', default=None)
    p.set_defaults(func=cmd_scan)

    p = sub.add_parser('tables', help='site table, evidence table, sample stats and cohort VCF')
    p.add_argument('--samples', required=True, help='cohort sample IDs, one per line (VCF column order)')
    p.add_argument('--ped', required=True)
    p.add_argument('--sites', required=True)
    p.add_argument('--calls', required=True)
    p.add_argument('--junctions', required=True)
    p.add_argument('--evidence', nargs='+', required=True, help='evidence files from the scan subcommand')
    p.add_argument('--cluster_stats', nargs='+', required=True)
    p.add_argument('--scan_stats', nargs='+', required=True)
    p.add_argument('--region_bed', action='append', default=[], help='BED of regions to flag (name from the file name, or NAME=PATH); can be given more than once')
    p.add_argument('--fai', default=None, help='reference .fai, for contig lengths in the VCF header')
    p.add_argument('--window', type=int, default=50, help='[default=50]')
    p.add_argument('--af_min_reads', type=int, default=2, help='count a sample as a carrier if it was called or has at least this many reads at the site [default=2]')
    add_mei_args(p)
    p.add_argument('--out_prefix', required=True)
    p.add_argument('--threads', type=int, default=2)
    p.set_defaults(func=cmd_tables)

    p = sub.add_parser('filter', help='rare breakends called in affected samples, with parental evidence')
    p.add_argument('--site_table', required=True)
    p.add_argument('--evidence', required=True, help='evidence table from the tables subcommand')
    p.add_argument('--sample_stats', required=True)
    p.add_argument('--ped', required=True)
    p.add_argument('--max_af', type=float, default=0.01, help='keep breakends with cohort AF at most this [default=0.01]')
    p.add_argument('--af_min_reads', type=int, default=2, help='count a sample as a carrier if it was called or has at least this many reads [default=2]')
    p.add_argument('--exclude_mei', action='store_true', help="drop candidates whose cluster is MEI-flagged in the candidate's sample")
    p.add_argument('--out_prefix', required=True)
    p.set_defaults(func=cmd_filter)

    args = parser.parse_args()
    args.func(args)


def add_mei_args(p):
    p.add_argument('--mei_vcf', default=None, help='cohort SV VCF/BCF with MEI-annotated calls')
    p.add_argument('--mei_info_key', default='MEI_type', help='[default=MEI_type]')
    p.add_argument('--mei_classes', default=None, help='comma-separated values of --mei_info_key to treat as MEIs [default: any]')
    p.add_argument('--mei_exclude_svtypes', default='DEL', help='[default=DEL]')
    p.add_argument('--mei_window', type=int, default=100, help='[default=100]')


def add_sort_args(p):
    p.add_argument('--tmpdir', default=None)
    p.add_argument('--sort_mem', default='2G', help='memory for GNU sort -S [default=2G]')
    p.add_argument('--threads', type=int, default=2)


##########
# helpers

def log(msg):
    sys.stderr.write("[%s] %s\n" % (time.strftime('%H:%M:%S'), msg))
    sys.stderr.flush()


def open_in(path):
    """Open a text file, plain or (b)gzipped."""
    if path == '-':
        return sys.stdin
    with open(path, 'rb') as f:
        magic = f.read(2)
    return gzip.open(path, 'rt') if magic == b'\x1f\x8b' else open(path)


class Out:
    """Text output; a .gz path is block-gzipped with bgzip if available (so tabix can index it), else gzipped."""

    def __init__(self, path, threads=1):
        self.path = path
        self.proc = None
        if path.endswith('.gz') and shutil.which('bgzip'):
            self.raw = open(path, 'wb')
            self.proc = subprocess.Popen(['bgzip', '-c', '-@', str(max(1, threads))], stdin=subprocess.PIPE, stdout=self.raw)
            self.handle = io.TextIOWrapper(self.proc.stdin, encoding='utf-8')
        elif path.endswith('.gz'):
            self.handle = gzip.open(path, 'wt', compresslevel=4)
        else:
            self.handle = open(path, 'w')
        self.write = self.handle.write

    def close(self):
        self.handle.close()
        if self.proc is not None:
            if self.proc.wait() != 0:
                sys.exit("Error: bgzip failed writing %s" % self.path)
            self.raw.close()


def read_rows(path):
    """Rows of a TSV, skipping '#' lines."""
    with open_in(path) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            yield line.rstrip('\n').split('\t')


def read_header(path, expected=None):
    with open_in(path) as f:
        first = f.readline()
    if not first.startswith('#'):
        sys.exit("Error: %s has no header line" % path)
    cols = first[1:].rstrip('\n').split('\t')
    if expected is not None and cols[:len(expected)] != expected:
        sys.exit("Error: %s does not have the expected columns (%s)" % (path, ', '.join(expected[:5]) + ' ...'))
    return {c: i for i, c in enumerate(cols)}


def read_list(path):
    with open_in(path) as f:
        return [line.rstrip('\n').rstrip('\r') for line in f]


def read_manifest(path):
    rows = []
    for fields in read_rows(path):
        if len(fields) < 2:
            sys.exit("Error: manifest %s needs at least 2 columns (sample, intermediate file)" % path)
        rows.append(fields)
    return rows


def chrom_key(chrom):
    """Sort key (a string, for GNU sort) putting chr1..chr22, X, Y, M first, then other contigs by name."""
    c = chrom[3:] if chrom[:3].lower() == 'chr' else chrom
    if c.isdigit():
        return "%04d" % int(c)
    special = {'X': 1001, 'Y': 1002, 'M': 1003, 'MT': 1003}
    if c.upper() in special:
        return "%04d" % special[c.upper()]
    return "9999:" + chrom


def local_bp(pos0, ori):
    """processbam's 0-based breakend position -> 1-based breakend base (last aligned base for '+', first for '-')."""
    return pos0 + 1 if ori == '-' else pos0


def partner_bp(pos0, pori):
    """processbam's partner position -> 1-based breakend base of the partner alignment (first aligned base for '+',
    last for '-')."""
    return pos0 + 1 if pori == '+' else pos0


def bnd_alt(ori, pchrom, ppos, pori, ref='N'):
    """VCF breakend ALT for a junction: ++ t[p[, +- t]p], -- ]p]t, -+ [p[t."""
    br = '[' if pori == '+' else ']'
    mate = "%s%s:%d%s" % (br, pchrom, ppos, br)
    return ref + mate if ori == '+' else mate + ref


def parse_bnd_alt(alt):
    """Inverse of bnd_alt: (ori, pchrom, ppos, pori)."""
    br = '[' if '[' in alt else ']'
    first, mate, last = alt.split(br)
    pchrom, ppos = mate.rsplit(':', 1)
    return ('+' if first else '-'), pchrom, int(ppos), ('+' if br == '[' else '-')


def median_int(values):
    values = sorted(values)
    return values[len(values) // 2]


def fmt_float(x):
    return "%.4g" % x


def gnu_sort():
    for exe in ('gsort', 'sort'):
        path = shutil.which(exe)
        if path:
            try:
                version = subprocess.run([path, '--version'], capture_output=True, text=True).stdout
            except OSError:
                continue
            if 'GNU' in version:
                return path
    sys.exit("Error: GNU sort is required")


def sort_rows(paths, keys, tmpdir, mem, threads):
    """Concatenate (gzipped) TSVs without headers and sort them with GNU sort into a temporary file."""
    out = tempfile.NamedTemporaryFile(dir=tmpdir, suffix='.tsv', delete=False)
    out.close()
    cmd = [gnu_sort(), '-t', '\t'] + keys + ['-S', mem, '--parallel=%d' % threads, '-T', tmpdir, '-o', out.name]
    proc = subprocess.Popen(cmd, stdin=subprocess.PIPE, env=dict(os.environ, LC_ALL='C'))
    for path in paths:
        with open(path, 'rb') as raw:
            magic = raw.read(2)
            raw.seek(0)
            src = gzip.GzipFile(fileobj=raw) if magic == b'\x1f\x8b' else raw
            shutil.copyfileobj(src, proc.stdin, 1 << 20)
    proc.stdin.close()
    if proc.wait() != 0:
        sys.exit("Error: '%s' failed" % ' '.join(cmd))
    return out.name


class Pedigree:
    """6-column PED: family, sample, father, mother, sex, phenotype (2 = affected)."""

    def __init__(self, path):
        self.family, self.father, self.mother, self.affected = {}, {}, {}, {}
        with open_in(path) as f:
            for n, line in enumerate(f, 1):
                if not line.strip() or line.startswith('#'):
                    continue
                fields = line.split()
                if len(fields) < 6:
                    sys.exit("Error: line %d of PED file %s has fewer than 6 columns" % (n, path))
                fam, sample, father, mother, _, phenotype = fields[:6]
                if not self.family and phenotype not in ('-9', '0', '1', '2'):
                    continue  # header line
                if sample in self.family:
                    sys.exit("Error: sample %s appears twice in PED file %s" % (sample, path))
                self.family[sample] = fam
                self.father[sample] = None if father in ('0', '.', '') else father
                self.mother[sample] = None if mother in ('0', '.', '') else mother
                self.affected[sample] = phenotype == '2'

    def members(self, samples):
        """family -> samples (restricted to `samples`)."""
        fam = defaultdict(list)
        for s in samples:
            if s in self.family:
                fam[self.family[s]].append(s)
        return fam


class Grouped:
    """Consecutive rows with the same site ID (integer column `col`), taken in increasing order."""

    def __init__(self, rows, name, col=0):
        self.it = itertools.groupby(rows, key=lambda r: int(r[col]))
        self.name = name
        self.last = -1
        self._advance()

    def _advance(self):
        kv = next(self.it, None)
        self.cur = None if kv is None else (kv[0], list(kv[1]))
        if self.cur is not None:
            if self.cur[0] <= self.last:
                sys.exit("Error: %s is not sorted by site ID (site %d after %d)" % (self.name, self.cur[0], self.last))
            self.last = self.cur[0]

    def take(self, site_id):
        if self.cur is not None and self.cur[0] < site_id:
            sys.exit("Error: %s has rows for site %d, which is not in the site list" % (self.name, self.cur[0]))
        if self.cur is not None and self.cur[0] == site_id:
            rows = self.cur[1]
            self._advance()
            return rows
        return []

    def check_done(self):
        if self.cur is not None:
            sys.exit("Error: %s has rows for site %d, which is not in the site list" % (self.name, self.cur[0]))


class SiteIndex:
    """Sites by (chrom, ori), for finding the sites whose range [start, end], widened by `margin` on each side,
    contains a position."""

    def __init__(self, path, margin):
        col = read_header(path, SITES_COLS)
        c_chrom, c_ori, c_start, c_end = col['chrom'], col['ori'], col['start'], col['end']
        c_id, c_pos, c_called = col['site_id'], col['pos'], col['n_called']
        self.lo, self.hi, self.ids, self.pos, self.n_called, self.maxw = {}, {}, {}, {}, {}, {}
        self.chroms = set()
        for r in read_rows(path):
            key = (r[c_chrom], r[c_ori])
            lo = self.lo.get(key)
            if lo is None:
                lo = self.lo[key] = array('q')
                self.hi[key], self.ids[key], self.pos[key], self.n_called[key] = array('q'), array('q'), array('q'), array('q')
                self.maxw[key] = 0
                self.chroms.add(key[0])
            start, end = int(r[c_start]) - margin, int(r[c_end]) + margin
            if lo and start < lo[-1]:
                sys.exit("Error: %s is not sorted by position (site %s)" % (path, r[c_id]))
            lo.append(start)
            self.hi[key].append(end)
            self.ids[key].append(int(r[c_id]))
            self.pos[key].append(int(r[c_pos]))
            self.n_called[key].append(int(r[c_called]))
            if end - start > self.maxw[key]:
                self.maxw[key] = end - start
        self.n_sites = sum(len(v) for v in self.lo.values())

    def _overlapping(self, key, start, end):
        """Indices of sites in key whose widened range overlaps [start, end]."""
        lo = self.lo.get(key)
        if lo is None:
            return []
        hi = self.hi[key]
        i = bisect.bisect_right(lo, end)
        floor = start - self.maxw[key]
        out = []
        while i > 0 and lo[i - 1] >= floor:
            i -= 1
            if hi[i] >= start:
                out.append(i)
        return out

    def find(self, chrom, ori, pos):
        key = (chrom, ori)
        return [self.ids[key][i] for i in self._overlapping(key, pos, pos)]

    def nearest(self, chrom, ori, pos):
        key = (chrom, ori)
        hits = self._overlapping(key, pos, pos)
        if not hits:
            return None
        return min((abs(self.pos[key][i] - pos), self.ids[key][i]) for i in hits)[1]

    def called_nearby(self, chrom, ori, start, end, site_id, extra):
        """Called samples (summed over sites) at other sites of the same orientation whose widened range is
        within `extra` bp of [start, end]."""
        key = (chrom, ori)
        return sum(self.n_called[key][i] for i in self._overlapping(key, start - extra, end + extra) if self.ids[key][i] != site_id)


class RegionSet:
    """Merged BED intervals, for overlap queries with 1-based closed ranges."""

    def __init__(self, spec):
        if '=' in spec and not os.path.exists(spec):
            self.name, path = spec.split('=', 1)
        else:
            path = spec
            name = os.path.basename(path)
            for ext in ('.gz', '.bgz', '.bed'):
                if name.endswith(ext):
                    name = name[:-len(ext)]
            self.name = name
        ivs = defaultdict(list)
        with open_in(path) as f:
            for line in f:
                if not line.strip() or line.startswith(('#', 'track', 'browser')):
                    continue
                fields = line.split('\t')
                ivs[fields[0]].append((int(fields[1]), int(fields[2])))
        self.starts, self.ends = {}, {}
        for chrom, lst in ivs.items():
            lst.sort()
            merged = [list(lst[0])]
            for s, e in lst[1:]:
                if s <= merged[-1][1]:
                    merged[-1][1] = max(merged[-1][1], e)
                else:
                    merged.append([s, e])
            self.starts[chrom] = array('q', (m[0] for m in merged))
            self.ends[chrom] = array('q', (m[1] for m in merged))

    def overlaps(self, chrom, start, end):
        starts = self.starts.get(chrom)
        if starts is None:
            return False
        i = bisect.bisect_left(starts, end)  # intervals starting before base `end` (0-based end - 1)
        return i > 0 and self.ends[chrom][i - 1] > start - 1


def load_mei_batch(path, samples, info_key, classes, exclude_svtypes, window):
    """Per-sample MEI indexes (calls where the sample is non-reference) from one pass over a cohort VCF."""
    import pysam
    start_time = time.time()
    vcf = pysam.VariantFile(path)
    if info_key not in vcf.header.info:
        sys.exit("Error: INFO/%s is not defined in the header of %s (set --mei_info_key)" % (info_key, path))
    in_vcf = set(vcf.header.samples)
    missing = [s for s in samples if s not in in_vcf]
    if missing:
        sys.exit("Error: %d sample(s) are not in the MEI VCF %s, e.g. %s" % (len(missing), path, ', '.join(missing[:5])))
    vcf.subset_samples(samples)
    excluded = set() if exclude_svtypes.upper() in ('NONE', '') else set(s.strip().upper() for s in exclude_svtypes.split(','))
    wanted = None if not classes else set(c.strip().upper() for c in classes.split(','))
    index = {s: gc.MEIIndex(window) for s in samples}
    n_used = 0
    for rec in vcf:
        value = rec.info.get(info_key)
        values = value if isinstance(value, tuple) else (value,)
        values = [str(v) for v in values if v is not None and str(v) not in gc.MISSING_VALUES]
        if not values:
            continue
        if wanted is not None and not any(v.upper() in wanted for v in values):
            continue
        svtype = gc.get_svtype(rec)
        if svtype in excluded:
            continue
        mei_class = ','.join(values)
        if svtype in ('INS', 'BND', 'TRA') or rec.stop <= rec.pos:
            start, end = rec.pos, rec.pos
        else:
            start, end = rec.pos, rec.stop
        n_used += 1
        for s in samples:
            gt = rec.samples[s].get('GT')
            if gt and any(a is not None and a > 0 for a in gt):
                index[s].add(rec.chrom, start, end, mei_class)
    vcf.close()
    for idx in index.values():
        idx.finalize()
    log("MEI calls from %s: %d used, %d sample-calls for %d samples [%.1fs]" % (
        path, n_used, sum(i.n_calls for i in index.values()), len(samples), time.time() - start_time))
    return index


def clean_class(c):
    return c.replace(',', '|').replace(';', '|').replace(' ', '_')


##########
# check

def file_prefix(path, suffixes):
    name = os.path.basename(path.rstrip('/'))
    for suffix in suffixes:
        if name.endswith(suffix):
            return name[:-len(suffix)]
    return None


def cmd_check(args):
    ids = read_list(args.samples)
    while ids and ids[-1] == '':
        ids.pop()
    problems = []
    inter = [p for p in read_list(args.intermediate_list) if p] if args.intermediate_list else None
    clips = [p for p in read_list(args.clips_list) if p] if args.clips_list else None
    n_inter = len(inter) if inter is not None else args.n_intermediate
    n_clips = len(clips) if clips is not None else args.n_clips
    for name, n in (('intermediate files', n_inter), ('unaligned-clips files', n_clips)):
        if n is not None and n != len(ids):
            problems.append("%d sample IDs but %d %s" % (len(ids), n, name))
    if inter is not None and clips is not None and len(inter) == len(clips):
        # catch shifted arrays (one sample's clips paired with another's intermediate file)
        unpaired = []
        for i, (a, b) in enumerate(zip(inter, clips)):
            pa = file_prefix(a, ('_intermediate.bed.gz', '_intermediate.bed'))
            pb = file_prefix(b, ('_unaligned_clips.bed.gz', '_unaligned_clips.bed'))
            if pa is None or pb is None or pa != pb:
                unpaired.append("%s: %s / %s" % (ids[i] if i < len(ids) else i + 1, os.path.basename(a), os.path.basename(b)))
        if unpaired:
            problems.append("%d intermediate and unaligned-clips files do not share a <prefix>_intermediate.bed / "
                            "<prefix>_unaligned_clips.bed name, e.g. %s" % (len(unpaired), '; '.join(unpaired[:5])))
    bad = [s for s in ids if not s or s != s.strip() or '\t' in s]
    if bad:
        problems.append("%d empty or badly formed sample IDs" % len(bad))
    dups = sorted(s for s, c in Counter(ids).items() if c > 1)
    if dups:
        problems.append("duplicate sample IDs: %s" % ', '.join(dups[:10]))
    ped = Pedigree(args.ped)
    missing = [s for s in ids if s not in ped.family]
    if missing:
        problems.append("%d samples are not in the PED file, e.g. %s" % (len(missing), ', '.join(missing[:10])))
    if problems:
        sys.exit("Error: " + "; ".join(problems))
    cohort = set(ids)
    n_aff = sum(ped.affected[s] for s in ids)
    trios = [s for s in ids if ped.father[s] in cohort and ped.mother[s] in cohort]
    log("%d samples in %d families; %d affected; %d with both parents in the cohort (%d of them affected)" % (
        len(ids), len(set(ped.family[s] for s in ids)), n_aff, len(trios), sum(ped.affected[s] for s in trios)))


##########
# clusters

def cmd_clusters(args):
    samples = read_manifest(args.manifest)
    names = [s[0] for s in samples]
    mei = load_mei_batch(args.mei_vcf, names, args.mei_info_key, args.mei_classes, args.mei_exclude_svtypes,
                         args.mei_window) if args.mei_vcf else None
    fc = Out(args.out_prefix + '.clusters.tsv.gz')
    fj = Out(args.out_prefix + '.junctions.tsv.gz')
    with open(args.out_prefix + '.cluster_stats.tsv', 'w') as fs:
        fs.write('#sample\tn_clusters\tn_clusters_mei\tn_junctions\n')
        for fields in samples:
            sample, path = fields[0], fields[1]
            t0 = time.time()
            index = mei[sample] if mei is not None else None
            n_cl = n_mei = n_j = 0
            for cluster in gc.read_clusters(path, args.cluster_dist):
                if len(set(line['read'] for line in cluster)) < args.min_cluster_count:
                    continue
                rec = gc.process_cluster(cluster, args.min_dist, args.min_cluster_count, args.max_unique_breakends,
                                         args.cluster_dist, index, args.mei_partner_frac)
                if rec is None:
                    continue
                chrom, ori = rec['ref'], rec['ori']
                key = chrom_key(chrom)
                passing = [line for line in cluster
                           if line['partner_ref'] != chrom or abs(line['pos'] - line['partner_pos']) >= args.min_dist]
                points = sorted(local_bp(line['pos'], ori) for line in passing)
                filt = ';'.join(rec['filters']) or 'PASS'
                mei_classes = ','.join(sorted(clean_class(c) for c in rec['mei_classes'])) or '.'
                top = rec['cluster_info'][0]
                top_str = "%s|%d|%s|%d" % (top['ref'], partner_bp(top['pos_med'], top['ori']), top['ori'], top['count'])
                fc.write("%s\t%s\t%d\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%s\t%s\t%s\t%d\n" % (
                    key, chrom, points[len(points) // 2], ori, sample, rec['ID'], points[0], points[-1], rec['nreads'],
                    rec['qual'], filt, mei_classes, top_str, len(rec['cluster_info'])))
                n_cl += 1
                n_mei += bool(rec['filters'])
                for g in rec['cluster_info']:
                    if g['count'] < args.junction_min_reads:
                        continue
                    fj.write("%s\t%s\t%d\t%s\t%s\t%d\t%s\t%s\t%d\t%d\t%s\n" % (
                        key, chrom, local_bp(g['local_pos'], ori), ori, g['ref'], partner_bp(g['pos_med'], g['ori']),
                        g['ori'], sample, rec['ID'], g['count'], filt))
                    n_j += 1
            fs.write("%s\t%d\t%d\t%d\n" % (sample, n_cl, n_mei, n_j))
            log("%s: %d clusters (%d MEI-flagged), %d junctions [%.1fs]" % (sample, n_cl, n_mei, n_j, time.time() - t0))
    fc.close()
    fj.close()


##########
# sites

def make_sites(rows, window):
    """Sites from one chromosome's cluster rows, sorted by position. A site's range spans its clusters' reads.
    Returns [(start, end, ori, rows)] sorted by start."""
    sites = []
    for ori in ('+', '-'):
        first = None
        for r in rows:
            if r[3] != ori:
                continue
            p = int(r[2])
            if first is None or p - first >= window:
                first = p
                cur = [int(r[6]), int(r[7]), ori, []]
                sites.append(cur)
            cur[3].append(r)
            cur[0] = min(cur[0], int(r[6]))
            cur[1] = max(cur[1], int(r[7]))
    sites.sort(key=lambda s: (s[0], s[2], int(s[3][0][2])))
    return sites


def top_partner(by_sample, window):
    """Most common top partner of a site's samples: (chrom, pos, ori, n_samples), grouping partner positions of
    the same chrom and ori less than `window` bp to the right of a group's first position."""
    groups = defaultdict(list)
    for sample, d in by_sample.items():
        pchrom, ppos, pori, _ = d['top'].rsplit('|', 3)
        groups[(pchrom, pori)].append(int(ppos))
    best = None
    for (pchrom, pori), positions in sorted(groups.items()):
        positions.sort()
        i = 0
        while i < len(positions):
            j = i
            while j < len(positions) and positions[j] - positions[i] < window:
                j += 1
            cand = (j - i, pchrom, positions[(i + j - 1) // 2], pori)
            if best is None or cand[0] > best[0]:
                best = cand
            i = j
    return best[1], best[2], best[3], best[0]


def match_junctions(junctions, window):
    """Group one chromosome's junction records into cohort junctions. Representatives are taken by reads (most
    first), then position, so the result doesn't depend on input order. Returns [(representative, members)]."""
    by_key = defaultdict(list)
    for j in junctions:
        by_key[(j.ori, j.pchrom, j.pori)].append(j)
    out = []
    w = window
    for key in sorted(by_key):
        recs = by_key[key]
        # grid of w x w cells: records in one cell always match, so dense pileups stay linear
        cells = defaultdict(list)
        for i, r in enumerate(recs):
            cells[(r.pos // w, r.ppos // w)].append(i)
        order = sorted(range(len(recs)), key=lambda i: (-recs[i].reads, recs[i].pos, recs[i].ppos, recs[i].sample, recs[i].cluster))
        assigned = bytearray(len(recs))
        for i in order:
            if assigned[i]:
                continue
            assigned[i] = 1
            rep = recs[i]
            members = [rep]
            cx, cy = rep.pos // w, rep.ppos // w
            for dx in (-1, 0, 1):
                for dy in (-1, 0, 1):
                    cell_key = (cx + dx, cy + dy)
                    cell = cells.get(cell_key)
                    if not cell:
                        continue
                    keep = []
                    for k in cell:
                        if assigned[k]:
                            continue
                        r = recs[k]
                        if abs(r.pos - rep.pos) < w and abs(r.ppos - rep.ppos) < w:
                            assigned[k] = 1
                            members.append(r)
                        else:
                            keep.append(k)
                    cells[cell_key] = keep
            out.append((rep, members))
    return out


def find_mates(index, window):
    """index: {(chrom, ori, pchrom, pori): [(pos, ppos, id)]}. Returns {id: mate id}: the closest junction seen from
    the other side (ends swapped, orientations flipped) with both positions less than `window` bp away."""
    for lst in index.values():
        lst.sort()
    mates = {}
    for (chrom, ori, pchrom, pori), lst in index.items():
        other = index.get((pchrom, FLIP[pori], chrom, FLIP[ori]))
        if not other:
            continue
        other_pos = [o[0] for o in other]
        for pos, ppos, jid in lst:
            lo = bisect.bisect_right(other_pos, ppos - window)
            hi = bisect.bisect_left(other_pos, ppos + window)
            best = None
            for k in range(lo, hi):
                opos, oppos, oid = other[k]
                if oid != jid and abs(oppos - pos) < window:
                    d = (abs(opos - ppos) + abs(oppos - pos), oid)
                    if best is None or d < best:
                        best = d
            if best is not None:
                mates[jid] = best[1]
    return mates


def cmd_sites(args):
    w = args.window
    tmpdir = tempfile.mkdtemp(prefix='bigclipper_sites_', dir=args.tmpdir)
    t0 = time.time()
    clusters = sort_rows(args.clusters, ['-k1,1', '-k3,3n', '-k4,4', '-k5,5', '-k6,6n'], tmpdir, args.sort_mem, args.threads)
    junctions = sort_rows(args.junctions, ['-k1,1', '-k3,3n'], tmpdir, args.sort_mem, args.threads)
    log("sorted cluster and junction rows [%.1fs]" % (time.time() - t0))

    fs = Out(args.out_prefix + '.sites.tsv.gz', args.threads)
    fs.write('#' + '\t'.join(SITES_COLS) + '\n')
    fc = Out(args.out_prefix + '.calls.tsv.gz', args.threads)
    fc.write('#' + '\t'.join(CALLS_COLS) + '\n')
    jtmp = os.path.join(tmpdir, 'junctions_nomates.tsv')
    fj = open(jtmp, 'w')
    mate_index = defaultdict(list)

    def by_chrom(path):
        with open(path) as f:
            for key, group in itertools.groupby((line.rstrip('\n').split('\t') for line in f), key=lambda r: r[0]):
                yield key, group

    jit = by_chrom(junctions)
    jcur = next(jit, None)
    n_sites = n_calls = n_junctions = n_records = 0
    next_site = next_junction = 0
    for key, group in by_chrom(clusters):
        rows = list(group)
        chrom = rows[0][1]
        site_of = {}
        for start, end, ori, members in make_sites(rows, w):
            site_id = next_site
            next_site += 1
            by_sample = {}
            for r in members:
                d = by_sample.get(r[4])
                if d is None:
                    d = by_sample[r[4]] = {'reads': 0, 'n': 0, 'filters': set(), 'mei': set(), 'pos': int(r[2]), 'top': r[12], 'best': -1}
                n = int(r[8])
                d['reads'] += n
                d['n'] += 1
                if r[10] != 'PASS':
                    d['filters'].update(r[10].split(';'))
                if r[11] != '.':
                    d['mei'].update(r[11].split(','))
                if n > d['best']:
                    d['best'], d['pos'], d['top'] = n, int(r[2]), r[12]
                site_of[(r[4], r[5])] = site_id
            points = sorted(int(r[2]) for r in members)
            pchrom, ppos, pori, psamples = top_partner(by_sample, w)
            mei_classes = set()
            for d in by_sample.values():
                mei_classes |= d['mei']
            fs.write("%d\t%s\t%s\t%d\t%s\t%d\t%d\t%d\t%d\t%d\t%s\t%s\t%d\t%s\t%d\n" % (
                site_id, key, chrom, points[len(points) // 2], ori, start, end, len(members), len(by_sample),
                sum(1 for d in by_sample.values() if d['filters']), ','.join(sorted(mei_classes)) or '.',
                pchrom, ppos, pori, psamples))
            for sample in sorted(by_sample):
                d = by_sample[sample]
                fc.write("%d\t%s\t%d\t%d\t%s\t%s\t%d\t%s\n" % (
                    site_id, sample, d['reads'], d['n'], ';'.join(sorted(d['filters'])) or 'PASS',
                    ','.join(sorted(d['mei'])) or '.', d['pos'], d['top']))
                n_calls += 1
            n_sites += 1

        jrows = []
        while jcur is not None and jcur[0] < key:  # cannot happen: junctions come from clusters
            jcur = next(jit, None)
        if jcur is not None and jcur[0] == key:
            jrows = [Junction(r[1], int(r[2]), r[3], r[4], int(r[5]), r[6], r[7], r[8], int(r[9]), r[10]) for r in jcur[1]]
            jcur = next(jit, None)
        n_records += len(jrows)
        matched = match_junctions(jrows, w)
        matched.sort(key=lambda m: (m[0].pos, m[0].ori, chrom_key(m[0].pchrom), m[0].pchrom, m[0].ppos, m[0].pori))
        out = []
        for rep, members in matched:
            jid = next_junction
            next_junction += 1
            reads = Counter()
            mei_samples = set()
            for m in members:
                reads[m.sample] += m.reads
                if m.filter != 'PASS':
                    mei_samples.add(m.sample)
            samples_str = ','.join("%s:%d" % (s, reads[s]) for s in sorted(reads))
            site_id = site_of[(rep.sample, rep.cluster)]
            out.append((site_id, jid, "%d\t%d\t%s\t%s\t%d\t%s\t%s\t%d\t%s\t%d\t%d\t%d\t%s" % (
                jid, site_id, key, chrom, rep.pos, rep.ori, rep.pchrom, rep.ppos, rep.pori, len(reads),
                sum(reads.values()), len(mei_samples), samples_str)))
            mate_index[(chrom, rep.ori, rep.pchrom, rep.pori)].append((rep.pos, rep.ppos, jid))
        out.sort()
        for _, _, line in out:
            fj.write(line + '\n')
        n_junctions += len(out)
        log("%s: %d sites, %d junctions so far [%.1fs]" % (chrom, n_sites, n_junctions, time.time() - t0))
    fj.close()
    fs.close()
    fc.close()
    log("%d sites, %d calls; %d junction records -> %d junctions [%.1fs]" % (
        n_sites, n_calls, n_records, n_junctions, time.time() - t0))

    mates = find_mates(mate_index, w)
    del mate_index
    fo = Out(args.out_prefix + '.junctions.tsv.gz', args.threads)
    fo.write('#' + '\t'.join(JUNCTION_COLS) + '\n')
    with open(jtmp) as f:
        for line in f:
            jid = int(line.split('\t', 1)[0])
            fo.write("%s\t%s\n" % (line.rstrip('\n'), mates.get(jid, '.')))
    fo.close()
    log("%d junctions have a mate [%.1fs]" % (len(mates), time.time() - t0))
    shutil.rmtree(tmpdir, ignore_errors=True)


##########
# scan

def cmd_scan(args):
    t0 = time.time()
    index = SiteIndex(args.sites, args.window - 1)
    log("indexed %d sites [%.1fs]" % (index.n_sites, time.time() - t0))
    samples = read_manifest(args.manifest)
    tmpdir = tempfile.mkdtemp(prefix='bigclipper_scan_', dir=args.tmpdir)
    part_files = []
    with open(args.out_prefix + '.scan_stats.tsv', 'w') as fstats:
        fstats.write('#sample\tn_split_lines\tn_split_reads\tn_clip_lines\tphased_frac\n')
        for fields in samples:
            if len(fields) < 3:
                sys.exit("Error: scan manifest needs 3 columns (sample, intermediate file, unaligned-clips file)")
            sample, inter, clips = fields[:3]
            t1 = time.time()
            split = defaultdict(set)
            clipped = defaultdict(set)
            hp = {}
            n_lines = n_phased = n_clip_lines = 0
            reads = set()
            with open_in(inter) as f:
                for line in f:
                    x = line.rstrip('\n').split('\t')
                    n_lines += 1
                    pos0, ori, read = int(x[1]), x[5], x[8]
                    h = x[12] if len(x) >= 14 else '.'
                    if h != '.':
                        n_phased += 1
                    reads.add(hash(read))
                    for site_id in index.find(x[0], ori, pos0 + 1 if ori == '-' else pos0):
                        split[site_id].add(read)
                        hp.setdefault(read, h)
            with open_in(clips) as f:
                for line in f:
                    x = line.rstrip('\n').split('\t')
                    n_clip_lines += 1
                    pos0, ori, read = int(x[1]), x[5], x[8]
                    h = x[10] if len(x) >= 12 else '.'
                    if h != '.':
                        n_phased += 1
                    for site_id in index.find(x[0], ori, pos0 + 1 if ori == '-' else pos0):
                        clipped[site_id].add(read)
                        hp.setdefault(read, h)
            part = os.path.join(tmpdir, '%d.tsv' % len(part_files))
            with open(part, 'w') as fp:
                for site_id in sorted(set(split) | set(clipped)):
                    sr, cl = split.get(site_id, ()), clipped.get(site_id, ())
                    union = set(sr) | set(cl)
                    h1 = sum(1 for r in union if hp[r] == '1')
                    h2 = sum(1 for r in union if hp[r] == '2')
                    fp.write("%d\t%s\t%d\t%d\t%d\t%d\t%d\t%d\n" % (site_id, sample, len(sr), len(cl), len(union), h1, h2, len(union) - h1 - h2))
            part_files.append(part)
            total = n_lines + n_clip_lines
            fstats.write("%s\t%d\t%d\t%d\t%s\n" % (sample, n_lines, len(reads), n_clip_lines, fmt_float(n_phased / total) if total else '0'))
            log("%s: %d split lines, %d clip lines, reads at %d sites [%.1fs]" % (
                sample, n_lines, n_clip_lines, len(set(split) | set(clipped)), time.time() - t1))

    out = Out(args.out_prefix + '.evidence.tsv.gz')
    out.write('#' + '\t'.join(SCAN_COLS) + '\n')
    handles = [open(p) for p in part_files]
    for line in heapq.merge(*handles, key=lambda line: int(line.split('\t', 1)[0])):
        out.write(line)
    for h in handles:
        h.close()
    out.close()
    shutil.rmtree(tmpdir, ignore_errors=True)
    log("done [%.1fs]" % (time.time() - t0))


##########
# tables

def sparse_row(entries, n, dots):
    """n tab-separated columns, '.' except for entries {index: value}; dots is '.\\t' * n. O(len(entries))."""
    parts = []
    prev = -1
    for i in sorted(entries):
        if i - prev > 1:
            parts.append(dots[:2 * (i - prev - 1) - 1])
        parts.append(entries[i])
        prev = i
    if n - prev > 1:
        parts.append(dots[:2 * (n - prev - 1) - 1])
    return '\t'.join(parts)


def read_stats(paths, cols):
    stats = {}
    for path in paths:
        for r in read_rows(path):
            stats[r[0]] = dict(zip(cols, r[1:]))
    return stats


def write_vcf_header(out, samples, contigs, lengths, args):
    out.write("##fileformat=VCFv4.2\n")
    out.write("##source=BigClipper_cohort_v%s\n" % VERSION)
    out.write("##bigclipperCommand=%s\n" % ' '.join(sys.argv))
    out.write('##FILTER=<ID=PASS,Description="All filters passed">\n')
    out.write('##FILTER=<ID=MEI,Description="More than half of the samples with this junction have it in an MEI-flagged cluster (MEI_SITE or MEI_PARTNER)">\n')
    out.write('##ALT=<ID=BND,Description="Breakend">\n')
    infos = [
        ('SVTYPE', '1', 'String', 'Type of structural variant'),
        ('ORI', '1', 'String', 'BigClipper orientation of the clipped reads: + if clipped on the right (sequence left of POS kept), - if clipped on the left'),
        ('PARTNER_ORI', '1', 'String', 'BigClipper orientation of the partner: equal to ORI if the partner piece is in the same orientation, opposite if inverted'),
        ('MATEID', '1', 'String', 'ID of the junction seen from the partner breakend, if any'),
        ('SITE', '1', 'Integer', 'Breakend site (site_id in the site table) of POS'),
        ('NSAMPLES', '1', 'Integer', 'Number of samples with this junction (a partner group of at least junction_min_reads reads)'),
        ('NREADS', '1', 'Integer', 'Total supporting reads over those samples'),
        ('NMEI', '1', 'Integer', 'Number of those samples in which the cluster is MEI-flagged'),
        ('SITE_CALLED', '1', 'Integer', 'Number of samples with a cluster called at the breakend site (any partner)'),
        ('SITE_CARRIERS', '1', 'Integer', 'Number of samples called at the site or with at least af_min_reads reads there'),
        ('SITE_AF', '1', 'Float', 'SITE_CARRIERS / (2 x number of samples)'),
        ('REGIONS', '.', 'String', 'Region sets (region_beds) overlapping the breakend site'),
    ]
    for i in infos:
        out.write('##INFO=<ID=%s,Number=%s,Type=%s,Description="%s">\n' % i)
    out.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype: 0/1 if the junction was called in the sample, ./. otherwise (phase 1 has no read depth, so no reference calls)">\n')
    out.write('##FORMAT=<ID=NR,Number=1,Type=Integer,Description="Reads supporting the junction in the sample">\n')
    out.write('##FORMAT=<ID=SR,Number=1,Type=Integer,Description="Reads with a split-read breakend at the site (same orientation, any partner)">\n')
    out.write('##FORMAT=<ID=CL,Number=1,Type=Integer,Description="Reads with an unaligned clip at the site (same orientation)">\n')
    for c in contigs:
        out.write("##contig=<ID=%s%s>\n" % (c, ",length=%d" % lengths[c] if c in lengths else ''))
    out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + '\t'.join(samples) + "\n")


def cmd_tables(args):
    t0 = time.time()
    w = args.window
    samples = [s for s in read_list(args.samples) if s]
    n_samples = len(samples)
    order = {s: i for i, s in enumerate(samples)}
    ped = Pedigree(args.ped)

    # sample stats
    cstats = read_stats(args.cluster_stats, ['n_clusters', 'n_clusters_mei', 'n_junctions'])
    sstats = read_stats(args.scan_stats, ['n_split_lines', 'n_split_reads', 'n_clip_lines', 'phased_frac'])
    missing = [s for s in samples if s not in cstats or s not in sstats]
    if missing:
        sys.exit("Error: no stats for %d samples, e.g. %s" % (len(missing), ', '.join(missing[:5])))
    med_split = median_int([int(sstats[s]['n_split_lines']) for s in samples]) or 1
    med_clip = median_int([int(sstats[s]['n_clip_lines']) for s in samples]) or 1
    with open(args.out_prefix + '.sample_stats.tsv', 'w') as f:
        f.write('#' + '\t'.join(SAMPLE_STATS_COLS) + '\n')
        for s in samples:
            c, sc = cstats[s], sstats[s]
            f.write('\t'.join([s, ped.family.get(s, '.'), ped.father.get(s) or '.', ped.mother.get(s) or '.',
                               '1' if ped.affected.get(s) else '0', c['n_clusters'], c['n_clusters_mei'], c['n_junctions'],
                               sc['n_split_lines'], sc['n_split_reads'], sc['n_clip_lines'], sc['phased_frac'],
                               fmt_float(int(sc['n_split_lines']) / med_split), fmt_float(int(sc['n_clip_lines']) / med_clip)]) + '\n')

    index = SiteIndex(args.sites, w - 1)
    regions = [RegionSet(spec) for spec in args.region_bed]
    mei = None
    if args.mei_vcf:
        mei = gc.MEIIndex(args.mei_window)
        mei.load_vcf(args.mei_vcf, args.mei_info_key, args.mei_classes, args.mei_exclude_svtypes, None, True)
    lengths = {}
    if args.fai:
        for r in read_rows(args.fai):
            lengths[r[0]] = int(r[1])
    log("loaded %d sites, %d region sets [%.1fs]" % (index.n_sites, len(regions), time.time() - t0))

    col = read_header(args.sites, SITES_COLS)
    read_header(args.calls, CALLS_COLS)
    read_header(args.junctions, JUNCTION_COLS)
    calls = Grouped(read_rows(args.calls), 'calls file')
    junctions = Grouped(read_rows(args.junctions), 'junctions file', col=JUNCTION_COLS.index('site_id'))
    for path in args.evidence:
        read_header(path, SCAN_COLS)
    evidence = Grouped(heapq.merge(*[read_rows(p) for p in args.evidence], key=lambda r: int(r[0])), 'evidence files')

    fsite = Out(args.out_prefix + '.site_table.tsv.gz', args.threads)
    fsite.write('#' + '\t'.join(SITE_TABLE_COLS) + '\n')
    fev = Out(args.out_prefix + '.evidence.tsv.gz', args.threads)
    fev.write('#' + '\t'.join(EVIDENCE_COLS) + '\n')
    vcf_path = args.out_prefix + '.junctions.vcf.gz'
    fvcf = Out(vcf_path, args.threads)
    # records come in chrom_key order, so the header lists contigs in that order too
    contigs = sorted(set(lengths) | index.chroms, key=lambda c: (chrom_key(c), c))
    write_vcf_header(fvcf, samples, contigs, lengths, args)
    dots = '.\t' * n_samples
    pending = []  # heap of (pos, junction_id, line) for the current chromosome
    cur_chrom = None
    last_written = (None, 0)
    n_unsorted = 0
    n_sites = n_ev_rows = n_vcf = 0

    def flush(below=None):
        nonlocal last_written, n_unsorted, n_vcf
        while pending and (below is None or pending[0][0] < below):
            pos, _, line = heapq.heappop(pending)
            if last_written[0] == cur_chrom and pos < last_written[1]:
                n_unsorted += 1
            last_written = (cur_chrom, pos)
            fvcf.write(line)
            n_vcf += 1

    for r in read_rows(args.sites):
        site_id = int(r[col['site_id']])
        chrom, ori = r[col['chrom']], r[col['ori']]
        start, end = int(r[col['start']]), int(r[col['end']])
        if chrom != cur_chrom:
            flush()
            if cur_chrom is not None:
                log("%s: %d sites so far [%.1fs]" % (cur_chrom, n_sites, time.time() - t0))
            cur_chrom = chrom
        else:
            flush(start - 10000)

        called = {c[1]: c for c in calls.take(site_id)}
        ev = {e[1]: e for e in evidence.take(site_id)}
        jrows = junctions.take(site_id)

        carriers = set(called) | set(s for s, e in ev.items() if int(e[4]) >= args.af_min_reads)
        n_reads1 = len(ev)
        n_reads2 = sum(1 for e in ev.values() if int(e[4]) >= 2)
        n_split1 = sum(1 for e in ev.values() if int(e[2]) >= 1)
        n_clip1 = sum(1 for e in ev.values() if int(e[3]) >= 1)
        called_fams = set(ped.family.get(s, s) for s in called)
        pchrom, ppos, pori = r[col['top_partner_chrom']], int(r[col['top_partner_pos']]), r[col['top_partner_ori']]
        mate = index.nearest(pchrom, FLIP[pori], ppos)
        region_names = [rs.name for rs in regions if rs.overlaps(chrom, start, end)]
        partner_regions = [rs.name for rs in regions if rs.overlaps(pchrom, ppos, ppos)]
        mei_near = ','.join(sorted(clean_class(c) for c in mei.query(chrom, start - 1, end - 1))) or '.' if mei is not None else '.'
        af = len(carriers) / (2.0 * n_samples)
        site_af = fmt_float(af)
        fsite.write('\t'.join(map(str, (
            site_id, chrom, r[col['pos']], ori, start, end, r[col['n_clusters']], len(called), len(called_fams),
            sum(1 for s in called if ped.affected.get(s)), r[col['n_called_mei']],
            index.called_nearby(chrom, ori, start, end, site_id, w), len(carriers), site_af,
            fmt_float(len(called) / (2.0 * n_samples)), n_reads1, n_reads2, n_split1, n_clip1,
            len(set(ped.family.get(s, s) for s in ev)), len(jrows), pchrom, ppos, pori, r[col['top_partner_samples']],
            '.' if mate is None else mate, r[col['mei_classes']], mei_near, ','.join(region_names) or '.',
            ','.join(partner_regions) or '.'))) + '\n')

        for s in sorted(set(called) | set(ev), key=lambda s: order.get(s, n_samples)):
            c = called.get(s)
            e = ev.get(s)
            partner = mate = '.'
            if c:
                # the sample's own partner (of its largest cluster at the site) and the site of the other breakend
                c_chrom, c_pos, c_ori, _ = c[7].rsplit('|', 3)
                partner = "%s:%s%s" % (c_chrom, c_pos, c_ori)
                m = index.nearest(c_chrom, FLIP[c_ori], int(c_pos))
                mate = '.' if m is None else str(m)
            fev.write("%d\t%s\t%s\t%s\t%s\t%s\t%s\n" % (site_id, s, c[2] if c else '0', c[4] if c else '.',
                                                        '\t'.join(e[2:8]) if e else '0\t0\t0\t0\t0\t0', partner, mate))
            n_ev_rows += 1

        if jrows:
            ev_cols = {order[s]: "./.:.:%s:%s" % (e[2], e[3]) for s, e in ev.items() if s in order}
        for j in jrows:
            jid = int(j[0])
            cols = dict(ev_cols)
            for item in j[12].split(','):
                s, n = item.rsplit(':', 1)
                i = order.get(s)
                if i is None:
                    sys.exit("Error: junction %d has sample %s, which is not in the sample list" % (jid, s))
                e = ev.get(s)
                cols[i] = "0/1:%s:%s:%s" % (n, e[2] if e else '0', e[3] if e else '0')
            n_j, n_mei = int(j[9]), int(j[11])
            info = ["SVTYPE=BND", "ORI=%s" % j[5], "PARTNER_ORI=%s" % j[8]]
            if j[13] != '.':
                info.append("MATEID=J%s" % j[13])
            info += ["SITE=%d" % site_id, "NSAMPLES=%d" % n_j, "NREADS=%s" % j[10], "NMEI=%d" % n_mei,
                     "SITE_CALLED=%d" % len(called), "SITE_CARRIERS=%d" % len(carriers), "SITE_AF=%s" % site_af]
            if region_names:
                info.append("REGIONS=" + ','.join(region_names))
            line = "%s\t%s\tJ%d\tN\t%s\t.\t%s\t%s\tGT:NR:SR:CL\t%s\n" % (
                chrom, j[4], jid, bnd_alt(j[5], j[6], int(j[7]), j[8]), 'MEI' if 2 * n_mei > n_j else 'PASS',
                ';'.join(info), sparse_row(cols, n_samples, dots))
            heapq.heappush(pending, (int(j[4]), jid, line))
        n_sites += 1
    flush()
    calls.check_done()
    junctions.check_done()
    evidence.check_done()
    fsite.close()
    fev.close()
    fvcf.close()
    log("%d sites, %d evidence rows, %d VCF records [%.1fs]" % (n_sites, n_ev_rows, n_vcf, time.time() - t0))
    if n_unsorted:
        log("%d VCF records were out of order; sorting with bcftools" % n_unsorted)
        tmp = vcf_path + '.unsorted.vcf.gz'
        os.rename(vcf_path, tmp)
        subprocess.run(['bcftools', 'sort', '-Oz', '-o', vcf_path, tmp], check=True)
        os.remove(tmp)
    if shutil.which('tabix'):
        subprocess.run(['tabix', '-f', '-p', 'vcf', vcf_path], check=True)


##########
# filter

def cmd_filter(args):
    """Three passes over the (sorted) tables: find the rare candidates and the mate sites whose parental evidence
    they need, collect that evidence, then write the candidates."""
    t0 = time.time()
    ped = Pedigree(args.ped)
    scol = read_header(args.sample_stats, SAMPLE_STATS_COLS)
    stats = {r[0]: r for r in read_rows(args.sample_stats)}
    n_samples = len(stats)
    if n_samples == 0:
        sys.exit("Error: no samples in %s" % args.sample_stats)
    families = ped.members(stats)
    col = read_header(args.site_table, SITE_TABLE_COLS)
    ecol = read_header(args.evidence, EVIDENCE_COLS)
    E_CALLED, E_FILTER, E_SR, E_CL, E_R = (ecol[c] for c in ('called_reads', 'filter', 'SR', 'CL', 'R'))
    E_PARTNER, E_MATE = ecol['partner'], ecol['mate_site']

    def candidates():
        """(site row, AF, carriers, evidence by sample, sample) for each affected sample called at a rare site."""
        evidence = Grouped(read_rows(args.evidence), 'evidence table')
        for r in read_rows(args.site_table):
            site_id = int(r[col['site_id']])
            ev = {e[1]: e for e in evidence.take(site_id)}
            carriers = [s for s, e in ev.items() if int(e[E_CALLED]) > 0 or int(e[E_R]) >= args.af_min_reads]
            af = len(carriers) / (2.0 * n_samples)
            if af > args.max_af:
                continue
            for s in sorted(ev):
                e = ev[s]
                if int(e[E_CALLED]) == 0 or not ped.affected.get(s):
                    continue
                if args.exclude_mei and e[E_FILTER] != 'PASS':
                    continue
                yield r, af, carriers, ev, s
        evidence.check_done()

    def parents(s):
        return [p for p in (ped.father.get(s), ped.mother.get(s)) if p is not None and p in stats]

    # pass 1: which parents' evidence is needed at which mate sites
    needed = defaultdict(set)
    n_sites_rare = set()
    for r, af, carriers, ev, s in candidates():
        n_sites_rare.add(r[col['site_id']])
        mate = ev[s][E_MATE]
        if mate != '.':
            needed[int(mate)].update(parents(s))
    # pass 2: that evidence (any read, or a call)
    at_mate = {}
    for e in read_rows(args.evidence):
        want = needed.get(int(e[0]))
        if want and e[1] in want:
            at_mate[(int(e[0]), e[1])] = (int(e[E_R]), int(e[E_CALLED]) > 0)
    log("%d sites with AF <= %g; parental evidence at %d mate sites [%.1fs]" % (
        len(n_sites_rare), args.max_af, len(needed), time.time() - t0))

    # pass 3: write
    out_all = open(args.out_prefix + '.rare_candidates.tsv', 'w')
    out_dn = open(args.out_prefix + '.denovo_candidates.tsv', 'w')
    header = '\t'.join(CANDIDATE_COLS) + '\n'
    out_all.write(header)
    out_dn.write(header)
    labels = Counter()
    for r, af, carriers, ev, s in candidates():
        e = ev[s]
        mate = None if e[E_MATE] == '.' else int(e[E_MATE])
        father, mother = ped.father.get(s), ped.mother.get(s)
        has = []
        pfields = []
        for parent in (father, mother):
            if parent is None or parent not in stats:
                pfields += [parent or '.', '.', '.', '.', '.', '.']
                has.append(None)
                continue
            pe = ev.get(parent)
            here = pe is not None and (int(pe[E_R]) > 0 or int(pe[E_CALLED]) > 0)
            m_reads, m_called = at_mate.get((mate, parent), (0, False))
            pfields += [parent, pe[E_R] if pe else '0', pe[E_SR] if pe else '0', pe[E_CL] if pe else '0',
                        '.' if mate is None else str(m_reads), stats[parent][scol['rel_split_lines']]]
            has.append(here or m_reads > 0 or m_called)
        fh, mh = has
        if fh and mh:
            label = 'both'
        elif fh:
            label = 'paternal'
        elif mh:
            label = 'maternal'
        elif fh is None or mh is None:
            label = 'parents_incomplete'
        else:
            label = 'no_parental_evidence'
        labels[label] += 1
        fam = ped.family.get(s)
        relatives = []
        for rel in families.get(fam, []):
            if rel in (s, father, mother) or rel not in ev:
                continue
            re_ = ev[rel]
            if int(re_[E_R]) == 0 and int(re_[E_CALLED]) == 0:
                continue
            relatives.append("%s(%s):%s%s" % (rel, 'A' if ped.affected.get(rel) else 'U', re_[E_R],
                                              ':called' if int(re_[E_CALLED]) > 0 else ''))
        row = ([r[col[c]] for c in ('site_id', 'chrom', 'pos', 'ori', 'start', 'end')] + [s, fam or '.'] +
               [e[ecol[c]] for c in ('called_reads', 'filter', 'SR', 'CL', 'R', 'H1', 'H2', 'HU', 'partner', 'mate_site')] +
               [fmt_float(af), str(len(carriers)), r[col['AF_called']], r[col['n_called']], r[col['n_called_families']],
                r[col['n_called_nearby']], r[col['n_samples_reads1']], label] + pfields +
               [','.join(relatives) or '.'] +
               [r[col[c]] for c in ('top_partner_chrom', 'top_partner_pos', 'top_partner_ori', 'top_partner_samples',
                                    'n_junctions', 'mei_classes', 'mei_near', 'regions', 'partner_regions')])
        line = '\t'.join(row) + '\n'
        out_all.write(line)
        if label == 'no_parental_evidence':
            out_dn.write(line)
    out_all.close()
    out_dn.close()
    log("%d candidates in affected samples: %s [%.1fs]" % (
        sum(labels.values()), ', '.join("%s=%d" % kv for kv in sorted(labels.items())) or 'none', time.time() - t0))


if __name__ == "__main__":
    main()
