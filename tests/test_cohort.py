"""Tests for bigclipper_cohort.py on a small synthetic cohort (all sample and family IDs are made up).

Run from the repository root:  python -m unittest discover -s tests -v
Requires pysam, bedtools, bgzip/tabix and GNU sort (gsort on macOS) on PATH.
"""
import csv
import gzip
import os
import random
import shutil
import subprocess
import sys
import tempfile
import unittest

import pysam

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
COHORT = os.path.join(REPO, 'scripts', 'bigclipper_cohort.py')
sys.path.insert(0, os.path.join(REPO, 'scripts'))
import bigclipper_cohort as bc  # noqa: E402

FLIP = {'+': '-', '-': '+'}
TRIOS = [('FAM1', 'F1_father', 'F1_mother', 'F1_child'), ('FAM2', 'F2_father', 'F2_mother', 'F2_child')]
SINGLETONS = ['S%02d' % i for i in range(1, 14)]
SAMPLES = [s for t in TRIOS for s in t[1:]] + ['F3_child'] + SINGLETONS  # 20 samples, so AF = carriers / 40
MAX_AF = 0.075  # 3 carriers kept, 4 excluded
OLD_FORMAT = 'S02'  # written as a 2.0 (12-column) intermediate file


def run(args, expect_fail=False):
    proc = subprocess.run([sys.executable, COHORT] + args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    if expect_fail:
        assert proc.returncode != 0, "expected failure: %s" % ' '.join(args)
    elif proc.returncode != 0:
        raise AssertionError("command failed: %s\n%s" % (' '.join(args), proc.stderr))
    return proc


class Cohort:
    """Per-sample intermediate and unaligned-clips lines in the bigclipper_processbam.py format."""

    def __init__(self):
        self.inter = {s: [] for s in SAMPLES}
        self.clips = {s: [] for s in SAMPLES}
        self.n = 0

    def junction(self, sample, chrom, pos0, ori, pchrom, ppos0, pori, n, hp=None, jitter=0):
        """n split reads joining chrom:pos0 (ori) to pchrom:ppos0 (pori), both sides written as processbam would."""
        for i in range(n):
            self.n += 1
            read = 'r%d' % self.n
            d = (i % 3 - 1) * jitter
            h = hp[i] if hp else '.'
            self.inter[sample].append((chrom, pos0 + d, ori, read, pchrom, ppos0 + d, pori, h))
            self.inter[sample].append((pchrom, ppos0 + d, FLIP[pori], read, chrom, pos0 + d, FLIP[ori], h))

    def clip(self, sample, chrom, pos0, ori, n):
        for i in range(n):
            self.n += 1
            self.clips[sample].append((chrom, pos0, ori, 'c%d' % self.n))

    def write(self, d):
        paths = {}
        for s in SAMPLES:
            inter = os.path.join(d, '%s_intermediate.bed' % s)
            with open(inter, 'w') as f:
                for chrom, pos, ori, read, pchrom, ppos, pori, h in sorted(self.inter[s]):
                    fields = [chrom, pos, pos + 1, '.', '.', ori, 60, 5000, read, '%s:%d%s' % (pchrom, ppos, pori), 60, 4000]
                    if s != OLD_FORMAT:
                        fields += [h, '.' if h == '.' else 1000]
                    f.write('\t'.join(map(str, fields)) + '\n')
            clips = os.path.join(d, '%s_unaligned_clips.bed' % s)
            with open(clips, 'w') as f:
                for chrom, pos, ori, read in sorted(self.clips[s]):
                    f.write('\t'.join(map(str, [chrom, pos, pos + 1, '.', '.', ori, 60, 5000, read, 300, '.', '.'])) + '\n')
            paths[s] = (inter, clips)
        return paths


def build_cohort():
    c = Cohort()
    # de novo deletion in F1_child: chr1:1,000,000 (last base kept) to chr1:1,020,001; 4 reads on haplotype 1
    c.junction('F1_child', 'chr1', 1000000, '+', 'chr1', 1020000, '+', 6, hp='1111..')
    # translocation in F2_child; its father has 3 split reads at the same breakend with a short (< min_dist) partner
    c.junction('F2_child', 'chr1', 2000000, '+', 'chr5', 500000, '+', 6)
    c.junction('F2_father', 'chr1', 2000010, '+', 'chr1', 2001000, '+', 3)
    # F2_child breakend whose only parental evidence is 2 unaligned clips in the mother
    c.junction('F2_child', 'chr1', 3000000, '-', 'chr1', 3050000, '-', 6)
    c.clip('F2_mother', 'chr1', 2999995, '-', 2)
    # 3 carriers (AF 0.075): kept
    for s in ('F1_child', 'S01', OLD_FORMAT):
        c.junction(s, 'chr1', 4000000, '+', 'chr9', 100000, '-', 6)
    # 4 carriers (AF 0.1): excluded
    for s in ('F1_child', 'S03', 'S04', 'S05'):
        c.junction(s, 'chr1', 5000000, '+', 'chr9', 200000, '-', 6)
    # called in F1_child only, but 3 samples have 2 reads each (not called): 4 carriers, excluded
    c.junction('F1_child', 'chr1', 6000000, '+', 'chr9', 300000, '-', 6)
    for s in ('S06', 'S07', 'S08'):
        c.junction(s, 'chr1', 6000005, '+', 'chr9', 300000, '-', 2)
    # called in F1_child; 5 samples with a single read each (noise, not carriers): kept
    c.junction('F1_child', 'chr1', 7000000, '+', 'chr9', 400000, '-', 6)
    for s in ('S06', 'S07', 'S08', 'S09', 'S10'):
        c.clip(s, 'chr1', 7000003, '+', 1)
    # common breakend (10 samples) and a rare one 60 bp away in F2_child: separate sites
    for s in SINGLETONS[:10]:
        c.junction(s, 'chr1', 8000000, '+', 'chr9', 500000, '-', 6, jitter=3)
    c.junction('F2_child', 'chr1', 8000060, '+', 'chr9', 600000, '-', 6)
    # rare '+' breakend in F1_child; 5 samples have a '-' breakend at the same position
    c.junction('F1_child', 'chr1', 9000000, '+', 'chr9', 700000, '-', 6)
    for s in SINGLETONS[:5]:
        c.junction(s, 'chr1', 9000000, '-', 'chr9', 800000, '+', 6)
    # unique to an unaffected sample: never a candidate
    c.junction('S11', 'chr1', 10000000, '+', 'chr9', 900000, '-', 6)
    # affected sample whose parents are not in the cohort
    c.junction('F3_child', 'chr1', 11000000, '+', 'chr9', 950000, '-', 6)
    # F2_child breakend at an MEI insertion called in F2_child
    c.junction('F2_child', 'chr1', 12000000, '+', 'chr7', 300000, '+', 6)
    # a '-' cluster with two partner groups 40 bp apart and a '+' site between them: junction IDs (by position)
    # interleave with site IDs
    c.junction('F1_child', 'chr1', 13000000, '-', 'chr9', 1000, '+', 3)
    c.junction('F1_child', 'chr1', 13000040, '-', 'chr8', 1000, '+', 3)
    c.junction('S12', 'chr1', 13000020, '+', 'chr9', 2000, '-', 6)
    # inherited deletion whose only parental evidence is 2 clips in the mother at the other breakend
    c.junction('F1_child', 'chr1', 14000000, '+', 'chr1', 14100000, '+', 6)
    c.clip('F1_mother', 'chr1', 14100000, '-', 2)
    # a cluster whose reads spread 45 bp either side of its median; the father's clips are 85 bp from the median
    c.junction('F2_child', 'chr1', 15000000, '+', 'chr9', 1500000, '-', 6, jitter=45)
    c.clip('F2_father', 'chr1', 15000085, '+', 2)
    return c


def write_ped(path, samples=SAMPLES):
    with open(path, 'w') as f:
        f.write('#family\tsample\tfather\tmother\tsex\tphenotype\n')
        for fam, father, mother, child in TRIOS:
            f.write('%s\t%s\t0\t0\t1\t1\n' % (fam, father))
            f.write('%s\t%s\t0\t0\t2\t1\n' % (fam, mother))
            f.write('%s\t%s\t%s\t%s\t1\t2\n' % (fam, child, father, mother))
        f.write('FAM3\tF3_child\tF3_father\tF3_mother\t2\t2\n')
        for s in SINGLETONS:
            f.write('%s_fam\t%s\t0\t0\t1\t1\n' % (s, s))


def write_mei_vcf(path):
    lines = ['##fileformat=VCFv4.2',
             '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="type">',
             '##INFO=<ID=MEI_type,Number=1,Type=String,Description="MEI class">',
             '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
             '##contig=<ID=chr1,length=20000000>',
             '\t'.join(['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'] + SAMPLES)]
    lines.append('\t'.join(['chr1', '12000000', 'mei1', 'N', '<INS>', '.', 'PASS', 'SVTYPE=INS;MEI_type=ALU', 'GT'] +
                           ['0/1' if s == 'F2_child' else '0/0' for s in SAMPLES]))
    with open(path, 'w') as f:
        f.write('\n'.join(lines) + '\n')


def read_tsv(path):
    opener = gzip.open if path.endswith('.gz') else open
    with opener(path, 'rt') as f:
        header = f.readline().rstrip('\n').lstrip('#').split('\t')
        return [dict(zip(header, line.rstrip('\n').split('\t'))) for line in f]


class TestHelpers(unittest.TestCase):
    def test_bnd_alt_round_trip(self):
        # deletion-type junctions are ++ and -- (a tandem duplication's too, with the partner on the other side),
        # inversion-type +- and -+
        cases = {('+', '+'): 'N[chr2:500[', ('+', '-'): 'N]chr2:500]', ('-', '-'): ']chr2:500]N', ('-', '+'): '[chr2:500[N'}
        for (ori, pori), alt in cases.items():
            self.assertEqual(bc.bnd_alt(ori, 'chr2', 500, pori), alt)
            self.assertEqual(bc.parse_bnd_alt(alt), (ori, 'chr2', 500, pori))
        self.assertEqual(bc.parse_bnd_alt('N[HLA-A*01:01:01:01:7['), ('+', 'HLA-A*01:01:01:01', 7, '+'))

    def test_breakend_bases(self):
        # deletion: read aligned to [999000, 1000000) then [1020000, 1021000) (0-based half-open)
        # processbam: '+' breakend at 1000000 with partner rs=1020000 '+'; '-' breakend at 1020000 with partner re=1000000 '-'
        self.assertEqual(bc.local_bp(1000000, '+'), 1000000)   # last kept base, 1-based
        self.assertEqual(bc.partner_bp(1020000, '+'), 1020001)  # first base after the deletion, 1-based
        self.assertEqual(bc.local_bp(1020000, '-'), 1020001)
        self.assertEqual(bc.partner_bp(1000000, '-'), 1000000)
        # so the two sides of the deletion are N[chr1:1020001[ at 1000000 and ]chr1:1000000]N at 1020001
        self.assertEqual(bc.bnd_alt('+', 'chr1', bc.partner_bp(1020000, '+'), '+'), 'N[chr1:1020001[')
        self.assertEqual(bc.bnd_alt('-', 'chr1', bc.partner_bp(1000000, '-'), '-'), ']chr1:1000000]N')

    def test_sparse_row(self):
        n = 7
        for entries in [{}, {0: 'a'}, {6: 'z'}, {0: 'a', 6: 'z'}, {3: 'c', 2: 'b'}, {i: str(i) for i in range(n)}]:
            self.assertEqual(bc.sparse_row(entries, n, '.\t' * n), '\t'.join(entries.get(i, '.') for i in range(n)))

    def test_chrom_key(self):
        chroms = ['chr2', '3', 'chr10', 'chrX', 'chrY', 'chrM', 'chr1_KI270706v1_random']
        keys = [bc.chrom_key(c) for c in chroms]
        self.assertEqual(keys, ['0002', '0003', '0010', '1001', '1002', '1003', '9999:chr1_KI270706v1_random'])
        self.assertEqual(sorted(keys), keys)  # GNU sort (LC_ALL=C) and Python order these keys the same way


def J(pos, ppos, sample, reads, ori='+', pchrom='chr2', pori='+', cluster='1'):
    return bc.Junction('chr1', pos, ori, pchrom, ppos, pori, sample, cluster, reads, 'PASS')


class TestMatchJunctions(unittest.TestCase):
    def summary(self, recs, window=50):
        return sorted((rep.pos, rep.ppos, rep.sample, sorted((m.sample, m.pos) for m in members))
                      for rep, members in bc.match_junctions(recs, window))

    def test_same_junction_merges(self):
        recs = [J(1000, 5000, 'A', 5), J(1020, 4990, 'B', 3), J(990, 5030, 'C', 8)]
        self.assertEqual(self.summary(recs), [(990, 5030, 'C', [('A', 1000), ('B', 1020), ('C', 990)])])
        # 50 bp from the representative at either end is too far
        recs = [J(1000, 5000, 'A', 5), J(1020, 4980, 'B', 3), J(990, 5030, 'C', 8)]
        self.assertEqual(len(self.summary(recs)), 2)

    def test_window_is_strict_on_both_ends(self):
        recs = [J(1000, 5000, 'A', 5), J(1050, 5000, 'B', 3), J(1000, 5050, 'C', 3), J(1049, 5049, 'D', 3)]
        out = self.summary(recs)
        self.assertEqual(len(out), 3)
        self.assertIn((1000, 5000, 'A', [('A', 1000), ('D', 1049)]), out)

    def test_distinct_partners_and_orientations_stay_apart(self):
        recs = [J(1000, 5000, 'A', 5), J(1000, 6000, 'B', 5), J(1000, 5000, 'C', 5, pori='-'),
                J(1000, 5000, 'D', 5, ori='-'), J(1000, 5000, 'E', 5, pchrom='chr3')]
        self.assertEqual(len(self.summary(recs)), 5)

    def test_no_chaining(self):
        recs = [J(1000, 5000, 'A', 9), J(1040, 5000, 'B', 2), J(1080, 5000, 'C', 2)]
        self.assertEqual(self.summary(recs), [(1000, 5000, 'A', [('A', 1000), ('B', 1040)]),
                                              (1080, 5000, 'C', [('C', 1080)])])

    def test_deterministic(self):
        rng = random.Random(1)
        recs = [J(1000 + rng.randrange(200), 5000 + rng.randrange(200), 'S%d' % i, rng.randrange(2, 9)) for i in range(300)]
        expected = self.summary(recs)
        for _ in range(3):
            rng.shuffle(recs)
            self.assertEqual(self.summary(recs), expected)
        # each record is in exactly one junction, within the window of its representative
        out = bc.match_junctions(recs, 50)
        self.assertEqual(sum(len(m) for _, m in out), len(recs))
        for rep, members in out:
            for m in members:
                self.assertTrue(abs(m.pos - rep.pos) < 50 and abs(m.ppos - rep.ppos) < 50)
                self.assertLessEqual(m.reads, rep.reads)

    def test_dense_pileup_is_fast(self):
        recs = [J(1000 + (i % 7), 5000 + (i % 11), 'S%d' % i, 2) for i in range(20000)]
        out = bc.match_junctions(recs, 50)
        self.assertEqual(len(out), 1)


class TestCohortPipeline(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp(prefix='bigclipper_cohort_test_')
        d = cls.tmp
        cls.paths = build_cohort().write(d)
        cls.ped = os.path.join(d, 'cohort.ped')
        write_ped(cls.ped)
        cls.mei = os.path.join(d, 'mei.vcf')
        write_mei_vcf(cls.mei)
        cls.samples = os.path.join(d, 'samples.txt')
        with open(cls.samples, 'w') as f:
            f.write('\n'.join(SAMPLES) + '\n')
        run(['check', '--samples', cls.samples, '--ped', cls.ped, '--n_intermediate', str(len(SAMPLES)), '--n_clips', str(len(SAMPLES))])

        batches = [SAMPLES[:11], SAMPLES[11:]]
        for b, batch in enumerate(batches):
            with open(os.path.join(d, 'batch%d.tsv' % b), 'w') as f:
                for s in batch:
                    f.write('%s\t%s\t%s\n' % (s, cls.paths[s][0], cls.paths[s][1]))
            run(['clusters', '--manifest', os.path.join(d, 'batch%d.tsv' % b), '--out_prefix', os.path.join(d, 'b%d' % b),
                 '--mei_vcf', cls.mei])
        run(['sites', '--clusters'] + [os.path.join(d, 'b%d.clusters.tsv.gz' % b) for b in range(2)] +
            ['--junctions'] + [os.path.join(d, 'b%d.junctions.tsv.gz' % b) for b in range(2)] +
            ['--out_prefix', os.path.join(d, 'cohort'), '--tmpdir', d])
        for b in range(2):
            run(['scan', '--sites', os.path.join(d, 'cohort.sites.tsv.gz'), '--manifest', os.path.join(d, 'batch%d.tsv' % b),
                 '--out_prefix', os.path.join(d, 'b%d' % b), '--tmpdir', d])
        with open(os.path.join(d, 'regions.bed'), 'w') as f:
            f.write('chr1\t10999990\t11000010\tsat\nchr9\t949000\t951000\tsat\n')
        run(['tables', '--samples', cls.samples, '--ped', cls.ped, '--sites', os.path.join(d, 'cohort.sites.tsv.gz'),
             '--calls', os.path.join(d, 'cohort.calls.tsv.gz'), '--junctions', os.path.join(d, 'cohort.junctions.tsv.gz'),
             '--evidence'] + [os.path.join(d, 'b%d.evidence.tsv.gz' % b) for b in range(2)] +
            ['--cluster_stats'] + [os.path.join(d, 'b%d.cluster_stats.tsv' % b) for b in range(2)] +
            ['--scan_stats'] + [os.path.join(d, 'b%d.scan_stats.tsv' % b) for b in range(2)] +
            ['--region_bed', 'satellites=' + os.path.join(d, 'regions.bed'), '--mei_vcf', cls.mei,
             '--out_prefix', os.path.join(d, 'cohort')])
        cls.run_filter('cohort', ['--max_af', str(MAX_AF)])
        cls.sites = read_tsv(os.path.join(d, 'cohort.site_table.tsv.gz'))
        cls.evidence = read_tsv(os.path.join(d, 'cohort.evidence.tsv.gz'))
        cls.cands = read_tsv(os.path.join(d, 'cohort.rare_candidates.tsv'))
        cls.denovo = read_tsv(os.path.join(d, 'cohort.denovo_candidates.tsv'))

    @classmethod
    def run_filter(cls, prefix, extra):
        d = cls.tmp
        return run(['filter', '--site_table', os.path.join(d, 'cohort.site_table.tsv.gz'),
                    '--evidence', os.path.join(d, 'cohort.evidence.tsv.gz'), '--sample_stats', os.path.join(d, 'cohort.sample_stats.tsv'),
                    '--ped', cls.ped, '--out_prefix', os.path.join(d, prefix)] + extra)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp)

    def site(self, pos, ori='+', chrom='chr1'):
        hits = [s for s in self.sites if s['chrom'] == chrom and s['ori'] == ori and int(s['start']) - 50 < pos < int(s['end']) + 50]
        self.assertEqual(len(hits), 1, "expected one %s site near %s:%d, found %d" % (ori, chrom, pos, len(hits)))
        return hits[0]

    def cand(self, pos, sample, ori='+', rows=None):
        hits = [c for c in (self.cands if rows is None else rows)
                if c['chrom'] == 'chr1' and c['ori'] == ori and abs(int(c['pos']) - pos) < 50 and c['sample'] == sample]
        return hits[0] if hits else None

    def test_de_novo(self):
        c = self.cand(1000000, 'F1_child')
        self.assertIsNotNone(c)
        self.assertEqual(c['parental_evidence'], 'no_parental_evidence')
        self.assertEqual(int(c['pos']), 1000000)
        self.assertEqual((c['SR'], c['CL'], c['R'], c['H1'], c['H2'], c['HU']), ('6', '0', '6', '4', '0', '2'))
        self.assertEqual((c['father_R'], c['mother_R'], c['father_mate_R'], c['mother_mate_R']), ('0', '0', '0', '0'))
        self.assertEqual((c['partner'], c['mate_site']), ('chr1:1020001+', self.site(1020001, ori='-')['site_id']))
        self.assertEqual((c['top_partner_chrom'], c['top_partner_pos'], c['top_partner_ori']), ('chr1', '1020001', '+'))
        self.assertIsNotNone(self.cand(1000000, 'F1_child', rows=self.denovo))
        # the other side of the deletion is a site too, and each is the other's mate
        other = self.site(1020001, ori='-')
        self.assertEqual(int(other['pos']), 1020001)
        self.assertEqual(self.site(1000000)['mate_site'], other['site_id'])
        self.assertEqual(other['mate_site'], self.site(1000000)['site_id'])

    def test_paternal_split_reads_with_short_partner(self):
        c = self.cand(2000000, 'F2_child')
        self.assertEqual(c['parental_evidence'], 'paternal')
        self.assertEqual((c['father_SR'], c['father_CL'], c['father_R']), ('3', '0', '3'))
        self.assertIsNone(self.cand(2000000, 'F2_child', rows=self.denovo))
        self.assertEqual(self.site(2000000)['n_called'], '1')

    def test_maternal_clip_only(self):
        c = self.cand(3000001, 'F2_child', ori='-')
        self.assertEqual(c['parental_evidence'], 'maternal')
        self.assertEqual((c['mother_SR'], c['mother_CL'], c['mother_R']), ('0', '2', '2'))
        self.assertEqual(int(c['pos']), 3000001)

    def test_parental_evidence_at_mate_breakend(self):
        here, there = self.cand(14000000, 'F1_child'), self.cand(14100001, 'F1_child', ori='-')
        self.assertEqual((here['parental_evidence'], here['mother_R'], here['mother_mate_R']), ('maternal', '0', '2'))
        self.assertEqual((there['parental_evidence'], there['mother_R'], there['mother_mate_R']), ('maternal', '2', '0'))
        self.assertFalse([c for c in self.denovo if c['site_id'] in (here['site_id'], there['site_id'])])

    def test_site_range_covers_read_spread(self):
        s = self.site(15000000)
        self.assertEqual((s['pos'], s['start'], s['end']), ('15000000', '14999955', '15000045'))
        c = self.cand(15000000, 'F2_child')
        self.assertEqual((c['parental_evidence'], c['father_CL']), ('paternal', '2'))
        # every called sample has at least its called reads at the site
        short = [e for e in self.evidence if int(e['R']) < int(e['called_reads'])]
        self.assertEqual(short, [])

    def test_af_threshold(self):
        kept = self.cand(4000000, 'F1_child')
        self.assertIsNotNone(kept)
        self.assertEqual((kept['n_carriers'], kept['AF'], kept['n_called']), ('3', '0.075', '3'))
        self.assertIsNone(self.cand(5000000, 'F1_child'))
        self.assertEqual(self.site(5000000)['n_carriers'], '4')

    def test_sub_threshold_carriers_count(self):
        s = self.site(6000000)
        self.assertEqual((s['n_called'], s['n_carriers'], s['AF_called']), ('1', '4', '0.025'))
        self.assertIsNone(self.cand(6000000, 'F1_child'))

    def test_single_read_noise_not_counted(self):
        s = self.site(7000000)
        self.assertEqual((s['n_carriers'], s['n_samples_reads1'], s['n_samples_clip1']), ('1', '6', '5'))
        self.assertIsNotNone(self.cand(7000000, 'F1_child'))

    def test_bounded_sites(self):
        common, rare = self.site(8000000), self.site(8000060)
        self.assertNotEqual(common['site_id'], rare['site_id'])
        self.assertEqual(common['n_called'], '10')
        self.assertEqual(rare['n_carriers'], '1')
        self.assertEqual(rare['n_called_nearby'], '10')
        self.assertIsNotNone(self.cand(8000060, 'F2_child'))

    def test_orientation_respected(self):
        plus, minus = self.site(9000000, '+'), self.site(9000001, '-')
        self.assertEqual((plus['n_called'], plus['n_carriers']), ('1', '1'))
        self.assertEqual(minus['n_called'], '5')
        self.assertIsNotNone(self.cand(9000000, 'F1_child'))

    def test_unaffected_carrier_excluded(self):
        self.assertEqual(self.site(10000000)['n_called'], '1')
        self.assertFalse([c for c in self.cands if c['sample'] == 'S11'])
        self.assertFalse([c for c in self.cands if c['sample'] not in ('F1_child', 'F2_child', 'F3_child')])

    def test_parents_incomplete_and_regions(self):
        c = self.cand(11000000, 'F3_child')
        self.assertEqual(c['parental_evidence'], 'parents_incomplete')
        self.assertEqual((c['father'], c['father_R']), ('F3_father', '.'))
        self.assertEqual((c['regions'], c['partner_regions']), ('satellites', 'satellites'))
        self.assertIsNone(self.cand(11000000, 'F3_child', rows=self.denovo))

    def test_old_format_intermediate(self):
        ev = [e for e in self.evidence if e['sample'] == OLD_FORMAT and e['site_id'] == self.site(4000000)['site_id']]
        self.assertEqual((ev[0]['called_reads'], ev[0]['SR'], ev[0]['HU']), ('6', '6', '6'))

    def test_mei(self):
        c = self.cand(12000000, 'F2_child')
        self.assertEqual(c['filter'], 'MEI_SITE')
        self.assertEqual((c['mei_classes'], c['mei_near']), ('ALU', 'ALU'))
        self.run_filter('nomei', ['--max_af', str(MAX_AF), '--exclude_mei'])
        rows = read_tsv(os.path.join(self.tmp, 'nomei.rare_candidates.tsv'))
        self.assertIsNone(self.cand(12000000, 'F2_child', rows=rows))
        # the other side of the junction (chr7) is dropped too: its partner is at the MEI (MEI_PARTNER)
        self.assertEqual(len(rows), len(self.cands) - 2)
        other = [r for r in self.cands if r['chrom'] == 'chr7' and r['sample'] == 'F2_child']
        self.assertEqual((other[0]['pos'], other[0]['ori'], other[0]['filter']), ('300001', '-', 'MEI_PARTNER'))

    def test_refilter(self):
        self.run_filter('strict', ['--max_af', '0.025'])
        rows = read_tsv(os.path.join(self.tmp, 'strict.rare_candidates.tsv'))
        self.assertIsNone(self.cand(4000000, 'F1_child', rows=rows))
        self.assertTrue(rows)
        self.assertTrue(all(float(r['AF']) <= 0.025 for r in rows))
        self.assertEqual({(r['site_id'], r['sample']) for r in rows},
                         {(r['site_id'], r['sample']) for r in self.cands if int(r['n_carriers']) <= 1})

    def test_sample_stats(self):
        stats = read_tsv(os.path.join(self.tmp, 'cohort.sample_stats.tsv'))
        self.assertEqual([s['sample'] for s in stats], SAMPLES)
        f1 = [s for s in stats if s['sample'] == 'F1_child'][0]
        self.assertEqual((f1['family'], f1['father'], f1['affected']), ('FAM1', 'F1_father', '1'))
        self.assertEqual(f1['n_split_reads'], str(len({r[3] for r in build_cohort().inter['F1_child']})))

    def test_vcf(self):
        path = os.path.join(self.tmp, 'cohort.junctions.vcf.gz')
        self.assertTrue(os.path.exists(path + '.tbi'))
        with pysam.VariantFile(path) as vcf:
            self.assertEqual(list(vcf.header.samples), SAMPLES)
            recs = list(vcf)
        by_id = {r.id: r for r in recs}
        last = {}
        for r in recs:  # sorted
            self.assertGreaterEqual(r.pos, last.get(r.chrom, 0))
            last[r.chrom] = r.pos
        dn = [r for r in recs if r.chrom == 'chr1' and r.pos == 1000000][0]
        self.assertEqual(dn.alts, ('N[chr1:1020001[',))
        self.assertEqual(dn.samples['F1_child']['GT'], (0, 1))
        self.assertEqual(dn.samples['F1_child']['NR'], 6)
        self.assertEqual(dn.samples['F1_father']['GT'], (None,))
        mate = by_id[dn.info['MATEID']]
        self.assertEqual((mate.chrom, mate.pos, mate.alts), ('chr1', 1020001, (']chr1:1000000]N',)))
        self.assertEqual(mate.info['MATEID'], dn.id)
        common = [r for r in recs if r.chrom == 'chr1' and abs(r.pos - 8000000) < 10]
        self.assertEqual(len(common), 1)
        self.assertEqual(common[0].info['NSAMPLES'], 10)
        self.assertEqual(common[0].info['NREADS'], 60)
        # the father's split reads at the translocation breakend show up as SR without a call
        tra = [r for r in recs if r.chrom == 'chr1' and r.pos == 2000000][0]
        self.assertEqual(tra.alts, ('N[chr5:500001[',))
        self.assertEqual(tra.samples['F2_father']['GT'], (None, None))
        self.assertEqual(tra.samples['F2_father']['SR'], 3)
        mei = [r for r in recs if r.chrom == 'chr1' and r.pos == 12000000][0]
        self.assertEqual(list(mei.filter.keys()), ['MEI'])

    def test_junctions_interleaved_with_sites(self):
        minus, plus = self.site(13000041, '-'), self.site(13000020, '+')
        self.assertEqual((minus['n_junctions'], plus['n_junctions']), ('2', '1'))
        # in the junctions file (sorted by site), junction IDs (by position) are out of order here
        rows = [r for r in read_tsv(os.path.join(self.tmp, 'cohort.junctions.tsv.gz')) if r['site_id'] in (minus['site_id'], plus['site_id'])]
        ids = [int(r['junction_id']) for r in rows]
        self.assertNotEqual(ids, sorted(ids))
        with pysam.VariantFile(os.path.join(self.tmp, 'cohort.junctions.vcf.gz')) as vcf:
            recs = [(r.pos, r.alts[0], r.info['SITE']) for r in vcf.fetch('chr1', 12999990, 13000060)]
        self.assertEqual(recs, [(13000001, '[chr9:1001[N', int(minus['site_id'])),
                                (13000020, 'N]chr9:2000]', int(plus['site_id'])),
                                (13000041, '[chr8:1001[N', int(minus['site_id']))])

    def test_check_fails(self):
        bad = os.path.join(self.tmp, 'bad_samples.txt')
        with open(bad, 'w') as f:
            f.write('\n'.join(SAMPLES + ['NOT_IN_PED']) + '\n')
        proc = run(['check', '--samples', bad, '--ped', self.ped], expect_fail=True)
        self.assertIn('NOT_IN_PED', proc.stderr)
        proc = run(['check', '--samples', self.samples, '--ped', self.ped, '--n_intermediate', '19'], expect_fail=True)
        self.assertIn('19 intermediate files', proc.stderr)
        # intermediate and clips arrays must pair up by file prefix
        inter, clips = os.path.join(self.tmp, 'inter.txt'), os.path.join(self.tmp, 'clips.txt')
        with open(inter, 'w') as f:
            f.write('\n'.join('gs://bucket/%s/%s' % (s, os.path.basename(self.paths[s][0])) for s in SAMPLES) + '\n')
        with open(clips, 'w') as f:
            f.write('\n'.join('gs://bucket/%s/%s' % (s, os.path.basename(self.paths[s][1])) for s in SAMPLES) + '\n')
        run(['check', '--samples', self.samples, '--ped', self.ped, '--intermediate_list', inter, '--clips_list', clips])
        with open(clips, 'w') as f:
            f.write('\n'.join(os.path.basename(self.paths[s][1]) for s in SAMPLES[1:] + SAMPLES[:1]) + '\n')
        proc = run(['check', '--samples', self.samples, '--ped', self.ped, '--intermediate_list', inter, '--clips_list', clips], expect_fail=True)
        self.assertIn('do not share', proc.stderr)
        with open(bad, 'w') as f:
            f.write('\n'.join(SAMPLES + ['S01']) + '\n')
        proc = run(['check', '--samples', bad, '--ped', self.ped], expect_fail=True)
        self.assertIn('duplicate', proc.stderr)

    def test_missing_mei_sample_fails(self):
        manifest = os.path.join(self.tmp, 'extra.tsv')
        with open(manifest, 'w') as f:
            f.write('NOT_IN_VCF\t%s\n' % self.paths['S01'][0])
        proc = run(['clusters', '--manifest', manifest, '--out_prefix', os.path.join(self.tmp, 'extra'), '--mei_vcf', self.mei], expect_fail=True)
        self.assertIn('NOT_IN_VCF', proc.stderr)


if __name__ == '__main__':
    unittest.main()
