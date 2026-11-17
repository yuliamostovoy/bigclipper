"""Tests for BigClipper on synthetic split reads with known breakpoints.

Run from the repository root:  python -m unittest discover -s tests -v
Requires pysam and bedtools on PATH.
"""
import base64
import os
import shlex
import shutil
import struct
import subprocess
import sys
import tempfile
import unittest

import pysam

try:
    import google_crc32c
    HAVE_CRC32C = google_crc32c.implementation == 'c'
except ImportError:
    HAVE_CRC32C = False

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PROCESSBAM = os.path.join(REPO, 'scripts', 'bigclipper_processbam.py')
GETCLUSTERS = os.path.join(REPO, 'scripts', 'bigclipper_getclusters.py')
VERIFYSTREAM = os.path.join(REPO, 'scripts', 'bigclipper_verifystream.py')

HLA = 'HLA-A*01:01:01:01'  # contig name containing ':'
CONTIGS = [('chr1', 10000000), ('chr2', 10000000), (HLA, 100000)]

# Each scenario: read segments in read order as (contig, 0-based ref start, aligned length, strand)
SCENARIOS = {
    # 4 kb deletion chr1:102000-106000, read on + strand
    'del_f': [('chr1', 100000, 2000, '+'), ('chr1', 106000, 2000, '+')],
    # same deletion seen by a read on the - strand
    'del_r': [('chr1', 106000, 2000, '-'), ('chr1', 100000, 2000, '-')],
    # inversion junction: chr1:202000 (+) joined to chr1:207000 read backwards
    'inv': [('chr1', 200000, 2000, '+'), ('chr1', 205000, 2000, '-')],
    # translocation to a contig whose name contains ':'
    'tra': [('chr1', 300000, 2000, '+'), (HLA, 1000, 2000, '+')],
    # MEI-like insertion at chr1:403000 (15 bp TSD) whose inserted sequence aligns to chr2:500000-506000
    'mei': [('chr1', 400000, 3000, '+'), ('chr2', 500000, 6000, '+'), ('chr1', 402985, 3000, '+')],
}

# Expected intermediate lines per scenario: (contig, pos, ori, partner) with 0-based breakend positions.
# Partner ori follows BigClipper's convention: same as the index ori unless the partner is inverted relative to it.
EXPECTED = {
    'del_f': {('chr1', 102000, '+', 'chr1:106000+'), ('chr1', 106000, '-', 'chr1:102000-')},
    'del_r': {('chr1', 102000, '+', 'chr1:106000+'), ('chr1', 106000, '-', 'chr1:102000-')},
    'inv': {('chr1', 202000, '+', 'chr1:207000-'), ('chr1', 207000, '+', 'chr1:202000-')},
    'tra': {('chr1', 302000, '+', HLA + ':1000+'), (HLA, 1000, '-', 'chr1:302000-')},
    'mei': {('chr1', 403000, '+', 'chr2:500000+'), ('chr2', 500000, '-', 'chr1:403000-'),
            ('chr2', 506000, '+', 'chr1:402985+'), ('chr1', 402985, '-', 'chr2:506000-')},
}


def cigar_string(cigar):
    return ''.join('%d%s' % (n, 'MIDNSHP=XB'[op]) for op, n in cigar)


def make_records(header, name, segments, mapq=60, hard_clip_supp=False, shift=0, lead=0, trail=0, tags=None):
    """Primary (segment 0) + supplementary records for one read, with SA tags (none for a single segment).
    lead/trail: unaligned bases at the start/end of the read (in read order); tags: extra tags for every record."""
    read_len = lead + sum(seg[2] for seg in segments) + trail
    qstarts = [lead + sum(seg[2] for seg in segments[:k]) for k in range(len(segments))]
    cigars = []
    for (contig, start, length, strand), q in zip(segments, qstarts):
        front, end = (q, read_len - q - length) if strand == '+' else (read_len - q - length, q)
        cigars.append([(4, front)] * bool(front) + [(0, length)] + [(4, end)] * bool(end))
    records = []
    for k, (contig, start, length, strand) in enumerate(segments):
        rec = pysam.AlignedSegment(header)
        rec.query_name = name
        rec.flag = (0x10 if strand == '-' else 0) | (0x800 if k > 0 else 0)
        rec.reference_name = contig
        rec.reference_start = start + shift
        rec.mapping_quality = mapq
        cigar = cigars[k]
        if k > 0 and hard_clip_supp:
            cigar = [(5 if op == 4 else op, n) for op, n in cigar]
        rec.cigartuples = cigar
        rec.query_sequence = 'A' * sum(n for op, n in cigar if op in (0, 1, 4))
        sa = ''.join('%s,%d,%s,%s,%d,0;' % (s[0], s[1] + shift + 1, s[3], cigar_string(cigars[j]), mapq)
                     for j, s in enumerate(segments) if j != k)
        if sa:
            rec.set_tag('SA', sa)
        for tag, value in (tags or {}).items():
            rec.set_tag(tag, value)
        records.append(rec)
    return records


def write_bam(path, reads):
    """reads: list of (name, segments, kwargs)"""
    header = pysam.AlignmentHeader.from_dict({'HD': {'VN': '1.6', 'SO': 'coordinate'},
                                              'SQ': [{'SN': c, 'LN': l} for c, l in CONTIGS]})
    records = []
    for name, segments, kwargs in reads:
        records += make_records(header, name, segments, **kwargs)
    records.sort(key=lambda r: (r.reference_id, r.reference_start))
    with pysam.AlignmentFile(path, 'wb', header=header) as out:
        for rec in records:
            out.write(rec)
    pysam.index(path)


def run(cmd, expect_fail=False):
    proc = subprocess.run([sys.executable] + cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    if expect_fail:
        assert proc.returncode != 0, "expected failure: %s" % ' '.join(cmd)
    elif proc.returncode != 0:
        raise AssertionError("command failed: %s\n%s" % (' '.join(cmd), proc.stderr))
    return proc


def read_intermediate(path):
    with open(path) as f:
        return [line.rstrip('\n').split('\t') for line in f]


# Reads for the unaligned-clips file: (name, segments, kwargs), and the clipped ends expected with the
# default --min_clip 20, as (contig, 0-based breakend, ori, clip_length, read name)
CLIP_READS = [
    ('solo_f', [('chr1', 150000, 2000, '+')], {'lead': 50, 'trail': 15}),   # 15 bp trailing clip is below 20
    ('solo_r', [('chr2', 150000, 2000, '-')], {'lead': 30, 'trail': 100}),  # - strand: read end is on the left
    ('noclip', [('chr2', 160000, 2000, '+')], {}),
    ('split_lead', SCENARIOS['del_f'], {'lead': 40}),     # split read whose start isn't aligned anywhere
    ('split_trail', SCENARIOS['del_r'], {'trail': 60, 'hard_clip_supp': True}),
]
EXPECTED_CLIPS = {
    ('chr1', 150000, '-', 50, 'solo_f'),
    ('chr2', 150000, '-', 100, 'solo_r'), ('chr2', 152000, '+', 30, 'solo_r'),
    ('chr1', 100000, '-', 40, 'split_lead'),
    ('chr1', 100000, '-', 60, 'split_trail'),
}


def read_clips(path):
    lines = read_intermediate(path)
    return {(l[0], int(l[1]), l[5], int(l[9]), l[8]) for l in lines}, lines


def read_vcf(path):
    with pysam.VariantFile(path) as vcf:  # also checks that the VCF parses
        return [{'chrom': r.chrom, 'pos': r.pos, 'ori': r.info['ORI'], 'alns': list(r.info['ALNS']),
                 'nreads': r.info['NREADS'], 'filter': list(r.filter.keys()), 'info': dict(r.info)} for r in vcf]


def write_mei_vcf(path, records, samples):
    lines = ['##fileformat=VCFv4.2',
             '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="type">',
             '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="len">',
             '##INFO=<ID=MEI_type,Number=1,Type=String,Description="MEI class">',
             '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">']
    lines += ['##contig=<ID=%s,length=%d>' % c for c in CONTIGS]
    lines.append('\t'.join(['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'] + samples))
    for chrom, pos, info, gts in records:
        lines.append('\t'.join([chrom, str(pos), '.', 'N', '<INS>', '.', 'PASS', info, 'GT'] + gts))
    with open(path, 'w') as f:
        f.write('\n'.join(lines) + '\n')


class TestProcessBam(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp()

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp)

    def check_scenario(self, name, **kwargs):
        bam = os.path.join(self.tmp, name + '.bam')
        write_bam(bam, [(name, SCENARIOS[name], kwargs)])
        run([PROCESSBAM, bam, '-o', name, '-d', self.tmp, '-t', '2'])
        lines = read_intermediate(os.path.join(self.tmp, name + '_intermediate.bed'))
        got = {(l[0], int(l[1]), l[5], l[9]) for l in lines}
        self.assertEqual(got, EXPECTED[name])
        self.assertEqual(len(lines), len(EXPECTED[name]))
        for l in lines:
            self.assertEqual(len(l), 14)
            self.assertEqual(int(l[2]), int(l[1]) + 1)
            self.assertEqual(l[12:], ['.', '.'])  # no HP/PS tags
        return lines

    def test_deletion_forward(self):
        self.check_scenario('del_f')

    def test_deletion_reverse_read(self):
        self.check_scenario('del_r')

    def test_inversion(self):
        self.check_scenario('inv')

    def test_translocation_contig_with_colon(self):
        self.check_scenario('tra')

    def test_mei_like_three_segments(self):
        lines = self.check_scenario('mei')
        # partner aligned length/MAPQ columns: the chr1 flanks' partner is the 6 kb chr2 segment
        for l in lines:
            if l[0] == 'chr1':
                self.assertEqual((l[10], l[11]), ('60', '6000'))

    def test_hard_clipped_supplementary(self):
        self.check_scenario('del_f', hard_clip_supp=True)

    def test_phase_tags(self):
        bam = os.path.join(self.tmp, 'phased.bam')
        write_bam(bam, [('del_f', SCENARIOS['del_f'], {'tags': {'HP': 2, 'PS': 98765}}),
                        ('solo_f', CLIP_READS[0][1], dict(CLIP_READS[0][2], tags={'HP': 1, 'PS': 140000}))])
        run([PROCESSBAM, bam, '-o', 'phased', '-d', self.tmp])
        lines = read_intermediate(os.path.join(self.tmp, 'phased_intermediate.bed'))
        self.assertEqual(len(lines), 2)
        self.assertTrue(all(l[12:] == ['2', '98765'] for l in lines))
        _, clips = read_clips(os.path.join(self.tmp, 'phased_unaligned_clips.bed'))
        self.assertEqual([l[10:] for l in clips], [['1', '140000']])

    def test_unaligned_clips(self):
        bam = os.path.join(self.tmp, 'clips.bam')
        write_bam(bam, CLIP_READS)
        proc = run([PROCESSBAM, bam, '-o', 'clips', '-d', self.tmp])
        got, lines = read_clips(os.path.join(self.tmp, 'clips_unaligned_clips.bed'))
        self.assertEqual(got, EXPECTED_CLIPS)
        self.assertEqual(len(lines), len(EXPECTED_CLIPS))
        for l in lines:
            self.assertEqual(len(l), 12)
            self.assertEqual((int(l[2]), l[3], l[4], l[6], l[7], l[10], l[11]), (int(l[1]) + 1, '.', '.', '60', '2000', '.', '.'))
        keys = [(l[0], int(l[1])) for l in lines]
        self.assertEqual(keys, sorted(keys))
        self.assertIn('5 unaligned clips (>= 20 bp)', proc.stderr)
        # the split reads' partnered ends are unchanged in the intermediate file
        inter = read_intermediate(os.path.join(self.tmp, 'clips_intermediate.bed'))
        self.assertEqual({(l[0], int(l[1]), l[5], l[9]) for l in inter}, EXPECTED['del_f'])
        self.assertEqual(len(inter), 4)
        # a lower threshold adds the 15 bp clip
        run([PROCESSBAM, bam, '-o', 'clips10', '-d', self.tmp, '--min_clip', '10'])
        got10, _ = read_clips(os.path.join(self.tmp, 'clips10_unaligned_clips.bed'))
        self.assertEqual(got10, EXPECTED_CLIPS | {('chr1', 152000, '+', 15, 'solo_f')})
        run([PROCESSBAM, bam, '-o', 'clips0', '-d', self.tmp, '--min_clip', '0'], expect_fail=True)

    def test_unmapped_read_with_sa_tag(self):
        """Unmapped records can carry an SA tag; they have no CIGAR and are skipped."""
        bam = os.path.join(self.tmp, 'unmapped_sa.bam')
        write_bam(bam, [('del_f', SCENARIOS['del_f'], {})])
        with pysam.AlignmentFile(bam) as src:
            header, records = src.header, list(src)
        rec = pysam.AlignedSegment(header)
        rec.query_name, rec.flag = 'unmapped', 0x1 | 0x4 | 0x8 | 0x40
        rec.reference_id, rec.reference_start, rec.mapping_quality = 0, 150000, 0
        rec.query_sequence = 'A' * 250
        rec.set_tag('SA', 'chr2,1000,+,250M,13,17;')
        with pysam.AlignmentFile(bam, 'wb', header=header) as out:
            for r in sorted(records + [rec], key=lambda r: (r.reference_id, r.reference_start)):
                out.write(r)
        run([PROCESSBAM, bam, '-o', 'unmapped_sa', '-d', self.tmp])
        lines = read_intermediate(os.path.join(self.tmp, 'unmapped_sa_intermediate.bed'))
        self.assertEqual({(l[0], int(l[1]), l[5], l[9]) for l in lines}, EXPECTED['del_f'])

    def test_url_input_warns(self):
        # nothing listens on port 9, so the open fails fast; the warning comes before it
        proc = run([PROCESSBAM, 'https://127.0.0.1:9/x.bam', '-o', 'url', '-d', self.tmp], expect_fail=True)
        self.assertIn('read over the network by htslib', proc.stderr)
        bam = os.path.join(self.tmp, 'local_nowarn.bam')
        write_bam(bam, [('del_f', SCENARIOS['del_f'], {})])
        self.assertNotIn('Warning', run([PROCESSBAM, bam, '-o', 'nowarn', '-d', self.tmp]).stderr)

    def test_stdin(self):
        bam = os.path.join(self.tmp, 'stdin.bam')
        write_bam(bam, [('del_r', SCENARIOS['del_r'], {})])
        with open(bam, 'rb') as f:
            subprocess.run([sys.executable, PROCESSBAM, '-', '-o', 'stdin', '-d', self.tmp], stdin=f, check=True)
        lines = read_intermediate(os.path.join(self.tmp, 'stdin_intermediate.bed'))
        self.assertEqual({(l[0], int(l[1]), l[5], l[9]) for l in lines}, EXPECTED['del_r'])

    def test_cram(self):
        """CRAM input (decoded without SEQ/QUAL), from a file and from stdin, matches BAM input."""
        ref = os.path.join(self.tmp, 'ref.fa')
        with open(ref, 'w') as f:
            for contig, length in CONTIGS:
                f.write('>%s\n' % contig)
                for i in range(0, length, 1000000):
                    f.write('A' * min(1000000, length - i) + '\n')
        pysam.faidx(ref)
        bam = os.path.join(self.tmp, 'cram_src.bam')
        write_bam(bam, [(n, s, {}) for n, s in SCENARIOS.items()] + CLIP_READS)
        cram = os.path.join(self.tmp, 'test.cram')
        with pysam.AlignmentFile(bam) as src, pysam.AlignmentFile(cram, 'wc', template=src, reference_filename=ref) as out:
            for rec in src:
                out.write(rec)
        run([PROCESSBAM, bam, '-o', 'cram_bam', '-d', self.tmp])
        run([PROCESSBAM, cram, '-o', 'cram_file', '-d', self.tmp, '-r', ref])
        with open(cram, 'rb') as f:
            subprocess.run([sys.executable, PROCESSBAM, '-', '-o', 'cram_stdin', '-d', self.tmp, '-r', ref],
                           stdin=f, check=True, stderr=subprocess.DEVNULL)
        expected = read_intermediate(os.path.join(self.tmp, 'cram_bam_intermediate.bed'))
        self.assertEqual(len(expected), sum(len(v) for v in EXPECTED.values()) + 4)  # + the two split CLIP_READS
        expected_clips = read_intermediate(os.path.join(self.tmp, 'cram_bam_unaligned_clips.bed'))
        self.assertEqual(read_clips(os.path.join(self.tmp, 'cram_bam_unaligned_clips.bed'))[0], EXPECTED_CLIPS)
        for name in ('cram_file', 'cram_stdin'):
            self.assertEqual(read_intermediate(os.path.join(self.tmp, name + '_intermediate.bed')), expected)
            # CRAM is decoded without SEQ, so clip lengths come from the CIGAR alone
            self.assertEqual(read_intermediate(os.path.join(self.tmp, name + '_unaligned_clips.bed')), expected_clips)

    def test_output_sorted(self):
        bam = os.path.join(self.tmp, 'all.bam')
        write_bam(bam, [(n, s, {}) for n, s in SCENARIOS.items()])
        run([PROCESSBAM, bam, '-o', 'all', '-d', self.tmp])
        lines = read_intermediate(os.path.join(self.tmp, 'all_intermediate.bed'))
        keys = [(l[0], int(l[1])) for l in lines]
        self.assertEqual(keys, sorted(keys))


def bgzf_block_offsets(data):
    offsets, p = [], 0
    while p < len(data):
        offsets.append(p)
        p += struct.unpack('<H', data[p + 16:p + 18])[0] + 1  # BSIZE = block size - 1
    return offsets


def crc32c_b64(data):
    return base64.b64encode(struct.pack('>I', google_crc32c.value(data))).decode()


@unittest.skipUnless(HAVE_CRC32C, 'google-crc32c (with its C extension) is not installed')
class TestVerifyStream(unittest.TestCase):
    """A stream that is cut short must fail the pipeline instead of producing a partial intermediate file."""

    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp()
        bam = os.path.join(cls.tmp, 'many.bam')
        write_bam(bam, [('r%d' % i, SCENARIOS['del_f'], {'shift': i * 10}) for i in range(1000)])
        with open(bam, 'rb') as f:
            cls.data = f.read()
        cls.size, cls.crc = len(cls.data), crc32c_b64(cls.data)
        run([PROCESSBAM, bam, '-o', 'direct', '-d', cls.tmp])
        cls.expected = read_intermediate(os.path.join(cls.tmp, 'direct_intermediate.bed'))

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp)

    def verify(self, data, size, crc):
        return subprocess.run([sys.executable, VERIFYSTREAM, '--size', str(size), '--crc32c', crc],
                              input=data, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    def stream(self, name, data, verify=True):
        """`cat data | [verify |] processbam -` under pipefail: (exit code, intermediate lines, stderr)"""
        src = os.path.join(self.tmp, name + '.input')
        with open(src, 'wb') as f:
            f.write(data)
        py = shlex.quote(sys.executable)
        cmd = 'cat %s' % shlex.quote(src)
        if verify:
            cmd += ' | %s %s --size %d --crc32c %s' % (py, shlex.quote(VERIFYSTREAM), self.size, shlex.quote(self.crc))
        cmd += ' | %s %s - -o %s -d %s' % (py, shlex.quote(PROCESSBAM), name, shlex.quote(self.tmp))
        proc = subprocess.run(['bash', '-o', 'pipefail', '-c', cmd],
                              stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        out = os.path.join(self.tmp, name + '_intermediate.bed')
        return proc.returncode, read_intermediate(out) if os.path.exists(out) else None, proc.stderr

    def test_crc32c_matches_gcs_encoding(self):
        proc = self.verify(b'hello', 5, 'mnG7TA==')  # value GCS reports for an object containing "hello"
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertEqual(proc.stdout, b'hello')

    def test_intact_stream(self):
        rc, lines, err = self.stream('intact', self.data)
        self.assertEqual(rc, 0, err)
        self.assertEqual(lines, self.expected)
        self.assertIn('verified %d bytes' % self.size, err)

    def test_cut_between_bgzf_blocks(self):
        offsets = bgzf_block_offsets(self.data)
        self.assertGreater(len(offsets), 20)
        cut = self.data[:offsets[len(offsets) // 2]]
        # the hazard: processbam alone exits 0 with about half the reads
        rc, lines, err = self.stream('cut_unverified', cut, verify=False)
        self.assertEqual(rc, 0, err)
        self.assertLess(len(lines), len(self.expected))
        rc, lines, err = self.stream('cut_verified', cut)
        self.assertNotEqual(rc, 0)
        self.assertIn('received %d bytes, expected %d' % (len(cut), self.size), err)

    def test_repeated_bytes(self):
        """e.g. a retry that resends data it already sent"""
        half = bgzf_block_offsets(self.data)[10]
        rc, lines, err = self.stream('repeated', self.data[:half] + self.data)
        self.assertNotEqual(rc, 0)
        self.assertIn('expected %d' % self.size, err)

    def test_corrupt_same_size(self):
        data = bytearray(self.data)
        data[len(data) // 2] ^= 1
        proc = self.verify(bytes(data), self.size, self.crc)
        self.assertNotEqual(proc.returncode, 0)
        self.assertIn(b'CRC32C is', proc.stderr)
        self.assertEqual(proc.stdout, bytes(data))  # passed through unchanged; only the exit code flags it

    def test_rejects_missing_or_malformed_expectations(self):
        for size, crc in ((self.size, ''), (self.size, 'abc'), (-1, self.crc)):
            self.assertNotEqual(self.verify(b'x', size, crc).returncode, 0, (size, crc))
        proc = subprocess.run([sys.executable, VERIFYSTREAM, '--size', '', '--crc32c', self.crc], input=b'x',
                              stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        self.assertNotEqual(proc.returncode, 0)

    def test_downstream_stops_reading(self):
        data = os.urandom(8 << 20)  # much more than a pipe buffer holds
        src = os.path.join(self.tmp, 'random.bin')
        with open(src, 'wb') as f:
            f.write(data)
        cmd = '%s %s --size %d --crc32c %s < %s | head -c 10 > /dev/null' % (
            shlex.quote(sys.executable), shlex.quote(VERIFYSTREAM), len(data), shlex.quote(crc32c_b64(data)),
            shlex.quote(src))
        proc = subprocess.run(['bash', '-o', 'pipefail', '-c', cmd], stderr=subprocess.PIPE, universal_newlines=True)
        self.assertNotEqual(proc.returncode, 0)
        self.assertIn('stopped reading', proc.stderr)
        self.assertNotIn('Traceback', proc.stderr)


class TestGetClusters(unittest.TestCase):
    N = 6  # reads per scenario

    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp()
        reads = []
        for name, segments in SCENARIOS.items():
            for i in range(cls.N):
                # jitter alignment positions by a few bp, alternate hard/soft clipping
                reads.append(('%s_%d' % (name, i), segments, {'shift': i % 3, 'hard_clip_supp': bool(i % 2)}))
        cls.bam = os.path.join(cls.tmp, 'sample.bam')
        write_bam(cls.bam, reads)
        run([PROCESSBAM, cls.bam, '-o', 'sample', '-d', cls.tmp])
        cls.intermediate = os.path.join(cls.tmp, 'sample_intermediate.bed')

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp)

    def getclusters(self, *extra, **kw):
        out = os.path.join(self.tmp, kw.get('out', 'out.vcf'))
        proc = run([GETCLUSTERS, kw.get('intermediate', self.intermediate), '-c', '5', '-o', out] + list(extra),
                   expect_fail=kw.get('expect_fail', False))
        return proc if kw.get('expect_fail') else (read_vcf(out), proc)

    def by_site(self, records):
        return {(r['chrom'], r['pos'], r['ori']): r for r in records}

    def test_clusters_and_partners(self):
        records, _ = self.getclusters()
        sites = self.by_site(records)
        # one cluster per breakend side (del_f and del_r reads pile up together: 12 reads)
        self.assertEqual(set(sites), {
            ('chr1', 102001, '+'), ('chr1', 106001, '-'), ('chr1', 202001, '+'), ('chr1', 207001, '+'),
            ('chr1', 302001, '+'), (HLA, 1001, '-'), ('chr1', 403001, '+'), ('chr1', 402986, '-'),
            ('chr2', 500001, '-'), ('chr2', 506001, '+')})
        self.assertEqual(sites[('chr1', 102001, '+')]['alns'], ['chr1:106001-106003|+|12'])
        self.assertEqual(sites[('chr1', 102001, '+')]['nreads'], 12)
        self.assertEqual(sites[('chr1', 106001, '-')]['alns'], ['chr1:102001-102003|-|12'])
        self.assertEqual(sites[('chr1', 202001, '+')]['alns'], ['chr1:207001-207003|-|6'])
        self.assertEqual(sites[(HLA, 1001, '-')]['alns'], ['chr1:302001-302003|-|6'])
        self.assertEqual(sites[('chr1', 302001, '+')]['alns'], [HLA + ':1001-1003|+|6'])
        self.assertTrue(all(r['filter'] == ['PASS'] for r in records))
        self.assertEqual(sites[('chr1', 403001, '+')]['info']['PARTNER_ALNLEN'], (6000,))
        # records are position-sorted (so the VCF can be bgzipped and indexed)
        order = [(r['chrom'], r['pos']) for r in records]
        self.assertEqual([o for o in order if o[0] == 'chr1'], sorted(o for o in order if o[0] == 'chr1'))

    def test_v2_0_intermediate(self):
        """A 2.0 intermediate (12 columns, no HP/PS) gives the same clusters as the current 14-column format."""
        v20 = os.path.join(self.tmp, 'v20_intermediate.bed')
        with open(self.intermediate) as f, open(v20, 'w') as out:
            for line in f:
                fields = line.rstrip('\n').split('\t')
                self.assertEqual(len(fields), 14)
                out.write('\t'.join(fields[:12]) + '\n')
        new, _ = self.getclusters()
        old, proc = self.getclusters(intermediate=v20, out='v20.vcf')
        self.assertNotIn('pre-2.0', proc.stderr)
        self.assertEqual(new, old)

    def test_legacy_intermediate(self):
        """A pre-2.0 intermediate (10 columns, partner position +1) gives the same clusters."""
        legacy = os.path.join(self.tmp, 'legacy_intermediate.bed')
        with open(self.intermediate) as f, open(legacy, 'w') as out:
            for line in f:
                fields = line.rstrip('\n').split('\t')[:10]
                ref, pos = fields[9][:-1].rsplit(':', 1)
                fields[9] = '%s:%d%s' % (ref, int(pos) + 1, fields[9][-1])
                out.write('\t'.join(fields) + '\n')
        new, _ = self.getclusters()
        old, proc = self.getclusters(intermediate=legacy, out='legacy.vcf')
        self.assertIn('pre-2.0 format', proc.stderr)
        strip = lambda rs: [(r['chrom'], r['pos'], r['ori'], r['alns'], r['nreads']) for r in rs]
        self.assertEqual(strip(new), strip(old))
        self.assertNotIn('PARTNER_MAPQ', old[0]['info'])

    def test_mei_filter_sample_genotype(self):
        mei_vcf = os.path.join(self.tmp, 'mei.vcf')
        write_mei_vcf(mei_vcf, [
            ('chr1', 403000, 'SVTYPE=INS;SVLEN=6000;MEI_type=L1', ['0/1', '0/0']),   # the MEI, carried by S1
            ('chr1', 102000, 'SVTYPE=INS;SVLEN=300;MEI_type=Alu', ['0/0', '0/1']),   # MEI near the deletion, not in S1
            ('chr1', 302000, 'SVTYPE=INS;SVLEN=500', ['0/1', '0/1']),                 # non-MEI insertion
        ], ['S1', 'S2'])
        records, proc = self.getclusters('--mei_vcf', mei_vcf, '--mei_sample', 'S1')
        flagged = {k: r['filter'] for k, r in self.by_site(records).items() if r['filter'] != ['PASS']}
        self.assertEqual(flagged, {('chr1', 403001, '+'): ['MEI_SITE'], ('chr1', 402986, '-'): ['MEI_SITE'],
                                   ('chr2', 500001, '-'): ['MEI_PARTNER'], ('chr2', 506001, '+'): ['MEI_PARTNER']})
        self.assertEqual(self.by_site(records)[('chr2', 500001, '-')]['info']['MEI_CLASS'], ('L1',))
        self.assertEqual(self.by_site(records)[('chr2', 500001, '-')]['info']['MEI_PARTNER_FRAC'], 1.0)

        # --pass_only drops them
        records, _ = self.getclusters('--mei_vcf', mei_vcf, '--mei_sample', 'S1', '--pass_only', out='pass.vcf')
        self.assertEqual(len(records), 6)

        # a sample that doesn't carry the L1 only gets the Alu near the deletion
        records, _ = self.getclusters('--mei_vcf', mei_vcf, '--mei_sample', 'S2', out='s2.vcf')
        flagged = {k for k, r in self.by_site(records).items() if r['filter'] != ['PASS']}
        self.assertEqual(flagged, {('chr1', 102001, '+'), ('chr1', 106001, '-')})  # site + partner

        # multi-sample VCF needs --mei_sample or --mei_all_sites
        proc = self.getclusters('--mei_vcf', mei_vcf, expect_fail=True)
        self.assertIn('--mei_sample', proc.stderr)
        records, _ = self.getclusters('--mei_vcf', mei_vcf, '--mei_all_sites', out='all.vcf')
        self.assertEqual(sum(r['filter'] != ['PASS'] for r in records), 6)

    def test_mei_filter_window_and_errors(self):
        mei_vcf = os.path.join(self.tmp, 'mei1.vcf')
        write_mei_vcf(mei_vcf, [('chr1', 403250, 'SVTYPE=INS;MEI_type=SVA', ['0/1'])], ['S1'])
        records, _ = self.getclusters('--mei_vcf', mei_vcf, out='w100.vcf')  # single sample: used automatically
        self.assertTrue(all(r['filter'] == ['PASS'] for r in records))  # 250 bp away > default window
        records, _ = self.getclusters('--mei_vcf', mei_vcf, '--mei_window', '300', out='w300.vcf')
        self.assertEqual(sum(r['filter'] != ['PASS'] for r in records), 4)

        proc = self.getclusters('--mei_vcf', mei_vcf, '--mei_info_key', 'NOPE', expect_fail=True)
        self.assertIn('INFO/NOPE', proc.stderr)

        other_ref = os.path.join(self.tmp, 'mei_chm13.bed')
        with open(other_ref, 'w') as f:
            f.write('NC_060925.1\t403000\t403001\tL1\n')
        _, proc = self.getclusters('--mei_bed', other_ref, out='bed.vcf')
        self.assertIn('none of the MEI call contigs', proc.stderr)

    def test_mei_svtype_and_class_filters(self):
        mei_vcf = os.path.join(self.tmp, 'mei_types.vcf')
        write_mei_vcf(mei_vcf, [
            ('chr1', 102000, 'SVTYPE=DEL;SVLEN=-300;MEI_type=Alu', ['0/1']),    # deleted reference Alu at the deletion
            ('chr1', 202000, 'SVTYPE=INS;SVLEN=80;MEI_type=TANDEM', ['0/1']),   # non-MEI value at the inversion
        ], ['S1'])
        flagged = lambda records: {(r['chrom'], r['pos']) for r in records if r['filter'] != ['PASS']}
        records, proc = self.getclusters('--mei_vcf', mei_vcf, out='types.vcf')
        self.assertEqual(flagged(records), {('chr1', 202001), ('chr1', 207001)})  # DEL ignored by default
        self.assertIn('skipped, by SVTYPE: DEL=1', proc.stderr)
        self.assertIn('used, by class: TANDEM=1', proc.stderr)
        records, _ = self.getclusters('--mei_vcf', mei_vcf, '--mei_classes', 'alu,L1,SVA', out='types2.vcf')
        self.assertEqual(flagged(records), set())
        records, _ = self.getclusters('--mei_vcf', mei_vcf, '--mei_classes', 'alu,L1,SVA', '--mei_exclude_svtypes', 'NONE', out='types3.vcf')
        self.assertEqual(flagged(records), {('chr1', 102001), ('chr1', 106001)})

    def test_mei_bed(self):
        bed = os.path.join(self.tmp, 'mei.bed')
        with open(bed, 'w') as f:
            f.write('chr1\t403000\t403001\tL1\n')
        records, _ = self.getclusters('--mei_bed', bed, out='bed2.vcf')
        self.assertEqual(sum(r['filter'] != ['PASS'] for r in records), 4)

    def write_intermediate(self, name, lines):
        path = os.path.join(self.tmp, name)
        with open(path, 'w') as f:
            for i, (pos, ori, read, partner) in enumerate(lines):
                f.write('chr5\t%d\t%d\t.\t.\t%s\t60\t5000\t%s\t%s\t60\t5000\n' % (pos, pos + 1, ori, read, partner))
        return path

    def test_min_dist_rechecks_cluster_size(self):
        """Reads whose partner fails --min_dist don't count toward --min_cluster_count."""
        lines = [(1000 + i, '+', 'local%d' % i, 'chr5:%d+' % (1100 + i)) for i in range(6)]
        lines.append((1003, '+', 'distant', 'chr9:5000+'))
        path = self.write_intermediate('mindist_intermediate.bed', lines)
        records, _ = self.getclusters('-d', '1', intermediate=path, out='md1.vcf')
        self.assertEqual(records[0]['nreads'], 7)
        records, _ = self.getclusters('-d', '1000', intermediate=path, out='md1000.vcf')
        self.assertEqual(records, [])

    def test_partner_cluster_distance(self):
        """Partner breakends up to -s bp apart are merged."""
        lines = [(1000, '+', 'r%d' % i, 'chr9:%d+' % (5000 + 45 * i)) for i in range(5)]
        path = self.write_intermediate('pdist_intermediate.bed', lines)
        records, _ = self.getclusters(intermediate=path, out='pd.vcf')
        self.assertEqual(records[0]['alns'], ['chr9:5001-5181|+|5'])

    def test_counts_unique_reads(self):
        """A read contributing two records to one cluster counts once."""
        lines = [(1000 + i, '+', 'r%d' % (i // 2), 'chr9:%d+' % (5000 + i)) for i in range(8)]
        path = self.write_intermediate('dupread_intermediate.bed', lines)
        records, _ = self.getclusters(intermediate=path, out='dup.vcf')
        self.assertEqual(records, [])  # only 4 distinct reads < -c 5


if __name__ == '__main__':
    unittest.main()
