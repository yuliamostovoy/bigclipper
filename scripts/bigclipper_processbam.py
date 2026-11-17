import sys
import os
import re
import time
import subprocess
import argparse
import pysam

# Output files (tab-separated, sorted by chrom and start; see "Intermediate files" in the README):
#   intermediate:     chrom start start+1 . . ori mapq aligned_length read_name partner("chr:pos<ori>")
#                     partner_mapq partner_aligned_length HP PS
#   unaligned clips:  chrom start start+1 . . ori mapq aligned_length read_name clip_length HP PS
# start is the 0-based breakend: the alignment end for a right clip (ori '+'), its start for a left clip ('-').
# A partner's ori equals the index ori unless the partner is inverted relative to it (deletion ++ or --,
# inversion +- or -+).

# skip secondary and unmapped records (unmapped ones can still carry an SA tag)
READ_FILTER_FLAGS = 0x100 | 0x4

# htslib required_fields bits (only matters for CRAM: lets htslib skip decoding SEQ/QUAL)
SAM_QNAME, SAM_FLAG, SAM_RNAME, SAM_POS, SAM_MAPQ, SAM_CIGAR, SAM_AUX = 0x1, 0x2, 0x4, 0x8, 0x10, 0x20, 0x800
CRAM_REQUIRED_FIELDS = SAM_QNAME | SAM_FLAG | SAM_RNAME | SAM_POS | SAM_MAPQ | SAM_CIGAR | SAM_AUX


def main():
    parser = argparse.ArgumentParser(description="Process a BAM/CRAM file to find soft-clipped reads with supplementary alignments, and output an intermediate file for input to bigclipper_getclusters.py. Also writes clipped read ends whose clipped sequence is not aligned anywhere, to a second file.\n")
    parser.add_argument('bamfile', help="BAM/CRAM file to read: a local path, '-' for stdin (e.g. streamed with `gcloud storage cat` through bigclipper_verifystream.py), or a URL htslib can read (gs://, https://, s3://; not integrity-checked, see README)")
    parser.add_argument('-o', '--output_prefix', dest="output_prefix", required=True, help="Prefix for output file")
    parser.add_argument('-d', '--output_dir', dest='output_dir', default='.', help="Directory for output file [default=current directory]")
    parser.add_argument('-t', '--threads', dest='threads', type=int, default=4, help="Number of threads for BAM/CRAM decompression [default=4]")
    parser.add_argument('-r', '--reference', dest='reference', default=None, help="Reference FASTA (needed for CRAM input)")
    parser.add_argument('-m', '--min_clip', dest='min_clip', type=int, default=20, help="Minimum length (bp) of a clipped read end without a split-read partner to write to the unaligned-clips file [default=20]")

    args = parser.parse_args()
    if args.min_clip < 1:
        parser.error("--min_clip must be at least 1")
    if '://' in args.bamfile and not args.bamfile.startswith('file://'):
        sys.stderr.write(
            "Warning: %s will be read over the network by htslib, with no check that the whole file arrives "
            "intact. htslib has silently dropped data from GCS before (https://github.com/samtools/htslib/issues/1037). "
            "For large files, download the file first or stream it through bigclipper_verifystream.py (see README).\n"
            % args.bamfile)

    outdir = args.output_dir.rstrip('/') or '/'
    os.makedirs(outdir, exist_ok=True)
    outpath = os.path.join(outdir, "%s_intermediate.bed" % args.output_prefix)
    clips_path = os.path.join(outdir, "%s_unaligned_clips.bed" % args.output_prefix)

    bamfile = open_alignment_file(args.bamfile, args.threads, args.reference)

    start_time = time.time()
    n_records = 0
    n_sa = 0
    n_lines = 0
    n_clips = 0
    n_skipped = 0
    min_clip = args.min_clip
    with open(outpath, 'w') as outfile, open(clips_path, 'w') as clipsfile:
        for read in bamfile:
            n_records += 1
            if read.flag & READ_FILTER_FLAGS:
                continue
            if read.has_tag('SA'):
                n_sa += 1
                result = process_read(read, min_clip)
                if result is None:
                    n_skipped += 1
                    continue
                lines, clips = result
                for line in lines:
                    outfile.write(line)
                n_lines += len(lines)
            else:
                clips = unaligned_clips(read, min_clip)
            for line in clips:
                clipsfile.write(line)
            n_clips += len(clips)
    bamfile.close()

    sort_bed(outpath, outdir)
    sort_bed(clips_path, outdir)

    sys.stderr.write("Read %d records (%d with SA tags) in %.1fs; wrote %d breakends to %s and %d unaligned clips (>= %d bp) to %s%s\n" % (
        n_records, n_sa, time.time() - start_time, n_lines, outpath, n_clips, min_clip, clips_path,
        "; skipped %d records with inconsistent SA/CIGAR read lengths" % n_skipped if n_skipped else ""))


def open_alignment_file(path, threads, reference=None):
    kwargs = {'threads': max(1, threads)}
    if reference:
        kwargs['reference_filename'] = reference
    if path == '-' or path.endswith('.cram') or reference:
        # CRAM: decode only the fields used (skips SEQ/QUAL); ignored for BAM
        kwargs['format_options'] = [("required_fields=0x%x" % CRAM_REQUIRED_FIELDS).encode()]
    verbosity = pysam.set_verbosity(0)  # suppress htslib's "could not retrieve index" noise for stdin/URLs
    try:
        bamfile = pysam.AlignmentFile(path, 'r', **kwargs)  # 'r' auto-detects BAM/CRAM/SAM
    finally:
        pysam.set_verbosity(verbosity)
    return bamfile


def sort_bed(path, tmpdir):
    env = dict(os.environ, LC_ALL='C')
    subprocess.run(['sort', '-k1,1', '-k2,2n', '-T', tmpdir, '-o', path, path], check=True, env=env)


def phase_tags(read):
    """HP and PS tags of a record (e.g. from HiPhase or WhatsHap haplotag), '.' if absent."""
    try:
        hp = read.get_tag('HP')
    except KeyError:
        hp = '.'
    try:
        ps = read.get_tag('PS')
    except KeyError:
        ps = '.'
    return hp, ps


def process_read(read, min_clip):
    """Return (intermediate-file lines, unaligned-clips lines) for one primary/supplementary record with
    an SA tag, or None if its SA tag and CIGAR are inconsistent."""
    alignments = parse_bam_record(read)  # list of dicts, entries: r, rs, re, qs, qe, mq, read_length, aligned_length, strand
    if alignments is None:
        return None
    index_aln = alignments[-1]

    if read.is_reverse:
        alignments = aln_reverse(alignments)  # now the index alignment is on '+' and the others are relative to it
    alns_sorted = sorted(alignments, key=lambda k: min(k['qs'], k['qe']))
    i = next(j for j, aln in enumerate(alns_sorted) if aln is index_aln)
    aln = index_aln
    hp, ps = phase_tags(read)

    lines = []
    clips = []
    if i > 0:  # a previous alignment exists: left clip, breakend at the alignment start
        aln_SA = alns_sorted[i - 1]
        if aln_SA['strand'] == '+':
            SA = "%s:%d%s" % (aln_SA['r'], aln_SA['re'], '-')
        else:
            SA = "%s:%d%s" % (aln_SA['r'], aln_SA['rs'], '+')
        lines.append(format_line(aln, aln['rs'], '-', read.query_name, SA, aln_SA, hp, ps))
    elif aln['qs'] >= min_clip:  # the read's first alignment, with unaligned sequence before it
        clips.append(format_clip(aln['r'], aln['rs'], '-', aln['mq'], aln['aligned_length'], read.query_name, aln['qs'], hp, ps))
    if i < len(alns_sorted) - 1:  # a subsequent alignment exists: right clip, breakend at the alignment end
        aln_SA = alns_sorted[i + 1]
        if aln_SA['strand'] == '+':
            SA = "%s:%d%s" % (aln_SA['r'], aln_SA['rs'], aln_SA['strand'])
        else:
            SA = "%s:%d%s" % (aln_SA['r'], aln_SA['re'], aln_SA['strand'])
        lines.append(format_line(aln, aln['re'], '+', read.query_name, SA, aln_SA, hp, ps))
    elif aln['read_length'] - aln['qe'] >= min_clip:  # the read's last alignment, with unaligned sequence after it
        clips.append(format_clip(aln['r'], aln['re'], '+', aln['mq'], aln['aligned_length'], read.query_name,
                                 aln['read_length'] - aln['qe'], hp, ps))
    return lines, clips


def unaligned_clips(read, min_clip):
    """Unaligned-clips lines for a record without an SA tag (soft clips only)."""
    left = read.query_alignment_start
    # query_length is 0 when SEQ isn't decoded (CRAM read with required_fields); infer it from the CIGAR then
    right = (read.query_length or read.infer_query_length()) - read.query_alignment_end
    if left < min_clip and right < min_clip:
        return []
    hp, ps = phase_tags(read)
    clips = []
    if left >= min_clip:
        clips.append(format_clip(read.reference_name, read.reference_start, '-', read.mapping_quality,
                                 read.query_alignment_length, read.query_name, left, hp, ps))
    if right >= min_clip:
        clips.append(format_clip(read.reference_name, read.reference_end, '+', read.mapping_quality,
                                 read.query_alignment_length, read.query_name, right, hp, ps))
    return clips


def format_line(aln, pos, ori, read_name, SA, aln_SA, hp, ps):
    return "%s\t%d\t%d\t.\t.\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%s\n" % (
        aln['r'], pos, pos + 1, ori, aln['mq'], aln['aligned_length'], read_name, SA,
        aln_SA['mq'], aln_SA['aligned_length'], hp, ps)


def format_clip(chrom, pos, ori, mapq, aligned_length, read_name, clip_length, hp, ps):
    return "%s\t%d\t%d\t.\t.\t%s\t%d\t%d\t%s\t%d\t%s\t%s\n" % (
        chrom, pos, pos + 1, ori, mapq, aligned_length, read_name, clip_length, hp, ps)


def aln_reverse(alignments):
    # flip the orientation of all alignments for a read
    ori_flip = {'+': '-', '-': '+'}
    for aln in alignments:
        aln['qs'] = aln['read_length'] - aln['qs']
        aln['qe'] = aln['read_length'] - aln['qe']
        aln['strand'] = ori_flip[aln['strand']]
    return alignments


def parse_bam_record(record):
    chrom = record.reference_name
    rstart = record.reference_start
    mq = record.mapping_quality
    cigar = record.cigartuples

    if record.is_reverse:
        strand = "-"
    else:
        strand = "+"

    alignments = []  # list of dicts

    alignments += parse_SA_field(record.get_tag("SA"))
    alignments.append(read_cigar(cigar, chrom, rstart, strand, mq))

    read_length = alignments[-1]['read_length']

    for aln in alignments:
        if aln['read_length'] != read_length:
            sys.stderr.write("Warning: read length of primary and supplementary alignments do not match for read %s (calculated using cigar strings); skipping this record\n" % record.query_name)
            return None

    return alignments


def parse_SA_field(sa):
    alignments = []
    aligns = sa.rstrip(';').split(";")
    for align in aligns:
        fields = align.rsplit(",", 5)  # split from the right: contig names may contain commas
        if len(fields) >= 6:
            chrom = fields[0]
            rstart = int(fields[1]) - 1  # SA-tag POS is 1-based; convert to 0-based like pysam's reference_start
            raw_cigar = fields[3]
            strand = fields[2]
            mq = int(fields[4])
            alignments.append(read_cigar(parse_cigar(raw_cigar), chrom, rstart, strand, mq))
        else:
            sys.stderr.write("ignoring alternate alignment because it doesn't have all 6 columns: %s\n" % align)

    return alignments  # returns an alignment list of dicts


CIGAR_OPS = {'M': 0, 'I': 1, 'D': 2, 'N': 3, 'S': 4, 'H': 5, 'P': 6, '=': 7, 'X': 8, 'B': 9}
CIGAR_RE = re.compile(r'(\d+)(\D)')


def parse_cigar(raw_cigar):
    # returns cigar as list of tuples just like in pysam
    cigar = []
    for num, element in CIGAR_RE.findall(raw_cigar):
        if element not in CIGAR_OPS:
            sys.stderr.write("CIGAR element not found: %s\n" % element)
            sys.exit(1)
        cigar.append((CIGAR_OPS[element], int(num)))
    return cigar


def read_cigar(cigar, chrom, rstart, strand, mq):
    # returns an alignment dict
    coordinates = cigar_coords(cigar)

    alignment = {'r': chrom, 'rs': rstart, 're': rstart + coordinates['ref_alignment_length']}

    if strand == "+":
        alignment['qs'] = coordinates['front_padding_length']
        alignment['qe'] = coordinates['front_padding_length'] + coordinates['read_alignment_length']
    else:
        alignment['qe'] = coordinates['end_padding_length']
        alignment['qs'] = coordinates['end_padding_length'] + coordinates['read_alignment_length']

    alignment['read_length'] = coordinates['front_padding_length'] + coordinates['read_alignment_length'] + coordinates['end_padding_length']
    alignment['mq'] = mq
    alignment['aligned_length'] = coordinates['read_alignment_length']
    alignment['strand'] = strand

    # skipping the part where we add a 'path' field to alignment because we don't need it here
    return alignment


def cigar_coords(cigar):
    coords = {'read_alignment_length': 0, 'ref_alignment_length': 0, 'front_padding_length': 0, 'end_padding_length': 0}

    no_matches_yet = True
    for element in cigar:
        num = element[1]
        if element[0] == 5:  # H
            if no_matches_yet:
                coords['front_padding_length'] += num
            else:
                coords['end_padding_length'] += num
        elif element[0] == 4:  # S
            if no_matches_yet:
                coords['front_padding_length'] += num
            else:
                coords['end_padding_length'] += num
        elif element[0] in [0, 7, 8]:  # M, =, X
            no_matches_yet = False
            coords['read_alignment_length'] += num
            coords['ref_alignment_length'] += num
        elif element[0] == 1:  # I
            no_matches_yet = False
            coords['read_alignment_length'] += num
        elif element[0] in [2, 3, 6]:  # D, N, P
            no_matches_yet = False
            coords['ref_alignment_length'] += num
        else:
            sys.stderr.write("Don't recognize cigar character: %s, assuming it advances both query and reference, like a match or mismatch\n" % element[0])
            coords['read_alignment_length'] += num
            coords['ref_alignment_length'] += num
    return coords


if __name__ == "__main__": main()
