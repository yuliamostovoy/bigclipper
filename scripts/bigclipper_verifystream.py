#!/usr/bin/env python3
"""Pass a stream from stdin to stdout unchanged, then check it against an expected size and CRC32C.

Use it between `gcloud storage cat` and bigclipper_processbam.py when streaming a BAM/CRAM
from Google Cloud Storage:

    set -o pipefail
    read -r SIZE CRC < <(gcloud storage objects describe gs://bucket/sample.bam --format='value(size,crc32c_hash)')
    gcloud storage cat gs://bucket/sample.bam \\
        | python bigclipper_verifystream.py --size "$SIZE" --crc32c "$CRC" \\
        | python bigclipper_processbam.py - -o sample

Exits non-zero if the stream differs from the stored object, so that with pipefail a truncated
stream (which processbam would read as complete) fails the pipeline.
"""
import argparse
import base64
import binascii
import os
import struct
import sys

CHUNK_SIZE = 1 << 20
PROG = 'bigclipper_verifystream'


def encode_crc32c(crc):
    """CRC32C as GCS reports it: big-endian, base64-encoded"""
    return base64.b64encode(struct.pack('>I', crc)).decode()


def decode_crc32c(value):
    try:
        raw = base64.b64decode(value, validate=True)
    except (binascii.Error, ValueError):
        raw = b''
    if len(raw) != 4:
        sys.exit("%s: invalid --crc32c %r; expected a base64-encoded 4-byte CRC32C, "
                 "as reported by `gcloud storage objects describe --format='value(crc32c_hash)'`" % (PROG, value))
    return struct.unpack('>I', raw)[0]


def load_crc32c():
    try:
        import google_crc32c
    except ImportError:
        sys.exit("%s: the google-crc32c Python package is required" % PROG)
    if google_crc32c.implementation != 'c':
        sys.exit("%s: google-crc32c has no compiled extension here; its pure-Python fallback "
                 "is too slow for BAM-sized streams" % PROG)
    return google_crc32c.extend


def copy_and_checksum(infd, outfd, extend):
    crc, nbytes = 0, 0
    while True:
        chunk = os.read(infd, CHUNK_SIZE)
        if not chunk:
            return crc, nbytes
        crc = extend(crc, chunk)
        nbytes += len(chunk)
        view = memoryview(chunk)
        while view:
            view = view[os.write(outfd, view):]


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--size", type=int, required=True,
                        help="Expected size of the stream in bytes")
    parser.add_argument("--crc32c", required=True,
                        help="Expected CRC32C of the stream, base64-encoded big-endian as GCS reports it")
    args = parser.parse_args()
    if args.size < 0:
        sys.exit("%s: invalid --size %d" % (PROG, args.size))
    expected_crc = decode_crc32c(args.crc32c)
    extend = load_crc32c()

    try:
        crc, nbytes = copy_and_checksum(sys.stdin.fileno(), sys.stdout.fileno(), extend)
    except BrokenPipeError:
        sys.exit("%s: the downstream command stopped reading the stream" % PROG)

    problems = []
    if nbytes != args.size:
        problems.append("received %d bytes, expected %d" % (nbytes, args.size))
    if crc != expected_crc:
        problems.append("CRC32C is %s, expected %s" % (encode_crc32c(crc), args.crc32c))
    if problems:
        sys.exit("%s: ERROR: the stream does not match the source object (%s). "
                 "Output computed from this stream is incomplete or corrupt." % (PROG, "; ".join(problems)))
    sys.stderr.write("%s: verified %d bytes, CRC32C %s\n" % (PROG, nbytes, args.crc32c))


if __name__ == "__main__":
    main()
