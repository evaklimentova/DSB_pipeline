#!/usr/bin/env python

import pysam
import argparse
import sys


def not_match(read):
    # return True if read contains D (deletion), I (insertion) or N (long deletion)
    cigar_string = read.cigartuples

    for c in cigar_string:
        if c[0] == 1 or c[0] == 2 or c[0] == 3:
            return True

    return False


def parse_input():
    # function for parsing input parameters, returns it as dictionary
    parser = argparse.ArgumentParser(description='Filter reads with 100% match')
    parser.add_argument('--input', required=True, metavar='<filename>')
    parser.add_argument('--output', required=True, metavar='<filename2>')
    args = parser.parse_args()
    return vars(args)


def browse_file():
    # function goes through the input sam file and filters reads that doesn't match with the original sequence

    print("Filtering reads")

    arguments = parse_input()

    try:
        samfile_original = pysam.AlignmentFile(arguments["input"], "r")
    except IOError:
        sys.stderr.write("'IOError': No such file: '" + arguments["input"] + "'\n")
        return

    samfile_new = pysam.AlignmentFile(arguments["output"], "w", template=samfile_original)

    for read in samfile_original.fetch():
        if read.cigarstring is None:
            continue
        if not_match(read):
            samfile_new.write(read)

    samfile_original.close()
    samfile_new.close()


browse_file()
