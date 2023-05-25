#!/usr/bin/env python

import pysam
import sys
import argparse


def add_one(start, num, array):
    # adds one to each element of array from position start to start + num - 1
    for i in range(num):
        array[start + i] += 1


def parse_cigar(cigar):
    # returns cigar string parsed as an array
    # every element is number + identifier
    start = 0
    result = []
    for i in range(len(cigar)):
        if cigar[i].isalpha():
            result.append(cigar[start:i + 1])
            start = i + 1
    return result


def count_one_cigar(cigar, start, array):
    # adds sequence coverage to 'array' from one cigar string
    parts = parse_cigar(cigar)
    for elem in parts:
        if elem[-1] == "M":
            add_one(start, int(elem[:-1]), array)
            start += int(elem[:-1])
        elif elem[-1] == "I":
            continue
        elif elem[-1] == "D":
            start += int(elem[:-1])
        elif elem[-1] == "S":
            continue
        elif elem[-1] == "N":
            start += int(elem[:-1])
        elif elem[-1] == "H":
            continue
        elif elem[-1] == "P":
            continue


def count_read_indel(read, file_ins, file_del):
    # writes to file_ins information about insertions:
    # position, sequence and from which end the insertion starts
    # writes to file_del information about deletion : position and length
    start_ref = read.get_reference_positions()[0]

    cigar_string = parse_cigar(read.cigarstring)
    md_index = 0

    insertion_end = "S5"
    current_pos = 0

    for elem in cigar_string:

        if elem[-1] == "I":
            # insertion
            seq = read.query_sequence[current_pos:(current_pos + int(elem[:-1]))]
            line = '{:^25} {:^10} {:^10} {:^30} {:^5}'.format(read.reference_name, start_ref + 1, int(elem[:-1]), seq, "I")
            file_ins.write(line+"\n")
            current_pos += int(elem[:-1])

        elif elem[-1] == "S":
            # soft clip
            seq = read.query_sequence[current_pos:(current_pos + int(elem[:-1]))]
            line = '{:^25} {:^10} {:^10} {:^30} {:^5}'.format(read.reference_name, start_ref + 1, int(elem[:-1]), seq, insertion_end)
            file_ins.write(line+"\n")
            current_pos += int(elem[:-1])

        elif elem[-1] == "M":
            # match or mismatch
            start_ref += int(elem[:-1])
            current_pos += int(elem[:-1])

        elif elem[-1] == "D":
            # deletion
            line = '{:^25} {:^10} {:^10} {:^5}'.format(read.reference_name, start_ref + 1, int(elem[:-1]), "D")
            file_del.write(line+"\n")
            start_ref += int(elem[:-1])

        elif elem[-1] == "N":
            # For mRNA-to-genome alignment, an N operation represents an intron.
            # For other types of alignments, the interpretation of N is not defined.
            line = '{:^25} {:^10} {:^10} {:^5}'.format(read.reference_name, start_ref + 1, int(elem[:-1]), "N")
            file_del.write(line + "\n")
            start_ref += int(elem[:-1])

        # "H" and "P" doesn't consume any of the two sequences

        insertion_end = "S3"


def print_dictionary(dic, output_name):
    # prints 'dic' dictionary to the file output_name
    file = open(output_name, "w")

    line = '{:^25} {:^10} {:^10}'.format("Name", "Position", "Score")
    file.write(line+"\n")

    for name, array in dic.items():
        for i in range(len(array)):
            line2 = '{:^25} {:^10} {:^10}'.format(name, str(i + 1), str(array[i]))
            file.write(line2+"\n")
        file.write("\n")


def parse_input():
    # function for parsing input parameters, returns it as dictionary
    parser = argparse.ArgumentParser(description='Count coverage and print deletions and insertions')
    parser.add_argument('--input', required=True, metavar='<filename>')
    parser.add_argument('--outputCov', default="coverage.txt", metavar='<filename2>')
    parser.add_argument('--outputIns', default="insertions.txt", metavar='<filename3>')
    parser.add_argument('--outputDel', default="deletions.txt", metavar='<filename4>')
    args = parser.parse_args()
    return vars(args)


def browse_file():
    # main function of the program
    # goes through SAM file, through all genes and sequences
    # counts it's coverage and prints the insertions + deletions
    
    print("Running 'count reference coverage' script.")

    arguments = parse_input()

    try:
        samfile = pysam.AlignmentFile(arguments["input"], "r")
    except IOError:
        sys.stderr.write("'IOError': No such file: '" + arguments["input"] + "'\n")
        return

    names = samfile.references
    lengths = samfile.lengths
    array_dict = {}

    for i in range(len(names)):
        array_dict[names[i]] = [0] * lengths[i]

    file_ins = open(arguments["outputIns"], "w")
    file_del = open(arguments["outputDel"], "w")
    line = '{:^25} {:^10} {:^10} {:^30} {:^5}'.format("Name", "Position", "Length", "Sequence", "End")
    file_ins.write(line+"\n")
    line = '{:^25} {:^10} {:^10} {:^5}'.format("Name", "Position", "Length", "Type")
    file_del.write(line+"\n")

    for read in samfile.fetch():
        ref = read.reference_name
        if ref is None:
            continue
        start = read.get_reference_positions()[0]
        cigar = read.cigarstring
        if cigar == "=":
            add_one(start, read.query_length, array_dict[ref])
            continue
        elif cigar == "X":
            continue

        count_one_cigar(cigar, start, array_dict[ref])
        count_read_indel(read, file_ins, file_del)

    print_dictionary(array_dict, arguments["outputCov"])

    file_ins.close()
    file_del.close()
    samfile.close()


browse_file()

