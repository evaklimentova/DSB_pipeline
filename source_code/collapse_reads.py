#!/usr/bin/env python

import argparse
import pandas as pd
from cigar import Cigar


def parse_input():
    # function for parsing input parameters, returns it as dictionary
    parser = argparse.ArgumentParser(description='Collapse reads by read profile')
    parser.add_argument('--input', required=True, metavar='<filename>')
    parser.add_argument('--output', required=True, metavar='<filename2>')
    args = parser.parse_args()
    return vars(args)


def add_profile(sam):
    sam['profile'] = ''
    sam['ref_cov'] = 0

    for index, row in sam.iterrows():
        cigar = Cigar(row['CIGAR'])
        profile = []
        tmp = ""
        for item in list(cigar.items()):
            if item[1] in ['D', 'N']:
                profile += [str(Cigar(tmp).reference_length() + row['POS']) + '-' + item[1] + str(item[0])]
            elif item[1] == 'I':
                if row['SEQ'] != '*':
                    profile += [str(Cigar(tmp).reference_length() + row['POS']) + '-' + item[1] + str(item[0]) +
                                row['SEQ'][len(Cigar(tmp)):(len(Cigar(tmp)) + item[0])]]
                else:
                    profile += [str(Cigar(tmp).reference_length() + row['POS']) + '-' + item[1] + str(item[0]) + '*']
            tmp += str(item[0]) + item[1]
        sam.at[index, 'profile'] = 'XP:Z:' + '_'.join(profile)
        sam.at[index, 'ref_cov'] = cigar.reference_length()
    sam = sam.sort_values(by=['ref_cov'], ascending=False, ignore_index=True)
    sam = sam.drop(columns=['ref_cov'])
    return sam


def browse_file():
    print("Collapsing reads")

    arguments = parse_input()

    sam_file = pd.read_csv(arguments["input"], sep='\t',
                           names=['QNAME', 'FLAG', 'RNAME', 'POS', 'MAPQ', 'CIGAR', 'RNEXT', 'PNEXT', 'TLEN', 'SEQ',
                                  'QUAL', 'NH', 'HI', 'AS', 'nM', 'MD'], skiprows=4)
    if sam_file.empty:
        print("There are no reads with insertions or deletions.")
        return

    sam_file = add_profile(sam_file)
    sam_file["N"] = sam_file.groupby(['profile'])['QNAME'].transform('count')
    sam_file = sam_file.drop_duplicates(subset=['profile'], ignore_index=True)
    sam_file['N'] = "XN:i:" + sam_file['N'].astype(str)
    sam_file.to_csv(arguments['output'], sep='\t', index=False, header=False)


browse_file()
