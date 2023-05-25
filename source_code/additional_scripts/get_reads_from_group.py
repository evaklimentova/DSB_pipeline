import pandas as pd
from cigar import Cigar

num = '8'
XP = '116-D8'

sam = pd.read_csv("../Results_rerun_11_2/" + num + ".Aligned.out.sam", sep='\t',
                       names=['QNAME', 'FLAG', 'RNAME', 'POS', 'MAPQ', 'CIGAR', 'RNEXT', 'PNEXT', 'TLEN', 'SEQ', 'QUAL', 'NH', 'HI', 'AS', 'nM', 'MD'], skiprows=4)

print(sam)

sam['profile'] = ''

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
    sam.at[index, 'profile'] = '_'.join(profile)

sam = sam[sam['profile'] == XP]

print(sam)

sam.to_csv(num + "_" + XP + ".sam", index=False, header=False, sep='\t')

