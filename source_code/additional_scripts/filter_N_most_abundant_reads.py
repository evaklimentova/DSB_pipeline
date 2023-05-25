import pandas as pd

num = '12'
N = 100

sam_file = pd.read_csv("../Results_rerun_8_2/" + num + "-collapsed.sam", sep='\t',
                       names=['QNAME', 'FLAG', 'RNAME', 'POS', 'MAPQ', 'CIGAR', 'RNEXT', 'PNEXT', 'TLEN', 'SEQ', 'QUAL', 'NH', 'HI', 'AS', 'nM', 'MD', 'XN', 'XP'])


sam_file['iXN'] = sam_file['XN'].apply(lambda x: int(str(x).split(':')[-1]))

sam_file = sam_file.sort_values(by='iXN', ascending=False).reset_index()
sam_file = sam_file[:N]
print(sam_file)
sam_file = sam_file.drop(columns=['iXN', 'index'])

sam_file.to_csv(num + "_collapsed_top_" + str(N) + ".sam", index=False, header=False, sep='\t')

