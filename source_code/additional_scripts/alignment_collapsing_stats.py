import pandas as pd

num = '12'

sam_file = pd.read_csv("../Results_2020/" + num + "-collapsed.sorted.sam", sep='\t',
                       names=['QNAME', 'FLAG', 'RNAME', 'POS', 'MAPQ', 'CIGAR', 'RNEXT', 'PNEXT', 'TLEN', 'SEQ', 'QUAL', 'NH', 'HI', 'AS', 'nM', 'MD', 'XP', 'XN'])


sam_file['XN'] = sam_file['XN'].apply(lambda x: int(str(x).split(':')[-1]))
sam_file['XP'] = sam_file['XP'].apply(lambda x: str(x).split(':')[-1])

result = sam_file[['XP', 'XN']]
result = result.sort_values(by='XN', ascending=False).reset_index()

total = result['XN'].sum()
result['cumsum'] = result['XN'].cumsum()
result['cumsum'] = round(100 * result['cumsum']/total, 2)
result['percentage'] = round(100 * result['XN']/total, 2)

result[['XP', 'XN', 'percentage', 'cumsum']].to_csv(num + "_rerun_collapsed_table.csv", index=False,
    header=['group', 'N of reads', 'percentage of total', 'cumulative percentage'])

print(result, total)
