import pandas as pd


def filter_insertions(num):
    sam = pd.read_csv("../Results_rerun_11_2/" + num + "-collapsed.sorted.sam", sep='\t',
                           names=['QNAME', 'FLAG', 'RNAME', 'POS', 'MAPQ', 'CIGAR', 'RNEXT', 'PNEXT', 'TLEN', 'SEQ', 'QUAL', 'NH', 'HI', 'AS', 'nM', 'MD', 'XP', 'XN'])

    csv = pd.read_csv("../Results_rerun_11_2/Reads_stats/" + num + "_rerun_collapsed_table.csv")

    sam_new = sam[sam['XP'].str.contains("I")]
    csv_new = csv[csv['group'].str.contains("I")]
    csv_new = csv_new.drop(columns=['cumulative percentage'])
    total = csv_new['N of reads'].sum()
    csv_new['percentage from insertions'] = round(100 * csv_new['N of reads'] / total, 2)

    sam_new.to_csv(num + "-onlyinsertions-collapsed.sam", index=False, header=False, sep='\t')
    csv_new.to_csv(num + "-onlyinsertions_collapsed_table.csv", index=False)


for n in range(1, 13):
    filter_insertions(str(n))
