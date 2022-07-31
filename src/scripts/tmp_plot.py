import re
import pandas as pd
from scipy import stats

from parse_gtf import get_gene_df, get_gene_merge_exon_dic, gene_df2dic
from plot import boxplot, scatterplot
gtf = '../../../data/Homo_sapiens.GRCh38.93.gtf'
gene_df = get_gene_df( gtf )
filter_gene_df = gene_df[ gene_df['genetype'].isin(['protein_coding', 'lincRNA', 'lncRNA']) ]
print( filter_gene_df[filter_gene_df['gene_name']=='MIR1302-2HG'])
print( filter_gene_df.shape )
print( set(filter_gene_df['genetype']) )
gene_range_dic = gene_df2dic(filter_gene_df)

chr_list = list( set( filter_gene_df['chrom'] ) )

sort_chr_list = [chrom for chrom in chr_list if re.match( r'chr\d+', chrom)]
sort_chr_list.sort(key=lambda arr: (arr[:3], int(arr[3:])))
if 'chrM' in chr_list:  sort_chr_list.append('chrM')
if 'chrMT' in chr_list:  sort_chr_list.append('chrMT')
if 'chrX' in chr_list:  sort_chr_list.append('chrX')
if 'chrY' in chr_list:  sort_chr_list.append('chrY')

output_root = '../../../nasap_GSM3618143/'
attr_df = pd.read_csv(output_root + 'csv/feature_attrs.csv', index_col=0)
rpkm_series = attr_df['rpkm']
gene_rpkm_dic = rpkm_series.to_dict()
print( dict(list(gene_rpkm_dic.items())[:3]) )

rpkm_list = []
for chrom in sort_chr_list:
  gene_list = filter_gene_df[ filter_gene_df['chrom'] == chrom ]['gene_name']

  cur_rpkm_list = []
  for gene in gene_list:
    try:
      cur_rpkm_list.append(gene_rpkm_dic[gene])
    except:
      continue
  rpkm_list.append( cur_rpkm_list )


boxplot(rpkm_list, 'chromsome genes RPKM', 'chromsome', 'RPKM', sort_chr_list, output_root + 'imgs/chr_rpkm.png' )


exon_intron_df = pd.read_csv(output_root + 'csv/exon_intron_ratio.csv', index_col=0)
filter_exon_intron_df = exon_intron_df[(exon_intron_df['exon'] > 0) & (exon_intron_df['intron'] >0) ]

slope, intercept, r_value, p_value, std_err = stats.linregress(filter_exon_intron_df['exon'],filter_exon_intron_df['intron'])
text= 'y = {:4.2e} x + {:4.2e}; \nR^2= {:2.2f}'.format(slope, intercept, r_value*r_value)
scatterplot( filter_exon_intron_df['exon'], filter_exon_intron_df['intron'], 'Exon/Intron ratio', 'Exon', 'Intron', text, output_root+'imgs/exon_intron_ratio.png' )
