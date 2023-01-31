import os, fire, re, sys, pickle
# sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))+ '/libs')
script_dir = os.path.abspath(os.path.join(os.path.abspath(__file__), ".."))
lib_dir = os.path.abspath(os.path.join(os.path.abspath(__file__), '../../libs') )
sys.path.append(lib_dir)
from collections import Counter

import numpy as np
import pandas as pd
from scipy import stats

import pyBigWig

from parse_gtf import get_gene_df, get_gene_merge_exon_dic, gene_df2dic
from plot import chr_rpkm_boxplot, expressed_pie, expressed_distribution, exon_intron_scatter, proximal_pausing_sites_distribution

## 这个模块可以单独使用
# nasap feature_assign --gtf=./a.gtf --forward=./a.bw --reverse=./b.bw --outpt=./tmp_out

def get_total_bases(forward_bw, reverse_bw):
  return pyBigWig.open(forward_bw).header()['sumData'] + pyBigWig.open(reverse_bw).header()['sumData']

def density2rpkm(density, total_bases):
  return density / total_bases * 1e6

def find_true_tss(chr, left, right, bw, strand):
  # ref: Defining the Status of RNA Polymerase at Promoters
  # Reannotation of TSSs
  # annotated TSS ± 150, highest number of reads within this window
  # If two bases both the same maximum, the most 5′-base was used as the TSS
  v_list = bw.values(chr, left, right)
  max_index_list = np.argwhere(v_list == np.amax(v_list)).flatten().tolist()
  if strand == '+':
    max_index = max_index_list[0]
  else:
    max_index = max_index_list[-1]
  return left + max_index

def pausing_index( chr, start, end, strand,  bw, tss_up, tss_down ):
  # 获取 表达的蛋白， lincRNA
  if strand == '+':
    true_tss = find_true_tss(chr, start-tss_up, start + tss_down, bw, strand)
    pp = bw.stats( chr, true_tss-tss_up, true_tss + tss_down)[0]
    gb = bw.stats( chr, true_tss +tss_down, end)[0]
  else:
    true_tss = find_true_tss(chr, end -tss_down, end + tss_up, bw, strand)
    pp = bw.stats( chr, true_tss -tss_down, true_tss + tss_up)[0]
    gb = bw.stats( chr, start, true_tss - tss_down)[0]
  ppl = tss_up + tss_down
  gbl = end - start - ppl
  pp_count, gb_count = int(pp * ppl ), int(gb* gbl )
  if gb == 0: gb = 1e-6
  # pi = float('nan')
  # pi= 'gene body count zero'
  pi = pp/gb
  return pi, pp_count, gb_count

def elongation_index( chr, start, end, strand,  bw ):
  # 获取 表达的蛋白， lincRNA
  if strand == '+':
    pp = bw.stats( chr, start - 100, start + 300)[0]
    gb = bw.stats( chr, start + 300, start + 2000)[0]
  else:
    pp = bw.stats( chr, end - 300, end + 100)[0]
    gb = bw.stats( chr, end - 2000, end - 300)[0]
  if pp == 0:
    #
    prr = 'proximal promoter count zero'
  else:
    prr = gb/pp
  return prr

def get_attrs(gene_range_dic, forward_bw, reverse_bw, total_bases, tss_up, tss_down):
  gene_rpkm_dic, gene_baseCount_dic, gene_rpm_dic =  {},{}, {}
  gene_pi_dic, gene_ei_dic = {}, {}
  gene_pp_count_dic, gene_gb_count_dic = {}, {}
  for strand in ['+', '-']:
    strand_gene_range_dic = gene_range_dic[strand]
    bw = pyBigWig.open( {'+': forward_bw, '-': reverse_bw}[strand])
    for gene_name, gene_range in strand_gene_range_dic.items():
      c, s, e = gene_range[0], gene_range[1], gene_range[2]
      try:
        density = bw.stats(c, s, e)[0]
        rpkm = density2rpkm(density, total_bases)
        base_count = int( density * (e - s) )
        rpm = base_count * 1e6 /total_bases
        if e - s < 2300:
          pi, prr='gene too short', 'gene too short'
          continue
        else:
          pi, pp_count, gb_count = pausing_index(c,s, e, strand, bw, tss_up, tss_down)
          prr = elongation_index(c,s, e, strand, bw)

      except:
        continue
      gene_rpkm_dic[gene_name] = rpkm
      gene_baseCount_dic[gene_name] = base_count
      gene_rpm_dic[gene_name] = rpm
      gene_pi_dic[gene_name] = pi
      gene_ei_dic[gene_name] = prr
      gene_pp_count_dic[gene_name] = pp_count
      gene_gb_count_dic[gene_name] = gb_count
    # dic2json(gene_rpkm_dic, rpkm_json_output)
  return gene_rpkm_dic, gene_baseCount_dic, gene_rpm_dic, \
  gene_pi_dic, gene_ei_dic, \
  gene_pp_count_dic, gene_gb_count_dic

def coding_vs_linc(protein_data, linc_data, output_root):
  protein_num, linc_num = len(protein_data), len(linc_data)
  big_labels = [f'Coding\n({protein_num})', f'LncRNA\n({linc_num})']
  expressed_protein_num, expressed_linc_num = len(protein_data[protein_data>0]), len(linc_data[linc_data>0])
  expressed_protein_ratio = np.round( expressed_protein_num/protein_num, 2) * 100
  expressed_linc_ratio = np.round( expressed_linc_num/linc_num, 2) * 100
  big_sizes = [protein_num, linc_num]

  small_labels = [f'{expressed_protein_ratio}% Protein-coding genes\nExpressed',
    f'{100-expressed_protein_ratio}% Protein-coding genes\nNot Expressed',
    f'{expressed_linc_ratio}% LncRNA genes\nExpressed',
    f'{100-expressed_linc_ratio}% LncRNA genes\nNot Expressed'
  ]
  small_sizes = [expressed_protein_num, protein_num-expressed_protein_num, expressed_linc_num, linc_num - expressed_linc_num]

  # bigger
  expressed_pie_data = { 'big_sizes': big_sizes,
  'big_labels': big_labels,
  'small_sizes': small_sizes,
  'small_labels': small_labels}
  expressed_pie(expressed_pie_data, output_root)
  # with open('expressed_pie_data.pickle', 'wb') as f:
  #   pickle.dump(expressed_pie_data, f)

  expressed_dist_data = {'protein_data': protein_data, 'linc_data': linc_data }
  expressed_distribution(expressed_dist_data, output_root)
  # with open('expressed_dist_data.pickle', 'wb') as f:
  #   pickle.dump(expressed_dist_data, f)

# count, RPKM, PI, EI, Exon/intron density
def main(gtf, forward_bw, reverse_bw, output_root, tss_up=150, tss_down=150 ):
  # 这里 多搞两个 参数用 tss_up, tss_down 去定义 pp区域，
  # tes的部分暂时不提，就用tss_right 到 tes 定义gb区域。
  total_bases = get_total_bases(forward_bw, reverse_bw)
  gene_df = get_gene_df( gtf )
  filter_gene_df = gene_df[ gene_df['genetype'].isin(['protein_coding', 'lincRNA', 'lncRNA']) ]
  gene_range_dic = gene_df2dic(filter_gene_df)
  # count, RPKM, promoter, geneBody count, PI, EI

  gene_rpkm_dic, gene_baseCount_dic, gene_rpm_dic, \
  gene_pi_dic, gene_ei_dic, \
  gene_pp_count_dic, gene_gb_count_dic = get_attrs(gene_range_dic, forward_bw, reverse_bw, total_bases, tss_up, tss_down)
  # print( list(gene_rpkm_dic.items())[:5] )

  attr_list = ['baseCount', 'rpkm', 'rpm', 'pp_count', 'gb_count', 'pi', 'ei']
  filter_gene_type = set(filter_gene_df['genetype'])
  for gene_type in filter_gene_type:
    gene_list = filter_gene_df[filter_gene_df['genetype'] == gene_type]['gene_name']
    for attr in attr_list:
      # type_attr_dic = {g: globals()['gene_' + attr +'_dic'][g] for g in gene_list}
      tmp_dic = eval('gene_' + attr +'_dic')
      tmp_list = list( tmp_dic.keys() )
      type_attr_dic = {g: tmp_dic[g] for g in gene_list if g in tmp_list}
      pd.DataFrame( {attr: pd.Series(type_attr_dic)} ).to_csv(output_root+'csv/' + gene_type + '_'+ attr + '.csv')


  gene_attrs_df = pd.concat({'rpkm':pd.Series(gene_rpkm_dic),
    'count': pd.Series(gene_baseCount_dic),
    'rpm': pd.Series(gene_rpm_dic),
    'pp_count': pd.Series(gene_pp_count_dic),
    'gb_count': pd.Series(gene_gb_count_dic),
    'pi': pd.Series(gene_pi_dic),
    'ei': pd.Series(gene_ei_dic)}, axis=1)

  gene_attrs_df.to_csv(output_root+'csv/%s_feature_attrs.csv'%gene_type)

  chr_list = list( set( filter_gene_df['chrom'] ) )

  sort_chr_list = [chrom for chrom in chr_list if re.match( r'chr\d+', chrom)]
  sort_chr_list.sort(key=lambda arr: (arr[:3], int(arr[3:])))
  # if 'chrM' in chr_list:  sort_chr_list.append('chrM')
  # if 'chrMT' in chr_list:  sort_chr_list.append('chrMT')
  if 'chrX' in chr_list:  sort_chr_list.append('chrX')
  if 'chrY' in chr_list:  sort_chr_list.append('chrY')
  rpkm_list = []
  for chrom in sort_chr_list:
    gene_list = filter_gene_df[ filter_gene_df['chrom'] == chrom ]['gene_name']

    cur_rpkm_list = []
    for gene in gene_list:
      try:
        cur_rpkm_list.append( np.log10(gene_rpkm_dic[gene]) )
      except:
        continue
    rpkm_list.append( cur_rpkm_list )
  chr_rpkm_data = {'rpkm_list': rpkm_list, 'sort_chr_list': sort_chr_list}
  chr_rpkm_boxplot( chr_rpkm_data, output_root )
  # with open('chr_rpkm_data.pickle', 'wb') as f:
  #   pickle.dump(chr_rpkm_data, f)

  # protein vs linc
  if ('protein_coding' in filter_gene_type) and ('lincRNA' in filter_gene_type):
    protein_data = pd.read_csv(output_root + 'csv/protein_coding_rpkm.csv').rpkm
    linc_data = pd.read_csv( output_root + 'csv/lincRNA_rpkm.csv').rpkm
    coding_vs_linc(protein_data, linc_data, output_root)
  if ('protein_coding' in filter_gene_type) and ('lncRNA' in filter_gene_type):
    protein_data = pd.read_csv(output_root + 'csv/protein_coding_rpkm.csv').rpkm
    linc_data = pd.read_csv( output_root + 'csv/lncRNA_rpkm.csv').rpkm
    coding_vs_linc(protein_data, linc_data, output_root)


  # Exon/Intron
  protein_exon_range_dic = get_gene_merge_exon_dic(gtf)
  forward_bw_signal = pyBigWig.open(  forward_bw )
  reverse_bw_signal = pyBigWig.open(  reverse_bw )
  exon_intron_ratio_dic = {}
  for gene_name, gr_dic in protein_exon_range_dic.items():
    chrom, strand, exon_range_list = gr_dic['chrom'], gr_dic['strand'], gr_dic['exon']
    try:
      gene_range = gene_range_dic[strand][gene_name]
    except:
      # print('no gene', gene_name)
      os.sys.exit(0)
      continue
    gene_start, gene_end = gene_range[1], gene_range[2]
    if strand == '+':
      bw = forward_bw_signal
    else:
      bw = reverse_bw_signal
    gene_length = gene_end - gene_start
    try:
      gene_total = bw.stats(chrom, gene_start, gene_end)[0] * gene_length
    except:
      # print(gene_name, chrom, gene_start, gene_end)
      continue
    exon_length, exon_total = 0, 0
    for exon_range in exon_range_list:
      cur_length = (exon_range[1] - exon_range[0])
      exon_length+=cur_length
      try:
        exon_total +=( cur_length* bw.stats(chrom, exon_range[0], exon_range[1])[0] )
      except:
        continue
    intron_total = gene_total - exon_total
    intron_length = gene_length - exon_length

    if intron_total <= 0:
      # print( "warning: intron density <=0")
      exon_density, intron_density,  ratio = exon_total/exon_length, 0, 0
    else:
      exon_density, intron_density = exon_total/exon_length, intron_total/intron_length
      ratio = exon_density/intron_density
    exon_intron_ratio_dic[gene_name] = {'exon': exon_density, 'intron': intron_density, 'ratio': ratio}

  exon_intron_df = pd.DataFrame(exon_intron_ratio_dic)
  exon_intron_df = exon_intron_df.T
  exon_intron_df.to_csv(output_root+'csv/exon_intron_ratio.csv')

  filter_exon_intron_df = exon_intron_df[(exon_intron_df['exon'] > 0) & (exon_intron_df['intron'] >0) ]

  '''
  slope, intercept, r_value, p_value, std_err = stats.linregress(filter_exon_intron_df['exon'],filter_exon_intron_df['intron'])
  text= 'y = {:4.2e} x + {:4.2e}; \nR^2= {:2.2f}; p= {:2.2f}'.format(slope, intercept, r_value*r_value, p_value)
  scatterplot( filter_exon_intron_df['exon'], filter_exon_intron_df['intron'], 'Exon/Intron ratio', 'Exon', 'Intron', text, output_root+'imgs/exon_intron_ratio.png' )
  '''

  # remove 两端的异常值
  exon_min, exon_max = filter_exon_intron_df['exon'].quantile( [0.1, 0.9] )
  intron_min, intron_max = filter_exon_intron_df['intron'].quantile( [0.1, 0.9] )
  filter_exon_intron_df1 = filter_exon_intron_df[(filter_exon_intron_df['exon'] > exon_min) & (filter_exon_intron_df['intron'] >intron_min) & (filter_exon_intron_df['exon'] < exon_max) & (filter_exon_intron_df['intron'] <intron_max) ]
  # 计算拟合值
  slope, intercept, r_value, p_value, std_err = stats.linregress( filter_exon_intron_df1['exon'], filter_exon_intron_df1['intron'] )
  text= 'y = {:4.2e} x + {:4.2e}; \nR^2= {:2.2f}; p= {:2.2e}'.format(slope, intercept, r_value*r_value, p_value)

  tmp_dir = output_root + 'csv/preprocess_report.csv'
  if os.path.exists(tmp_dir):
    with open(tmp_dir, encoding="utf-8",mode="a") as file:
        file.write('mRNA_contamination' + ',' + str(slope) + '\n')
  else:
    with open(output_root + 'csv/preprocess_report.csv', 'w') as file:
      file.write( 'mRNA_contamination' + ',' + str(slope) + '\n' )

  exon_intron_data = { 'filter_exon_intron_df1': filter_exon_intron_df1, 'text': text}
  exon_intron_scatter( exon_intron_data, output_root )
  # with open('exon_intron_data.pickle', 'wb') as f:
  #   pickle.dump(exon_intron_data, f)


  # proximal pause site
  pp_site_file = open(output_root + 'bed/proximal_pause_sites.bed', 'w')
  if ('protein_coding' in filter_gene_type):
    protein_data = pd.read_csv(output_root + 'csv/protein_coding_pi.csv', index_col=0)
    filter_genes = protein_data[pd.to_numeric( protein_data.pi, errors='coerce' ) > 2].index
    site_index_list = []
    for i, row in filter_gene_df[filter_gene_df['gene_name'].isin(filter_genes)].iterrows():
      chr = row['chrom']
      gene_name = row['gene_name']
      if row['strand'] == '+':
        tss = row['start']
        # density_list = [forward_bw_signal.stats(chr, tss+j*50, tss+(j+1)*50 )[0] for j in range(7)]
        signals = forward_bw_signal.values(chr, tss, tss+350 )
        signal_median = np.median( signals )
        if signal_median < 0.5: continue
        max_idx = np.argmax( signals )
        if signals[max_idx]/signal_median < 5: continue
        ln = chr + '\t' + str( tss + max_idx ) + '\t' +str( tss + max_idx+1 ) + '\t' + gene_name + '\t' + str(max_idx ) + '\n'
      else:
        tss = row['end']
        signals = forward_bw_signal.values(chr, tss-350, tss )
        signal_median = np.median( signals )
        if signal_median < 0.5: continue
        max_idx = np.argmax( signals )
        if signals[max_idx]/signal_median < 5: continue
        # density_list = [reverse_bw_signal.stats(chr, tss-(j+1)*50, tss-j*50 )[0] for j in range(7)]
        ln = chr + '\t' + str( tss - (350-max_idx+1) ) + '\t' + str( tss -  (350-max_idx) ) + '\t' + gene_name + '\t' +  str(max_idx ) + '\n'
      # print(ln)
      pp_site_file.write(ln)
      site_index_list.append( max_idx )
  # print( site_index_list )
  pp_site_file.close()
  if len(site_index_list) > 0:
    proximal_pausing_sites_distribution(site_index_list, output_root)
  # with open('site_index_list.pickle', 'wb') as f:
  #   pickle.dump(site_index_list, f)


if __name__ == '__main__':
  fire.Fire( main )