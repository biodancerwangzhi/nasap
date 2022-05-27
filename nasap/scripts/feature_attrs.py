import os, fire
import pandas as pd
import pyBigWig

import sys

# sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))+ '/libs')
script_dir = os.path.abspath(os.path.join(os.path.abspath(__file__), ".."))
lib_dir = os.path.abspath(os.path.join(os.path.abspath(__file__), '../../libs') )
sys.path.append(lib_dir)

## 这个模块可以单独使用
# nasap feature_assign --gtf=./a.gtf --forward=./a.bw --reverse=./b.bw --outpt=./tmp_out

from py_ext import json2dic, dic2json, read_tmp_variable
from parse_gtf import get_gene_df, get_gene_merge_exon_dic, gene_df2dic

def get_total_bases(forward_bw, reverse_bw):
  return pyBigWig.open(forward_bw).header()['sumData'] + pyBigWig.open(reverse_bw).header()['sumData']

def density2rpkm(density, total_bases):
  return density / total_bases * 1e6

def pausing_index( chr, start, end, strand,  bw ):
  # 获取 表达的蛋白， lincRNA
  if strand == '+':
    pp = bw.stats( chr, start, start + 1000)[0]
    gb = bw.stats( chr, start +1000, end)[0]
  else:
    pp = bw.stats( chr, end -1000, end)[0]
    gb = bw.stats( chr, start, end -1000)[0]
  pp_count, gb_count = int(pp* 1000), int(gb* (end-start +1000) )
  if gb == 0:
    #pi = float('nan')
    pi= 'gene body count zero'
  else:
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

def get_attrs(gene_range_dic, forward_bw, reverse_bw, total_bases):
  gene_rpkm_dic, gene_baseCount_dic, gene_pi_dic, gene_ei_dic, gene_pp_count_dic, gene_gb_count_dic = {}, {}, {}, {},{}, {}
  for strand in ['+', '-']:
    strand_gene_range_dic = gene_range_dic[strand]
    bw = pyBigWig.open( {'+': forward_bw, '-': reverse_bw}[strand])
    for gene_name, gene_range in strand_gene_range_dic.items():
      c, s, e = gene_range[0], gene_range[1], gene_range[2]
      try:
        density = bw.stats(c, s, e)[0]
        rpkm = density2rpkm(density, total_bases)
        base_count = int( density * (e - s) )
        if e - s < 2300:
          pi, prr='gene too short', 'gene too short'
          continue
        else:
          pi, pp_count, gb_count = pausing_index(c,s, e, strand, bw)
          prr = elongation_index(c,s, e, strand, bw)

      except:
        continue
      gene_rpkm_dic[gene_name] = rpkm
      gene_baseCount_dic[gene_name] = base_count
      gene_pi_dic[gene_name] = pi
      gene_ei_dic[gene_name] = prr
      gene_pp_count_dic[gene_name] = pp_count
      gene_gb_count_dic[gene_name] = gb_count
    # dic2json(gene_rpkm_dic, rpkm_json_output)
  return gene_rpkm_dic, gene_baseCount_dic, gene_pi_dic, gene_ei_dic, gene_pp_count_dic, gene_gb_count_dic

# count, RPKM, PI, EI, Exon/intron density
def main(gtf, forward_bw, reverse_bw, output_root='./tmp_output/' ):
  total_bases = get_total_bases(forward_bw, reverse_bw)

  gene_df = get_gene_df( gtf )
  filter_gene_df = gene_df[ gene_df['genetype'].isin(['protein_coding', 'lincRNA', 'lncRNA']) ]
  gene_range_dic = gene_df2dic(filter_gene_df)
  # count, RPKM, promoter, geneBody count, PI, EI

  gene_rpkm_dic, gene_baseCount_dic, gene_pi_dic, gene_ei_dic, gene_pp_count_dic, gene_gb_count_dic = get_attrs(gene_range_dic, forward_bw, reverse_bw, total_bases)
  # print( list(gene_rpkm_dic.items())[:5] )
  gene_attrs_df = pd.concat({'rpkm':pd.Series(gene_rpkm_dic), 'count': pd.Series(gene_baseCount_dic),
   'pp_count': pd.Series(gene_pp_count_dic), 'gb_count': pd.Series(gene_gb_count_dic),
   'pi': pd.Series(gene_pi_dic), 'ei': pd.Series(gene_ei_dic)}, axis=1)

  gene_attrs_df.to_csv(output_root+'csv/feature_attrs.csv')

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
      print('no gene', gene_name)
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
  exon_intron_df.T.to_csv(output_root+'csv/exon_intron_ratio.csv')

if __name__ == '__main__':
  fire.Fire( main ) # pipeline )