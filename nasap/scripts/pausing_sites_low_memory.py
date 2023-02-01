import os, fire, pickle

import numpy as np
import pandas as pd
import pyBigWig
from collections import defaultdict

from parse_gtf import sort_chr_list
from plot import global_distribution_lollipopplot

def get_bw_pause_sites(para_list):
  [bw_dir, chr, size, isForward] = para_list
  win_side= 100 # 左右100，总共201
  bw = pyBigWig.open(bw_dir)

  start, end = 0, size
  gene_depth_array = bw.values(chr, start, end, numpy=True)
  # 获取小窗的window np arrays
  suit_index = np.where(gene_depth_array>=4)[0]
  suit_index = suit_index[suit_index >win_side]
  suit_index = suit_index[suit_index < end-win_side]

  site_list = []
  count_list = []

  for i in suit_index:
    window=gene_depth_array[i-win_side:i+win_side+1]
    win_count_max=window.max()
    if win_count_max <4: continue
    # 201bp小窗 中点
    if window[win_side] <4: continue
    if window[win_side] != win_count_max: continue

    win_count_mean=window.mean()
    # win_count_median=np.median( window )
    win_count_std=window.std()
    if win_count_max>=(win_count_mean + win_count_std*3):
      # return 1
      site_list.append(i)
      count_list.append( win_count_max)
    # return 0
  # suit_value = np.array( list( map(judge_window, suit_index) ) )
  # chr_sites = suit_index[suit_value >0]
  if isForward:
    chr_strand = chr + '+'
  else:
    chr_strand = chr + '-'
  # print( {chr_strand:  {'site': site_list, 'count': count_list} } )
  return {chr_strand:  {'site': site_list, 'count': count_list} }

# 思路:
# 1 获取所有染色体
# 2 单个染色体运行并生成结果
# 3 写入识别的暂停位点
# 4 画图
def main(forward_bw, reverse_bw, output_root='./tmp_output/', cores=1):
  filter200_chr_site_dic = {}
  for file_type, bw_dir in {'forward': forward_bw, 'reverse': reverse_bw}.items():
    isForward = True
    if file_type == 'reverse': isForward = False
    bw_value = pyBigWig.open(bw_dir)
    chr_sizes = {chr:size for chr, size in bw_value.chroms().items() }
    chr_size_list = list(chr_sizes.items())
    for chr_size_tuple in chr_size_list:
      # print( chr_size_tuple)
      # print( chr_size_tuple )
      chrom, size = chr_size_tuple[0], chr_size_tuple[1]
      chr_strand_site_dic = get_bw_pause_sites([bw_dir, chrom, size, isForward])
      # print( chr_strand_site_dic)

      for chr_strand, site_count_dic in chr_strand_site_dic.items():
        sort_site_list = site_count_dic['site']
        count_list = site_count_dic['count']

        final_index = len(sort_site_list) - 1
        if final_index == -1: continue
        if final_index == 0:
          filter200_chr_site_dic[chr_strand] = {'site': sort_site_list}
          filter200_chr_site_dic[chr_strand]['count'] =  count_list
          continue
        if final_index == 1:
          if (sort_site_list[1] - sort_site_list[0]) > 200:
            filter200_chr_site_dic[chr_strand] = {'site': sort_site_list}
            filter200_chr_site_dic[chr_strand]['count'] =  count_list
          continue

        tmp_list = []
        tmp_count_list = []
        last_index, cur_index, next_index =0, 1, 2
        if (sort_site_list[1] - sort_site_list[0]) >200:
          tmp_list.append( sort_site_list[0])
          tmp_count_list.append( count_list[0])

        while cur_index < final_index:
          last_value = sort_site_list[last_index]
          cur_value = sort_site_list[cur_index]
          next_value = sort_site_list[next_index]
          if (cur_value - last_value) > 200 and (next_value - cur_value) > 200:
            tmp_list.append(cur_value)
            tmp_count_list.append( count_list[cur_index] )
          last_index+=1
          cur_index+=1
          next_index+=1
        if (sort_site_list[final_index] - sort_site_list[final_index-1]) >200:
          tmp_list.append( sort_site_list[final_index])
          tmp_count_list.append( count_list[final_index])
        if len(tmp_list) == 0: continue
        filter200_chr_site_dic[chr_strand] = {'site': tmp_list, 'count': tmp_count_list}

  with open(output_root + 'bed/pausing_sites.bed', 'a+') as f:
    chr_strands = list(filter200_chr_site_dic.keys() )
    # print( list(set([c[:-1] for c in chr_strands]) ) )
    sort_chrs = sort_chr_list( list(set([c[:-1] for c in chr_strands]) ) )
    # print( sort_chrs )
    if not sort_chrs:
      print( 'Warning: no pausing sites to draw lollipop plot')
      os.sys.exit(1)
    filter_chr_site_count_dic  = defaultdict(lambda: {})
    for chr in sort_chrs:
      if chr+'+' in chr_strands:
        filter_chr_site_count_dic[chr]['forward'] = filter200_chr_site_dic[chr + '+']
      else:
        filter_chr_site_count_dic[chr]['forward'] = {'site': [], 'count': []}

      if chr+'-' in chr_strands:
        filter_chr_site_count_dic[chr]['reverse'] = filter200_chr_site_dic[chr + '-']
      else:
        filter_chr_site_count_dic[chr]['reverse'] = {'site': [], 'count': []}

    global_pause_sites_data = {'chrs': sort_chrs, 'lollipop_dic': filter_chr_site_count_dic}
    global_distribution_lollipopplot(global_pause_sites_data, output_root)
    # with open('global_pause_sites_data.pickle', 'wb') as f:
    #   pickle.dump(global_pause_sites_data, f)

    n =0
    for chr_strand, site_count_dic in filter200_chr_site_dic.items():
      site_list, count_list = site_count_dic['site'], site_count_dic['count']
      if len(site_list)==0: continue
      chr, strand = chr_strand[:-1], chr_strand[-1]
      for site, count in zip(site_list, count_list):
        f.write(chr+'\t'+ str(site) + '\t' + str(site+1) + '\t' + 'ps_'+str(n) + '\t'+ str(count) + '\t' + strand + '\n')
        n+=1

if __name__ == '__main__':
  fire.Fire( main )