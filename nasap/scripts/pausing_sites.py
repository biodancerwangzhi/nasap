import fire

import multiprocessing as mp

import numpy as np
import pandas as pd
import pyBigWig

def get_bw_pause_sites(para_list):
  [bw_dir, chr, size] = para_list
  win_side= 100 # 左右100，总共201
  bw = pyBigWig.open(bw_dir)

  start, end = 0, size
  gene_depth_array = bw.values(chr, start, end, numpy=True)
  # 获取小窗的window np arrays
  suit_index = np.where(gene_depth_array>=4)[0]
  suit_index = suit_index[suit_index >win_side]
  suit_index = suit_index[suit_index < end-win_side]

  def judge_window(i):
    window=gene_depth_array[i-win_side:i+win_side+1]
    win_count_max=window.max()
    if win_count_max <4: return 0
    # 201bp小窗 中点
    if window[win_side] <4: return 0
    if window[win_side] != win_count_max: return 0

    win_count_mean=window.mean()
    # win_count_median=np.median( window )
    win_count_std=window.std()
    if win_count_max>=(win_count_mean + win_count_std*3):
      return 1
    return 0
  suit_value = np.array( list( map(judge_window, suit_index) ) )
  chr_sites = suit_index[suit_value >0]
  if 'forward' in bw_dir:
    chr_strand = chr + '+'
  else:
    chr_strand = chr + '-'
  return {chr_strand:  chr_sites}

def main(forward_bw, reverse_bw, output_root='./tmp_output/', cores=1):
  filter200_chr_site_dic = {}
  for bw_dir in [forward_bw, reverse_bw]:
    bw_value = pyBigWig.open(bw_dir)
    chr_sizes = {chr:size for chr, size in bw_value.chroms().items() }
    chr_size_list = list(chr_sizes.items())

    pool = mp.Pool(cores)
    results = []
    for chr_size_tuple in chr_size_list:
      chrom, size = chr_size_tuple[0], chr_size_tuple[1]
      process = pool.apply_async(get_bw_pause_sites, ([bw_dir, chrom, size],))
      results.append(process)
    pool.close()
    pool.join()

    for process in results:
      res_dic = process.get()
      for chr_strand, sort_site_list in res_dic.items():
        final_index = len(sort_site_list) - 1
        if final_index == -1: continue
        if final_index == 0:
            filter200_chr_site_dic[chr_strand] = sort_site_list
            continue
        if final_index == 1:
            if (sort_site_list[1] - sort_site_list[0]) > 200:
                filter200_chr_site_dic[chr_strand] = sort_site_list
            continue

        tmp_list = []
        last_index, cur_index, next_index =0, 1, 2
        if (sort_site_list[1] - sort_site_list[0]) >200: tmp_list.append( sort_site_list[0])

        while cur_index < final_index:
            last_value = sort_site_list[last_index]
            cur_value = sort_site_list[cur_index]
            next_value = sort_site_list[next_index]
            if (cur_value - last_value) > 200 and (next_value - cur_value) > 200:
                tmp_list.append(cur_value)
            last_index+=1
            cur_index+=1
            next_index+=1
        if (sort_site_list[final_index] - sort_site_list[final_index-1]) >200: tmp_list.append( sort_site_list[final_index])
        # print(chr, len(site_list), len(tmp_list) )
        if len(tmp_list) == 0: continue
        filter200_chr_site_dic[chr_strand] = tmp_list

  with open(output_root + 'bed/pausing_sites.bed', 'w') as f:
    n =0
    for chr_site, site_list in filter200_chr_site_dic.items():
      if len(site_list)==0: continue
      chr, strand = chr_site[:-1], chr_site[-1]
      for site in site_list:
        f.write(chr+'\t'+ str(site) + '\t' + str(site+1) + '\t' + 'ps_'+str(n) + '\t0\t' + strand + '\n')
        n+=1

if __name__ == '__main__':
  fire.Fire( main )
