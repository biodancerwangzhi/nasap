import os, sys, pickle
script_dir = os.path.abspath(os.path.join(os.path.abspath(__file__), ".."))
lib_dir = os.path.abspath(os.path.join(os.path.abspath(__file__), '../../libs') )
sys.path.append(lib_dir)
output_root = sys.argv[1]
if not output_root.endswith('/'): output_root = output_root +'/'
from collections import Counter

import pandas as pd
import numpy as np

from py_ext import json2dic, get_file_size
from plot import insertion_distribution, readsLength_distribution, q20_q30_line, preprocess_summary_stack

# 目标:
# 1 生成 preprocess_report.txt
# 2 生成 3张图片

global stats_list
stats_list = []

def parse_fastp_json(json_dir):
  fastp_dic = json2dic(json_dir)
  summary = fastp_dic['summary']
  before = summary['before_filtering']
  after = summary['after_filtering']
  res = fastp_dic['filtering_result']

  para_dic = {
   'before_reads': int(before['total_reads']),
   'before_bases': int(before['total_bases']),
   'after_reads': int(after['total_reads']),
   'after_bases': int(after['total_bases']),

   'before_q20_rate': float(before['q20_rate']),
   'before_q30_rate': float(before['q30_rate']),
   'before_read1_mean_length': float(before['read1_mean_length']),
   'before_gc_rate': float(before['gc_content']),
   'after_q20_rate': float(after['q20_rate']),
   'after_q30_rate': float(after['q30_rate']),
   'after_reads_mean_length': float(after['read1_mean_length']),
   'after_gc_rate': float(after['gc_content']),

   'passed_filter_reads': int(res['passed_filter_reads']),
   'low_quality_reads': int(res['low_quality_reads']),
   'too_short_reads': int(res['too_short_reads']),
   'too_many_N_reads': int(res['too_many_N_reads'])
  }
  if 'read2_mean_length' in before.keys():
    para_dic['before_read2_mean_length'] = float(before['read2_mean_length'])
  if 'read2_mean_length' in after.keys():
    para_dic['after_read2_mean_length'] = float(after['read2_mean_length'])
  if 'adapter_cutting' in fastp_dic.keys():
    para_dic['adapter_trimmed_reads'] = fastp_dic['adapter_cutting']['adapter_trimmed_reads']
  if 'polyx_trimming' in fastp_dic.keys():
    para_dic['polyx_trimmed_reads'] = fastp_dic['polyx_trimming']['total_polyx_trimmed_reads']
  return para_dic


# 0 get variable
tmp_variable_dic = {ln.split('--')[0]: ln.split('--')[1].strip() for ln in open(output_root+"tmp_variable.txt")}
# stats_list.append( ['Read1_name', os.path.basename(tmp_variable_dic['read1_name'])] )
# stats_list.append( ['Read1_size',get_file_size(tmp_variable_dic['read1_name'])] )
stats_list.append( ['Read1_num', tmp_variable_dic['read1_num']] )

try:
  # stats_list.append( ['Read2_name', os.path.basename(tmp_variable_dic['read2_name'])] )
  # stats_list.append( ['Read2_size',get_file_size(tmp_variable_dic['read2_name'])] )
  stats_list.append( ['Read2_num', tmp_variable_dic['read2_num']] )
except:
  pass

stats_list.append( ['Reads_with_adapter', tmp_variable_dic['reads_with_adapter']] )
stats_list.append( ['Uninformative_adapter_reads', tmp_variable_dic['uninformative_adapter_reads']] )
stats_list.append( ['Pct_uninformative_adapter_reads', tmp_variable_dic['pct_uninformative_adapter_reads']] )
stats_list.append( ['Peak_adapter_insertion_size', tmp_variable_dic['peak_adapter_insertion_size'] ] )

remove_adapter_dic = parse_fastp_json(output_root+"json/remove_adapter.json")
cut_twoEnd_dic = parse_fastp_json(output_root+"json/trim.json")
filter_quality_dic = parse_fastp_json(output_root+"json/filter_quality.json")
remove_polyX_dic = parse_fastp_json(output_root+"json/remove_polyX.json")

q20_list, q30_list = [], []
# 1.1 raw reads
total_reads = remove_adapter_dic['before_reads']
s1_1 = total_reads
s1_2, s1_3 = 0, 0

table_data = []

annotate_text1 = 'Total reads: {reads}\nTotal bases: {bases}\nMean length: {mean_length}\nGC rate: {gc}\nQ20: {q20}\nQ30: {q30}'.format(
  reads=total_reads, bases=remove_adapter_dic['before_bases'], mean_length=remove_adapter_dic['before_read1_mean_length'],
  q20=remove_adapter_dic['before_q20_rate'], q30=remove_adapter_dic['before_q30_rate'], gc=remove_adapter_dic['before_gc_rate']
)
q20_list.append( remove_adapter_dic['before_q20_rate'] )
q30_list.append( remove_adapter_dic['before_q30_rate'] )
table_data.append( [total_reads, remove_adapter_dic['before_bases'], remove_adapter_dic['before_read1_mean_length'],
  remove_adapter_dic['before_gc_rate'], remove_adapter_dic['before_q20_rate'], remove_adapter_dic['before_q30_rate']]
)

# 1.2 remove adapter
pass_adapter = remove_adapter_dic['passed_filter_reads']
fail_adapter = remove_adapter_dic['too_short_reads']
with_adapter = remove_adapter_dic['adapter_trimmed_reads']
s2_1 = with_adapter
s2_2 = pass_adapter - with_adapter
s2_3 = fail_adapter
annotate_text2 = 'Total reads: {reads}\nTotal bases: {bases}\nReads mean length: {mean_length}\nQ20 rate: {q20}\nQ30 rate: {q30}\nGC rate: {gc}'.format(
  reads=pass_adapter, bases=remove_adapter_dic['after_bases'], mean_length=remove_adapter_dic['after_reads_mean_length'],
  q20=remove_adapter_dic['after_q20_rate'], q30=remove_adapter_dic['after_q30_rate'], gc=remove_adapter_dic['after_gc_rate']
)
q20_list.append(remove_adapter_dic['after_q20_rate'] )
q30_list.append(remove_adapter_dic['after_q30_rate'] )
table_data.append([pass_adapter, remove_adapter_dic['after_bases'], remove_adapter_dic['after_reads_mean_length'],
  remove_adapter_dic['after_gc_rate'], remove_adapter_dic['after_q20_rate'], remove_adapter_dic['after_q30_rate']
])
stats_list.append( ['Reads_with_adapter', with_adapter] )
stats_list.append( ['Uninformative_adapter_reads', fail_adapter] )
stats_list.append( ['Pct_uninformative_adapter_reads', round(fail_adapter/(pass_adapter+fail_adapter)*100, 3)] )
# stats_list.append( ['Peak_adapter_insertion_size', peak_adapter_insertion_size] )
before_adapter_bases = remove_adapter_dic['before_bases']
after_adapter_bases = remove_adapter_dic['after_bases']
adapter_loss_rate = (before_adapter_bases - after_adapter_bases)/before_adapter_bases
stats_list.append( ['Adapter_loss_rate', adapter_loss_rate] )

insertion_size_list = [int( ln.split('\t')[1].strip() ) for ln in open(output_root + 'txt/adapter_read_len1.txt')]
# 计算 y->insertion size的频率(Counter就行)
insert_size_counter = Counter(insertion_size_list)
x_list = list( insert_size_counter.keys() )
y_list = list( insert_size_counter.values() )
x_max = max(x_list)  # x_max对应reads 没有去掉 adapter
insertion_length, length_count = [], []
degraded_reads, intact_reads = 0, 0
for x, y in zip(x_list, y_list):
  if x==x_max:
    continue
  l = x_max - x
  insertion_length.append( l )
  length_count.append(y)
  if x >= 10 and x < 20:
    degraded_reads+=y
  if x >= 30 and x < 40:
    intact_reads+=y

degradation_ratio = degraded_reads/float(intact_reads)
stats_list.append( ['Degradation_ratio', degradation_ratio] )

insertion_dist_data =  {
  'insertion_length': insertion_length,
  'length_count': length_count,
  'degradation_ratio': degradation_ratio,
}

insertion_distribution( insertion_dist_data, output_root  )
# with open('insertion_dist_data.pickle', 'wb') as f:
#   pickle.dump(insertion_dist_data, f)

# 1.3 trim two ends
pass_twoEnd = cut_twoEnd_dic['passed_filter_reads']
fail_twoEnd = cut_twoEnd_dic['too_short_reads']
with_twoEnd = int(tmp_variable_dic['reads_with_cutTwoEnd'])
s3_1 = with_twoEnd
s3_2 = pass_twoEnd - with_twoEnd
s3_3 = fail_twoEnd
annotate_text3 = 'Total reads: {reads}\nTotal bases: {bases}\nReads mean length: {mean_length}\nQ20 rate: {q20}\nQ30 rate: {q30}\nGC rate: {gc}'.format(
  reads=pass_twoEnd, bases=cut_twoEnd_dic['after_bases'], mean_length=cut_twoEnd_dic['after_reads_mean_length'],
  q20=cut_twoEnd_dic['after_q20_rate'], q30=cut_twoEnd_dic['after_q30_rate'], gc=cut_twoEnd_dic['after_gc_rate']
)
q20_list.append( cut_twoEnd_dic['after_q20_rate'] )
q30_list.append( cut_twoEnd_dic['after_q30_rate'] )
table_data.append([ pass_twoEnd, cut_twoEnd_dic['after_bases'], cut_twoEnd_dic['after_reads_mean_length'],
  cut_twoEnd_dic['after_gc_rate'], cut_twoEnd_dic['after_q20_rate'], cut_twoEnd_dic['after_q30_rate']
])

stats_list.append( ['Trimmed_reads', with_twoEnd ])
before_trim_bases = cut_twoEnd_dic['before_bases']
after_trim_bases = cut_twoEnd_dic['after_bases']
trim_loss_rate= round( (before_trim_bases-after_trim_bases)/ before_trim_bases*100, 2)
stats_list.append( ['Trim_loss_rate', trim_loss_rate] )

# 1.4 remove polyX
pass_polyX = remove_polyX_dic['passed_filter_reads']
fail_polyX = remove_polyX_dic['too_short_reads']
with_polyX = remove_polyX_dic['polyx_trimmed_reads']
s4_1 = with_polyX
s4_2 = pass_polyX - with_polyX
s4_3 = fail_polyX
annotate_text4 = 'Total reads: {reads}\nTotal bases: {bases}\nReads mean length: {mean_length}\nQ20 rate: {q20}\nQ30 rate: {q30}\nGC rate: {gc}'.format(
  reads=pass_polyX, bases=remove_polyX_dic['after_bases'], mean_length=remove_polyX_dic['after_reads_mean_length'],
  q20=remove_polyX_dic['after_q20_rate'], q30=remove_polyX_dic['after_q30_rate'], gc=remove_polyX_dic['after_gc_rate']
)
q20_list.append( remove_polyX_dic['after_q20_rate'] )
q30_list.append( remove_polyX_dic['after_q30_rate'] )
table_data.append([pass_polyX, remove_polyX_dic['after_bases'], remove_polyX_dic['after_reads_mean_length'],
  remove_polyX_dic['after_gc_rate'], remove_polyX_dic['after_q20_rate'], remove_polyX_dic['after_q30_rate']
])
stats_list.append( ['Reads_with_polyX', with_polyX] )
stats_list.append( ['Uninformative_polyX_reads', fail_polyX] )

# 1.5 mean quality
pass_quality = filter_quality_dic['passed_filter_reads']
fail_quality = filter_quality_dic['low_quality_reads'] + filter_quality_dic['too_many_N_reads']
s5_1 = pass_quality
s5_2 = 0
s5_3 = fail_quality
annotate_text5 = 'Total reads: {reads}\nTotal bases: {bases}\nReads mean length: {mean_length}\nQ20 rate: {q20}\nQ30 rate: {q30}\nGC rate: {gc}'.format(
  reads=pass_quality, bases=filter_quality_dic['after_bases'], mean_length=filter_quality_dic['after_reads_mean_length'],
  q20=filter_quality_dic['after_q20_rate'], q30=filter_quality_dic['after_q30_rate'], gc=filter_quality_dic['after_gc_rate']
)
q20_list.append( filter_quality_dic['after_q20_rate'] )
q30_list.append( filter_quality_dic['after_q30_rate'] )
table_data.append([ pass_quality, filter_quality_dic['after_bases'], filter_quality_dic['after_reads_mean_length'],
  filter_quality_dic['after_gc_rate'], filter_quality_dic['after_q20_rate'], filter_quality_dic['after_q30_rate']
])

stackbar_list = [
  # have not failed
  [s1_1, s2_1, s3_1, s4_1, s5_1],
  [s1_2, s2_2, s3_2, s4_2, s5_2],
  [s1_3, s2_3, s3_3, s4_3, s5_3]
]
steps=['raw reads', 'adapter removed', 'trim', 'polyX removed', 'quality filtered']
step_stackbar = {'steps': steps, 'stackbar_list': stackbar_list}
preprocess_summary_stack(step_stackbar, output_root)
# with open('step_stackbar.pickle', 'wb') as f:
#   pickle.dump(step_stackbar, f)

q20_q30_data = {'q20_list': q20_list,
'q30_list': q30_list,
'steps': steps}
q20_q30_line( q20_q30_data, output_root)
# with open('q20_q30_data.pickle', 'wb') as f:
#   pickle.dump(q20_q30_data, f)

def reads_length_dist(reads_len_file, min, max):
  reads_length = [int(ln.split('\t')[1]) for ln in open(reads_len_file) ]
  length_counter = Counter(reads_length)
  x = list( range(min, max+1) )
  y = [length_counter[i] if i in length_counter.keys() else 0 for i in x ]
  return [x, y]

read_max_length = int( remove_adapter_dic['before_read1_mean_length'] )

remove_adapter_list = reads_length_dist(output_root + "/txt/adapter_read_len1.txt", 1, read_max_length)
cut_twoEnd_list = reads_length_dist(output_root + "/txt/filter_trim_len1.txt", 1, read_max_length)
remove_polyX_list = reads_length_dist(output_root + "/txt/filter_polyX_len1.txt", 1, read_max_length)
filter_quality_list = reads_length_dist(output_root + "/txt/filter_quality_len1.txt", 1, read_max_length)

bar_list = [[list(range(1, read_max_length+1)), [0]*(read_max_length-1) + [total_reads]], filter_quality_list, remove_adapter_list, remove_polyX_list, cut_twoEnd_list]
table_col = ['Total reads', 'Total bases', 'Mean length', 'GC rate', 'Q20', 'Q30']
table_row = ['raw', 'adapter removed', 'trimed', 'polyX removed', 'quality filtered']

reads_length_data = {'steps': steps, 'bar_list': bar_list,
'annotate_list': [annotate_text1, annotate_text2, annotate_text3, annotate_text4, annotate_text5]}
readsLength_distribution(reads_length_data, output_root)
# with open('reads_length_data.pickle', 'wb') as f:
#   pickle.dump(reads_length_data, f)

c = open(output_root + 'csv/preprocess_report.csv', 'w')
for stat_item in stats_list:
  c.write(stat_item[0] + ',' + str(stat_item[1]) + '\n')
c.close()
