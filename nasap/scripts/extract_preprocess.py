import os, sys
from collections import Counter

import matplotlib.pyplot as plt
from matplotlib.pyplot import Polygon

script_dir = os.path.abspath(os.path.join(os.path.abspath(__file__), ".."))
lib_dir = os.path.abspath(os.path.join(os.path.abspath(__file__), '../../libs') )
sys.path.append(lib_dir)

from py_ext import json2dic

output_root = sys.argv[1]
if not output_root.endswith('/'): output_root = output_root +'/'
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
tmp_variable_dir = output_root+"tmp_variable.txt"
variable_dic = {ln.split('--')[0]: ln.split('--')[1].strip() for ln in open(tmp_variable_dir)}
stats_list.append( ['Read1_num', variable_dic['read1_num']] )
stats_list.append( ['Reads_with_adapter', variable_dic['reads_with_adapter']] )
stats_list.append( ['Uninformative_adapter_reads', variable_dic['uninformative_adapter_reads']] )
stats_list.append( ['Pct_uninformative_adapter_reads', variable_dic['pct_uninformative_adapter_reads']] )
stats_list.append( ['Peak_adapter_insertion_size', variable_dic['peak_adapter_insertion_size'] ] )

remove_adapter_dic = parse_fastp_json(output_root+"json/remove_adapter.json")
cut_twoEnd_dic = parse_fastp_json(output_root+"json/trim.json")
filter_quality_dic = parse_fastp_json(output_root+"json/filter_quality.json")
remove_polyX_dic = parse_fastp_json(output_root+"json/remove_polyX.json")

# 1.1 raw reads
total_reads = remove_adapter_dic['before_reads']
s1_1 = total_reads
s1_2, s1_3 = 0, 0
annotate_text1 = 'Total reads: {reads}\nTotal bases: {bases}\nReads mean length: {mean_length}\nQ20 rate: {q20}\nQ30 rate: {q30}\nGC rate: {gc}'.format(
  reads=total_reads, bases=remove_adapter_dic['before_bases'], mean_length=remove_adapter_dic['before_read1_mean_length'],
  q20=remove_adapter_dic['before_q20_rate'], q30=remove_adapter_dic['before_q30_rate'], gc=remove_adapter_dic['before_gc_rate']
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
stats_list.append( ['Reads_with_adapter', with_adapter] )
stats_list.append( ['Uninformative_adapter_reads', fail_adapter] )
stats_list.append( ['Pct_uninformative_adapter_reads', round(fail_adapter/(pass_adapter+fail_adapter)*100, 3)] )
# stats_list.append( ['Peak_adapter_insertion_size', peak_adapter_insertion_size] )
before_adapter_bases = remove_adapter_dic['before_bases']
after_adapter_bases = remove_adapter_dic['after_bases']
adapter_loss_rate = (before_adapter_bases - after_adapter_bases)/before_adapter_bases
stats_list.append( ['Adapter_loss_rate', adapter_loss_rate] )

def draw_adapter_distribution(output_root):
  # x->insertion size,
  # y->number of reads,
  # text->degradation ratio, x-label, y-label
  # scatter
  # 两段阴影
  insertion_size_list = [int( ln.split('\t')[1].strip() ) for ln in open(output_root + 'txt/adapter_read_len1.txt')]
  # 计算 y->insertion size的频率(Counter就行)
  insert_size_counter = Counter(insertion_size_list)
  x_list = list( insert_size_counter.keys() )
  y_list = list( insert_size_counter.values() )
  x_max = max(x_list)  # x_max对应reads 没有去掉 adapter
  adapter_length, length_count = [], []
  degraded_reads, intact_reads = 0, 0
  for x, y in zip(x_list, y_list):
    if x==x_max:
      continue
    l = x_max - x
    adapter_length.append( l )
    length_count.append(y)
    if x >= 10 and x < 20:
      degraded_reads+=y
    if x >= 30 and x < 40:
      intact_reads+=y

  degradation_ratio = degraded_reads/float(intact_reads)
  stats_list.append( ['Degradation ratio', degradation_ratio] )

  fig,ax = plt.subplots()
  max_length =max(length_count)
  ax.set_ylim(ymin=0, ymax= max_length )

  # 散点
  # plt.scatter(x_list,y_list )
  ax.scatter(adapter_length,length_count, color='black', alpha=0.6 )
  # 阴影 用polygon
  poly1=Polygon([(0, 0),(0, max_length*1.1), (20, max_length*1.1), (20, 0)],alpha=0.3, facecolor='#FFE6D0', edgecolor="black", linestyle="-." )
  poly2=Polygon([(20, 0), (20, max_length*1.1), (30, max_length*1.1), (30, 0)],alpha=0.3, facecolor='#FFFDE5', edgecolor="black", linestyle="-."  )
  ax.add_patch(poly1)
  ax.add_patch(poly2)

  # text
  # plt.text(0.5, 0.5, 'degradation rate' + str(round(degradation_ratio, 2) ) )
  plt.text(10, int(max_length/3),  'high degradation', fontsize=12, rotation='vertical', zorder=99)
  plt.text(25, int(max_length/3),  'partial degradation', fontsize=12, rotation='vertical', zorder=99)
  plt.text(0.7, 0.9, 'degradation ratio ' + str(round(degradation_ratio, 2) ) , ha='center', va='center', transform=ax.transAxes, zorder=99, backgroundcolor='white')
  plt.xlim([0, max(adapter_length)])
  plt.ylim([0, max_length * 1.1])
  plt.xlabel('Size of insertion')
  plt.ylabel('Number of reads')
  plt.savefig(output_root + 'imgs/adapter_insertion_distribution.png' )
  plt.savefig(output_root + 'imgs/adapter_insertion_distribution.pdf' )

# adapter_dist_list = defaultdict(lambda: [])
# for ln in open(output_root + '/adapter_distribution.txt'):
#   ls = ln.split('\t')
#   adapter_dist_list[ls[0]].append( int(ls[1]) )

# adapter_dist_dic = {k: abs(v[0]-v[1]) for k, v in adapter_dist_list.items()  if len(v) == 2}
# # 这里的peak，我认为是 数量最多的 insertion 长度
# peak_adapter_insertion_size = max( adapter_dist_dic.values() )
# print('peak_adapter_insertion_size', peak_adapter_insertion_size )
# adapter_length_list = list(adapter_dist_dic.values())
draw_adapter_distribution( output_root  )

# 1.3 trim two ends
pass_twoEnd = cut_twoEnd_dic['passed_filter_reads']
fail_twoEnd = cut_twoEnd_dic['too_short_reads']
with_twoEnd = int(variable_dic['reads_with_cutTwoEnd'])
s3_1 = with_twoEnd
s3_2 = pass_twoEnd - with_twoEnd
s3_3 = fail_twoEnd
annotate_text3 = 'Total reads: {reads}\nTotal bases: {bases}\nReads mean length: {mean_length}\nQ20 rate: {q20}\nQ30 rate: {q30}\nGC rate: {gc}'.format(
  reads=pass_twoEnd, bases=cut_twoEnd_dic['after_bases'], mean_length=cut_twoEnd_dic['after_reads_mean_length'],
  q20=cut_twoEnd_dic['after_q20_rate'], q30=cut_twoEnd_dic['after_q30_rate'], gc=cut_twoEnd_dic['after_gc_rate']
)

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

stackbar_list = [
  # have not failed
  [s1_1, s2_1, s3_1, s4_1, s5_1],
  [s1_2, s2_2, s3_2, s4_2, s5_2],
  [s1_3, s2_3, s3_3, s4_3, s5_3]
]
def draw_readsNum(steps, stackbar_list):
  # 多堆叠图 [total-> keep(with or not), failed]
  # 同样的y ['original', 'quantity filter', 'adapter', 'trim two side', 'polyX']
  # 堆叠图 有多个x, x_have, x_not, x_failed
  x_have, x_not, x_failed = stackbar_list
  plt.figure(figsize=(10,7))#设置画布的尺寸
  plt.bar(steps, x_have, label="reads with feature",edgecolor = 'black')
  plt.bar(steps, x_not, label="reads without feature",edgecolor = 'black', bottom = x_have)
  plt.bar(steps, x_failed, label="failed reads",edgecolor = 'black', bottom = [i+j for i, j in zip(x_have,x_not)])
  plt.legend( loc=3,fontsize=14)  # 设置图例位置
  plt.ylabel('Reads number',fontsize=16)
  plt.xlabel('Preprocess steps',fontsize=16)
  # plt.title("stack plot",fontsize=18 )

  for x1, y1, y2, y3 in zip(steps, x_have, x_not, x_failed):
    p1 = y1/(y1+y2+y3)
    p2 = y2/(y1+y2+y3)
    p3 = y3/(y1+y2+y3)
    if p1 >0.05:
      plt.text(x1, y1 * 0.4, '{:.0%}'.format(p1), ha='center',fontsize = 15)
    if p2>0.05:
      plt.text(x1, y1 + (y2)* 0.4,  '{:.0%}'.format(p2), ha='center',fontsize = 15)
    if p3>0.05:
      plt.text(x1, y1 + y2 + (y3)* 0.4, '{:.0%}'.format(p3), ha='center',fontsize = 15)

  # plt.show()
  plt.savefig(output_root +'imgs/reads_distribution.png')
  plt.savefig(output_root +'imgs/reads_distribution.pdf')

steps=['raw reads', 'remove adapter', 'trim two ends', 'remove polyX', 'filter quality']
draw_readsNum(steps, stackbar_list)


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

def readsLength_dist_subplots(steps, bar_list, text_list): #
  fig, axes = plt.subplots(len(steps), 1, sharex=True, sharey=True, figsize=(14, 20))
  for i, step in enumerate(steps):
    x, y = bar_list[i]
    axes[i].bar(x, y)
    axes[i].set_ylabel('log10(Reads number)' )
    axes[i].set_title('Preprocess: ' + step, fontsize = 14, fontweight ='bold', verticalalignment='bottom', loc ='left')
    axes[i].text(0.01, 0.9, text_list[i],
        horizontalalignment='left',
        verticalalignment='top',
        transform=axes[i].transAxes)
  plt.xlabel('Reads length')
  plt.yscale("log")
  # plt.show()
  plt.savefig(output_root + 'imgs/reads_distribution_after_preprocess.png')
  plt.savefig(output_root + 'imgs/reads_distribution_after_preprocess.pdf')

readsLength_dist_subplots(steps, bar_list, [annotate_text1, annotate_text2, annotate_text3, annotate_text4, annotate_text5])


c = open(output_root + 'csv/preprocess_report.csv', 'w')
for stat_item in stats_list:
  c.write(stat_item[0] + ',' + str(stat_item[1]) + '\n')
c.close()


