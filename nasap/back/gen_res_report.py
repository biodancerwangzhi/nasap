import json, sys
from collections import Counter
import matplotlib.pyplot as plt
from matplotlib.pyplot import Polygon

output_root = sys.argv[1]
if not output_root.endswith('/'): output_root = output_root +'/'

try:
  isPair=sys.argv[2]=='pair'
except:
  isPair=False

stats_list = []

def json2dic(json_dir):
    with open(json_dir, 'r') as f:
        dict = json.load(fp=f)
    return dict

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

output_variable_dir = output_root+ "variable.txt"
variable_dic = {ln.split('--')[0]: ln.split('--')[1] for ln in open(output_variable_dir)}

filter_quality_dic = parse_fastp_json(output_root+"/json/filter_quality.json")
remove_adapter_dic = parse_fastp_json(output_root+"/json/remove_adapter.json")
remove_polyX_dic = parse_fastp_json(output_root+"/json/remove_polyX.json")
cut_twoEnd_dic = parse_fastp_json(output_root+"/json/cut_twoEnd.json")

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
  plt.title("stack plot",fontsize=18 )

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

steps=['raw reads', 'filter quality', 'remove adapter', 'remove polyX', 'trim low quality ends']

total_reads = filter_quality_dic['before_reads']
s1_1 = total_reads
s1_2, s1_3 = 0, 0
annotate_text1 = 'Total reads: {reads}\nTotal bases: {bases}\nReads mean length: {mean_length}\nQ20 rate: {q20}\nQ30 rate: {q30}\nGC rate: {gc}'.format(
  reads=total_reads, bases=filter_quality_dic['before_bases'], mean_length=filter_quality_dic['before_read1_mean_length'],
  q20=filter_quality_dic['before_q20_rate'], q30=filter_quality_dic['before_q30_rate'], gc=filter_quality_dic['before_gc_rate']
)
stats_list.append( ['Raw_reads', total_reads] )

pass_quality = filter_quality_dic['passed_filter_reads']
fail_quality = filter_quality_dic['low_quality_reads'] + filter_quality_dic['too_many_N_reads']
s2_1 = pass_quality
s2_2 = 0
s2_3 = fail_quality
annotate_text2 = 'Total reads: {reads}\nTotal bases: {bases}\nReads mean length: {mean_length}\nQ20 rate: {q20}\nQ30 rate: {q30}\nGC rate: {gc}'.format(
  reads=pass_quality, bases=filter_quality_dic['after_bases'], mean_length=filter_quality_dic['after_reads_mean_length'],
  q20=filter_quality_dic['after_q20_rate'], q30=filter_quality_dic['after_q30_rate'], gc=filter_quality_dic['after_gc_rate']
)
quality_loss_rate= round( (filter_quality_dic['before_bases']-filter_quality_dic['after_bases'])/ filter_quality_dic['before_bases']*100, 2)
stats_list.append( ['Trimmed_reads', pass_quality ])
stats_list.append( ['Trim_loss_rate', quality_loss_rate] )


pass_adapter = remove_adapter_dic['passed_filter_reads']
fail_adapter = remove_adapter_dic['too_short_reads']
with_adapter = remove_adapter_dic['adapter_trimmed_reads']
s3_1 = with_adapter
s3_2 = pass_adapter - with_adapter
s3_3 = fail_adapter
annotate_text3 = 'Total reads: {reads}\nTotal bases: {bases}\nReads mean length: {mean_length}\nQ20 rate: {q20}\nQ30 rate: {q30}\nGC rate: {gc}'.format(
  reads=pass_adapter, bases=remove_adapter_dic['after_bases'], mean_length=remove_adapter_dic['after_reads_mean_length'],
  q20=remove_adapter_dic['after_q20_rate'], q30=remove_adapter_dic['after_q30_rate'], gc=remove_adapter_dic['after_gc_rate']
)
stats_list.append( ['Reads_with_adapter', with_adapter] )
stats_list.append( ['Uninformative_adapter_reads', fail_adapter] )
stats_list.append( ['Pct_uninformative_adapter_reads', round(fail_adapter/(pass_adapter+fail_adapter)*100, 3)] )
# stats_list.append( ['Peak_adapter_insertion_size', peak_adapter_insertion_size] )
num_fragment_10_20, num_fragment_20_30 = 0, 0

for ln in open(output_root+'txt/filter_adapter_fq_len.txt'):
  fragment_length = int( ln.split('\t')[1] )
  if fragment_length >=10 and fragment_length < 20:
    num_fragment_10_20+=1
  if fragment_length >= 20 and fragment_length <30:
    num_fragment_20_30+=1
degradation_ratio = round(num_fragment_10_20/num_fragment_20_30, 3)
stats_list.append( ['Degradation_ratio', degradation_ratio] )


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

pass_twoEnd = cut_twoEnd_dic['passed_filter_reads']
fail_twoEnd = cut_twoEnd_dic['too_short_reads']
with_twoEnd = int(variable_dic['reads_with_cutTwoEnd'])
s5_1 = with_twoEnd
s5_2 = pass_twoEnd - with_twoEnd
s5_3 = fail_twoEnd
annotate_text5 = 'Total reads: {reads}\nTotal bases: {bases}\nReads mean length: {mean_length}\nQ20 rate: {q20}\nQ30 rate: {q30}\nGC rate: {gc}'.format(
  reads=pass_twoEnd, bases=cut_twoEnd_dic['after_bases'], mean_length=cut_twoEnd_dic['after_reads_mean_length'],
  q20=cut_twoEnd_dic['after_q20_rate'], q30=cut_twoEnd_dic['after_q30_rate'], gc=cut_twoEnd_dic['after_gc_rate']
)

stackbar_list = [
  # have not failed
  [s1_1, s2_1, s3_1, s4_1, s5_1],
  [s1_2, s2_2, s3_2, s4_2, s5_2],
  [s1_3, s2_3, s3_3, s4_3, s5_3]
]
draw_readsNum(steps, stackbar_list)


def reads_length_dist(reads_len_file, min, max):
  reads_length = [int(ln.split('\t')[1]) for ln in open(reads_len_file) ]
  length_counter = Counter(reads_length)
  x = list( range(min, max+1) )
  y = [length_counter[i] if i in length_counter.keys() else 0 for i in x ]
  return [x, y]

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


filter_quality_list = reads_length_dist(output_root + "/txt/filter_quality_fq_len.txt", 1, 100)
remove_adapter_list = reads_length_dist(output_root + "/txt/filter_adapter_fq_len.txt", 1, 100)
remove_polyX_list = reads_length_dist(output_root + "/txt/filter_polyX_fq_len.txt", 1, 100)
cut_twoEnd_list = reads_length_dist(output_root + "/txt/filter_cutTwoEnd_fq_len.txt", 1, 100)

raw_mean_length = json2dic( output_root + "/json/filter_quality.json" )['read1_before_filtering']['total_cycles']
bar_list = [[list(range(1, raw_mean_length+1)), [0]*(raw_mean_length-1) + [total_reads]], filter_quality_list, remove_adapter_list, remove_polyX_list, cut_twoEnd_list]
readsLength_dist_subplots(steps, bar_list, [annotate_text1, annotate_text2, annotate_text3, annotate_text4, annotate_text5])


def draw_adapter_distribution(insertion_size_list):
  # x->insertion size, y->number of reads, text->degradation ratio, x-label, y-label
  # scatter
  # 两段阴影

  # 计算 y->insertion size的频率(Counter就行)
  x_list, y_list, y_max = [], [], 0
  for x, y in Counter(insertion_size_list).items():
    x_list.append(x)
    y_list.append(y)
    if y > y_max: y_max=y
  # degradation ration -> 20-30/ 30-40
  degradation_ratio = len( [x for x in insertion_size_list if x >=10 and x<20] ) / \
    len( [x for x in insertion_size_list if x >=30 and x<40] )

  fig,ax = plt.subplots()
  ax.set_ylim(ymin=0, ymax=y_max)
  # 阴影 用polygon
  poly1=Polygon([(0, 0),(0, y_max), (20, y_max), (20, 0)],alpha=0.3, color=(0.89, 0.1, 0.5) )
  poly2=Polygon([(20, 0), (20, y_max), (30, y_max), (30, 0)],alpha=0.3, color=(0.2, 0.5, 0.4) )
  ax.add_patch(poly1)
  ax.add_patch(poly2)

  # 散点
  # plt.scatter(x_list,y_list )
  ax.scatter(x_list,y_list )

  # text
  # plt.text(0.5, 0.5, 'degradation rate' + str(round(degradation_ratio, 2) ) )
  plt.text(10, int(y_max/3),  'high degradation', fontsize=12, rotation='vertical' )
  plt.text(25, int(y_max/3),  'partial degradation', fontsize=12, rotation='vertical' )
  plt.text(0.7, 0.9, 'degradation rate ' + str(round(degradation_ratio, 2) ) , ha='center', va='center', transform=ax.transAxes)
  plt.savefig(output_root + 'imgs/adapter_insertion_distribution.png')
  plt.savefig(output_root + 'imgs/adapter_insertion_distribution.pdf')


insertion_size_list = [int(ln.split('\t')[1].strip() ) for ln in open(output_root + 'txt/filter_cutTwoEnd_fq_len.txt')]
draw_adapter_distribution(insertion_size_list)

stats_list.append(['mapped_reads', pass_twoEnd])
stats_list.append(['QC_filtered_reads', int(variable_dic['fail_qc_num'].strip())])
align_num = pass_twoEnd - int(variable_dic['fail_qc_num'].strip()) - int(variable_dic['unmap_num'].strip())
stats_list.append(['Unaligned_reads', int(variable_dic['unmap_num'].strip())])
stats_list.append(['Aligned_reads', align_num])
stats_list.append(['Unique_aligned_reads', int(variable_dic['unique_num'].strip())])


print(stats_list)
c = open(output_root + 'stats.csv', 'w')
for stat_item in stats_list:
  c.write(stat_item[0] + '\t' + str(stat_item[1]) + '\n')
c.close()
sys.exit(1)


def get_file_feature(file, feature):
  feature_list, length_list = [], []
  for ln in open(file):
    feature_list.append(feature)
    length_list.append( int( ln.split('\t')[1].strip('\n') ) )
  return feature_list, length_list

feature_list = ['original']*12500
length_list = [100]*12500
def update_list(file, feature, feature_list, length_list):
  a2_1, a2_2=get_file_feature(file, feature)
  feature_list.extend(a2_1)
  length_list.extend( a2_2 )
  return feature_list, length_list

# feature_list, length_list = update_list("./tmp_output/txt/filter_quality_fq_len.txt", 'filter_quality', feature_list, length_list)
# feature_list, length_list = update_list("./tmp_output/txt/filter_adapter_fq_len.txt", 'filter_adapter', feature_list, length_list)
# feature_list, length_list = update_list("./tmp_output/txt/filter_polyX_fq_len.txt", 'filter_polyX', feature_list, length_list)
# feature_list, length_list = update_list("./tmp_output/txt/filter_cutTwoEnd_fq_len.txt", 'cut TwoEnd', feature_list, length_list)

# import pandas as pd
# import matplotlib.pyplot as plt
# import seaborn as sns
# df = pd.DataFrame({'length': length_list, 'feature': feature_list})
# sns.rugplot(
#   data=df, x="length", y="feature"
# )
# plt.show()