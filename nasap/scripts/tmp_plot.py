import re, math, pickle
from collections import Counter, defaultdict

import numpy as np
import pandas as pd
from scipy import stats

import matplotlib.pyplot as plt
plt.style.use('seaborn-whitegrid')
from matplotlib.pyplot import Polygon
from matplotlib import ticker
import seaborn as sns
sns.set(style='whitegrid')

output_root = '/home/nasap_performance/GSM1480325/'
output_root = './'

'''
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
stats_list = [] # delete later
stats_list.append( ['Degradation ratio', degradation_ratio] )
max_length =max(length_count)

data = {
  'adapter_length': adapter_length,
  'length_count': length_count,
  'max_length': max_length,
  'x_max': x_max,
  'degradation_ratio': degradation_ratio,
}

with open(output_root + 'adapter_dist.pickle', 'wb') as f:
  pickle.dump(data, f)
'''


# x->insertion size,
# y->number of reads,
# text->degradation ratio
# i scatter
# i 两段阴影
# todo draw_adapter_distribution->insertion distribution
def insertion_distribution(data, output_root):
  max_x, max_y = max(data['adapter_length']), max(data['length_count'])
  fig, ax = plt.subplots(figsize=(8, 8))
  # 散点
  ax.scatter( data['adapter_length'], data['length_count'], color='black', alpha=0.6 )
  # 阴影 用polygon
  poly1=Polygon([(0, 0),(0, max_y), (20, max_y), (20, 0)],alpha=0.3, facecolor='#FFE6D0', edgecolor="black", linestyle="-." )
  poly2=Polygon([(20, 0), (20, max_y), (30, max_y), (30, 0)],alpha=0.3, facecolor='#FFFDE5', edgecolor="black", linestyle="-."  )
  ax.add_patch(poly1)
  ax.add_patch(poly2)
  # text
  plt.text(max_x *0.5, 1, 'degradation rate:' + str(round(data['degradation_ratio'], 2) ), zorder=99 )
  plt.text(10, max_y *0.1, 'high degradation', ha="left", va="bottom", rotation=90, zorder=99)
  plt.text(25, max_y *0.1, 'partial degradation', ha="left", va="bottom", rotation=90, zorder=99)
  # ax.text(0.7, 0.9, 'degradation ratio ' + str(round(degradation_ratio, 2) ) , ha='center', va='center', transform=ax.transAxes, zorder=99, backgroundcolor='white')
  plt.xlim([0, max_x*1.05])
  plt.xlabel('Size of insertion')
  plt.ylabel('Number of reads')
  # ax.set_yscale('log')
  # plt.tight_layout()
  formatter = ticker.ScalarFormatter(useMathText=True)
  formatter.set_scientific(True)
  formatter.set_powerlimits((-1,1))
  ax.yaxis.set_major_formatter(formatter)

  plt.savefig(output_root + 'imgs/adapter_insertion_distribution.png')
  plt.savefig(output_root + 'imgs/adapter_insertion_distribution.pdf' )
  plt.close(fig)

#{'adapter_length': adapter_length,
# 'length_count': length_count,
# 'degradation_ratio': degradation_ratio}
fr = open(output_root + 'adapter_dist.pickle', 'rb')
adapter_dist_data = pickle.load(fr)
insertion_distribution( adapter_dist_data, output_root  )


# sharex -> reads_length
# y -> reads_count
# axs -> steps
# annotate -> qc
def readsLength_distribution(reads_length_data, output_root ): # table_data, table_col, table_row
  fig, axes = plt.subplots(len(reads_length_data['steps']), 1, sharex=True, sharey=True, figsize=(14, 20))

  for i, step in enumerate(reads_length_data['steps']):
    x, y = reads_length_data['bar_list'][i]
    axes[i].bar(x, y)
    axes[i].set_ylabel('' )
    axes[i].grid( axis='y')
    axes[i].set_title('Preprocess: ' + step, fontsize = 14, fontweight ='bold', verticalalignment='bottom', loc ='left')
    axes[i].text(0.01, 0.9, reads_length_data['annotate_list'][i],
        horizontalalignment='left',
        verticalalignment='top',
        transform=axes[i].transAxes)
    axes[i].spines['top'].set_visible(False)
    axes[i].spines['right'].set_visible(False)
    axes[i].spines['bottom'].set_visible(False)
    # axes[i].spines['left'].set_visible(False)
  plt.xlabel('Reads length')
  plt.yscale("log")
  sns.despine(left=True)
  plt.savefig(output_root + 'imgs/reads_distribution.png')
  plt.savefig(output_root + 'imgs/reads_distribution.pdf')
  plt.close(fig)

# {'steps': steps,
# 'bar_list': bar_list,
# 'annotate_list': [annotate_text1, annotate_text2, annotate_text3, annotate_text4, annotate_text5]}
fr = open(output_root + 'readsLength.pickle', 'rb')
reads_length_data = pickle.load(fr)
readsLength_distribution(reads_length_data, output_root)

# x-> steps
# y->q20, q30
def q20_q30_line( data, output_root):
  fig, ax = plt.subplots(figsize=(10, 4))
  x = [1, 2, 3, 4, 5]
  l1, l2 = ax.plot(x, data['q20_list'], x, data['q30_list'] , lw=2)
  l1.set_label('Q20')
  l2.set_label('Q30')
  ax.legend()
  ax.set_xticks([1, 2, 3, 4, 5])
  ax.set_xticklabels(data['steps'], fontsize=18) #tex code
  plt.xticks(rotation=30)

  plt.tight_layout()
  plt.savefig(output_root + 'imgs/quality_trend.png')
  plt.savefig(output_root + 'imgs/quality_trend.pdf')
  plt.close(fig)

q20_list = [0.975902, 0.975577, 0.993754, 0.993774, 0.99425, ]
q30_list = [0.958546, 0.959483, 0.978441, 0.9785, 0.97924 ]
q20_q30_data = {
  'q20_list': q20_list,
  'q30_list': q30_list,
  'steps': ['raw reads', 'adapter removed', 'trim', 'polyX removed', 'quality filtered']
}
# {'q20_list': q20_list,
# 'q30_list': q30_list,
# 'steps': steps}
q20_q30_line( q20_q30_data, output_root)


# x->steps
# y->x_have, x_not, x_failed
def preprocess_summary_stack(data, output_root ):
  x_have, x_not, x_failed = data['stackbar_list']
  fig, ax = plt.subplots(figsize=(8,8))
  plt.bar(data['steps'], x_have, label="reads with feature",edgecolor = 'black')
  plt.bar(data['steps'], x_not, label="reads without feature",edgecolor = 'black', bottom = x_have)
  plt.bar(data['steps'], x_failed, label="failed reads",edgecolor = 'black', bottom = [i+j for i, j in zip(x_have,x_not)])
  plt.legend( bbox_to_anchor=(0.5, 1.2), loc='upper center', ncol=3, fontsize='large' )
  plt.xticks(rotation=30)
  plt.ylabel('Reads number',fontsize=12)

  for x1, y1, y2, y3 in zip(data['steps'], x_have, x_not, x_failed):
    p1 = y1/(y1+y2+y3)
    p2 = y2/(y1+y2+y3)
    p3 = y3/(y1+y2+y3)
    # print( p1, p2, p3)
    if p1 >0.05:
      plt.text(x1, y1 * 0.8, '{:.0%}'.format(p1), ha='center', va="top")
    if p2>0.05:
      plt.text(x1, y1 + (y2)* 0.4,  '{:.0%}'.format(p2), ha='center')
    if p3>0.05:
      plt.text(x1, y1 + y2 + (y3)* 0.2, '{:.0%}'.format(p3), ha='center', va="bottom")

  # plt.yscale('log')
  formatter = ticker.ScalarFormatter(useMathText=True)
  formatter.set_scientific(True)
  formatter.set_powerlimits((-1,1))
  ax.yaxis.set_major_formatter(formatter)
  plt.tight_layout()
  plt.savefig(output_root +'imgs/reads_ratio.png')
  plt.savefig(output_root +'imgs/reads_ratio.pdf')
  plt.close(fig)
# {'steps': steps
# 'stackbar_list': stackbar_list}
fr = open(output_root + 'step_stackbar.pickle', 'rb')
step_stackbar = pickle.load(fr)
preprocess_summary_stack( step_stackbar, output_root)

# x->mapping tools
# y->uni, multi, unmap
def mapping_ratio_stack(data, output_root):
  x_un, x_uni, x_mul = data['map_ratio_list']
  fig = plt.figure(figsize=(8,8)) #设置画布的尺寸
  plt.bar(data['soft'], x_un, label="unmapped reads",edgecolor = 'black')
  plt.bar(data['soft'], x_uni, label="unique mapped reads",edgecolor = 'black', bottom = x_un)
  plt.bar(data['soft'], x_mul, label="non-unique mapped reads",edgecolor = 'black', bottom = [i+j for i, j in zip(x_un,x_uni)])
  # plt.legend( loc=3,fontsize=14)  # 设置图例位置
  plt.legend( bbox_to_anchor=(0.5, 1.2), loc='upper center', ncol=3, fontsize='large' )
  plt.ylabel('Mapping reads number',fontsize=16)
  plt.xlabel('Mapping reads split method',fontsize=16)

  for x1, y1, y2, y3 in zip(data['soft'], x_un, x_uni, x_mul):
    p1 = y1/(y1+y2+y3)
    p2 = y2/(y1+y2+y3)
    p3 = y3/(y1+y2+y3)
    if p1 >0.05:
      plt.text(x1, y1 * 0.4, '{:.0%}'.format(p1), ha='center',fontsize = 15)
    if p2>0.05:
      plt.text(x1, y1 + (y2)* 0.4,  '{:.0%}'.format(p2), ha='center',fontsize = 15)
    if p3>0.05:
      plt.text(x1, y1 + y2 + (y3)* 0.4, '{:.0%}'.format(p3), ha='center',fontsize = 15)

  plt.tight_layout()
  plt.savefig(output_root +'imgs/mapping_split.png')
  plt.savefig(output_root +'imgs/mapping_split.pdf')
  plt.close(fig)

mapping_ratio_data = {'map_ratio_list': [
    [0.15, 0.28 ],
    [0.83, 0.7],
    [0.02, 0.02]
],
'soft': ['bwa-mem', 'bowtie2']
}
# {'map_ratio_list': map_ratio_list
# 'soft': soft_list}
mapping_ratio_stack(mapping_ratio_data, output_root)

# todo replace boxplot to chr_rpkm_boxplot
# x-> chrs
# y-> boxs
def chr_rpkm_boxplot( data, output_root ):
  fig, ax = plt.subplots(figsize=(4,3))
  bp = ax.boxplot(data['rpkm_list'], patch_artist=True) #, boxprops={'facecolor': fixed_color_list[0], 'color':fixed_color_list[0] })

  # for box in bp['boxes']:
  #   box.set(facecolor='#087E8B', alpha=0.6, linewidth=2)

  for whisker in bp['whiskers']:
    whisker.set(linewidth=2)

  for median in bp['medians']:
    median.set(color='black', linewidth=3)

  # ax.set_title( data['title'] )
  # ax.set_xlabel( x_label)
  ax.set_ylabel( 'log10(RPKM)' )

  ax.set_xticklabels( [x.replace('chr', '') for x in data['sort_chr_list']], fontsize=7)
  fig.savefig(output_root + 'imgs/chr_rpkm.png')
  fig.savefig(output_root + 'imgs/chr_rpkm.pdf')
  # plt.yscale('log')
  plt.tight_layout()
  plt.close(fig)

# {'rpkm_list': rpkm_list,
# 'title': title,
# 'x_label': x_label,
# 'y_label': y_label,
# 'xtick_labels': xtick_labels}

fr = open(output_root + 'chr_rpkm.pickle', 'rb')
chr_rpkm= pickle.load(fr)
chr_rpkm_boxplot( chr_rpkm, output_root )

# big_size, big_label
# small_size, small_label
def expressed_pie(data, output_root):
  fig = plt.figure(figsize=(12, 7))
  # bigger
  bigger = plt.pie(data['big_sizes'], labels=data['big_labels'], labeldistance=0.75, startangle=90, frame=True, textprops={'fontsize': 12, 'color': 'w'})
  # [autotext.set_color('white') for autotext in autotexts]
  '''
  smaller = plt.pie(sizes_vegefruit, labels=labels_vegefruit,
                    colors=colors_vegefruit, radius=0.7,
                    startangle=90, labeldistance=0.7)
  '''

  wedges, texts = plt.pie(data['small_sizes'], wedgeprops=dict(width=0.5),
    radius=0.7, labeldistance=0.7,  startangle=90)
  bbox_props = dict(boxstyle="square,pad=0.3", fc="white", ec="k", lw=0.72)
  kw = dict(arrowprops=dict(arrowstyle="-", color='black'),
    bbox=bbox_props, zorder=100, va="center")

  for i, p in enumerate(wedges):
    ang = (p.theta2 - p.theta1)/2. + p.theta1
    y = np.sin(np.deg2rad(ang))
    x = np.cos(np.deg2rad(ang))

    horizontalalignment = {-1: "right", 1: "left"}[int(np.sign(x))]
    connectionstyle = "angle,angleA=0,angleB={}".format(ang)
    kw["arrowprops"].update({"connectionstyle": connectionstyle})
    plt.annotate(data['small_labels'][i],
    xy=(x*0.6, y*0.6), xytext=(1.35*np.sign(x), 1.4*y),
      horizontalalignment=horizontalalignment, **kw)

  centre_circle = plt.Circle((0, 0), 0.4, color='white', linewidth=0)
  # fig = plt.gcf()
  fig.gca().add_artist(centre_circle)

  plt.axis('equal')
  plt.tight_layout()
  plt.savefig(output_root + 'imgs/expressed_genes.png')
  plt.close(fig)


protein_num, linc_num = 18628, 5194
big_labels = [f'Coding\n({protein_num})', f'LncRNA\n({linc_num})']
expressed_protein_num, expressed_linc_num = 12000, 1500
expressed_protein_ratio = np.round( expressed_protein_num/protein_num, 2) * 100
expressed_linc_ratio = np.round( expressed_linc_num/linc_num, 2) * 100
big_sizes = [protein_num, linc_num]

small_labels = [f'{expressed_protein_ratio}% Protein-coding genes\nExpressed',
  f'{100-expressed_protein_ratio}% Protein-coding genes\nNot Expressed',
  f'{expressed_linc_ratio}% LncRNA genes\nExpressed',
  f'{100-expressed_linc_ratio}% LncRNA genes\nNot Expressed'
]
small_sizes = [expressed_protein_num, protein_num-expressed_protein_num, expressed_linc_num, linc_num - expressed_linc_num]
# { 'big_sizes': big_sizes,
# 'big_labels': big_labels,
# 'small_sizes': small_sizes,
# 'small_labels': small_labels}
expressed_pie_data = { 'big_sizes': big_sizes,
'big_labels': big_labels,
'small_sizes': small_sizes,
'small_labels': small_labels}
expressed_pie(expressed_pie_data, output_root)

def expressed_distribution(protein_x,protein_density,linc_x, linc_density, output_root):
  fig, ax = plt.subplots(figsize=(8,8))
  plt.plot(protein_x,protein_density(protein_x), label='Coding')
  plt.plot(linc_x,linc_density(linc_x), label='LncRNA')
  ax.set_xlabel('log10(RPKM)')
  ax.set_ylabel('density')
  ax.legend(loc='best')
  plt.tight_layout()
  plt.savefig(output_root + 'imgs/expression_distribution.png')
  plt.savefig(output_root + 'imgs/expression_distribution.pdf')
  plt.close(fig)

# fr = open(output_root + 'protein_linc_distribution.pickle', 'rb')
# data3 = pickle.load(fr)
# expressed_distribution( data3['protein_x'], data3['protein_density'],
# data3['linc_x'], data3['linc_density'], output_root)
# expressed_distribution()

# x->exon_intron_df.exon
# y->exon_intron_df.intron
# text
def exon_intron_scatter( data, output_root ):
  fig, ax = plt.subplots(figsize=(8,8))
  sns.regplot( data=data['filter_exon_intron_df1'], x='exon', y='intron', marker="x",
  scatter_kws = {'alpha' : 0.2 },
  line_kws={"color": 'black'}, x_jitter = 0.5, y_jitter = 0.5,  ax=ax, truncate = True )
  plt.text(0.2, 0.75, data['text'], transform = ax.transAxes)
  # plt.title( 'Exon/Intron ratio',  )
  plt.xlabel('Exon density')
  plt.ylabel('Intron density')
  plt.xlim([0, None])
  plt.ylim([0, None])

  fig.savefig(output_root + 'imgs/exon_intron_ratio.png')
  fig.savefig(output_root + 'imgs/exon_intron_ratio.pdf')
  fig.tight_layout()
  plt.close(fig)

fr = open(output_root + 'exonIntron.pickle', 'rb')
exon_intron_data = pickle.load(fr)
exon_intron_scatter( exon_intron_data, output_root )

def proximal_pausing_sites_distribution(site_index_list, output_root):
  fig, ax = plt.subplots(figsize=(12,7))
  sns.displot(data=pd.DataFrame({'pause site': site_index_list}), x='pause site', kde=True, ax=ax)
  plt.xlabel('Pause site')
  plt.ylabel('Frequency')
  plt.text(325, int(np.max(site_index_list)*0.9), s='Mean=%s\nMedian=%s'%(int(np.mean(site_index_list)), int(np.median(site_index_list))), verticalalignment="top",horizontalalignment="right")
  # plt.show()
  plt.savefig(output_root+'imgs/proximal_pause_site.png')
  plt.close(fig)

fr = open(output_root + 'proximal_pausing_sites.pickle', 'rb')
data5 = pickle.load(fr)
proximal_pausing_sites_distribution(data5['site_index_list'], output_root)



chr_strands = set([ln.strip().split('\t')[-1] for ln in open(output_root + 'pausing_sites.bed', 'a+') ] )
def sort_chr_list( chr_list ):
  sort_chr_list = [chrom for chrom in chr_list if re.match( r'chr\d+$', chrom)]
  sort_chr_list.sort(key=lambda arr: (arr[:3], int(arr[3:])))
  if 'chrM' in chr_list:  sort_chr_list.append('chrM')
  if 'chrMT' in chr_list:  sort_chr_list.append('chrMT')
  if 'chrX' in chr_list:  sort_chr_list.append('chrX')
  if 'chrY' in chr_list:  sort_chr_list.append('chrY')
  return sort_chr_list

sort_chrs = sort_chr_list( list(set([c[:-1] for c in chr_strands]) ) )
# chr_site_count_dic = defaultdict( {'site': [], 'count': []})
filter_chr_site_count_dic = {}
for ln in open(output_root + 'pausing_sites.bed', 'a+'):
  ls = ln.split('\t')
  chrom = ls[0]
  if chrom not in filter_chr_site_count_dic.keys():
    filter_chr_site_count_dic[chrom] = {'forward': {'site': [], 'count': []}, 'reverse': {'site': [], 'count': []}}
  if ls[-1].strip() == '+':
    filter_chr_site_count_dic[chrom]['forward']['site'].append(int(ls[1]) )
  else:
    filter_chr_site_count_dic[chrom]['reverse']['site'].append(int(ls[2]) )
  filter_chr_site_count_dic[chrom]['forward']['count'].append( int(ls[4]) )

print( filter_chr_site_count_dic )
def lollipopplot(chrs, lollipop_dic, output_file):
  fig, axes = plt.subplots(len(chrs)*2, 1, figsize=(  30, 2*len(chrs) ),  sharex=True, sharey=False)
  plt.xticks([])

  i = 0
  for chr in chrs:
    try:
      axes[i].stem( lollipop_dic[chr]['forward']['site'], lollipop_dic[chr]['forward']['count'], linefmt=color, markerfmt='.', basefmt=color  )
    except:
      pass
    axes[i].set_ylabel( chr+':+') #, color=color )
    axes[i].yaxis.tick_right()
    axes[i].yaxis.get_major_ticks()[1].set_visible(False)
    axes[i].yaxis.get_major_ticks()[-2].set_visible(False)
    axes[i].spines['top'].set_visible(False)
    axes[i].spines['left'].set_visible(False)
    # axes[i].spines['right'].set_color( color )
    i+=1
    try:
      axes[i].stem( lollipop_dic[chr]['reverse']['site'], [-1*x for x in lollipop_dic[chr]['reverse']['count']], linefmt=color, markerfmt='.', basefmt=color )
    except:
      pass
    axes[i].set_ylabel( chr+':-') # , color=color)
    axes[i].yaxis.tick_right()
    axes[i].yaxis.get_major_ticks()[1].set_visible(False)
    axes[i].yaxis.get_major_ticks()[-2].set_visible(False)
    axes[i].spines['bottom'].set_visible(False)
    axes[i].spines['left'].set_visible(False)
    # axes[i].spines['right'].set_color(color)
    i+=1

  # plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)
  # plt.title('Genome wide pausing sites distribution', y=1.1, fontsize=18)
  fig.tight_layout()
  fig.savefig(output_file, bbox_inches='tight')
  # plt.show()
  plt.close(fig)
lollipopplot(sort_chrs, filter_chr_site_count_dic, output_root + 'imgs/pause_sites.png')
