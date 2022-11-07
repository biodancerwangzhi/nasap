import re, math, pickle
from collections import Counter, defaultdict

import numpy as np
import pandas as pd
from scipy import stats
import networkx as nx

import matplotlib.pyplot as plt
plt.style.use('seaborn-whitegrid')
from matplotlib.pyplot import Polygon
from matplotlib import ticker
import seaborn as sns
sns.set(style='whitegrid')

# fixed_color_list = ["#E64B35B2", "#4DBBD5B2", "#00A087B2", "#3C5488B2", "#F39B7FB2", "#8491B4B2",
# "#91D1C2B2", "#DC0000B2", "#7E6148B2"]
# chr_color_list = list(colors._colors_full_map.values())

def insertion_distribution(data, output_root):
  max_x, max_y = max(data['insertion_length']), max(data['length_count'])
  fig, ax = plt.subplots(figsize=(8, 8))
  # 散点
  ax.scatter( data['insertion_length'], data['length_count'], color='black', alpha=0.6 )
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

def q20_q30_line( data, output_root):
  fig, ax = plt.subplots(figsize=(8, 8))
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

def chr_rpkm_boxplot( data, output_root ):
  fig, ax = plt.subplots(figsize=(8, 8))
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

def expressed_pie(data, output_root):
  fig = plt.figure(figsize=(12, 12))
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
    xy=(x*0.6, y*0.6), xytext=(0.8*np.sign(x), 0.8*y),
      horizontalalignment=horizontalalignment, **kw)

  centre_circle = plt.Circle((0, 0), 0.4, color='white', linewidth=0)
  # fig = plt.gcf()
  fig.gca().add_artist(centre_circle)

  plt.axis('equal')
  plt.tight_layout()
  plt.savefig(output_root + 'imgs/expressed_genes.png')
  plt.savefig(output_root + 'imgs/expressed_genes.pdf')
  plt.close(fig)

def expressed_distribution(data, output_root):
  protein_data, linc_data = data['protein_data'], data['linc_data']
  protein_data1 = np.log10( protein_data[protein_data >0] )
  linc_data1 = np.log10( linc_data[linc_data >0] )
  protein_density = stats.gaussian_kde( protein_data1 )
  linc_density = stats.gaussian_kde( linc_data1 )
  # print( min(protein_data1), max(protein_data1), min(linc_data1), max(linc_data1) )
  protein_x = np.linspace(min(protein_data1), max(protein_data1) )
  linc_x = np.linspace(min(linc_data1), max(linc_data1) )

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

def proximal_pausing_sites_distribution(site_index_list, output_root):
  fig, ax = plt.subplots(figsize=(8,8))
  ax = sns.histplot(data=pd.DataFrame({'pause site': site_index_list}), x='pause site', kde=True, ax=ax,  alpha=0.75)
  plt.xlabel('Pause site')
  plt.ylabel('Frequency')
  # plt.text(325, int(np.max(site_index_list)*0.9), s='Mean=%s\nMedian=%s'%(int(np.mean(site_index_list)), int(np.median(site_index_list))), verticalalignment="top",horizontalalignment="right")
  ax.annotate('Mean=%s\nMedian=%s'%(int(np.mean(site_index_list)), int(np.median(site_index_list)) ), xy=(0.05, 0.85), xycoords='axes fraction')
  sns.despine(fig=None, ax=None, top=False, right=False, left=False, bottom=False, offset=None, trim=False)
  plt.savefig(output_root+'imgs/proximal_pause_site.png')
  plt.savefig(output_root+'imgs/proximal_pause_site.pdf')
  plt.close(fig)

def global_distribution_lollipopplot(data, output_root):
  fig, axes = plt.subplots(len(data['chrs'])*2, 1, figsize=(  30, 2*len(data['chrs']) ),  sharex=True, sharey=False)
  plt.xticks([])

  i = 0
  for chr in data['chrs']:
    # color = chr_color_list[i]
    try:
      axes[i].stem( data['lollipop_dic'][chr]['forward']['site'], data['lollipop_dic'][chr]['forward']['count'], markerfmt='.' ) #, linefmt=color,  basefmt=color  )
    except:
      pass
    axes[i].set_ylabel( chr+':+') # , color=color )
    axes[i].yaxis.tick_right()
    axes[i].yaxis.get_major_ticks()[1].set_visible(False)
    axes[i].yaxis.get_major_ticks()[-2].set_visible(False)
    axes[i].spines['top'].set_visible(False)
    axes[i].spines['left'].set_visible(False)
    # axes[i].spines['right'].set_color( color )
    i+=1
    try:
      axes[i].stem( data['lollipop_dic'][chr]['reverse']['site'], [-1*x for x in data['lollipop_dic'][chr]['reverse']['count']], markerfmt='.') #, linefmt=color,  basefmt=color )
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
  fig.savefig(output_root + 'imgs/pause_sites.png', bbox_inches='tight')
  fig.savefig(output_root + 'imgs/pause_sites.pdf', bbox_inches='tight')
  # plt.show()
  plt.close(fig)


def in_out_degree(data, output_root):
  plt.figure(figsize=(8, 8))
  plt.loglog(range(len(data['in_degree'])), data['in_degree'], marker='o', label='in-degree')
  plt.loglog(range(len(data['out_degree'])), data['out_degree'], marker='o', label='out-degree')
  plt.xlabel('Degree')
  plt.ylabel('Frequency')
  plt.legend()
  plt.savefig(output_root + 'imgs/network_degree.png')
  plt.savefig(output_root + 'imgs/network_degree.pdf')
  plt.close()

def network_motif(data, output_root):
  fig = plt.figure(figsize=(16, 9))
  axgrid = fig.add_gridspec(5, 14)

  ax0 = fig.add_subplot(axgrid[0, 0], ylabel='Triads' )
  ax0.set_ylabel('motif', fontdict={'size': 16}, rotation=0)

  # ax0.set_axis_off()
  ax0.spines['right'].set_visible(False)
  ax0.spines['top'].set_visible(False)
  ax0.spines['bottom'].set_visible(False)
  ax0.spines['left'].set_visible(False)
  ax0.set_xticks([])
  ax0.set_yticks([])

  triad_list = [
    "021D", "021U", "021C", "111D", "111U", "030T", "030C",
    "201", "120D", "120U", "120C", "210", "300"
  ]
  triad_graph_list = [nx.triad_graph( triad ) for triad in triad_list]
  for i, graph in enumerate(triad_graph_list):
      tmp_ax = fig.add_subplot(axgrid[0, i+1] )
      nx.draw_spectral(graph, ax =tmp_ax, node_size=12 )
      plt.title( i+1 )

  ax1 = fig.add_subplot(axgrid[1:5, :])
  #ax0.grid(True, which='minor')
  #ax0.axhline(y=0, color='k')
  plt.xticks(np.arange(min(data['xAxis']), max(data['xAxis'])+1, 1.0))
  ax1.bar(data['xAxis'], data['triad_count_list']  )
  ax1.set_ylabel('Quantity', fontdict={'size': 16})
  ax1.set_xlabel('Motif', fontdict={'size': 16})
  ax1.set_yscale('log')
  plt.tight_layout()
  plt.savefig(output_root + 'imgs/network_motif.png')
  plt.savefig(output_root + 'imgs/network_motif.pdf')
  plt.close()

