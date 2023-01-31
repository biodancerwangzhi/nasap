import os, sys, pickle, fire
script_dir = os.path.abspath(os.path.join(os.path.abspath(__file__), ".."))
lib_dir = os.path.abspath(os.path.join(os.path.abspath(__file__), '../../libs') )
sys.path.append(lib_dir)
from collections import defaultdict, Counter

import numpy as np

import matplotlib.pyplot as plt
plt.style.use('seaborn-whitegrid')
from matplotlib.pyplot import Polygon
from matplotlib import ticker
import seaborn as sns
sns.set(style='whitegrid')

import networkx as nx
import hvplot.networkx as hvnx
from community import community_louvain

from py_ext import dic2csv, two_list_2_dict
from plot import in_out_degree

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



def get_source_target(source_file):
  source_set, target_set, source_target_edge_set = set(), set(), set()
  for ln in open(source_file):
    # if ln.startswith('hsa-'): continue
    ls = ln.split(',')
    source, target =ls[0], ls[1].strip()
    source_set.add(source)
    target_set.add(target)
    source_target_edge_set.add( ( source, target ) )

  target_set = target_set - source_set
  return  list( source_set ), list( target_set ), source_target_edge_set

def generate_graph(tf_nodes, protein_nodes, edges):
  G = nx.DiGraph()
  G.add_nodes_from( tf_nodes, shape = 'D')
  G.add_nodes_from( protein_nodes, shape='o' )
  G.add_edges_from( edges )
  return G

def filter_nodes(G, attrs_dic, over_condition=0):
  filter_nodes = []
  for node, attr in attrs_dic.items():
    try:
      float(attr)
    except:
      continue
    if float(attr) > over_condition:
      filter_nodes.append(node)
  sub_G = G.subgraph(filter_nodes)
  return sub_G

def network_summary_info(G, output_root):
  graph_info = {}
  graph_info['nodes_num'] = nx.number_of_nodes(G)
  graph_info['edges_num'] = nx.number_of_edges(G)
  graph_info['mean_degree'] = 2* graph_info['edges_num']/graph_info['nodes_num']
  graph_info['assortativity_coefficient'] = nx.degree_assortativity_coefficient(G) # 匹配系数
  graph_info['correlation_coefficient'] = nx.degree_pearson_correlation_coefficient(G) # 相关性系数
  # graph_info['closeness_centrality'] = nx.closeness_centrality(G) # 节点距离中心系数
  # graph_info['betweenness_centrality'] = nx.betweenness_centrality(G) # 节点介数中心系数
  graph_info['transitivity'] = nx.transitivity(G) # 图或网络的传递性。
  # graph_info['clustering'] = nx.clustering(G) # 图或网络中节点的聚类系数,用来度量连接的密度。
  # graph_info['average_clustering'] = nx.average_clustering(G) # nx.clustering(G) 的均值
  # graph_info['square_clustering'] = nx.square_clustering(G)
  graph_info['density'] = nx.density(G)
  dic2csv(graph_info, output_root + 'csv/graph_info.csv')


def degree_analysis(G, output_root):
  def degree_histogram_directed(G, in_degree=False, out_degree=False):
    """Return a list of the frequency of each degree value.

    Parameters
    ----------
    G : Networkx graph
        A graph
    in_degree : bool
    out_degree : bool

    Returns
    -------
    hist : list
        A list of frequencies of degrees.
        The degree values are the index in the list.

    Notes
    -----
    Note: the bins are width one, hence len(list) can be large
    (Order(number_of_edges))
    """
    nodes = G.nodes()
    if in_degree:
      in_degree = dict(G.in_degree())
      degseq=[in_degree.get(k,0) for k in nodes]
    elif out_degree:
      out_degree = dict(G.out_degree())
      degseq=[out_degree.get(k,0) for k in nodes]
    else:
      degseq=[v for k, v in G.degree()]
    dmax=max(degseq)+1
    freq= [ 0 for d in range(dmax) ]
    for d in degseq:
      freq[d] += 1
    return freq

  in_degree_freq = degree_histogram_directed(G, in_degree=True)
  out_degree_freq = degree_histogram_directed(G, out_degree=True)
  # degrees = range(len(in_degree_freq))
  in_out_degree_dic = {'in_degree': in_degree_freq, 'out_degree': out_degree_freq}
  in_out_degree( in_out_degree_dic, output_root)
  # with open('in_out_degree_dic.pickle', 'wb') as f:
  #   pickle.dump(in_out_degree_dic, f)


def computeTriads(graph):
    triads_16 = nx.triadic_census(graph)
    triads_13 = {x: y for x,y in triads_16.items() if(x != '003' and x != '012' and x != '102')}
    return triads_13

def motif_analysis( G, output_root ):
  outdeg = G.out_degree()
  # 这步要 缩减一下几点 否则 triadic 跑不动
  to_keep = [n[0] for n in outdeg if n[1] > 3]
  g = G.subgraph(to_keep)

  # triads_16 = nx.triadic_census( g )
  triads_13 = computeTriads(g)

  xAxis = [1,2,3,4,5,6,7,8,9,10,11,12,13]
  triad_list = [
    "021D", "021U", "021C", "111D", "111U", "030T", "030C",
    "201", "120D", "120U", "120C", "210", "300"
  ]
  triad_count_list =[triads_13[t] for t in triad_list]

  motif_dic = {'xAxis': xAxis, 'triad_count_list': triad_count_list}
  network_motif(motif_dic, output_root)
  with open('motif_dic.pickle', 'wb') as f:
    pickle.dump(motif_dic, f)


def community_analysis( G, node_type_dic, output_root ):
  type_list = list( set( node_type_dic.values() ) )
  shape_list = 'od^ps>v<h8'
  if len(type_list) > len( shape_list):
    print('Error, node type is over', len(shape_list) )
  type_shape_dic = two_list_2_dict( type_list, shape_list[: len(type_list)] )
  partition = community_louvain.best_partition( nx.to_undirected( G ),resolution=1)
  # dic2json(partition, output_root + 'json/partition.json' )

  community_dic = defaultdict( lambda: [])
  for n, c in partition.items():
    community_dic[c].append(n)

  dic2csv({'community_' + str(c): len(n_list) for c, n_list in community_dic.items()}, output_root + 'csv/community_node_num.csv')
  for c, n_list in community_dic.items():
    g = G.subgraph( n_list )
    keep_nodes = [n for n,d in g.degree if d >= 2 ]
    g = g.subgraph( keep_nodes )
    with open( output_root + 'txt/community_%s.txt'%c, 'w' ) as f:
      for n in g:
        f.write(n + '\n')

    if len(g) <=10:
      print( 'community %s nodes less than 10'%str(c) )
      continue
    if len( g ) > 1000:
      print( 'community %s nodes more than 1000'%str(c) )
      continue

    pos = nx.spring_layout(g, k=0.4, iterations=10)
    type_nodes_dic = {type: [] for type in type_list}
    type_pos_dic = {type: {} for type in type_list}
    for n in g:
      try:
        type = node_type_dic[n]
        type_nodes_dic[type].append(n)
        type_pos_dic[type][n] = pos[n]
        # node_communitity_dic[n] = c
      except:
        continue
    fig = plt.figure(figsize=(12,12))
    # nx.draw( g, pos, node_shape=[type_shape_dic[node_type_dic[n] ] for n in g] )

    hv_final = hvnx.draw_networkx_edges(g, pos, connectionstyle="arc3,rad=0.1", width=950, height=950, edge_line_width=0.1 )

    for type in type_nodes_dic.keys():
      # if len(g) >= 400:
        # hv_final = hv_final * hvnx.draw_networkx_nodes(g.subgraph(type_nodes_dic[type]), type_pos_dic[type], node_shape=type_shape_dic[type], node_color='#FCCC25', alpha=0.65, node_size=5)

      # if 10 < len( g ) < 400:
      hv_final = hv_final * hvnx.draw_networkx_nodes(g.subgraph(type_nodes_dic[type]), type_pos_dic[type], node_shape=type_shape_dic[type], node_size = [g.degree[node] * 10 for node in type_nodes_dic[type]], label=type  )

      nx.draw_networkx_nodes(g, pos, node_shape=type_shape_dic[type], nodelist = [x for x in g.subgraph(type_nodes_dic[type])], label='%s'%type )

    hvnx.save( hv_final, output_root + 'html/community_' +str(c)+'.html')

    nx.draw_networkx_edges(g,pos)
    plt.legend(scatterpoints = 1)
    plt.savefig(output_root + "imgs/community_%s.png"%str(c), format="PNG")
    plt.close(fig)



def main(tf_source=None, tf_filter_nodes=None, enhancer_source=None, enhancer_filter_nodes=None, output_root='./tmp_output/' ):
  filter_source, filter_target, edge_set =[], [], set()
  node_type_dic = {}

  # 1 使用source 构建网络
  if tf_source:
    tf_list, gene_list, tf_target_edge_set = get_source_target(tf_source)
    print( 'tf gene number is %s, target gene number is %s'%( str(len(tf_list) ), str(len(gene_list)) ))
    if tf_filter_nodes:
      filter_tf_nodes = [ln.strip() for ln in open(tf_filter_nodes)]
      for tf in tf_list:
        if tf in filter_tf_nodes:
          filter_source.append( tf )

      for gene in gene_list:
        if gene in filter_tf_nodes:
          filter_target.append( gene )

    else:
      filter_source.extend(tf_list)
      filter_target.extend(gene_list)
    edge_set = tf_target_edge_set
    for gene in gene_list:
      node_type_dic[gene] = 'gene'
    for tf in tf_list:
      node_type_dic[tf] = 'tf'

  if enhancer_source:
    enhancer_list, gene_list, enhancer_gene_edge_set = get_source_target(enhancer_source)
    print( 'enhancer gene number is %s, target gene number is %s'%( str(len(enhancer_list) ), str(len(gene_list)) ))
    if enhancer_filter_nodes:
      filter_enhancer_nodes = [ln.strip() for ln in open(enhancer_filter_nodes)]
      for enhancer in enhancer_list:
        if enhancer in filter_enhancer_nodes:
          filter_source.append( enhancer )
      for gene in gene_list:
        if gene in filter_enhancer_nodes:
          filter_target.append( gene )
    else:
      filter_source.extend( enhancer_list  )
      filter_target.extend( gene_list )

    for edge in enhancer_gene_edge_set:
      edge_set.add(edge)
    for gene in gene_list:
      if gene not in node_type_dic.keys():
        node_type_dic[gene] = 'gene'
    for enhancer in enhancer_list:
      node_type_dic[enhancer] = 'enhancer'

  filter_source, filter_target = list(set(filter_source)), list(set(filter_target))
  print( 'after filtering, source gene number is %s, target gene number is %s'%( str(len(filter_source) ), str(len(filter_target) )) )

  G = generate_graph(filter_target, filter_source, edge_set)
  network_summary_info(G, output_root)
  # attr_dic = {ln.split(',')[0]: ln.split(',')[1] for ln in open(filter_nodes)}
  # print( len(protein_list), len(tf_list), len(tf_protein_edge_set))
  # G = generate_graph(protein_list, tf_list, tf_protein_edge_set)
  # filter_G = filter_nodes(G, attr_dic)
  # 2 degree 分析
  degree_analysis(G, output_root)

  # 3 motif 分析
  motif_analysis(G, output_root)

  # 4 社区发现
  community_analysis( G, node_type_dic, output_root )

if __name__ == '__main__':
  fire.Fire(main)