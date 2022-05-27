import os, sys
import matplotlib.pyplot as plt
import numpy as np
import networkx as nx
import hvplot.networkx as hvnx
from community import community_louvain
import fire

script_dir = os.path.abspath(os.path.join(os.path.abspath(__file__), ".."))
lib_dir = os.path.abspath(os.path.join(os.path.abspath(__file__), '../../libs') )
sys.path.append(lib_dir)

from py_ext import dic2json, two_list_2_dict


def get_source_protein_tf(source_file):
  tf_set, protein_set, tf_protein_edge_set = set(), set(), set()
  for ln in open(source_file):
    if ln.startswith('hsa-'): continue
    ls = ln.split('\t')
    tf, gene =ls[0], ls[2]
    tf_set.add(tf)
    protein_set.add(gene)
    tf_protein_edge_set.add( ( tf, gene ) )

  protein_set = protein_set - tf_set
  return list( protein_set ), list( tf_set ), tf_protein_edge_set

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

def degree_plot(G, output_root):
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
  plt.figure(figsize=(12, 8))
  plt.loglog(range(len(in_degree_freq)), in_degree_freq, 'bo-', label='in-degree')
  plt.loglog(range(len(out_degree_freq)), out_degree_freq, color='orange', marker='o', label='out-degree')
  plt.xlabel('Degree')
  plt.ylabel('Frequency')
  plt.legend()
  plt.savefig(output_root + 'imgs/network_degree.png')
  plt.savefig(output_root + 'imgs/network_degree.pdf')
  plt.close()


def computeTriads(graph):
    triads_16 = nx.triadic_census(graph)
    triads_13 = {x: y for x,y in triads_16.items() if(x != '003' and x != '012' and x != '102')}
    return triads_13

def motif_plot( G, output_root ):
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

  triad_graph_list = [nx.triad_graph( triad ) for triad in triad_list]

  # fig = plt.figure()
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

  for i, graph in enumerate(triad_graph_list):
      tmp_ax = fig.add_subplot(axgrid[0, i+1] )
      nx.draw_spectral(graph, ax =tmp_ax, node_size=12 )
      plt.title( i+1 )

  ax1 = fig.add_subplot(axgrid[1:5, :])
  #ax0.grid(True, which='minor')
  #ax0.axhline(y=0, color='k')
  plt.xticks(np.arange(min(xAxis), max(xAxis)+1, 1.0))
  ax1.bar(xAxis, triad_count_list )
  ax1.set_ylabel('Quantity', fontdict={'size': 16})
  ax1.set_xlabel('Motif', fontdict={'size': 16})
  plt.tight_layout()
  plt.savefig(output_root + 'imgs/network_motif.png')
  plt.savefig(output_root + 'imgs/network_motif.pdf')
  plt.close()



def analysis_community( G, node_type_dic, output_root ):
  type_list = list( set( node_type_dic.values() ) )
  shape_list = 'so^>v<dph8'
  if len(type_list) > len( shape_list):
    print('Error, node type is over', len(shape_list) )
  type_shape_dic = two_list_2_dict( type_list, shape_list[: len(type_list)] )

  partition = community_louvain.best_partition( nx.to_undirected( G ),resolution=1)
  communities = set( partition.values() )
  # print( Counter( list(partition.values()) ) )
  node_communitity_dic = {}
  for c in communities:
      part_value = []
      for k, v in partition.items():
          if v == c:
              part_value.append(k)

      g = G.subgraph( part_value)
      # keep_node = [n for n, d in g.degree if d >5]
      # filter_g = g.subgraph( keep_node )

      # print( len(g.nodes ))
      # 移去 度为 <2的节点， 度为2表示 入度
      keep_nodes = [n for n,d in g.degree if d >= 2 ]
      g = g.subgraph( keep_nodes )
      # print( len(g.nodes ))

      pos = nx.spring_layout(g, k=0.4, iterations=70)
      type_nodes_dic = {type: [] for type in type_list}
      type_pos_dic = {type: {} for type in type_list}
      for n in g:
        try:
          type = node_type_dic[n]
          type_nodes_dic[type].append(n)
          type_pos_dic[type][n] = pos[n]
          node_communitity_dic[n] = c
        except:
          continue

      if len(g) <=10: continue

      hv_final = hvnx.draw_networkx_edges(g, pos, connectionstyle="arc3,rad=0.1", width=950, height=950, edge_line_width=0.1 )

      for type in type_nodes_dic.keys():
        if len(g) >= 400:
          hv_final = hv_final * hvnx.draw_networkx_nodes(g.subgraph(type_nodes_dic[type]), type_pos_dic[type], node_shape=type_shape_dic[type], node_color='#FCCC25', alpha=0.65, node_size=5)

        if 10 < len( g ) < 400:
          hv_final = hv_final * hvnx.draw_networkx_nodes(g.subgraph(type_nodes_dic[type]), type_pos_dic[type], node_shape=type_shape_dic[type], node_size = [g.degree[node] * 10 for node in type_nodes_dic[type]], label=type  )

      hvnx.save( hv_final, output_root + 'html/community_' +str(c)+'.html')

  dic2json(node_communitity_dic, output_root + 'json/network_community.json' )


def main(regulatory_source, attribute_file, output_root='./tmp_output/' ):
  # 1 用source 构建网络
  protein_list, tf_list, tf_protein_edge_set = get_source_protein_tf(regulatory_source)
  attr_dic = {ln.split(',')[0]: ln.split(',')[1] for ln in open(attribute_file)}
  # print( len(protein_list), len(tf_list), len(tf_protein_edge_set))
  G = generate_graph(protein_list, tf_list, tf_protein_edge_set)
  filter_G = filter_nodes(G, attr_dic)
  # 2 degree 分析
  degree_plot(filter_G, output_root)

  # 3 motif 分析
  motif_plot(filter_G, output_root)

  # 4 社区发现
  node_type_dic = {}
  for p in protein_list:
    node_type_dic[p] = 'protein'
  for tf in tf_list:
    node_type_dic[tf] = 'tf'
  analysis_community( filter_G, node_type_dic, output_root )

if __name__ == '__main__':
  fire.Fire(main)