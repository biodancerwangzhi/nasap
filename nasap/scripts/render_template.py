import os, re
import fire


def filter_config(output_root, config_list):
  ## 用tuple(type, file, key) 表示
  ## a 元素 是否存在，('ele', file, key)
  ## b 文件存在。 ('file', file, key)
  para_dic = {}
  for config_tuple in config_list:
    if config_tuple[0] == 'para':
      file, para = config_tuple[1], config_tuple[2]
      value = ''
      if os.path.exists(output_root + file):
        for ln in open(output_root + file):
          if ln.startswith(para):
            try:
              value = ln.split(',')[1].strip()
            except:
              print('no', file, para)
              continue
      # if value:
        # filter_config_list.append( (para, value) )
      para_dic[para] =value
    else:
      file, para = config_tuple[1], config_tuple[2]
      if os.path.exists(output_root + file):
        # filter_config_list.append( (para , file) )
        para_dic[para] =file
      else:
        para_dic[para] =""
  return para_dic

def render_template(para_dic,  output_root):
  # 把statics拷贝进output文件夹，
  template_dir = os.path.split(os.path.realpath(__file__))[0]+'/templates/'
  os.system( 'cp -r %s %s'%(template_dir + 'static/', output_root+'static/'))

  # 修改html的data，把参数字典导入到vue的data中,然后写入output的index.html
  f = open(output_root + 'index.html', 'w')
  template = open( template_dir + 'template.html').read()
  res = template.replace("data:{'replace':'here'}", 'data: %s'%str(para_dic) )
  # print( res )
  f.write(res)
  f.close()

def dic2str(dic):
  s = '{ '
  for k, v in dic.items():
    if isinstance(v, str):
      s+= k
      s+=': "'
      s+= v
      s+='",'
    if isinstance(v, list):
      s+= k + ': ["' + '","'.join(str(x) for x in v) + '"]'
  s+= ' }'
  return s

def main(output_root='./tmp_output/', type='all', is_server=''):
  ## type参数 可选 ['all', 'server',
  ## 'accessment', 'quantification',
  ## 'pausing', 'network']

  # 思路：
  ## html中引入 vue, 从而实现模板, 用变量和 v-if=变量 决定是否渲染
  ## 对所有参数值一一进行判断。例如文件是否存在
  ## 所有参数有个默认值 '', 不会显示。检查成功，这个默认值变成具体值

  # 步骤：
  # 1 生成所有参数，并对每个参数进行判断
  # 2 根据传入的type参数，决定哪些模块需要被show
  # 3 对所有值进行筛选
  # 3 渲染模板

  type_dic = {
    ## 所有参数按照模块分类(现在看分类没什么意义)
    # tuple包含3个值
    ## 1 para或file或img
    ## 2 file地址
    ## 3 变量名
    'preprocess': [
      # 1.0 file info
      ('para', 'csv/preprocess_report.csv', 'Read1_name'),
      ('para', 'csv/preprocess_report.csv', 'Read2_name'),
      ('para', 'csv/preprocess_report.csv', 'Read1_num'),
      ('para', 'csv/preprocess_report.csv', 'Read2_num'),
      ('para', 'csv/preprocess_report.csv', 'Read1_size'),
      ('para', 'csv/preprocess_report.csv', 'Read2_size'),
      # 1.1 adapter ratio
      ('para', 'csv/preprocess_report.csv', 'Reads_with_adapter'),
      ('para', 'csv/preprocess_report.csv', 'Uninformative_adapter_reads'),
      ('para', 'csv/preprocess_report.csv', 'Pct_uninformative_adapter_reads'),
      ('para', 'csv/preprocess_report.csv', 'Peak_adapter_insertion_size'),
      ('para', 'csv/preprocess_report.csv', 'Adapter_loss_rate'),
      # 1.2 RNA intergrity
      ('para', 'csv/preprocess_report.csv', 'Degradation_ratio'),
      ('img', 'imgs/adapter_insertion_distribution.png', 'img_adapter_distribution'),
      ('img', 'imgs/reads_distribution.png', 'img_reads_distribution'),
      ('img', 'imgs/reads_ratio.png', 'img_reads_ratio'),
      # 1.3 library complexity
      ('para', 'csv/mapping_report.csv', 'NRF'),
      ('para', 'csv/mapping_report.csv', 'PBC1'),
      ('para', 'csv/mapping_report.csv', 'PBC2'),
      # 1.4 QC trend
      ('img', 'imgs/quality_trend.png', 'img_quality_trend'),
      # 1.5 Nascent RNA purity
      ('para', 'csv/mapping_report.csv', 'assign_mapped'),
      ('para', 'csv/mapping_report.csv', 'chrM_mapped'),

      # 1.6 Other preprocess metric
      ('para', 'csv/preprocess_report.csv', 'Trimmed_reads'),
      ('para', 'csv/preprocess_report.csv', 'Trim_loss_rate'),
      ('para', 'csv/preprocess_report.csv', 'Reads_with_polyX'),
      ('para', 'csv/preprocess_report.csv', 'Uninformative_polyX_reads'),
      ('img', 'imgs/mapping_split.png', 'img_mapping_split'),

      # 1.7 output files
      ('file', 'fastq/clean_read1.fq.gz', 'clean_read1'),
      ('file', 'fastq/clean_read2.fq.gz', 'clean_read2'),
      ('file', 'sam/original.sam', 'original_sam'),
      ('file', 'sam/unmapped.sam', 'unmapped_sam'),
      ('file', 'sam/assign.sam', 'assign_sam'),
      ('file', 'sam/unassign.sam', 'unassign_sam'),
      ('file', 'sam/uniquemapped.sam', 'uniquemapped_sam'),
      ('file', 'sam/multimapped.sam', 'multimapped_sam'),
      ('file', 'sam/redundant.sam', 'redundant_sam'),
      ('file', 'sam/nonredundant.sam', 'nonredundant_sam'),
      ('file', 'bw/forward_bw', 'forward_bw'),
      ('file', 'bw/reverse_bw', 'reverse_bw'),

    ],
    'featureAssign': [
      # 2.1 Quantification & Normalization
      ('file', 'csv/all_feature_attrs.csv', 'all_feature_attrs_csv'),
      ('file', 'csv/lincRNA_baseCount.csv', 'lincRNA_baseCount_csv'),
      ('file', 'csv/lincRNA_gb_count.csv', 'lincRNA_gb_count_csv'),
      ('file', 'csv/lincRNA_pp_count.csv', 'lincRNA_pp_count_csv'),
      ('file', 'csv/lincRNA_rpkm.csv', 'lincRNA_rpkm_csv'),
      ('file', 'csv/lincRNA_rpm.csv', 'lincRNA_rpm_csv'),
      ('file', 'csv/protein_coding_baseCount.csv', 'protein_coding_baseCount_csv'),
      ('file', 'csv/protein_coding_gb_count.csv', 'protein_coding_gb_count_csv'),
      ('file', 'csv/protein_coding_pp_count.csv', 'protein_coding_pp_count_csv'),
      ('file', 'csv/protein_coding_rpkm.csv', 'protein_coding_rpkm_csv'),
      ('file', 'csv/protein_coding_rpkm.csv', 'protein_coding_rpm_csv'),
      ('file', 'csv/erna_quant.csv', 'erna_quant_csv'),
      # 3.1 PI $EI
      ('file', 'csv/lincRNA_ei.csv', 'lincRNA_ei_csv'),
      ('file', 'csv/lincRNA_pi.csv', 'lincRNA_pi_csv'),
      ('file', 'csv/protein_coding_ei.csv', 'protein_coding_ei_csv'),
      ('file', 'csv/protein_coding_pi.csv', 'protein_coding_pi_csv'),
      # 1.5
      ('file', 'csv/exon_intron_ratio.csv', 'exon_intron_ratio_csv'),
      ('img', 'imgs/exon_intron_ratio.png', 'img_exon_intron_ratio'),
      ('para', 'csv/preprocess_report.csv', 'mRNA_contamination'),
      # 2.2 Abundance esitimate
      ('img', 'imgs/expressed_genes.png', 'img_expressed_genes'),
      ('img', 'imgs/expression_distribution.png', 'img_expression_distribution'),
      # 2.3 Gene expression distribution
      ('img', 'imgs/chr_rpkm.png', 'img_chr_rpkm'),
      # 3.2 proximal pausing sites
      ('file', 'bed/proximal_pause_site.bed', 'proximal_pausing_site_bed'),
      ('img', 'imgs/proximal_pause_site.png', 'img_proximal_pause_site')
    ],
    'pausing_sites': [
      # 3.3 Global pausing sites
      ('file', 'bed/pausing_sites.bed', 'pausing_sites_bed'),
      ('img', 'imgs/pause_sites.png', 'img_pause_sites')
    ],
    'network': [
      # 4.1 Network measurement
      ('para', 'csv/graph_info.csv', 'nodes_num'),
      ('para', 'csv/graph_info.csv', 'edges_num'),
      ('para', 'csv/graph_info.csv', 'mean_degree'),
      ('para', 'csv/graph_info.csv', 'assortativity_coefficient'),
      ('para', 'csv/graph_info.csv', 'correlation_coefficient'),
      ('para', 'csv/graph_info.csv', 'transitivity'),
      ('para', 'csv/graph_info.csv', 'density'),
      # 4.2 Degree distribution
      ('img', 'imgs/network_degree.png', 'img_network_degree'),
      # 4.3 motif distribution
      ('img', 'imgs/network_motif.png', 'img_network_motif')
    ]
  }
  if not output_root.endswith('/'): output_root = output_root +'/'
  config_list = type_dic['preprocess'] + type_dic['featureAssign'] + type_dic['pausing_sites'] + type_dic['network']
  para_dic = filter_config(output_root, config_list)

  ## 渲染模块的字典
  type_module_dic = {
    'server': {'show_basic': '', 'show_assess': 'true', 'show_quant': 'true', 'show_pause': 'true', 'show_network': 'true'},
    'all': {'show_basic':'true', 'show_assess': 'true', 'show_quant': 'true', 'show_pause': 'true', 'show_network': 'true'},
    'assessment': {'show_basic':'', 'show_assess': 'true', 'show_quant': '', 'show_pause': '', 'show_network': ''},
    'quantification': {'show_basic':'', 'show_assess': '', 'show_quant': 'true', 'show_pause': '', 'show_network': ''},
    'pausing': {'show_basic':'', 'show_assess': '', 'show_quant': '', 'show_pause': 'true', 'show_network': ''},
    'network': {'show_basic':'', 'show_assess': '', 'show_quant': '', 'show_pause': '', 'show_network': 'true'},
  }
  show_module_dic = type_module_dic[type]
  para_dic.update( show_module_dic )

  para_dic['is_server'] = is_server



  ## 单独 为 network模块 识别一下 communities
  if para_dic['show_network']:
    community_list = []
    for community_img in os.listdir(output_root +'imgs/'):
      if community_img.startswith('community_'):
        num = int( community_img.split('_')[1].replace('.png', '') )
        community_list.append( num )
    if community_list == []:
      para_dic['community_list'] = ''
    else:
      # 4.4 community detection
      para_dic['community_list'] = ['community_%s'%str(x) for x in sorted(community_list)]

  para_str = dic2str( para_dic )
  # print(para_str)
  # 定义好参数 给模板的默认值，大多数为None
  render_template(para_str, output_root)


if __name__ == '__main__':
  fire.Fire( main )