import os, re, json, requests, fire
import pandas as pd
from loguru import logger

# 思路:
# (提取参数，检验软件，检验参数 在nasap文件中做)
# 根据调控数据下载调控数据
# 生成模板

def parse_regulatory_data(chrom, start, end,  resp_data):
  source_target_gene_list, link_list = [], []

  for regulator_dic in resp_data["regulator"]:
    regulator_chr = regulator_dic['c']
    regulator_start = regulator_dic['s']
    regulator_end = regulator_dic['e']
    regulator_name = regulator_dic['n']

    for target_obj in regulator_dic['targets']:
      target_dic = target_obj["target_site"]
      target_chr = target_dic['c']
      target_start = target_dic['s']
      target_end = target_dic['e']
      target_name = target_dic['n']

      if target_start >= int(start) and target_end <= int(end):
        ln=str(regulator_chr) + '\t' + str(regulator_start) + '\t' + str(regulator_end) + '\t' + str(target_chr) + '\t' + str(target_start) + '\t' + str(target_end) + '\n'
        link_list.append(ln)
        source_target_gene_list.append( (regulator_name, target_name) )

  return link_list, source_target_gene_list


def generate_link_file( specie, regulation_region, regulation_type ):
  chrom, start, end = re.split(r'[:-]+', regulation_region)
  try:
    # 网路获取 enhancer调控数据
    # print( 'http://grobase.top/meta/api/rest/{regulation}/{specie}/{chr}/{start}/{end}'.format(regulation=regulation_type, specie=specie, chr=chrom, start=start, end=end) )
    resp = requests.get('http://grobase.top/meta/api/rest/{regulation}/{specie}/{chr}/{start}/{end}'.format(regulation=regulation_type, specie=specie, chr=chrom, start=start, end=end))
    resp_data = json.loads( resp.text )

    link_list, source_target_gene_list = parse_regulatory_data(chrom, start, end, resp_data )
    return [link_list, source_target_gene_list]
  except:
    print( 'request failed for %s data'%regulation_type)
    return ''


def render_track_template(para_dic, show_dic, gtf, output_root):
  template = open(os.path.split(os.path.realpath(__file__))[0] + '/templates/template_track.txt').read()

  f = open(output_root + 'txt/track.txt', 'w')
  # print( para_dic )
  for k, v in para_dic.items():
    # print('$'+k, v)
    template = template.replace( '$'+k, v )

  for k, v in show_dic.items():
    if v == False:
      regex = re.compile('#show\-%s[^"]*#'%k)
      template = re.sub(regex, '', template)

  f.write(template)
  f.close()

def main( specie, region, gtf, forward, reverse, output_root, rpkm_file='' ):
  show_dic = {
    'enhancer': False,
    'tf': False
  }
  if not output_root.endswith('/'): output_root = output_root +'/'

  specie_regulation_dic = {
    'arabidopsis': ['tf'],
    'chicken': ['enhancer'],
    'chimp': ['enhancer'],
    'fly': ['enhancer', 'tf'],
    'frog': ['enhancer'],
    'mouse': ['enhancer', 'tf'],
    'human': ['enhancer', 'tf'],
    'rat': ['enhancer', 'tf'],
    'rhesus': ['enhancer'],
    'sheep': ['enhancer'],
    'worm': ['enhancer', 'tf'],
    'yeast': ['tf'],
    'zebrafish': ['enhancer', 'tf']
  }

  if specie not in specie_regulation_dic.keys():
    print('no regulatory data for specie ' + specie)
    os.sys.exit(1)

  para_dic = { 'specie': specie, 'gtf': gtf, 'forward': forward, 'reverse': reverse }
  regulation_list = specie_regulation_dic[specie]

  # 表达为0的基因 被filter
  filter_gene_list = []
  # print( 'rpkm', rpkm_file)
  if rpkm_file:
    rpkm_df = pd.read_csv( rpkm_file, index_col = 0 )
    filter_gene_list.extend( list( rpkm_df[rpkm_df.rpkm ==0].index ) )
  print('filter genes', len(filter_gene_list))
  filter_n = 0
  for regulation_type in regulation_list:
    logger.info('genome_tracks_visualize--generate %s regulatory link file'%regulation_type)
    # 获取 调控区
    link_data = generate_link_file( specie, region, regulation_type )
    # print(link_data )
    if link_data != '':
      show_dic[regulation_type] = True
      link_list, source_target_gene_list = link_data

      with open(output_root + 'txt/%s_link.txt'%regulation_type, 'w' ) as f:
        for link_ln, g in zip(link_list, source_target_gene_list):
          if (g[0] not in filter_gene_list) and (g[1] not in filter_gene_list):
            f.write(link_ln)
          else:
            filter_n +=1

      para_dic[ regulation_type ] = os.path.abspath(output_root + 'txt/%s_link.txt'%regulation_type)

  print( 'filtered num is', filter_n)

  render_track_template(para_dic, show_dic, gtf, output_root)
  # 渲染图片
  logger.info('genome_tracks_visualize--start plotting genome tracks')
  os.system('pyGenomeTracks --tracks {track_dir} --region {region} -o {img_dir}'.format(track_dir=output_root + 'txt/track.txt', region=region, img_dir=output_root + 'imgs/genome_track.png' ) )
  os.system('pyGenomeTracks --tracks {track_dir} --region {region} -o {img_dir}'.format(track_dir=output_root + 'txt/track.txt', region=region, img_dir=output_root + 'imgs/genome_track.pdf' ) )
  logger.info('genome_tracks_visualize--your results can be found at %s'%output_root)




if __name__ == "__main__":
  fire.Fire( main )
