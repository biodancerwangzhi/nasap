import os, sys, fire, psutil, time

from loguru import logger

# from .libs.py_ext import
# 思路:
# 每个命令都有参数，参数是值或文件(文件夹), 每个参数都要被检验
# 参数可以是required 或者 非required
# 如果是文件则检查文件是否存在
# 文件夹如果不存在，并创建


self_scripts = os.path.abspath( os.path.join(os.path.dirname(__file__)) ) + '/scripts/'


def _check_soft(name_list, isPython=False):
  if isPython:
    for name in name_list:
      try:
        exec('import '+ name )
      except:
        logger.error( 'import %s error.'%name)
        if name == 'community':
          name = 'python-louvain'
        logger.error( 'you need to install %s '%name)
        os.sys.exit(1)
  else:
    for name in name_list:
      soft_dep = os.system(name + ' --version')
      if soft_dep != 0:
        logger.error( 'import %s error.'%name)
        logger.error( 'you need to install %s first'%name)
        os.sys.exit(1)

def _check_para(name, value, type, isFile=False, required=False, output_root='./test_output/'):
  # 先检查参数的类型，
  # 如果是文件就检查文件是否存在
  # 如果是文件夹就判断文件夹是否存在，不存在就创建
  if not output_root.endswith('/'): output_root = output_root + '/'
  if required:
    if value == None:
      logger.error( name + ' is a required parameter')
      with open( output_root + 'error.log', 'w') as f:
        f.write( name + ' is a required parameter' )
      os.sys.exit(0)
  if not isinstance(value, type):
    logger.error('Error'+ name + 'should be a', type, 'type.')
    with open( output_root + 'error.log', 'w') as f:
      f.write( 'Error', name, 'should be a', type, 'type.')
    os.sys.exit(0)
  if isFile:
    if not os.path.exists(value):
      logger.error('Error' + value + 'not exist!')
      with open( output_root + 'error.log', 'w') as f:
        f.write( 'Error'+ value + 'not exist!')
      os.sys.exit(0)
  return True

def _get_cores(value):
  total_cores = psutil.cpu_count()
  cur_cores = 1
  if isinstance(value, int):
    if value > total_cores:
      print("warning: your cpu cores number is", total_cores, "we set the parameter cores to", total_cores)
      cur_cores = total_cores
    else:
      cur_cores = value
    return cur_cores
  if isinstance(value, str):
    if str.lower() == 'max':
      cur_cores = total_cores
      return cur_cores
    if str.lower == 'max/2':
      cur_cores = total_cores/2
      return cur_cores
  print("warning: cores parameter error use default cores=1")
  return 1

def _cur_time():
  return time.asctime( time.localtime(time.time()) )

def _build_project_path(root_dir):
  if not root_dir.endswith('/'): root_dir = root_dir + '/'
  json_dir=root_dir + 'json/'
  fq_dir=root_dir + 'fastq/'
  txt_dir=root_dir + 'txt/'
  img_dir=root_dir + 'imgs/'
  sam_dir=root_dir + 'sam/'
  bed_dir=root_dir + 'bed/'
  bw_dir=root_dir + 'bw/'
  csv_dir=root_dir + 'csv/'
  html_dir=root_dir + 'html/'
  for d in [root_dir, json_dir, fq_dir, txt_dir, img_dir, sam_dir, bed_dir, bw_dir, csv_dir, html_dir]:
    if not os.path.exists(d):
      os.makedirs(d)

def _check_output(root_dir):
  output_subdir = ['json/','fastq/','txt/', 'imgs/','sam/','bed/','bw/', 'html/']

  for sub_dir in output_subdir:
    if not root_dir.endswith('/'): root_dir = root_dir + '/'
    sub_dir = root_dir + sub_dir
    if not os.listdir(sub_dir):
      os.removedirs(sub_dir)


def server(output_root='./test_output/', forward_bw=None, reverse_bw=None, gtf=None,  cores=1, tf_source=None, tf_filter_nodes=None, enhancer_source=None, enhancer_filter_nodes=None):
  logger.info('server--check installed software')
  _check_soft( ['pandas', 'pyBigWig', 'numpy', 'scipy', 'networkx','community', 'hvplot'], isPython=True)

  logger.info('server--construct output dir')
  if not output_root.endswith('/'): output_root = output_root + '/'
  _build_project_path(output_root)

  logger.info('server--check parameter and input files')

  feature_assign_cmd_list = ['python', self_scripts + 'feature_attrs.py']
  pausing_sites_cmd_list = ['python', self_scripts + 'pausing_sites_low_memory.py']
  network_cmd_list = ['python', self_scripts + 'network_analysis.py']
  render_cmd_list = ['python', self_scripts + 'render_template.py']
  render_cmd_list.extend(['--type', 'server'])

  feature_assign_cmd_list.extend(['--output_root', output_root])
  pausing_sites_cmd_list.extend(['--output_root', output_root])
  network_cmd_list.extend(['--output_root', output_root])
  render_cmd_list.extend(['--output_root', output_root])

  # 检查必须参数
  # 检查文件
  _check_para('--forward_bw', forward_bw, str, isFile=True, required=True, output_root=output_root)
  feature_assign_cmd_list.extend(['--forward_bw', forward_bw] )
  pausing_sites_cmd_list.extend(['--forward_bw', forward_bw] )
  _check_para('--reverse_bw', reverse_bw, str, isFile=True, required=True, output_root=output_root)
  feature_assign_cmd_list.extend(['--reverse_bw', reverse_bw] )
  pausing_sites_cmd_list.extend(['--reverse_bw', reverse_bw] )
  _check_para('--gtf', gtf, str, isFile=True, required=True, output_root=output_root)
  feature_assign_cmd_list.extend(['--gtf', gtf] )

  # 可选参数
  cur_cores = _get_cores(cores)
  pausing_sites_cmd_list.extend(['--cores', str(cur_cores)])

  if not (tf_source or enhancer_source):
    logger.error( 'Error, you have to provide at least one parameter for tf_source or enhancer_source')
    os.sys.exit(1)

  if tf_source:
    _check_para('--tf_source', tf_source, str, isFile=True, required=True, output_root=output_root)
    network_cmd_list.extend( ['--tf_source', tf_source])
  if tf_filter_nodes:
    _check_para('--tf_filter_nodes', tf_filter_nodes, str, isFile=True, required=True, output_root=output_root)
    network_cmd_list.extend( ['--tf_filter_nodes', tf_filter_nodes])
  if enhancer_source:
    _check_para('--enhancer_source', enhancer_source, str, isFile=True, required=True, output_root=output_root)
    network_cmd_list.extend( ['--enhancer_source', enhancer_source])
  if enhancer_filter_nodes:
    _check_para('--enhancer_filter_nodes', enhancer_filter_nodes, str, isFile=True, required=True, output_root=output_root)
    network_cmd_list.extend( ['--enhancer_filter_nodes', enhancer_filter_nodes])

  all_steps = [feature_assign_cmd_list, pausing_sites_cmd_list, network_cmd_list, render_cmd_list]
  steps_name = ['feature_assign', 'pausing_sites', 'network_analysis', 'render_output']
  for cmd_list, step in zip(all_steps, steps_name):
    try:
      logger.info(step + ' start')
      os.system( ' '.join(cmd_list) )
      logger.info(step + ' finished')
    except Exception as e:
      with open( output_root + 'error.log', 'w') as f:
        logger.error( e )
        f.write(e)
      os.sys.exit(1)

def assessment(output_root='./test_output/', read1=None,  cores=1, read2=None, adapter1=None, adapter2=None, umi=None, bowtie_index=None, gtf=None):
  logger.info('assessment--check installed software')
  _check_soft( ['fastp', 'bioawk', 'python', 'bowtie2', 'samtools', 'bedtools', 'deeptools'], isPython=False)

  logger.info('assessment--construct output dir')
  if not output_root.endswith('/'): output_root = output_root + '/'
  _build_project_path(output_root)

  # 这里不检查 文件是否存在，因为在_preprocess函数中有检查
  _preprocess(output_root=output_root, read1=read1, read2=read2, adapter1=adapter1, adapter2=adapter2, umi=umi )
  if read2:
    _alignment(output_root=output_root, read1=output_root + 'fastq/clean_read1.fq.gz', read2=output_root + 'fastq/clean_read2.fq.gz', bowtie_index=bowtie_index, gtf=gtf, cores=cores )
  else:
    _alignment(output_root=output_root, read1=output_root + 'fastq/clean_read1.fq.gz', bowtie_index=bowtie_index, gtf=gtf, cores=cores )

  _genome_tracks(output_root=output_root, bam=output_root + 'sam/uniquemapped_sort.bam')
  _render_template(output_root=output_root, type='assessment', is_server='False')

def all(output_root='./test_output/', read1=None, bowtie_index=None,  gtf=None,  cores=1, read2=None, adapter1=None, adapter2=None, umi=None, tf_source=None, tf_filter_nodes=None, enhancer_source=None, enhancer_filter_nodes=None):
  logger.info('all--check installed software')
  _check_soft( ['fastp', 'bioawk', 'python', 'bowtie2', 'samtools', 'bedtools', 'deeptools'], isPython=False)
  _check_soft( ['pandas', 'pyBigWig', 'numpy', 'scipy', 'networkx','community', 'hvplot'], isPython=True)
  if not output_root.endswith('/'): output_root = output_root + '/'
  _build_project_path(output_root)
  assessment(output_root=output_root, read1=read1,  cores=cores, read2=read2, adapter1=adapter1, adapter2=adapter2, umi=umi, bowtie_index=bowtie_index, gtf=gtf)
  server(output_root=output_root, forward_bw=output_root + 'bw/forward.bw', reverse_bw=output_root + 'bw/reverse.bw', gtf=gtf, cores=cores, tf_source=tf_source, tf_filter_nodes=tf_filter_nodes, enhancer_source=enhancer_source, enhancer_filter_nodes=enhancer_filter_nodes )

  # 把模块组合
  # preprocess_cmd_list = ['bash', self_scripts + 'preprocess.bash']
  # extract_preprocess_cmd_list = ['python', self_scripts + 'extract_preprocess.py']
  # map1_cmd_list = ['bash', self_scripts + 'map1.bash']
  # map2_cmd_list = ['python', self_scripts + 'map_split.py', '--sam_file=original.sam']
  # map3_cmd_list = ['bash', self_scripts + 'map2.bash']
  # genome_track_cmd_list = ['bamCoverage', '--binSize', '1']
  # feature_assign_cmd_list = ['python', self_scripts + 'feature_attrs.py']
  # pausing_sites_cmd_list = ['python', self_scripts + 'pausing_sites.py']

  # # 必须参数 参数什么都不加
  # _check_para('--read1', read1, str, isFile=True, required=True, output_root=output_root)
  # preprocess_cmd_list.extend(['--read1', read1])

  # _check_para('--bowtie_index', bowtie_index, str, required=True, output_root=output_root)
  # map1_cmd_list.extend(['--bowtie_index', bowtie_index])
  # map3_cmd_list.extend(['--bowtie_index', bowtie_index])

  # _check_para('--gtf', gtf, str, isFile=True, required=True, output_root=output_root)
  # map1_cmd_list.extend(['--gtf', gtf])
  # map3_cmd_list.extend(['--gtf', gtf])
  # feature_assign_cmd_list.extend(['--gtf', gtf])

  # if output_root:
  #   if not output_root.endswith('/'): output_root = output_root +'/'
  #   preprocess_cmd_list.extend(['--output_root', output_root])
  #   extract_preprocess_cmd_list.append(output_root)
  #   map1_cmd_list.extend(['--output_root', output_root])
  #   map1_cmd_list.extend(['--read1', output_root + 'fastq/clean_read1.fq.gz'])
  #   map2_cmd_list.append('--sam_dir='+output_root+'sam/')
  #   map3_cmd_list.extend(['--output_root', output_root])
  #   map3_cmd_list.extend(['--read1', output_root + 'fastq/clean_read1.fq.gz'])
  #   genome_track_cmd_list.extend(['--bam', output_root + 'sam/uniquemapped_sort.bam'])
  #   feature_assign_cmd_list.extend(['--forward_bw', output_root + 'bw/forward.bw'] )
  #   feature_assign_cmd_list.extend(['--reverse_bw', output_root + 'bw/reverse.bw'] )
  #   feature_assign_cmd_list.extend(['--output_root', output_root])
  #   pausing_sites_cmd_list.extend(['--forward_bw', output_root + 'bw/forward.bw'] )
  #   pausing_sites_cmd_list.extend(['--reverse_bw', output_root + 'bw/reverse.bw'] )
  #   pausing_sites_cmd_list.extend(['--output_root', output_root])
  # if cores:
  #   cur_cores = str( _get_cores(cores) )
  #   preprocess_cmd_list.extend(['--cores', cur_cores])
  #   map1_cmd_list.extend(['--cores', cur_cores])
  #   map3_cmd_list.extend(['--cores', cur_cores])
  #   genome_track_cmd_list.extend(['-p', cur_cores])
  #   pausing_sites_cmd_list.extend(['--cores', cur_cores])

  # # 可选参数
  # if read2:
  #   _check_para('--read2', read2, str, isFile=True, required=True, output_root=output_root)
  #   preprocess_cmd_list.extend(['--read2', read2])
  #   clean_read2=output_root + 'fastq/clean_read2.fq.gz'
  #   map1_cmd_list.extend(['--read2', clean_read2])
  #   map3_cmd_list.extend(['--read2', clean_read2])
  # if adapter1:
  #   _check_para('--adapter1', adapter1, str, required=True, output_root=output_root)
  #   preprocess_cmd_list.extend(['--adapter1', adapter1])
  # if adapter2:
  #   _check_para('--adapter2', adapter2, str, required=True, output_root=output_root)
  #   preprocess_cmd_list.extend(['--adapter2', adapter2])

  # foward_genome_track_cmd_list, reverse_genome_track_cmd_list = genome_track_cmd_list.copy(), genome_track_cmd_list.copy()
  # foward_genome_track_cmd_list.extend(['--filterRNAstrand', 'forward'])
  # reverse_genome_track_cmd_list.extend(['--filterRNAstrand', 'reverse'])
  # foward_genome_track_cmd_list.extend(['-o', output_root + 'bw/forward.bw'])
  # reverse_genome_track_cmd_list.extend(['-o', output_root + 'bw/reverse.bw'])

  # all_steps = [preprocess_cmd_list, extract_preprocess_cmd_list, map1_cmd_list, map2_cmd_list, map3_cmd_list,
  # foward_genome_track_cmd_list, reverse_genome_track_cmd_list, pausing_sites_cmd_list]
  # _build_project_path(output_root)
  # for cmd_list in all_steps:
  #   os.system( ' '.join(cmd_list) )
    # print( ' '.join(cmd_list) )

  # _check_output(output_root)



def _preprocess(output_root=os.getcwd()+'/tmp_output', read1=None,  cores=1, read2=None, adapter1=None, adapter2=None, umi=None):
  logger.info('preprocess--check installed software')
  _check_soft( ['fastp', 'bioawk', 'python'], isPython=False)
  if read2:
    _check_soft( ['flash'], isPython=False)

  logger.info('preprocess--check parameter and input files')
  if not output_root.endswith('/'):
    output_root = output_root + '/'
  logger.info('preprocess--construct output dir')
  _build_project_path(output_root)
  log_file = open(output_root + 'tmp.log', 'w')
  cmd_list = ['bash', self_scripts + 'preprocess.bash']
  cmd_list.extend(['--output_root', output_root])
  # 必须参数 参数什么都不加
  _check_para('--read1', read1, str, isFile=True, required=True, output_root=output_root)
  cmd_list.extend(['--read1', read1])
  log_file.write( 'read1--'+os.path.basename(read1)+'\n' )
  log_file.write('read1_size--' + str(os.stat(read1).st_size / (1024 * 1024)) + 'Mb'+'\n')

  # 可选参数
  cur_cores = _get_cores(cores)
  cmd_list.extend(['--cores', str(cur_cores)])

  if read2:
    _check_para('--read2', read2, str, isFile=True, required=True, output_root=output_root)
    cmd_list.extend(['--read2', read2])
    log_file.write( 'read2--'+os.path.basename(read2)+'\n' )
    log_file.write('read2_size--' + str(os.stat(read2).st_size / (1024 * 1024)) + 'Mb'+'\n')

  if adapter1:
    _check_para('--adapter1', adapter1, str, required=True, output_root=output_root)
    cmd_list.extend(['--adapter1', adapter1])
    log_file.write( 'adapter1--'+ adapter1 +'\n' )


  if adapter2:
    _check_para('--adapter2', adapter2, str, required=True, output_root=output_root)
    cmd_list.extend(['--adapter2', adapter2])
    log_file.write( 'adapter2--'+ adapter2 +'\n' )

  # print( ' '.join(cmd_list) )
  log_file.write(_cur_time()+ 'preprocess start'+'\n')
  logger.info('preprocess--running script')
  os.system( ' '.join(cmd_list) )
  # extract preprocess
  os.system( ' '.join(['python', self_scripts + 'extract_preprocess.py', output_root]) )

  log_file.close()
  logger.success(_cur_time() + 'preprocess--Finished. Find the results in ' + output_root)

def _alignment(read1=None, bowtie_index=None, gtf=None, output_root=None, cores=1, read2=None ):
  logger.info('alignment--check installed software')
  _check_soft( ['bowtie2', 'samtools', 'bedtools'], isPython=False)

  logger.info('alignment--check parameter and input files')
  if not output_root.endswith('/'): output_root = output_root + '/'
  logger.info('alignment--construct output dir')
  _build_project_path(output_root)
  log_file = open(output_root + 'tmp.log', 'a+' )
  cmd_list = ['bash', self_scripts + 'map1.bash']
  cmd_list.extend(['--output_root', output_root])

  _check_para('--read1', read1, str, isFile=True, required=True, output_root=output_root)
  cmd_list.extend(['--read1', read1])
  log_file.write('clean_read1--'+os.path.basename(read1) + '\n' )
  log_file.write('clean_read1_size--' + str(os.stat(read1).st_size / (1024 * 1024)) + 'Mb\n')

  _check_para('--bowtie_index', bowtie_index, str, required=True, output_root=output_root)
  cmd_list.extend(['--bowtie_index', bowtie_index])

  _check_para('--gtf', gtf, str, isFile=True, required=True, output_root=output_root)
  cmd_list.extend(['--gtf', gtf])


  cur_cores = _get_cores(cores)
  cmd_list.extend(['--cores', str(cur_cores)])
  if read2:
    _check_para('--read2', read2, str, isFile=True, required=True, output_root=output_root)
    cmd_list.extend(['--read2', read2])
    log_file.write('clean_read2--'+os.path.basename(read2) + '\n' )
    log_file.write('clean_read2_size--' + str(os.stat(read2).st_size / (1024 * 1024)) + 'Mb\n')

  log_file.write(_cur_time()+ 'alignment start\n')
  os.system( ' '.join(cmd_list) )
  os.system( ' '.join(['python', self_scripts + 'map_split.py',
  '--sam_dir='+output_root+'sam/', '--sam_file=original.sam',
  '--output_root='+output_root
  ]) )
  cmd_list[1] = self_scripts + 'map2.bash'
  os.system( ' '.join(cmd_list) )
  log_file.close()
  logger.success(_cur_time()+ 'alignment--Finished. Find the results in ' + output_root)

def _genome_tracks(bam=None, output_root=None, cores=1):
  logger.info('genome_tracks--check installed software')
  _check_soft( ['deeptools'], isPython=False)

  logger.info('genome_tracks--check parameter and input files')
  if not output_root.endswith('/'): output_root = output_root + '/'
  logger.info('genome_tracks--construct output dir')
  _build_project_path(output_root)
  log_file = open(output_root + 'tmp.log', 'a+' )
  cmd_list = ['bamCoverage', '--binSize', '1']
  _check_para('--bam', bam, str, isFile=True, required=True, output_root=output_root)
  cmd_list.extend(['--bam', bam])
  log_file.write('bam--'+os.path.basename(bam) + '\n' )
  log_file.write('bam_size--' + str(os.stat(bam).st_size / (1024 * 1024)) + 'Mb\n')
  # 这里运行时间太长 cores直接写死 40，40是线程
  # if cores:
  #   cur_cores = _get_cores(cores)
  #   cmd_list.extend(['-p', str(cur_cores)])
  cmd_list.extend(['-p', str(40)])
  cmd_list_forward, cmd_list_reverse = cmd_list.copy(), cmd_list.copy()
  cmd_list_forward.extend(['--filterRNAstrand', 'forward'])
  cmd_list_reverse.extend(['--filterRNAstrand', 'reverse'])
  cmd_list_forward.extend(['-o', output_root + 'bw/forward.bw'])
  cmd_list_reverse.extend(['-o', output_root + 'bw/reverse.bw'])
  log_file.write(_cur_time()+ 'generate forward genome track start\n')
  os.system( ' '.join(cmd_list_forward) )
  log_file.write(_cur_time()+ 'generate reverse genome track start\n')
  os.system( ' '.join(cmd_list_reverse) )
  log_file.close()
  logger.success(_cur_time()+'genome_tracks--Finished. Find the results in ' + output_root)

def feature_assign(forward_bw=None, reverse_bw=None, gtf=None, output_root=None):
  logger.info('feature_assign--check installed software')
  _check_soft( ['pandas', 'pyBigWig', 'scipy'], isPython=True)

  logger.info('feature_assign--check parameter and input files')
  if not output_root.endswith('/'): output_root = output_root + '/'
  logger.info('feature_assign--construct output dir')
  _build_project_path(output_root)
  log_file = open(output_root + 'tmp.log', 'w' )

  cmd_list = ['python', self_scripts + 'feature_attrs.py']
  _check_para('--forward_bw', forward_bw, str, isFile=True, required=True, output_root=output_root)
  cmd_list.extend(['--forward_bw', forward_bw])
  log_file.write('forward_bw--'+os.path.basename(forward_bw) + '\n' )
  log_file.write('forward_bw_size--' + str(os.stat(forward_bw).st_size / (1024 * 1024)) + 'Mb\n')

  _check_para('--reverse_bw', reverse_bw, str, isFile=True, required=True, output_root=output_root)
  cmd_list.extend(['--reverse_bw', reverse_bw])
  log_file.write('reverse_bw--'+os.path.basename(reverse_bw) + '\n' )
  log_file.write('reverse_bw_size--' + str(os.stat(reverse_bw).st_size / (1024 * 1024)) + 'Mb\n')

  _check_para('--gtf', gtf, str, isFile=True, required=True, output_root=output_root)
  cmd_list.extend(['--gtf', gtf])

  cmd_list.extend(['--output_root', output_root])
  log_file.write(_cur_time()+ 'feature assign start\n')

  try:
    os.system( ' '.join(cmd_list) )
    logger.success('feature_assign--Finished. Find the results in ' + output_root)
  except:
    logger.error('feature_assign--Failed')
  log_file.write(_cur_time()+ 'feature assign finished\n')
  log_file.close()
  _render_template(output_root=output_root, type='quantification', is_server='False')

def pausing_sites(forward_bw=None, reverse_bw=None, output_root=None, cores=1):
  logger.info('pausing_sites--check installed software')
  _check_soft( ['pandas', 'pyBigWig', 'numpy'], isPython=True)

  logger.info('pausing_sites--check parameter and input files')
  if not output_root.endswith('/'): output_root = output_root + '/'
  logger.info('pausing_sites--construct output dir')
  _build_project_path(output_root)
  log_file = open(output_root + 'tmp.log', 'w' )

  cmd_list = ['python', self_scripts + 'pausing_sites_low_memory.py']

  _check_para('--forward_bw', forward_bw, str, isFile=True, required=True, output_root=output_root)
  cmd_list.extend(['--forward_bw', forward_bw])
  log_file.write('forward_bw--'+os.path.basename(forward_bw) )
  log_file.write('forward_bw_size--' + str(os.stat(forward_bw).st_size / (1024 * 1024)) + 'Mb')

  _check_para('--reverse_bw', reverse_bw, str, isFile=True, required=True, output_root=output_root)
  cmd_list.extend(['--reverse_bw', reverse_bw])
  log_file.write('reverse_bw--'+os.path.basename(reverse_bw) )
  log_file.write('reverse_bw_size--' + str(os.stat(reverse_bw).st_size / (1024 * 1024)) + 'Mb')

  cmd_list.extend(['--output_root', output_root])
  if cores:
    cur_cores = _get_cores(cores)
    cmd_list.extend(['--cores', str(cur_cores)])

  log_file.close()
  try:
    os.system( ' '.join(cmd_list) )
    logger.success('pausing_sites--Finished. Find the results in ' + output_root)
  except:
    logger.error('pausing_sites--Failed')
  _render_template(output_root=output_root, type='pausing', is_server='False')

def network_analysis(tf_source=None, tf_filter_nodes=None, enhancer_source=None, enhancer_filter_nodes=None, output_root=None ):
  logger.info('network_analysis--check installed software')
  _check_soft(['networkx','community', 'hvplot', 'numpy'], isPython=True)

  cmd_list = ['python', self_scripts + 'network_analysis.py']

  logger.info('network_analysis--check parameter and input files')
  if not output_root.endswith('/'): output_root = output_root + '/'
  logger.info('network_analysis--construct output dir')
  cmd_list.extend(['--output_root', output_root])
  _build_project_path(output_root)
  log_file = open(output_root + 'tmp.log', 'w' )

  if not(tf_source or enhancer_source):
    logger.error( 'network_analysis--Failed ' + 'no tf or enhancer source file')
    os.sys.exit(1)

  if tf_source:
    _check_para('--tf_source', tf_source, str, isFile=True, required=True, output_root=output_root)
    cmd_list.extend(['--tf_source', tf_source])
  if tf_filter_nodes:
    _check_para('--tf_filter_nodes', tf_filter_nodes, str, isFile=True, required=True, output_root=output_root)
    cmd_list.extend(['--tf_filter_nodes', tf_filter_nodes])
  if enhancer_source:
    _check_para('--enhancer_source', enhancer_source, str, isFile=True, required=True, output_root=output_root)
    cmd_list.extend(['--enhancer_source', enhancer_source])
  if enhancer_filter_nodes:
    _check_para('--enhancer_filter_nodes', enhancer_filter_nodes, str, isFile=True, required=True, output_root=output_root)
    cmd_list.extend(['--enhancer_filter_nodes', enhancer_filter_nodes])

  os.system( ' '.join(cmd_list) )
  logger.success('network analysis--Finished. Find the results in ' + output_root)
  _render_template(output_root=output_root, type='network', is_server='False')

def genome_tracks_visualize(forward_bw=None, reverse_bw=None, region=None, gene=None, output_root=None, cores=1):
  logger.info('genome_tracks_visualize--check installed software')
  _check_soft('pyGenomeTrack', isPython=False)
  logger.info('genome_tracks_visualize--check parameter and input files')
  if not output_root.endswith('/'):
    output_root = output_root + '/'
  logger.info('genome_tracks_visualize--construct output dir')
  _build_project_path(output_root)
  log_file = open(output_root + 'tmp.log', 'w')
  cmd_list = ['bash', self_scripts + 'preprocess.bash']
  cmd_list.extend(['--output_root', output_root])
  # 必须参数 参数什么都不加
  _check_para('--forward_bw', forward_bw, str, isFile=True, required=True, output_root=output_root)
  _check_para('--reverse_bw', reverse_bw, str, isFile=True, required=True, output_root=output_root)

  # check region和gene
  # 生成 pyGenomeTrack ini文件


def _render_template(type='server', output_root=None, is_server=None):
  # 每个步骤对应了输出文件， 用来丰富报表的数据
  # 考虑到 要和论文的内容对应上，所以不能按照 流程来生成报告
  # 把报告分成 1 基础信息 2 qc 3 转录级别评价 4 暂停因子和暂停位点 5 网络分析
  # 只有两种情况 会用到render_template
  logger.info('render_template--check installed software')
  cmd_list = ['python', self_scripts + 'render_template.py']

  cmd_list.extend(['--type', type])
  if not output_root.endswith('/'): output_root = output_root + '/'
  _build_project_path(output_root)
  cmd_list.extend(['--output_root', output_root])
  if is_server:
    cmd_list.extend( ['--is_server', is_server])
  os.system( ' '.join(cmd_list) )


  # 基础信息 table
  # 输入文件名
  # Forward_bw, Reverse_bw,

def main():
  fire.core.Display = lambda lines, out: print(*lines, file=out)
  fire.Fire({
    #'preprocess': preprocess,
    #'alignment': alignment,
    #'genome_tracks': genome_tracks,
    'assessment': assessment,
    'feature_assign': feature_assign,
    'pausing_sites': pausing_sites,
    'network_analysis': network_analysis,
    'render_template': _render_template,
    'server': server,
    'all': all
  })

if __name__ == '__main__':
  main()