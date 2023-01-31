import os,  subprocess, fire, psutil, time

from loguru import logger

# from .libs.py_ext import
# 思路:
# 每个命令都有参数，参数是值或文件(文件夹), 每个参数都要被检验
# 参数可以是required 或者 非required
# 如果是文件则检查文件是否存在
# 文件夹如果不存在，并创建

self_src = os.path.abspath( os.path.join(os.path.dirname(__file__)) ) + '/scripts/'
global_root = './'

def _ver():
  with open(os.path.abspath(os.path.join(os.path.dirname(__file__) ) ) + '/__version__.py') as f:
    ver = f.read().split('version=')[1].replace('"', '').strip()

  print("nASAP version is " + ver)

def _write_err( content, output_root=None):
  logger.error( content )
  if not output_root: output_root = global_root
  with open( output_root + 'error.log', 'w') as f:
    f.write( content )
  os.sys.exit(0)

def _write_info(content, output_root):
  logger.info(content)
  with open( output_root + 'info.log', 'a') as f:
    f.write("\n")
    f.write( content )

def _check_soft(name_list, isPython=False):
  if isPython:
    for name in name_list:
      try:
        exec('import '+ name )
      except:
        if name == 'community':
          name = 'python-louvain'
        _write_err( 'import %s error.\nyou need to install %s first'%(name, name) )
  else:
    for name in name_list:
      soft_dep = os.system(name + ' --version')
      if soft_dep != 0:
        _write_err( 'import %s error.\nyou need to install %s first'%(name, name) )

def _check_para(name, value, type, isFile=False, required=False, output_root='./test_output/'):
  # 先检查参数的类型，
  # 如果是文件就检查文件是否存在
  # 如果是文件夹就判断文件夹是否存在，不存在就创建
  if not output_root.endswith('/'): output_root = output_root + '/'
  if required:
    if value == None:
      _write_err(name + ' is a required parameter')
  if not isinstance(value, type):
    _write_err('Error: '+ name + ' should be a ', type, ' type.')
  if isFile:
    if not os.path.exists(value):
      _write_err('Error: ' + value + ' not exist!')
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

def _create_project_dir(output_root):
  global global_root
  global_root = output_root
  if not output_root.endswith('/'): output_root = output_root + '/'
  json_dir=output_root + 'json/'
  fq_dir=output_root + 'fastq/'
  txt_dir=output_root + 'txt/'
  img_dir=output_root + 'imgs/'
  sam_dir=output_root + 'sam/'
  bed_dir=output_root + 'bed/'
  bw_dir=output_root + 'bw/'
  csv_dir=output_root + 'csv/'
  html_dir=output_root + 'html/'
  for d in [output_root, json_dir, fq_dir, txt_dir, img_dir, sam_dir, bed_dir, bw_dir, csv_dir, html_dir]:
    if not os.path.exists(d):
      os.makedirs(d)

def _check_output(output_root):
  output_subdir = ['json/','fastq/','txt/', 'imgs/','sam/','bed/','bw/', 'html/']

  for sub_dir in output_subdir:
    if not output_root.endswith('/'): output_root = output_root + '/'
    sub_dir = output_root + sub_dir
    if not os.listdir(sub_dir):
      os.removedirs(sub_dir)


def server(output_root='./test_output/', forward_bw=None, reverse_bw=None, gtf=None,  cores=1, tf_source=None, tf_filter_nodes=None, enhancer_source=None, enhancer_filter_nodes=None):
  global global_root
  global_root = output_root
  logger.info('server--check installed software')
  _check_soft( ['pandas', 'pyBigWig', 'numpy', 'scipy', 'networkx','community', 'hvplot'], isPython=True)

  logger.info('server--construct output dir')
  if not output_root.endswith('/'): output_root = output_root + '/'
  _create_project_dir(output_root)

  logger.info('server--check parameter and input files')

  feature_assign_cmd_list = ['python', self_src + 'feature_attrs.py']
  pausing_sites_cmd_list = ['python', self_src + 'pausing_sites_low_memory.py']
  network_cmd_list = ['python', self_src + 'network_analysis.py']
  render_cmd_list = ['python', self_src + 'render_template.py']
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
      _write_err(step + ' --Failed\n' + str(e) )

def assessment(output_root='./test_output/', read1=None,  cores=1, read2=None, adapter1=None, adapter2=None, umi_loc=None, umi_len=None, bowtie_index=None, gtf=None, genome=None, scale_factor=None):
  global global_root
  global_root = output_root
  logger.info('assessment--check installed software')
  _check_soft( ['fastp', 'bioawk', 'python', 'bowtie2', 'samtools', 'bedtools', 'deeptools'], isPython=False)

  logger.info('assessment--construct output dir')
  if not output_root.endswith('/'): output_root = output_root + '/'
  _create_project_dir(output_root)

  # 这里不检查 文件是否存在，因为在_pp函数中有检查
  _pp(output_root=output_root, read1=read1, read2=read2, adapter1=adapter1, adapter2=adapter2, umi_loc=umi_loc, umi_len=umi_len, genome=genome)
  if read2:
    _align(output_root=output_root, read1=output_root + 'fastq/clean_read1.fq.gz', read2=output_root + 'fastq/clean_read2.fq.gz', bowtie_index=bowtie_index, gtf=gtf, cores=cores )
  else:
    _align(output_root=output_root, read1=output_root + 'fastq/clean_read1.fq.gz', bowtie_index=bowtie_index, gtf=gtf, cores=cores )

  # remove umi with gencore
  if umi_loc:
    # 不用check 值, 因为上面check 过了
    try:
      subprocess.run(self_src + 'gencore/gencore', shell=True, check=True)
      gencore_dir = self_src + 'gencore/gencore'
    except:
      _check_soft( ['gencore'], isPython=False)
      gencore_dir = 'gencore'
    try:
      os.system(' '.join([gencore_dir, '-i', output_root + 'sam/uniquemapped_sort.bam', '-o',
      output_root + 'sam/uniquemapped_sort_umi.bam', '-r', genome]) )
      # 排序 index
      os.system('samtools sort -@ %s -o %s %s'%(cores, output_root + 'sam/uniquemapped_sort.bam', output_root + 'sam/uniquemapped_sort_umi.bam') )
      os.system('samtools index -@ %s %s'%(cores, output_root + 'sam/uniquemapped_sort.bam') )
    except Exception as e:
      _write_err('gencore remove umi--Failed\n' + str(e) )

  if scale_factor:
    try:
      scale_factor = float(scale_factor)
    except Exception as e:
      _write_err('--scale_factor parsed--Failed\n' + str(e) )

  _tracks(output_root=output_root, bam=output_root + 'sam/uniquemapped_sort.bam', scale_factor=scale_factor)
  _render_template(output_root=output_root, type='assessment', is_server='No')

def all(output_root='./test_output/', read1=None, bowtie_index=None,  gtf=None,  cores=1, read2=None, adapter1=None, adapter2=None, umi_loc=None, umi_len=None, tf_source=None, tf_filter_nodes=None, enhancer_source=None, enhancer_filter_nodes=None):
  global global_root
  global_root = output_root
  logger.info('all--check installed software')
  _check_soft( ['fastp', 'bioawk', 'python', 'bowtie2', 'samtools', 'bedtools', 'deeptools'], isPython=False)
  _check_soft( ['pandas', 'pyBigWig', 'numpy', 'scipy', 'networkx','community', 'hvplot'], isPython=True)
  if not output_root.endswith('/'): output_root = output_root + '/'
  _create_project_dir(output_root)
  assessment(output_root=output_root, read1=read1,  cores=cores, read2=read2, adapter1=adapter1, adapter2=adapter2, umi_loc=umi_loc, umi_len=umi_len, bowtie_index=bowtie_index, gtf=gtf)
  server(output_root=output_root, forward_bw=output_root + 'bw/forward.bw', reverse_bw=output_root + 'bw/reverse.bw', gtf=gtf, cores=cores, tf_source=tf_source, tf_filter_nodes=tf_filter_nodes, enhancer_source=enhancer_source, enhancer_filter_nodes=enhancer_filter_nodes )
  _render_template(output_root=output_root, type='all', is_server='No')


def _pp(output_root=os.getcwd()+'/tmp_output', read1=None,  cores=1, read2=None, adapter1=None, adapter2=None, umi_loc=None, umi_len=None, genome=None, scaleFactor=None):
  logger.info('preprocess--check installed software')
  _check_soft( ['fastp', 'bioawk', 'python'], isPython=False)
  if read2:
    _check_soft( ['flash'], isPython=False)

  logger.info('preprocess--check parameter and input files')
  if not output_root.endswith('/'):
    output_root = output_root + '/'
  logger.info('preprocess--construct output dir')
  _create_project_dir(output_root)
  log_file = open(output_root + 'tmp.log', 'w')
  cmd_list = ['bash', self_src + 'preprocess.bash']
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

  if umi_loc:
    if genome:
      _check_para('--genome', genome, str, isFile=True, required=True, output_root=output_root)
    else:
      logger.error('umi module must provide genome file with --genome')
      log_file.write( 'umi module must provide genome file with --genome' )
      os.sys.exit(1)

    # 不用check_para, 因为会检测 值必须是在 列表内的值
    if umi_loc not in ['read1', 'read2', 'per_read', 'index1', 'index2', 'per_index']:
      logger.error('umi_loc can be one of index1/index2/read1/read2/per_index/per_read')
      log_file.write( 'umi_loc can be index1/index2/read1/read2/per_index/per_read'+os.path.basename(read2)+'\n' )
      os.sys.exit(1)
    else:
      if umi_loc in ['read1', 'read2', 'per_read']:
        if isinstance(umi_len, int):
          # read 模式
          cmd_list.extend(['--umi_loc', umi_loc])
          cmd_list.extend(['--umi_len', str(umi_len)])
        else:
          logger.error('umi on ' + umi_loc + ' must provide --umi_len')
          os.sys.exit(1)
      else:
        # index模式
        cmd_list.extend(['--umi_loc', umi_loc])

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
  try:
    os.system( ' '.join(cmd_list) )
  except Exception as e:
    _write_err('preprocess--Failed\n' + str(e) )

  # extract preprocess
  try:
    os.system( ' '.join(['python', self_src + 'extract_preprocess.py', output_root]) )
  except Exception as e:
    _write_err('extract preprocess--Failed\n' + str(e) )

  log_file.close()
  logger.success(_cur_time() + 'preprocess--Finished. Find the results in ' + output_root)

def preprocess_fast(output_root=os.getcwd()+'/tmp_output', read1=None,  cores=1, read2=None, adapter1=None, adapter2=None, umi=None, bowtie_index=None):
  global global_root
  global_root = output_root
  logger.info('preprocess--check installed software')
  _check_soft( ['fastp', 'bioawk', 'python'], isPython=False)
  if read2:
    _check_soft( ['flash'], isPython=False)

  logger.info('preprocess--check parameter and input files')
  if not output_root.endswith('/'):
    output_root = output_root + '/'
  logger.info('preprocess--construct output dir')
  _create_project_dir(output_root)
  log_file = open(output_root + 'tmp.log', 'w')
  cmd_list = ['bash', self_src + 'preprocess_fast.bash']
  cmd_list.extend(['--output_root', output_root])
  # 必须参数 参数什么都不加
  _check_para('--read1', read1, str, isFile=True, required=True, output_root=output_root)
  cmd_list.extend(['--read1', read1])
  _check_para('--bowtie_index', bowtie_index, str, required=True, output_root=output_root)
  cmd_list.extend(['--bowtie_index', bowtie_index])

  # 可选参数
  # cur_cores = _get_cores(cores)
  cur_cores = cores
  cmd_list.extend(['--cores', str(cur_cores)])

  if read2:
    _check_para('--read2', read2, str, isFile=True, required=True, output_root=output_root)
    cmd_list.extend(['--read2', read2])

  if adapter1:
    _check_para('--adapter1', adapter1, str, required=True, output_root=output_root)
    cmd_list.extend(['--adapter1', adapter1])

  if adapter2:
    _check_para('--adapter2', adapter2, str, required=True, output_root=output_root)
    cmd_list.extend(['--adapter2', adapter2])

  # print( ' '.join(cmd_list) )
  log_file.write(_cur_time()+ 'preprocess start'+'\n')
  logger.info('preprocess--running script')
  print( ' '.join(cmd_list) )
  try:
    os.system( ' '.join(cmd_list) )
  except Exception as e:
    _write_err('fast preprocess--Failed\n' + str(e) )

  # todo 没有uniquemapped_sort.bam
  _tracks(output_root=output_root, bam=output_root + 'sam/uniquemapped_sort.bam')

  log_file.close()
  logger.success(_cur_time() + 'preprocess--Finished. Find the results in ' + output_root)


def _align(read1=None, bowtie_index=None, gtf=None, output_root=None, cores=1, read2=None ):
  logger.info('alignment--check installed software')
  _check_soft( ['bowtie2', 'samtools', 'bedtools'], isPython=False)

  logger.info('alignment--check parameter and input files')
  if not output_root.endswith('/'): output_root = output_root + '/'
  logger.info('alignment--construct output dir')
  _create_project_dir(output_root)
  log_file = open(output_root + 'tmp.log', 'a+' )
  cmd_list = ['bash', self_src + 'map1.bash']
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
  try:
    os.system( ' '.join(cmd_list) )
  except Exception as e:
    _write_err('original_mapping--Failed\n' + str(e) )

  try:
    os.system( ' '.join(['python', self_src + 'map_split.py',
    '--sam_dir='+output_root+'sam/', '--sam_file=original.sam',
    '--output_root='+output_root
    ]) )
  except Exception as e:
    _write_err('map_split--Failed\n' + str(e) )

  cmd_list[1] = self_src + 'map2.bash'
  try:
    os.system( ' '.join(cmd_list) )
  except Exception as e:
    _write_err('map_comutation--Failed\n' + str(e) )

  log_file.close()
  logger.success(_cur_time()+ 'alignment--Finished. Find the results in ' + output_root)

def _tracks(bam=None, output_root=None, scale_factor=None):
  logger.info('genome_tracks--check installed software')
  _check_soft( ['deeptools'], isPython=False)

  logger.info('genome_tracks--check parameter and input files')
  if not output_root.endswith('/'): output_root = output_root + '/'
  logger.info('genome_tracks--construct output dir')
  _create_project_dir(output_root)
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
  if scale_factor:
    cmd_list.extend(['--scaleFactor', str(scale_factor)])
  cmd_list_forward, cmd_list_reverse = cmd_list.copy(), cmd_list.copy()
  cmd_list_forward.extend(['--filterRNAstrand', 'forward'])
  cmd_list_reverse.extend(['--filterRNAstrand', 'reverse'])
  cmd_list_forward.extend(['-o', output_root + 'bw/forward.bw'])
  cmd_list_reverse.extend(['-o', output_root + 'bw/reverse.bw'])
  log_file.write(_cur_time()+ ': generate forward genome track start\n')
  try:
    subprocess.run(' '.join(cmd_list_forward),  shell=True,  check=True)
    cmd_list_forward[-1] = output_root + 'bw/forward_5_end.bw'
    cmd_list_forward.extend(['--Offset', '1'])
    subprocess.run(' '.join(cmd_list_forward),  shell=True,  check=True)
  except Exception as e:
    _write_err('forward bamCoverage--Failed\n' + str(e) )

  log_file.write(_cur_time()+ ': generate reverse genome track start\n')
  try:
    subprocess.run(' '.join(cmd_list_reverse),  shell=True,  check=True)
    cmd_list_reverse[-1] = output_root + 'bw/reverse_5_end.bw'
    cmd_list_reverse.extend(['--Offset', '1'])
    subprocess.run(' '.join(cmd_list_reverse),  shell=True,  check=True)
  except Exception as e:
    _write_err('reverse bamCoverage--Failed\n' + str(e) )

  log_file.close()
  logger.success(_cur_time()+'genome_tracks--Finished. Find the results in ' + output_root)

def feature_assign(forward_bw=None, reverse_bw=None, gtf=None, output_root=None):
  global global_root
  global_root = output_root
  logger.info('feature_assign--check installed software')
  _check_soft( ['pandas', 'pyBigWig', 'scipy'], isPython=True)

  logger.info('feature_assign--check parameter and input files')
  if not output_root.endswith('/'): output_root = output_root + '/'
  logger.info('feature_assign--construct output dir')
  _create_project_dir(output_root)
  log_file = open(output_root + 'tmp.log', 'w' )

  cmd_list = ['python', self_src + 'feature_attrs.py']
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
  except Exception as e:
    _write_err('feature_assign--Failed\n' + str(e) )

  log_file.write(_cur_time()+ 'feature assign finished\n')
  log_file.close()
  _render_template(output_root=output_root, type='quantification', is_server='No')

def pausing_sites(forward_bw=None, reverse_bw=None, output_root=None, cores=1):
  global global_root
  global_root = output_root
  logger.info('pausing_sites--check installed software')
  _check_soft( ['pandas', 'pyBigWig', 'numpy'], isPython=True)

  logger.info('pausing_sites--check parameter and input files')
  if not output_root.endswith('/'): output_root = output_root + '/'
  logger.info('pausing_sites--construct output dir')
  _create_project_dir(output_root)
  log_file = open(output_root + 'tmp.log', 'w' )

  cmd_list = ['python', self_src + 'pausing_sites_low_memory.py']

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
  except Exception as e:
    _write_err('pausing_sites--Failed\n' + str(e) )

  _render_template(output_root=output_root, type='pausing', is_server='No')

def network_analysis(tf_source=None, tf_filter_nodes=None, enhancer_source=None, enhancer_filter_nodes=None, output_root=None ):
  global global_root
  global_root = output_root
  logger.info('network_analysis--check installed software')
  _check_soft(['networkx','community', 'hvplot', 'numpy'], isPython=True)

  cmd_list = ['python', self_src + 'network_analysis.py']

  logger.info('network_analysis--check parameter and input files')
  if not output_root.endswith('/'): output_root = output_root + '/'
  logger.info('network_analysis--construct output dir')
  cmd_list.extend(['--output_root', output_root])
  _create_project_dir(output_root)
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

  try:
    os.system( ' '.join(cmd_list) )
    logger.success('network analysis--Finished. Find the results in ' + output_root)
    _render_template(output_root=output_root, type='network', is_server='No')
  except Exception as e:
    _write_err('network analysis--Failed\n' + str(e) )

def network_links(specie=None, region='', gtf=None, forward_bw=None, reverse_bw=None, rpkm_file=None, output_root=None):
  global global_root
  global_root = output_root
  logger.info('network_links--check installed software')
  _check_soft(['pyGenomeTracks'], isPython=False)
  logger.info('network_links--check parameter and input files')
  if not output_root.endswith('/'):
    output_root = output_root + '/'
  logger.info('network_links--construct output dir')
  _create_project_dir(output_root)
  log_file = open(output_root + 'tmp.log', 'w')
  cmd_list = ['python', self_src + 'genome_track_visualization.py']
  cmd_list.extend(['--output_root', output_root])
  # 必须参数 参数什么都不加
  if not specie:
    print( 'you must provide specie')
    os.sys.exit(1)
  if region == '':
    print( 'you must provide regulatory region with --region i.e --region chr1:1:5000000')
    os.sys.exit(1)
  if not gtf:
    print( 'you must provide gtf file')
    os.sys.exit(1)
  if not forward_bw:
    print( 'you must provide forward bw file')
    os.sys.exit(1)
  if not reverse_bw:
    print( 'you must provide reverse bw file')
    os.sys.exit(1)

  cmd_list.extend(['--specie', specie])
  cmd_list.extend(['--region', region])
  cmd_list.extend(['--gtf', os.path.abspath(gtf)])
  cmd_list.extend(['--forward', os.path.abspath(forward_bw)])
  cmd_list.extend(['--reverse', os.path.abspath(reverse_bw)])
  # 非必须参数
  if rpkm_file:
    _check_para('--rpkm_file', rpkm_file, str, isFile=True, required=True, output_root=output_root)
    cmd_list.extend(['--rpkm_file', rpkm_file])
  print( ' '.join(cmd_list) )
  try:
    os.system( ' '.join(cmd_list) )
  except Exception as e:
    _write_err('network links--Failed\n' + str(e) )

def _render_template(type='server', output_root=None, is_server=None):
  global global_root
  global_root = output_root
  # 每个步骤对应了输出文件， 用来丰富报表的数据
  # 考虑到 要和论文的内容对应上，所以不能按照 流程来生成报告
  # 把报告分成 1 基础信息 2 qc 3 转录级别评价 4 暂停因子和暂停位点 5 网络分析
  # 只有两种情况 会用到render_template
  logger.info('render_template--check installed software')
  cmd_list = ['python', self_src + 'render_template.py']

  cmd_list.extend(['--type', type])
  if not output_root.endswith('/'): output_root = output_root + '/'
  _create_project_dir(output_root)
  cmd_list.extend(['--output_root', output_root])
  if is_server:
    cmd_list.extend( ['--is_server', is_server])
  print( ' '.join(cmd_list) )
  try:
    os.system( ' '.join(cmd_list) )
  except Exception as e:
    _write_err('render template--Failed\n' + str(e) )


  # 基础信息 table
  # 输入文件名
  # Forward_bw, Reverse_bw,

def main():
  fire.core.Display = lambda lines, out: print(*lines, file=out)
  fire.Fire({
    #'preprocess': preprocess,
    #'alignment': alignment,
    #'genome_tracks': genome_tracks,
    '--version': _ver,
    '-V': _ver,
    'assessment': assessment,
    'feature_assign': feature_assign,
    'pausing_sites': pausing_sites,
    'network_analysis': network_analysis,
    'network_links': network_links,
    'render_template': _render_template,
    'server': server,
    'all': all,
    'preprocess_fast': preprocess_fast,
  })

if __name__ == '__main__':
  main()