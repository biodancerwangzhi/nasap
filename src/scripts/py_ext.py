import os, sys, json, shutil, subprocess, multiprocessing
from shutil import which

def get_file_size( filePath ):
  fsize = os.path.getsize(filePath)
  fsize = fsize/float(1024*1024)
  return str(round(fsize,2) ) + 'Mb'

def wc_count(file_name):
  out = subprocess.getoutput("wc -l %s" % file_name)
  return int(out.split()[0])

def json2dic(json_dir):
  with open(json_dir, 'r') as f:
    dict = json.load(fp=f)
  return dict

def dic2json(dict, json_dir):
  with open(json_dir, 'w') as f:
    json.dump(dict, f)

def two_list_2_dict(key_list, value_list):
  return dict(zip(key_list, value_list))

def read_tmp_variable(name=None):
  # 写死内部数据
  if not os.path.exists('./tmp_variable.json'):
    dic2json({}, './tmp_variable.json' )
  if name:
    return json2dic('./tmp_variable.json')[name]
  else:
    return json2dic('./tmp_variable.json')

def write_tmp_variable(name, value):
  tmp_variable_dic = read_tmp_variable()
  tmp_variable_dic[name] = value
  dic2json(tmp_variable_dic, './tmp_variable.json')

def write_string_to_file(string,  file):
  with open(file, 'w') as f:
    f.write(string)

def create_output_file( output_dir ):
  output_img_dir = os.path.join(output_dir, 'imgs/')
  output_express_dir = os.path.join(output_dir, 'expression/')
  output_pause_dir = os.path.join(output_dir, 'pausing/')

  setDir( output_dir )
  os.mkdir( output_img_dir )
  os.mkdir( output_express_dir )
  os.mkdir( output_pause_dir )
  write_tmp_variable("output_dir", output_dir)
  write_tmp_variable("output_img_dir", output_img_dir)
  write_tmp_variable("output_express_dir", output_express_dir)
  write_tmp_variable("output_pause_dir", output_pause_dir)

def remove_tmp():
  # 所有tmp开头的文件删掉
  for file in os.listdir('./'):
    if file.startswith('tmp_'):
      os.remove(file)

def add_string_to_file(string, file):
  with open(file, 'a+') as f:
    if not string.endswith('\n'): string = string + '\n'
    f.write(string )

def write_error(string):
  output = read_tmp_variable('output_dir')
  add_string_to_file(string, os.path.join( output,'./error.log'))
  print(string)
  sys.exit()

def write_log(string):
  output = read_tmp_variable('output_dir')
  add_string_to_file(string, os.path.join( output,'./process.log'))
  print(string)

def write_warning(string):
  output = read_tmp_variable('output_dir')
  add_string_to_file(string, os.path.join( output,'./warning.log'))
  print(string)

def print_error( error_str):
  print(error_str)
  sys.exit()

def filter_list_none(raw_list):
  return list(filter(None, raw_list ))

def copy_dir(source_dir, target_dir):
    source_path = os.path.abspath(source_dir)
    target_path = os.path.abspath( target_dir )
    if not os.path.exists(target_path):
      os.makedirs(target_path)

    if os.path.exists(source_path):
      for root, dirs, files in os.walk(source_path):
        for file in files:
          src_file = os.path.join(root, file)
          shutil.copy(src_file, target_path)

def make_zip(source_dir, output_filename):
    # zipf = zipfile.ZipFile(output_filename, 'w')
    # pre_len = len(os.path.dirname(source_dir))
    # for parent, dirnames, filenames in os.walk(source_dir):
    #     for filename in filenames:
    #         pathfile = os.path.join(parent, filename)
    #         arcname = pathfile[pre_len:].strip(os.path.sep)     #相对路径
    #         zipf.write(pathfile, arcname)
    # zipf.close()
    # 上面压缩，浏览器检测错误，无法下载
    os.system(' tar czvf %s %s'%(output_filename, source_dir) )

def setDir(path):
  folder = os.path.exists(path)
  if folder:
    shutil.rmtree(path)
  os.makedirs(path)

def zip_plus_move(output_dir, output_gz_file, mv_dir):
  make_zip(output_dir, output_gz_file)
  try:
    shutil.move(output_gz_file, mv_dir)
  except:
    pass

def run_bash_script(cmd, stdout_file=None, shell = True):
  '''
  Run a command
  如果返回值以Error: 开头就是运行失败, 反之就运行成功, 返回值就是output
  程序需要输出什么文本结果 根据指定的参数, 不在这个函数的控制范围内
  '''
  if not stdout_file:
    # 把 output 输出到控制台
    process = subprocess.Popen(cmd, shell = shell, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
  else:
    with open( stdout_file, 'w') as f:
      process = subprocess.Popen(cmd, shell = shell, stdout = f, stderr=subprocess.PIPE)

  stdout, stderr = process.communicate()
  return_code = process.poll()
  if return_code != 0:
    # 把 错误信息 打印到 控制台
    print("Process failed!")
    print(cmd)
    print(stderr)
    return('Error: ' + stderr)
  return(stdout)

def run_in_parallel(input_list, args, func, kwargs_dict = None, cpus = None, onebyone = False):
    '''
    Take an input list, divide into chunks and then apply a function to each of the chunks in parallel.
    '''
    if not cpus:
        #divide by two to get the number of physical cores
        #subtract one to leave one core free
        cpus = int(os.cpu_count()/2 - 1)  # 12超线程, 算出来cpus为5
    elif cpus == "all":
        cpus = os.cpu_count()
    # 一个无用参数的位置用于标记 第几轮并行
    arg_to_parallelize = args.index("foo")
    if not onebyone:
        # 列表中0-20，chunk_list是 [[0, 5, 10, 15], [1, 6, 11, 16], [2, 7, 12, 17], [3, 8, 13, 18], [4, 9, 14, 19]]
        chunk_list = [input_list[i::cpus] for i in range(cpus)]
    else:
        #each element in the input list will constitute a chunk of its own.
        chunk_list = input_list
    pool = multiprocessing.Pool(cpus)
    results = []
    for i in chunk_list:
        current_args = args.copy()
        # current_args = copy.deepcopy(args)
        # 把input数据放进参数中
        current_args[arg_to_parallelize] = i
        if kwargs_dict:
            process = pool.apply_async(func, tuple(current_args), kwargs_dict)
        else:
            process = pool.apply_async(func, tuple(current_args))
        results.append(process)
    pool.close()
    pool.join()
    return(results)

def is_tool(name):
  """Check whether `name` is on PATH and marked as executable."""
  # from whichcraft import which
  return which(name) is not None