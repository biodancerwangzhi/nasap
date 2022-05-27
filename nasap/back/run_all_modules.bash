# common functions
pkg_exist(){
  Tool=$1
  if ! [ -x "$(command -v $Tool)" ]; then
    echo "Error: $Tool is not installed." >&2
    exit 1
  fi
}

dir_exist(){
  FILE=$1
  if [ ! -d "$FILE" ]; then
    echo "$FILE directory not exist."
    exit 1
  fi
}

file_exist(){
  FILE=$1
  if [ ! -f "$FILE" ]; then
    echo "$FILE does not exist."
    exit 1
  fi
}
# 0 提取参数
ARGUMENT_LIST=(
  "bowtie_index"
  "read1"
  "read2"
  "output"
  "cores"
  "gtf"
)

opts=$(getopt \
  --longoptions "$(printf "%s:," "${ARGUMENT_LIST[@]}")" \
  --name "$(basename "$0")" \
  --options "" \
  -- "$@"
)

eval set -- $opts

# 参数列表
# 单端/双端
# 标记umi: 是/否
# 去接头: 提供adapter(序列文本或序列文件)/自己检测
# dedup: 是否用UMI

# (之后考虑是否参数了) 去polyX: 是/否 # 保留最短的reads长度: 推荐 16bp
while [[ $# -gt 0 ]]; do
  case "$1" in
    --bowtie_index)
      bowtie_index=$2
      shift 2
      ;;

    --gtf)
      gtf=$2
      shift 2
      ;;

    --read1)
      read1=$2
      shift 2
      ;;

    --read2)
      read2=$2
      shift 2
      ;;

    --output)
      output_dir=$2
      shift 2
      ;;

    --cores)
      cores=$2
      shift 2
      ;;

    *)
      break
      ;;
  esac
done

file_exist $read1

if [ $read2 ]; then
  file_exist $read2
fi

# check cpu para
cpu=1
total_cpu=`grep -c "model name" /proc/cpuinfo`
if [ $cores ]; then
  if [ "$cores" = 'max' ]; then
    cpu=$total_cpu
  fi
  if [ "$cores" = 'max/2' ]; then
    cpu=$(( $total_cpu/2 ))
  fi
  if [ -n "$(echo $cores| sed -n "/^[0-9]\+$/p")" ];then
    if (( $cores > $total_cpu )); then
      cpu=$total_cpu
    else
      cpu=$cores
    fi
  fi
fi

# software dependency:
pkg_exist bioawk
pkg_exist fastp
pkg_exist bowtie2
pkg_exist samtools
# pkg_exist umi_tools
pkg_exist python

# output parameter
if [ -z $output_dir ]; then
  # 给个默认值
  output_dir='./tmp_output/'
fi

# 给输出文件最后添加 /
echo "$output_dir" | grep -q -E '\/$' && output_dir=$output_dir || output_dir=$output_dir'/'

echo "all parameters here: "
echo "  bowtie_index: "$bowtie_index
echo "  read1: "$read1
echo "  read2: "$read2
echo "  output:"$output_dir
echo "  cores:"$cores
echo "  gtf:"$gtf

# 1 preprocess_map
echo "preprocess and mapping are running:"
if [ $read2 ]; then
bash /home/nasap_performance/nasap/preprocess_map.bash \
--bowtie_index $bowtie_index --cores $cores \
--read1 $read1 --read2 $read2 --output $output_dir
else
bash /home/nasap_performance/nasap/preprocess_map.bash \
--bowtie_index $bowtie_index --cores $cores \
--read1 $read1 --output $output_dir
fi
# 2 feature assign

# 3 network analysis

# 4 extract results
# if [ $read2 ]; then
# python /home/nasap_performance/nasap/gen_res_report.py $output_dir pair
# else
# python /home/nasap_performance/nasap/gen_res_report.py $output_dir
# fi

