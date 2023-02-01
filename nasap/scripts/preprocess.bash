diff_fq_length(){
  # 用 $1 $2 为 fq文件传参
  # 把 name 和 length 写进文件，从而比较差集
  # bioawk -c fastx '{print $name "\t" length($seq)}' $1 | sort > ./tmp_output/tmp_fq_len1.txt
  # bioawk -c fastx '{print $name "\t" length($seq)}' $2 | sort > ./tmp_output/tmp_fq_len2.txt
  # diff_reads_num=$(comm -23 ./tmp_output/tmp_fq_len1.txt ./tmp_output/tmp_fq_len2.txt | wc -l)
  # echo $diff_reads_num
  # return $diff_reads_num
  comm -23 ./tmp_output/tmp_fq_len1.txt ./tmp_output/tmp_fq_len2.txt | wc -l
}

get_fq_length(){
  # 指定 fq文件和 read长度文件
  bioawk -c fastx '{print $name "\t" length($seq)}' $1 | sort > $2
  num=$(more $2 | wc -l)
  echo $num
}

batch_fq_length(){
  # 参数是一个 array, 单数是测序数据 双数是text文件
  #   echo "start"
  # 获取数组个数
	arr=($(echo "$@"))
  newarr=()
  num=$[ $# - 1 ]
  # echo "arr" $arr
  j=0
  for ((i=0;i<=num;i=i+2)){
    # echo ${arr[$i]} ${arr[$i+1]}
    reads_num=$(get_fq_length ${arr[$i]} ${arr[$i+1]})
    newarr[$j]=$reads_num
    j=$(( $j+1 ))
  }
  echo ${newarr[@]}
}

read_length_change_num(){
  comm_num=$(comm -12 $1 $2 | wc -l)
  a_num=$(more $1 | wc -l )
  echo $(( a_num - comm_num ))
}

parse_pair_fastp(){
  # 参数 before_read1_len, before_read2_len,
  # after_read1_len, after_read2_len, change_read
  # failed_out_len, failed_pair_read1_len, failed_pair_read2_len
  echo "pass"

}

echo_log(){
# 参数: 1 pair_end: 1或2, 2
  if (( $1==2 )); then
  echo "    pair 1 success reads number is:"$2
  echo "    pair 1 failed reads number is:"$3
  echo "    pair 2 success reads number is:"$4
  echo "    pair 2 failed reads number is:"$5
  else
  echo "    success reads number is:"$2
  echo "    failed reads number is:"$3
  fi
}

# read arguments
ARGUMENT_LIST=(
  "read1"
  "read2"
  "output_root"
  "cores"
  "adapter1"
  "adapter2"
  "umi_loc"
  "umi_len"
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
    --read1)
      read1=$2
      shift 2
      ;;
    --output_root)
      output_dir=$2
      shift 2
      ;;
    --cores)
      cpu=$2
      shift 2
      ;;

    --read2)
      read2=$2
      shift 2
      ;;
    --adapter1)
      adapter1=$2
      shift 2
      ;;
    --adapter2)
      adapter2=$2
      shift 2
      ;;
    --umi_loc)
      umi_loc=$2
      shift 2
      ;;
    --umi_len)
      umi_len=$2
      shift 2
      ;;
    *)
      break
      ;;
  esac
done

# check meta data 依赖：
# file_exist $gtf

json_dir=$output_dir'json/'
fq_dir=$output_dir'fastq/'
txt_dir=$output_dir'txt/'
img_dir=$output_dir'imgs/'

tmp_txt=$output_dir'tmp.txt'
tmp_html=$output_dir'tmp.html'
tmp_fq=$output_dir'tmp.fq.gz'
tmp_variable=$output_dir'tmp_variable.txt'
if [ -f $tmp_variable ]; then
  rm $tmp_variable
fi

if ! [ -x "$(command -v fastp)" ]; then
  echo 'Error: fastp is not installed.' >&2
  exit 1
fi

if ! [ -x "$(command -v bioawk)" ]; then
  echo 'Error: bioawk is not installed.' >&2
  exit 1
fi

#### 1 fq level 思路
# preprocess每一步本质
## 1 计算从而去掉一些bp
## 2 去掉bp后产生一些失败的reads(失败原因是太短，质量太低)
## 3剩下的reads



echo "preprocess starts:"

# 1.0 原始数据 统计
if [ $read2 ]; then
raw_read_arr=($read1 $txt_dir'read_len1.txt' $read2 $txt_dir'read_len2.txt')
arg=$(echo ${raw_read_arr[*]})
res=($(batch_fq_length $arg))
read1_num=${res[0]}
read2_num=${res[1]}
echo "  Raw read1 number:" $read1_num
echo "  Raw read2 number:" $read2_num
echo read1_num--$read1_num >> $tmp_variable
echo read2_num--$read2_num >> $tmp_variable
else
read1_num=$(get_fq_length $read1 $txt_dir'read_len1.txt')
echo "  Raw read1 number:" $read1_num
echo read1_num--$read1_num >> $tmp_variable
fi

echo "  preprocess--1 remove UMI:"
if [ $read2 ]; then
  if [[ $umi_loc ]]; then
    if [[ $umi_len ]]; then
      fastp --thread $cpu -i $read1 -I $read2 \
      -o $fq_dir'remove_umi1.fq.gz' -O $fq_dir'remove_umi2.fq.gz' \
      -U --umi_loc $umi_loc --umi_len $umi_len \
      --json $json_dir'remove_umi.json' \
      --html $tmp_html &>$tmp_txt
    else
      fastp --thread $cpu -i $read1 -I $read2 \
      -o $fq_dir'remove_umi1.fq.gz' -O $fq_dir'remove_umi2.fq.gz' \
      -U --umi_loc $umi_loc \
      --json $json_dir'remove_umi.json' \
      --html $tmp_html &>$tmp_txt
    fi
  else
    cp $read1 $fq_dir'remove_umi1.fq.gz'
    cp $read2 $fq_dir'remove_umi2.fq.gz'
  fi
else
  if [[ $umi_loc ]]; then
    echo "umi_loc"
    if [[ $umi_len ]]; then
      echo "umi_len"
      fastp --thread $cpu -i $read1 \
      -o $fq_dir'remove_umi1.fq.gz' \
      -U --umi_loc $umi_loc --umi_len $umi_len \
      --json $json_dir'remove_umi.json' \
      --html $tmp_html &>$tmp_txt
    else
      fastp --thread $cpu -i $read1 \
      -o $fq_dir'remove_umi1.fq.gz' \
      -U --umi_loc $umi_loc \
      --json $json_dir'remove_umi.json' \
      --html $tmp_html &>$tmp_txt
    fi
  else
    echo "copy read"
    cp $read1 $fq_dir'remove_umi1.fq.gz'
  fi
fi


# 1.1 标记UMI(可选)
    # mark umi 如果有umi, umi会被加到 seq name后面 :UMI, 序列长度会从 100到100-umi长度
    # fastp -i $read1 -o ./tmp_output/umi.fq.gz -A -G -Q -L --umi --umi_loc read1 --umi_len 6 &>tmp.txt

# fastp -i testdata/R1.fq -o testdata/out.R1.fq -U --umi_loc=read1 --umi_len=8

# todo
# 2 查看umi例子，思考其逻辑 参数判断
#   umi_loc 指定 index(name的第二部分) 放进name的第一部分的尾部
#   umi_loc 指定 read 是需要 提供 umi_len 会把 5'端的地方 移到 name的第一部分
#   先判断上面两个条件是否成立 同时给个状态值 UMI, 后面gencore 处理bam 去掉 duplicates

# 3 添加fastp 代码
# 4 测试一下 test.fa




# 1.2 去接头  同时筛选短于 2bp的reads->(失败reads占比 + 剩下insert分布) 提示建库测序片段太短
# ( 可以给接头、也可以自动检测 )
echo "  preprocess--2 remove adapter:"
if [ $read2 ]; then
# 这步 具有adapter的 reads 统计量不同，因为 peppro设置 cutadapt -O 为1
# -O 1表示 adapter 在3‘端只有 1bp也被看做有adapter，这明显不合理，fastp没有这个选项
# 看了fastp的结果统计，
  if [[ $adapter1 && $adapter2 ]]; then
    fastp -G -Q -l 2 --adapter_sequence $adapter1 --thread $cpu \
      --adapter_sequence_r2 $adapter2 \
      -i $fq_dir'remove_umi1.fq.gz' \
      -o $fq_dir'filter_adapter1.fq.gz' \
      -I $fq_dir'remove_umi2.fq.gz' \
      -O $fq_dir'filter_adapter2.fq.gz' \
      --failed_out $fq_dir'failed_adapter.fq.gz' \
      --unpaired1 $fq_dir'failed_adapter_unpair1.fq.gz' \
      --unpaired2 $fq_dir'failed_adapter_unpair2.fq.gz' \
      --json $json_dir'remove_adapter.json' \
      --html $tmp_html &>$tmp_txt
  else
    fastp -G -Q -l 2 --detect_adapter_for_pe --thread $cpu \
      -i $fq_dir'remove_umi1.fq.gz' \
      -o $fq_dir'filter_adapter1.fq.gz' \
      -I $fq_dir'remove_umi2.fq.gz' \
      -O $fq_dir'filter_adapter2.fq.gz' \
      --failed_out $fq_dir'failed_adapter.fq.gz' \
      --unpaired1 $fq_dir'failed_adapter_unpair1.fq.gz' \
      --unpaired2 $fq_dir'failed_adapter_unpair2.fq.gz' \
      --json $json_dir'remove_adapter.json' \
      --html $tmp_html &>$tmp_txt
  fi
else
  if [[ $adapter1 ]]; then
    fastp -G -Q -l 2 --adapter_sequence $adapter1 --thread $cpu \
      -i $fq_dir'remove_umi1.fq.gz' \
      -o $fq_dir'filter_adapter1.fq.gz' \
      --failed_out $fq_dir'failed_adapter.fq.gz' \
      --json $json_dir'remove_adapter.json' \
      --html $tmp_html &>$tmp_txt
  else
    fastp -G -Q -l 2 --thread $cpu \
      -i $fq_dir'remove_umi1.fq.gz' \
      -o $fq_dir'filter_adapter1.fq.gz' \
      --failed_out $fq_dir'failed_adapter.fq.gz' \
      --json $json_dir'remove_adapter.json' \
      --html $tmp_html &>$tmp_txt
  fi
fi

# 问题:
# 1 Reads_with_adapter: 双端 read1总数-长度没变的reads-unpair1中长度没变的reads
# 2 Uninformative_adapter_reads: read1 failed reads数
# 3 Pct_uninformative_adapter_reads: read1 failed reads数/read1总数
# 4 Peak_adapter_insertion_size: 统计去掉adapter后的insert长度 次数
# 5 Degradation_ratio: 10-20 / 30-40

if [ $read2 ]; then
  adapter_arr=($fq_dir"filter_adapter1.fq.gz" $txt_dir'adapter_read_len1.txt'
  $fq_dir"filter_adapter2.fq.gz" $txt_dir'adapter_read_len2.txt'
  $fq_dir"failed_adapter_unpair1.fq.gz" $txt_dir'failed_adapter_unpair1.txt'
  )

  arg=$(echo ${adapter_arr[*]})
  res=($(batch_fq_length $arg))
  filter_adapter_read1_num=${res[0]}
  filter_adapter_read2_num=${res[1]}
  failed_adapter_num=$(( $read1_num-$filter_adapter_read1_num ))
  failed_adapter_read1_unpair_num=${res[2]}
  failed_adapter_read2_unpair_num=$failed_adapter_num-$failed_adapter_read1_unpair_num

  echo_log 2 $filter_adapter_read1_num $failed_adapter_num $filter_adapter_read2_num $failed_adapter_num

  echo "    Failed read1 unpaired number:" $failed_adapter_read1_unpair_num
  no_adapter1_num=$(comm -12 $txt_dir'read_len1.txt' $txt_dir'adapter_read_len1.txt' | wc -l)
  # no_adapter2_num=$(comm -12 $txt_dir'read_len2.txt' $txt_dir'adapter_read_len2.txt' | wc -l)
  # failed_adapter_unpair1是 在read2中错误，导致的read1 unpair
  no_adapter1_unpair_num=$(comm -12 $txt_dir'read_len1.txt' $txt_dir'failed_adapter_unpair1.txt' | wc -l)
  # no_adapter2_unpair_num=$(comm -12 $txt_dir'read_len2.txt' $txt_dir'failed_adapter_unpair2.txt' | wc -l)
  # echo "    No adapter read1 number:" $no_adapter1_num
  # echo "    No adapter read2 number:" $no_adapter2_num
  # echo "    No adapter unpair read1 number:" $no_adapter1_unpair_num
  # echo "    No adapter unpair read2 number:" $no_adapter2_unpair_num
  # reads_with_adapter=$(( $read1_num + $read2_num - $no_adapter1_num - $no_adapter2_num - $unpair_no_adapter1 - $unpair_no_adapter2))
  # echo 'reads with adapter' $reads_with_adapter

  # bioawk -t -c fastx 'END {print "    failed reads number is:"NR}' $fq_dir'failed_short_adapter.fq'
  # bioawk -c fastx 'BEGIN {short=0} {if(length($seq) < 4) short +=1} END {print "    short seq total:",short}' $fq_dir'failed_adapter1.fq'

  # 1 Reads_with_adapter: read1总数-长度没变的reads-unpair1中长度没变的reads
  reads_with_adapter=$(( $read1_num - $no_adapter1_num - $no_adapter1_unpair_num ))
  echo "    Reads_with_adapter:" $reads_with_adapter
else
  adapter_arr=($fq_dir"filter_adapter1.fq.gz" $txt_dir'adapter_read_len1.txt'
  $fq_dir"failed_adapter.fq.gz" $txt_dir'failed_adapter.txt'
  )
  arg=$(echo ${adapter_arr[*]})
  res=($(batch_fq_length $arg))
  filter_adapter_read1_num=${res[0]}
  failed_adapter_num=$(( $read1_num-$filter_adapter_read1_num ))
  echo_log 1 $filter_adapter_read1_num  $failed_adapter_num
  no_adapter1_num=$(comm -12 $txt_dir'read_len1.txt' $txt_dir'adapter_read_len1.txt' | wc -l)
  reads_with_adapter=$(( $read1_num - $no_adapter1_num  ))
  echo "    Reads_with_adapter:" $reads_with_adapter
fi
echo reads_with_adapter--$reads_with_adapter >> $tmp_variable

# 2 Uninformative_adapter_reads: read1 failed reads数
echo "    Uninformative_adapter_reads:" $failed_adapter_num
echo uninformative_adapter_reads--$failed_adapter_num >> $tmp_variable
# 3 Pct_uninformative_adapter_reads: read1 failed reads数/read1总数
pct_uninformative_adapter_reads=$(echo | awk "{print $failed_adapter_num/$read1_num*100}")
echo "    Pct_uninformative_adapter_reads:"$pct_uninformative_adapter_reads
echo pct_uninformative_adapter_reads--$pct_uninformative_adapter_reads >> $tmp_variable

# 4 Peak_adapter_insertion_size:
if [ $read2 ]; then
  flash -q -t 40 --compress-prog=pigz --suffix=gz \
  $fq_dir'filter_adapter1.fq.gz' $fq_dir'filter_adapter2.fq.gz' \
  -o flash -d $fq_dir
  peak_adapter_insertion_size=$(cat $fq_dir'flash.hist' | sort -k2nr | head -1 | awk '{print $1;}')
  # 5 Degradation_ratio
  # degradation_ratio=$(cat $fq_dir'flash.hist' | awk -F'\t' '{if (match($1, /1[0-9]/)) low+=$2} {if (match($1, /3[0-9]/)) high+=$2} END {print low/high}')

else
  # 这里的peak，我认为是 数量最多的 insertion 长度
  peak_adapter_insertion_size=$(cat $txt_dir'adapter_read_len1.txt' | awk -F'\t' '{print $2;}' | sort | uniq -c | sort -k1nr | head -1 | awk '{print $2;}')
  # degradation_ratio=$(cat $txt_dir'adapter_read_len1.txt' | awk -F'\t' '{print $2;}' | sort | uniq -c | awk '{if (match($2, /1[0-9]/)) low+=$1} {if (match($2, /3[0-9]/)) high+=$1} END {print low/high}')
fi
echo "    Peak_adapter_insertion_size:"$peak_adapter_insertion_size
echo peak_adapter_insertion_size--$peak_adapter_insertion_size >> $tmp_variable
# echo "    Degradation_ratio:"$degradation_ratio
# echo degradation_ratio--$degradation_ratio >> $tmp_variable

echo "  preprocess--3 trimming low-quality bases from both ends:"
if [ $read2 ]; then
fastp -A -G -Q -l 16 --cut_front --cut_tail --thread $cpu \
  -i $fq_dir'filter_adapter1.fq.gz' \
  -o $fq_dir'filter_trim1.fq.gz' \
  -I $fq_dir'filter_adapter2.fq.gz' \
  -O $fq_dir'filter_trim2.fq.gz' \
  --failed_out $fq_dir'failed_trim.fq.gz' \
  --json $json_dir'trim.json' \
  --html $tmp_html &>$tmp_txt
else
fastp -A -G -Q -l 16 --cut_front --cut_tail --thread $cpu \
  -i $fq_dir'filter_adapter1.fq.gz' \
  -o $fq_dir'filter_trim1.fq.gz' \
  --failed_out $fq_dir'failed_trim.fq.gz' \
  --json $json_dir'trim.json' \
  --html $tmp_html &>$tmp_txt
fi

if [ $read2 ]; then
  trim_arr=($fq_dir'filter_trim1.fq.gz' $txt_dir'filter_trim_len1.txt'
  $fq_dir'filter_trim2.fq.gz' $txt_dir'filter_trim_len2.txt'
  )
  arg=$(echo ${trim_arr[*]})
  res=($(batch_fq_length $arg))
  filter_trim1_num=${res[0]}
  failed_trim1_num=$(( $filter_adapter_read1_num-$filter_trim1_num ))
  echo_log 2 $filter_trim1_num $failed_trim1_num $filter_trim1_num $failed_trim1_num
else
  trim_arr=($fq_dir'filter_trim1.fq.gz' $txt_dir'filter_trim_len1.txt')
  arg=$(echo ${trim_arr[*]})
  res=($(batch_fq_length $arg))
  filter_trim1_num=${res[0]}
  failed_trim1_num=$(( $filter_adapter_read1_num-$filter_trim1_num ))
  echo_log 1 $filter_trim1_num $failed_trim1_num
fi
reads_with_cutTwoEnd=$(read_length_change_num $txt_dir'adapter_read_len1.txt' $txt_dir'filter_trim_len1.txt')
echo reads_with_cutTwoEnd--$reads_with_cutTwoEnd >> $tmp_variable

# 1.4 去polyX  同时筛选短于 16bp的reads ->(具有polyX的reads占比，) 提示TES暂停
echo "  preprocess--4 remove polyX:"
if [ $read2 ]; then
fastp -A -G -Q -l 16 --trim_poly_x --thread $cpu \
  -i $fq_dir'filter_trim1.fq.gz' \
  -o $fq_dir'filter_polyX1.fq.gz' \
  -I $fq_dir'filter_trim2.fq.gz' \
  -O $fq_dir'filter_polyX2.fq.gz' \
  --failed_out $fq_dir'failed_polyX.fq' \
  --json $json_dir'remove_polyX.json' \
  --html $tmp_html &>$tmp_txt
else
fastp -A -G -Q -l 16 --trim_poly_x --thread $cpu \
  -i $fq_dir'filter_trim1.fq.gz' \
  -o $fq_dir'filter_polyX1.fq.gz' \
  --failed_out $fq_dir'failed_polyX.fq' \
  --json $json_dir'remove_polyX.json' \
  --html $tmp_html &>$tmp_txt
fi

if [ $read2 ]; then
  polyX_arr=($fq_dir'filter_polyX1.fq.gz' $txt_dir'filter_polyX_len1.txt'
  $fq_dir'filter_polyX2.fq.gz' $txt_dir'filter_polyX_len2.txt'
  )
  arg=$(echo ${polyX_arr[*]})
  res=($(batch_fq_length $arg))
  filter_polyX1_num=${res[0]}
  failed_polyX1_num=$(( $filter_trim1_num-$filter_polyX1_num ))
  echo_log 2 $filter_polyX1_num $failed_polyX1_num $filter_polyX1_num $failed_polyX1_num
else
  polyX_arr=($fq_dir"filter_polyX1.fq.gz" $txt_dir'filter_polyX_len1.txt')
  arg=$(echo ${polyX_arr[*]})
  res=($(batch_fq_length $arg))
  filter_polyX1_num=${res[0]}
  failed_polyX1_num=$(( $filter_trim1_num-$filter_polyX1_num ))
  echo_log 1 $filter_polyX1_num $failed_polyX1_num
fi


# Adaptered_reads=$(diff_fq_length ./tmp_output/trim.fq.gz ./tmp_output/trim_adapter.fq.gz)
# # Reads_with_adapter
# echo "adaptered_reads--$Adaptered_reads"
# echo "adaptered_reads--$Adaptered_reads" >>$output
# comm -3 ./tmp_output/tmp_fq_len1.txt ./tmp_output/tmp_fq_len2.txt | sed 's/^\t//' > ./tmp_output/adapter_distribution.txt

# polyX_reads=$(diff_fq_length ./tmp_output/trim_adapter.fq.gz ./tmp_output/trim_adapter_polyX.fq.gz)
# echo "polyX_reads--$polyX_reads"
# echo "polyX_reads--$polyX_reads" >>$output

# 1.5 QC 去 平均质量低(q20) 的reads->(失败reads占比) 提示建库测序片段太短
    # 使用 序列的平均质控 筛选原始序列 phred 筛选以q20为标准
echo "  preprocess--5 drop low quality reads(mean phred quality less than 20):"
if [ $read2 ]; then
fastp -A -G -q 20 --thread $cpu \
  -i $fq_dir'filter_polyX1.fq.gz' \
  -o $fq_dir'filter_quality1.fq.gz' \
  -I $fq_dir'filter_polyX2.fq.gz' \
  -O $fq_dir'filter_quality2.fq.gz' \
  --failed_out $fq_dir'failed_quality.fq.gz' \
  --json $json_dir'filter_quality.json' \
  --html $tmp_html &>$tmp_txt
else
fastp -A -G -q 20 --thread $cpu \
  -i $fq_dir'filter_polyX1.fq.gz' \
  -o $fq_dir'filter_quality1.fq.gz' \
  --failed_out $fq_dir'failed_quality.fq' \
  --json $json_dir'filter_quality.json' \
  --html $tmp_html &>$tmp_txt
fi

if [ $read2 ]; then
  quality_arr=($fq_dir'filter_quality1.fq.gz' $txt_dir'filter_quality_len1.txt'
  $fq_dir'filter_quality2.fq.gz' $txt_dir'filter_quality_len2.txt'
  )
  arg=$(echo ${quality_arr[*]})
  res=($(batch_fq_length $arg))
  filter_quality1_num=${res[0]}
  failed_quality1_num=$(( $filter_polyX1_num-$filter_quality1_num ))
  echo_log 2 $filter_quality1_num $failed_quality1_num $filter_quality1_num $failed_quality1_num
else
  quality_arr=($fq_dir'filter_quality1.fq.gz' $txt_dir'filter_quality_len1.txt')
  arg=$(echo ${quality_arr[*]})
  res=($(batch_fq_length $arg))
  filter_quality1_num=${res[0]}
  failed_quality1_num=$(( $filter_polyX1_num-$filter_quality1_num ))
  echo_log 1 $filter_quality1_num $failed_quality1_num
fi


# 4 insert distribution plot
# mv ./tmp_output/tmp_fq_len2.txt ./tmp_output/fragment_distribution.txt

if [ $read2 ]; then
mv $fq_dir'filter_quality1.fq.gz' $fq_dir'clean_read1.fq.gz'
mv $fq_dir'filter_quality2.fq.gz' $fq_dir'clean_read2.fq.gz'
else
mv $fq_dir'filter_quality1.fq.gz' $fq_dir'clean_read1.fq.gz'
fi

# todo : delete later
# find $fq_dir* | grep -v clean_read* | xargs rm

echo "  preprocess--finished, the clean reads in "$fq_dir