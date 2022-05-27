diff_fq_length(){
  # 用 $1 $2 为 fq文件传参
  # 把 name 和 length 写进文件，从而比较差集
  bioawk -c fastx '{print $name "\t" length($seq)}' $1 | sort > ./tmp_output/tmp_fq_len1.txt
  bioawk -c fastx '{print $name "\t" length($seq)}' $2 | sort > ./tmp_output/tmp_fq_len2.txt
  # diff_reads_num=$(comm -23 ./tmp_output/tmp_fq_len1.txt ./tmp_output/tmp_fq_len2.txt | wc -l)
  # echo $diff_reads_num
  # return $diff_reads_num
  comm -23 ./tmp_output/tmp_fq_len1.txt ./tmp_output/tmp_fq_len2.txt | wc -l
}

get_fq_length(){
  bioawk -c fastx '{print $name "\t" length($seq)}' $1 | sort > $2
  num=$(more $2 | wc -l)
  echo $num
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
      cpu=$2
      shift 2
      ;;

    *)
      break
      ;;
  esac
done

# check meta data 依赖：
# file_exist $gtf

# check input data 依赖
# read1="SRR064994.fastq"
# read1="../share_data/test_data/test_r1.fq.gz"
# update running env
variable_output=$output_dir'variable.txt'

json_dir=$output_dir'json/'
fq_dir=$output_dir'fastq/'
txt_dir=$output_dir'txt/'
sam_dir=$output_dir'sam/'
img_dir=$output_dir'imgs/'
bed_dir=$output_dir'bed/'

# rm -rf $output_dir
# mkdir -p $output_dir $json_dir $fq_dir $txt_dir $sam_dir $img_dir $bed_dir

tmp_txt=$output_dir'tmp.txt'
tmp_html=$output_dir'tmp.html'
tmp_fq=$output_dir'tmp.fq.gz'

#### 1 fq level 思路
#preprocess每一步本质
## 1 计算从而去掉一些bp
## 2 去掉bp后产生一些失败的reads(失败原因是太短，质量太低)
## 3 剩下的reads
## 4 输出日志，输出失败的reads数目

# total reads
#:<<BLOCK
if [ $read2 ]; then
bioawk -t -c fastx 'END {print "  Raw read1 number:"NR}' $read1
bioawk -t -c fastx 'END {print "  Raw read2 number:"NR}' $read2
else
bioawk -t -c fastx 'END {print "  Raw read1 number:"NR}' $read1
fi
#BLOCK

# 1.1 标记UMI(可选)
    # mark umi 如果有umi, umi会被加到 seq name后面 :UMI, 序列长度会从 100到100-umi长度
    # fastp -i $read1 -o ./tmp_output/umi.fq.gz -A -G -Q -L --umi --umi_loc read1 --umi_len 6 &>tmp.txt

# 1.1 去接头  同时筛选短于 16bp的reads->(失败reads占比 + 剩下insert分布) 提示建库测序片段太短
# ( 可以给接头、也可以自动检测 )

echo "  preprocess--1 remove adapter:"
# :<<BLOCK
if [ $read2 ]; then
fastp -G -Q -l 16 --adapter_sequence TGGAATTCTCGGGTGCCAAGG \
  -i $read1 \
  -o $fq_dir'filter_adapter1.fq.gz' \
  --failed_out $fq_dir'failed_adapter1.fq.gz' \
  --json $json_dir'remove_adapter1.json' \
  --html $tmp_html &>$tmp_txt
fastp -G -Q -l 16 --adapter_sequence GATCGTCGGACTGTAGAACTCTGAAC \
  -i $read2 \
  -o $fq_dir'filter_adapter2.fq.gz' \
  --failed_out $fq_dir'failed_adapter2.fq.gz' \
  --json $json_dir'remove_adapter2.json' \
  --html $tmp_html &>$tmp_txt

success_adapter_num1=$(get_fq_length $fq_dir'filter_adapter1.fq.gz' $txt_dir'filter_adapter_fq_len1.txt')
success_adapter_num2=$(get_fq_length $fq_dir'filter_adapter2.fq.gz' $txt_dir'filter_adapter_fq_len2.txt')
failed_adapter_num1=$(get_fq_length $fq_dir'failed_adapter1.fq.gz' $txt_dir'failed_adapter_fq_len1.txt')
failed_adapter_num2=$(get_fq_length $fq_dir'failed_adapter2.fq.gz' $txt_dir'failed_adapter_fq_len2.txt')
echo_log 2 $success_adapter_num1 $failed_adapter_num1 $success_adapter_num2 $failed_adapter_num2

else
fastp -G -Q -l 16 --adapter_sequence TGGAATTCTCGGGTGCCAAGG \
  -i $read1 \
  -o $fq_dir'filter_adapter1.fq.gz' \
  --failed_out $fq_dir'failed_adapter1.fq.gz' \
  --json $json_dir'remove_adapter1.json' \
  --html $tmp_html &>$tmp_txt

success_adapter_num1=$(get_fq_length $fq_dir'filter_adapter1.fq.gz' $txt_dir'filter_adapter_fq_len1.txt')
failed_adapter_num1=$(get_fq_length $fq_dir'failed_adapter1.fq.gz' $txt_dir'failed_adapter_fq_len1.txt')
echo_log 1 $success_adapter_num1 $failed_adapter_num1

fi
# BLOCK
bioawk -c fastx 'BEGIN {short=0} {if(length($seq) < 3) short +=1} END {print "    short seq total:",short}' $fq_dir'failed_adapter1.fq.gz'

# 1.2 trim 两端 低质量bp 同时筛选短于 16bp的reads -> (失败reads占比 + )
echo "  preprocess--2 trimming low-quality bases from both ends:"
:<<BLOCK
BLOCK
if [ $read2 ]; then
fastp -A -G -Q -l 16 --cut_front --cut_tail \
  -i $fq_dir'filter_adapter1.fq.gz' \
  -o $fq_dir'filter_trim1.fq.gz' \
  --failed_out $fq_dir'failed_trim1.fq.gz' \
  --json $json_dir'cut_twoEnd1.json' \
  --html $tmp_html &>$tmp_txt

fastp -A -G -Q -l 16 --cut_front --cut_tail \
  -i $fq_dir'filter_adapter2.fq.gz' \
  -o $fq_dir'filter_trim2.fq.gz' \
  --failed_out $fq_dir'failed_trim2.fq.gz' \
  --json $json_dir'cut_twoEnd2.json' \
  --html $tmp_html &>$tmp_txt

success_trim_num1=$(get_fq_length $fq_dir'filter_trim1.fq.gz' $txt_dir'filter_trim_fq_len1.txt')
success_trim_num2=$(get_fq_length $fq_dir'filter_trim2.fq.gz' $txt_dir'filter_trim_fq_len2.txt')
failed_trim_num1=$(get_fq_length $fq_dir'failed_trim1.fq.gz' $txt_dir'failed_trim_fq_len1.txt')
failed_trim_num2=$(get_fq_length $fq_dir'failed_trim2.fq.gz' $txt_dir'failed_trim_fq_len2.txt')
echo_log 2 $success_trim_num1 $failed_trim_num1 $success_trim_num2 $failed_trim_num2

else
fastp -A -G -Q -l 16 --cut_front --cut_tail \
  -i $fq_dir'filter_adapter1.fq.gz' \
  -o $fq_dir'filter_trim1.fq.gz' \
  --failed_out $fq_dir'failed_trim1.fq.gz' \
  --json $json_dir'cut_twoEnd1.json' \
  --html $tmp_html &>$tmp_txt

success_trim_num1=$(get_fq_length $fq_dir'filter_trim1.fq.gz' $txt_dir'filter_trim_fq_len1.txt')
failed_trim_num1=$(get_fq_length $fq_dir'failed_trim1.fq.gz' $txt_dir'failed_trim_fq_len1.txt')
echo_log 1 $success_trim_num1 $failed_trim_num1
fi
# echo 'reads_with_cutTwoEnd--'$reads_with_cutTwoEnd>>$variable_output


# BLOCK
# Trimmed_reads=$(diff_fq_length $read1 ./tmp_output/trim.fq.gz)
# Trimmed_reads=$(diff_fq_length ./tmp_output/umi.fq.gz ./tmp_output/trim.fq.gz)
# echo "trimmed_reads--$Trimmed_reads"
# echo "trimmed_reads--$Trimmed_reads" >>$output

# 1.3 去polyX  同时筛选短于 16bp的reads ->(具有polyX的reads占比，) 提示TES暂停
echo "  preprocess--3 remove polyX:"

if [ $read2 ]; then
fastp -A -G -Q -l 16 --trim_poly_x \
  -i $fq_dir'filter_trim1.fq.gz' \
  -o $fq_dir'filter_polyX1.fq.gz' \
  --failed_out $fq_dir'failed_polyX1.fq.gz' \
  --json $json_dir'remove_polyX1.json' \
  --html $tmp_html &>$tmp_txt

fastp -A -G -Q -l 16 --trim_poly_x \
  -i $fq_dir'filter_trim2.fq.gz' \
  -o $fq_dir'filter_polyX2.fq.gz' \
  --failed_out $fq_dir'failed_polyX2.fq.gz' \
  --json $json_dir'remove_polyX2.json' \
  --html $tmp_html &>$tmp_txt

success_polyX_num1=$(get_fq_length $fq_dir'filter_polyX1.fq.gz' $txt_dir'filter_polyX_fq_len1.txt')
success_polyX_num2=$(get_fq_length $fq_dir'filter_polyX2.fq.gz' $txt_dir'filter_polyX_fq_len2.txt')
failed_polyX_num1=$(get_fq_length $fq_dir'failed_polyX1.fq.gz' $txt_dir'failed_polyX_fq_len1.txt')
failed_polyX_num2=$(get_fq_length $fq_dir'failed_polyX2.fq.gz' $txt_dir'failed_polyX_fq_len2.txt')
echo_log 2 $success_polyX_num1 $failed_polyX_num1 $success_polyX_num2 $failed_polyX_num2

else
fastp -A -G -Q -l 16 --trim_poly_x \
  -i $fq_dir'filter_trim1.fq.gz' \
  -o $fq_dir'filter_polyX1.fq.gz' \
  --failed_out $fq_dir'failed_polyX1.fq.gz' \
  --json $json_dir'remove_polyX1.json' \
  --html $tmp_html &>$tmp_txt

success_polyX_num1=$(get_fq_length $fq_dir'filter_polyX1.fq.gz' $txt_dir'filter_polyX_fq_len1.txt')
failed_polyX_num1=$(get_fq_length $fq_dir'failed_polyX1.fq.gz' $txt_dir'failed_polyX_fq_len1.txt')
echo_log 1 $success_polyX_num1 $failed_polyX_num1

fi

# 1.4 QC 去 平均质量低(q20) 的reads->(失败reads占比) 提示建库测序片段太短
    # 使用 序列的平均质控 筛选原始序列 phred 筛选以q20为标准

echo "  preprocess--4 drop low quality reads(mean phred quality less than 20):"
# :<<BLOCK
if [ $read2 ]; then
fastp -A -G -q 20 \
  -i $fq_dir'filter_polyX1.fq.gz' \
  -o $fq_dir'filter_quality1.fq' \
  --failed_out $fq_dir'failed_quality1.fq.gz' \
  --json $json_dir'filter_quality1.json' \
  --html $tmp_html &>$tmp_txt

fastp -A -G -q 20 \
  -i $fq_dir'filter_polyX2.fq.gz' \
  -o $fq_dir'filter_quality2.fq' \
  --failed_out $fq_dir'failed_quality2.fq.gz' \
  --json $json_dir'filter_quality2.json' \
  --html $tmp_html &>$tmp_txt
success_quality_num1=$(get_fq_length $fq_dir'filter_quality1.fq' $txt_dir'filter_quality_fq_len1.txt')
success_quality_num2=$(get_fq_length $fq_dir'filter_quality2.fq' $txt_dir'filter_quality_fq_len2.txt')
failed_quality_num1=$(get_fq_length $fq_dir'failed_quality1.fq.gz' $txt_dir'failed_quality_fq_len1.txt')
failed_quality_num2=$(get_fq_length $fq_dir'failed_quality2.fq.gz' $txt_dir'failed_quality_fq_len2.txt')
echo_log 2 $success_quality_num1 $failed_quality_num1 $success_quality_num2 $failed_quality_num2
else
fastp -A -G -q 20 \
  -i $fq_dir'filter_polyX1.fq.gz' \
  -o $fq_dir'filter_quality1.fq.gz' \
  --failed_out $fq_dir'failed_quality1.fq.gz' \
  --json $json_dir'filter_quality1.json' \
  --html $tmp_html &>$tmp_txt
success_quality_num1=$(get_fq_length $fq_dir'filter_quality1.fq.gz' $txt_dir'filter_quality_fq_len1.txt')
failed_quality_num1=$(get_fq_length $fq_dir'failed_quality1.fq.gz' $txt_dir'failed_quality_fq_len1.txt')
echo_log 1 $success_quality_num1 $failed_quality_num1
fi

if [ $read2 ]; then
echo "  preprocess--5 remove unpaired reads"
fastq_pair -t 5000000 $fq_dir'filter_quality1.fq' $fq_dir'filter_quality2.fq'
mv $fq_dir'filter_quality1.fq.paired.fq' $fq_dir'clean_read1.fq.gz'
mv $fq_dir'filter_quality2.fq.paired.fq' $fq_dir'clean_read2.fq.gz'
else
mv $fq_dir'filter_quality1.fq.gz' $fq_dir'clean_read1.fq.gz'
fi

# mv $fq_dir'clean_read1.fq.gz' $fq_dir'filter_quality1.fq.gz'
# mv $fq_dir'clean_read2.fq.gz' $fq_dir'filter_quality2.fq.gz'

# todo : delete later
# find $fq_dir* | grep -v filter_cutTwoEnd* | xargs rm

# 4 insert distribution plot
# mv ./tmp_output/tmp_fq_len2.txt ./tmp_output/fragment_distribution.txt

echo "  preprocess--finished, the clean reads in "$fq_dir

## bam
# 1 比对
echo "  mapping--1 bowtie2 mapping:"
if [ $read2 ]; then
bowtie2 -1 $fq_dir'clean_read1.fq.gz' -2 $fq_dir'clean_read2.fq.gz' -x $bowtie_index -S $sam_dir'original.sam' --threads $cpu &> $txt_dir'bowtie2_out.txt'
else
bowtie2 -U $fq_dir'clean_read1.fq.gz' -x $bowtie_index -S $sam_dir'original.sam' --threads $cpu &> $txt_dir'bowtie2_out.txt'
fi
cat $txt_dir'bowtie2_out.txt' | awk '{print"    "$0}'



# 2 质控 基于mapq
# samtools view -h -q 10 -U $sam_dir'fail_qc.sam' $sam_dir'original.sam' -@ $cpu -o $sam_dir'filter_qc.sam'
# # bioawk -c sam 'BEGIN{total=0} {total = total + 1 } END {print "mapped_qc_sam--"total}' ./tmp_output/filter_qc.sam >> $output
# fail_qc_sam_num=$(samtools view -@ $cpu -c $sam_dir'fail_qc.sam')
# echo 'fail_qc_num--'$fail_qc_sam_num>> $variable_output

## deduplication(这步放的位置 要考究)
# 2 sam 文件 拆分  original.sam -> 1 unmap.sam 2 map.sam
# map.sam -> 1 low_quality.sam 2 high_quality.sam
# high_quality.sam -> 1 unique_map.sam 2 multiple_map.sam
echo "  mapping--2 split mapping reads:"
samtools view -@ $cpu -h -f 4 -o $sam_dir'map.sam' -U $sam_dir'unmap.sam' $sam_dir'original.sam'
samtools view -@ $cpu -h -q 10 -o $sam_dir'high_quality.sam' -U $sam_dir'low_quality.sam' $sam_dir'map.sam'
grep -v "XS:" $sam_dir'high_quality.sam' >$sam_dir'unique_map.sam'
# sam 统计量
total_sam=$(samtools view -@ $cpu -c $sam_dir'original.sam')
echo 'total_num--'$total_sam >> $variable_output
unmap_sam=$(samtools view -@ $cpu -c $sam_dir'unmap.sam')
echo 'unmap_num--'$unmap_sam >> $variable_output
map_sam=$(samtools view -@ $cpu -c $sam_dir'map.sam')
echo 'map_num--'$map_sam >> $variable_output
low_quality_sam=$(samtools view -@ $cpu -c $sam_dir'low_quality.sam')
echo 'low_quality_num--'$low_quality_sam >> $variable_output
high_quality_sam=$(samtools view -@ $cpu -c $sam_dir'high_quality.sam')
echo 'high_quality_num--'$high_quality_sam >> $variable_output
unique_sam=$(samtools view -@ $cpu -c $sam_dir'unique_map.sam')
echo 'unique_num--'$unique_sam >> $variable_output

# 3 unique_map压缩 排序 index
samtools view -@ $cpu -Sb $sam_dir'unique_map.sam' > $sam_dir'unique_map.bam'
samtools sort -@ $cpu -o $sam_dir'unique_map_sort.bam' $sam_dir'unique_map.bam'
samtools index -@ $cpu $sam_dir'unique_map_sort.bam'


# 4 bam基础统计量 包括complexity depth coverage
# samtools depth ./tmp_output/sort.bam |  awk '{sum+=$3} END { print "average_depth--"sum/NR}' >>$output


bedtools bamtobed -i $sam_dir'unique_map_sort.bam' | \
awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$6}' |  \
grep -v 'chrM' | \
sort |  \
uniq -c |  \
awk 'BEGIN{mt=0;m0=0;m1=0;m2=0} ($1==1){m1=m1+1} ($1==2){m2=m2+1} {m0=m0+1} {mt=mt+$1} END{m1_m2=-1.0; if(m2>0) m1_m2=m1/m2; printf "%d\t%d\t%d\t%d\t%f\t%f\t%f\n",mt,m0,m1,m2,m0/mt,m1/m0,m1_m2}' > $output_dir'pbc_qc.txt'

# 输出文件 7列:
# 1 TotalReadPairs
# 2 DistinctReadPairs
# 3 OneReadPair
# 4 TwoReadPairs
# 5 NRF=Distinct/Total
# 6 PBC1=OnePair/Distinct
# 7 PBC2=OnePair/TwoPair

# 双端的
# bedtools bamtobed -bedpe -i $sam_dir'unique_map_sort.bam' | \
# awk 'BEGIN{OFS="\t"}{print $1,$2,$4,$6,$9,$10}' | \
# grep -v 'chrM' | \
# sort | \
# uniq -c | \
# awk 'BEGIN{mt=0;m0=0;m1=0;m2=0} ($1==1){m1=m1+1} ($1==2){m2=m2+1} {m0=m0+1} {mt=mt+$1} END{m1_m2=-1.0; if(m2>0) m1_m2=m1/m2; printf "%d\t%d\t%d\t%d\t%f\t%f\t%f\n",mt,m0,m1,m2,m0/mt,m1/m0,m1_m2}' > pbc_qc.txt

# 统计 线粒体 mt 含量
samtools idxstats $sam_dir'unique_map_sort.bam' | grep 'chrM' | cut -f 3

# bam coverage
# deeptools bamCoverage

# 5 feature assign
# 提取基于gtf的 bed(左右扩展3kb)附近的reads
cat $gtf |  awk 'OFS="\t" {if ($3=="gene") {print $1,$4-3001,$5+3000,$10,".",$7}}' | tr -d '";' | awk '{if ($2 >1) {print}}' > $bed_dir'annotation_extend_3k.bed'
samtools view -@ $cpu -L $bed_dir'annotation_extend_3k.bed' -o $sam_dir'assign.sam' -U $sam_dir'unassign.sam' $sam_dir'unique_map_sort.bam'
assign_sam=$(samtools view -@ $cpu -c $sam_dir'assign.sam')
echo 'assign_num--'$assign_sam >> $variable_output


python $(pwd)'/nasap_performance/nasap/x.py' $output_dir




# 三种 dedup方法
# 1 利用umi  umi-tools
# 2 利用序列完全相同  fastp --dedup
# 3 (自创)先比对 排除掉比对到promoter近端的序列 然后利用sam的unique去重
# fastp -i test_r1_filter1_r1.fq -o test_r1_filter2_r1.fq -A -G -Q -L --dedup &>tmp.txt
# echo 'dedup ' && python parse_json.py

#  with umi
# umi_tools dedup -I ./tmp_output/sort.bam --umi-separator=":" -S ./tmp_output/sort_dedup.bam -L ./tmp_output/dedup.log

#

