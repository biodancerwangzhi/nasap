ARGUMENT_LIST=(
  "bowtie_index"
  "gtf"
  "read1"
  "read2"
  "output_root"
  "cores"
)

opts=$(getopt \
  --longoptions "$(printf "%s:," "${ARGUMENT_LIST[@]}")" \
  --name "$(basename "$0")" \
  --options "" \
  -- "$@"
)

eval set -- $opts

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
    --output_root)
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


fq_dir=$output_dir'fastq/'
txt_dir=$output_dir'txt/'
sam_dir=$output_dir'sam/'
img_dir=$output_dir'imgs/'
bed_dir=$output_dir'bed/'


variable_report=$output_dir'csv/mapping_report.csv'
if [ -f $variable_report ]; then
  rm $variable_report
fi

# 接着map1.bash


# 2 质控 基于mapq(delete)
# samtools view -h -q 10 -U $sam_dir'fail_qc.sam' $sam_dir'original.sam' -@ $cpu -o $sam_dir'filter_qc.sam'
# # bioawk -c sam 'BEGIN{total=0} {total = total + 1 } END {print "mapped_qc_sam--"total}' ./tmp_output/filter_qc.sam >> $output
# fail_qc_sam_num=$(samtools view -@ $cpu -c $sam_dir'fail_qc.sam')
# echo 'fail_qc_num--'$fail_qc_sam_num>> $variable_output

## deduplication(这步放的位置 要考究)
# 2 sam 文件 拆分  original.sam -> 1 unmap.sam 2 map.sam
# map.sam -> 1 low_quality.sam 2 high_quality.sam
# high_quality.sam -> 1 unique_map.sam 2 multiple_map.sam
# 1 unmap= flag without AS and XS 或者 flag==4 这个值不等于 align 0 times
# more $sam_dir'original.sam' | grep -v "AS:" | grep -v "XS:" | wc -l
# 12,298,548


# samtools view -@ $cpu -h -f 4 $sam_dir'original.sam' | wc -l
# samtools view -@ $cpu -h -f 4 -o $sam_dir'unmap.sam' -U $sam_dir'map.sam' $sam_dir'original.sam'

# 2

# samtools view -@ $cpu -h -q 10 -o $sam_dir'high_quality.sam' -U $sam_dir'low_quality.sam' $sam_dir'map.sam'
# grep -v "XS:" $sam_dir'map.sam' >$sam_dir'unique_map.sam'
# grep -v "AS:|XS:" $sam_dir'original.sam' >$sam_dir'unmap.sam'
# # sam 统计量
# total_sam=$(samtools view -@ $cpu -c $sam_dir'original.sam')
# echo 'total_num--'$total_sam >> $variable_output
# unmap_sam=$(samtools view -@ $cpu -c $sam_dir'unmap.sam')
# echo 'unmap_num--'$unmap_sam >> $variable_output
# map_sam=$(samtools view -@ $cpu -c $sam_dir'map.sam')
# echo 'map_num--'$map_sam >> $variable_output
# low_quality_sam=$(samtools view -@ $cpu -c $sam_dir'low_quality.sam')
# echo 'low_quality_num--'$low_quality_sam >> $variable_output
# high_quality_sam=$(samtools view -@ $cpu -c $sam_dir'high_quality.sam')
# echo 'high_quality_num--'$high_quality_sam >> $variable_output
# unique_sam=$(samtools view -@ $cpu -c $sam_dir'unique_map.sam')
# echo 'unique_num--'$unique_sam >> $variable_output
# exit 1


# 4 bam基础统计量 包括complexity
# 这里sort为了使用bam2bed可以用的格式
# echo "  mapping--3 compute complexity:"
samtools view -@ $cpu -Sb $sam_dir'uniquemapped.sam' > $sam_dir'uniquemapped.bam'
samtools sort -n -@ $cpu -o $sam_dir'uniquemapped_sort.bam' $sam_dir'uniquemapped.bam'

cat $gtf |  awk 'OFS="\t" {if ($3=="gene") {print "chr"$1,$4-3001,$5+3000,$10,".",$7}}' | tr -d '";' | awk '{if ($2 >1) {print}}' > $bed_dir'annotation_extend_3k.bed'
samtools view -h -@ $cpu -L $bed_dir'annotation_extend_3k.bed' -o $sam_dir'assign.sam' -U $sam_dir'unassign.sam' $sam_dir'uniquemapped_sort.bam'
assign_sam=$(samtools view -@ $cpu -c $sam_dir'assign.sam')
echo 'assign_mapped,'$assign_sam >> $variable_report

if [ $read2 ]; then
# 输出文件 7列:
# 1 TotalReadPairs
# 2 DistinctReadPairs
# 3 OneReadPair
# 4 TwoReadPairs
# 5 NRF=Distinct/Total
# 6 PBC1=OnePair/Distinct
# 7 PBC2=OnePair/TwoPair
bedtools bamtobed -bedpe -i $sam_dir'uniquemapped_sort.bam' | \
awk 'BEGIN{OFS="\t"}{print $1,$2,$4,$6,$9,$10}' | \
grep -v 'chrM' | \
sort | \
uniq -c | \
awk 'BEGIN{mt=0;m0=0;m1=0;m2=0} ($1==1){m1=m1+1} ($1==2){m2=m2+1} {m0=m0+1} {mt=mt+$1} END{m1_m2=-1.0; if(m2>0) m1_m2=m1/m2; printf "%d\t%d\t%d\t%d\t%f\t%f\t%f\n",mt,m0,m1,m2,m0/mt,m1/m0,m1_m2}' > $txt_dir'pbc_qc.txt'

else
bedtools bamtobed -i $sam_dir'uniquemapped_sort.bam' | \
awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$6}' |  \
grep -v 'chrM' | \
sort |  \
uniq -c |  \
awk 'BEGIN{mt=0;m0=0;m1=0;m2=0} ($1==1){m1=m1+1} ($1==2){m2=m2+1} {m0=m0+1} {mt=mt+$1} END{m1_m2=-1.0; if(m2>0) m1_m2=m1/m2; printf "%d\t%d\t%d\t%d\t%f\t%f\t%f\n",mt,m0,m1,m2,m0/mt,m1/m0,m1_m2}' > $txt_dir'pbc_qc.txt'

fi
more $txt_dir'pbc_qc.txt' | awk '{print "NRF,"$5}' >> $variable_report
more $txt_dir'pbc_qc.txt' | awk '{print "PBC1,"$6}' >> $variable_report
more $txt_dir'pbc_qc.txt' | awk '{print "PBC2,"$7}' >> $variable_report

# 5 sort and  index
# echo "  mapping--4 sort and index:"
samtools sort -@ $cpu -o $sam_dir'uniquemapped_sort.bam' $sam_dir'uniquemapped.bam'
samtools index -@ $cpu $sam_dir'uniquemapped_sort.bam'
# depth coverage
# samtools depth ./tmp_output/sort.bam |  awk '{sum+=$3} END { print "average_depth--"sum/NR}' >>$output

# 统计 线粒体 mt 含量
# echo "  mapping--5 mapping reads stat:"
chrM_reads_num=$(samtools idxstats $sam_dir'uniquemapped_sort.bam' | grep 'chrM' | cut -f 3)
echo 'chrM_mapped,'$chrM_reads_num >> $variable_report


# 三种 dedup方法
# 1 利用umi  umi-tools
# 2 利用序列完全相同  fastp --dedup
# 3 (自创)先比对 排除掉比对到promoter近端的序列 然后利用sam的unique去重
# fastp -i test_r1_filter1_r1.fq -o test_r1_filter2_r1.fq -A -G -Q -L --dedup &>tmp.txt
# echo 'dedup ' && python parse_json.py

#  with umi
# umi_tools dedup -I ./tmp_output/sort.bam --umi-separator=":" -S ./tmp_output/sort_dedup.bam -L ./tmp_output/dedup.log
