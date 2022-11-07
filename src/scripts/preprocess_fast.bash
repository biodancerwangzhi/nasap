ARGUMENT_LIST=(
  "read1"
  "read2"
  "output_root"
  "cores"
  "adapter1"
  "adapter2"
  "bowtie_index"
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
    --bowtie_index)
      bowtie_index=$2
      shift 2
      ;;
    *)
      break
      ;;
  esac
done

json_dir=$output_dir'json/'
fq_dir=$output_dir'fastq/'
sam_dir=$output_dir'sam/'
txt_dir=$output_dir'txt/'
img_dir=$output_dir'imgs/'

# todo
# rm -rf $output_dir
# mkdir -p $output_dir $json_dir $fq_dir $txt_dir $sam_dir $img_dir $bed_dir

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

if [ $read2 ]; then
# 这步 具有adapter的 reads 统计量不同，因为 peppro设置 cutadapt -O 为1
# -O 1表示 adapter 在3‘端只有 1bp也被看做有adapter，这明显不合理，fastp没有这个选项
# 看了fastp的结果统计，
  if [[ $adapter1 && $adapter2 ]]; then
    fastp -G -Q --adapter_sequence $adapter1 --thread $cpu \
      -A -l 16 --cut_front --cut_tail \
      --trim_poly_x \
      -q 20 \
      --adapter_sequence_r2 $adapter2 \
      -i $read1 \
      -o $fq_dir'clean1.fq.gz' \
      -I $read2 \
      -O $fq_dir'clean2.fq.gz' \
      --failed_out $fq_dir'failed_adapter.fq.gz' \
      --unpaired1 $fq_dir'failed_adapter_unpair1.fq.gz' \
      --unpaired2 $fq_dir'failed_adapter_unpair2.fq.gz' \
      --json $json_dir'remove_adapter.json' \
      --html $tmp_html &>$tmp_txt
  else
    fastp -G -Q --detect_adapter_for_pe --thread $cpu \
      -A -l 16 --cut_front --cut_tail \
      --trim_poly_x \
      -q 20 \
      -i $read1 \
      -o $fq_dir'clean1.fq.gz' \
      -I $read2 \
      -O $fq_dir'clean2.fq.gz' \
      --failed_out $fq_dir'failed_adapter.fq.gz' \
      --unpaired1 $fq_dir'failed_adapter_unpair1.fq.gz' \
      --unpaired2 $fq_dir'failed_adapter_unpair2.fq.gz' \
      --json $json_dir'remove_adapter.json' \
      --html $tmp_html &>$tmp_txt
  fi
else
  if [[ $adapter1 ]]; then
    fastp -G -Q --adapter_sequence $adapter1 --thread $cpu \
      -A -l 16 --cut_front --cut_tail \
      --trim_poly_x \
      -q 20 \
      -i $read1 \
      -o $fq_dir'clean1.fq.gz' \
      --failed_out $fq_dir'failed_adapter.fq.gz' \
      --json $json_dir'remove_adapter.json' \
      --html $tmp_html &>$tmp_txt
  else
    fastp -G -Q --thread $cpu \
      -A -l 16 --cut_front --cut_tail \
      --trim_poly_x \
      -q 20 \
      -i $read1 \
      -o $fq_dir'clean1.fq.gz' \
      --failed_out $fq_dir'failed_adapter.fq.gz' \
      --json $json_dir'remove_adapter.json' \
      --html $tmp_html &>$tmp_txt
  fi
fi

echo "mapping starts:"

if ! [ -x "$(command -v bowtie2)" ]; then
  echo 'Error: bowtie2 is not installed.' >&2
  exit 1
fi

if [ $read2 ]; then
  bowtie2 -1 $fq_dir'clean1.fq.gz' -2 $fq_dir'clean2.fq.gz' -x $bowtie_index -S $sam_dir'original.sam' --threads $cpu &> $txt_dir'bowtie2_out.txt'
else
  bowtie2 -U $fq_dir'clean1.fq.gz' -x $bowtie_index -S $sam_dir'original.sam' --threads $cpu &> $txt_dir'bowtie2_out.txt'
fi

# extract unqiue mapped reads
echo 'unique mapping'
grep -E "@|NM:" $sam_dir'original.sam' | grep -v "XS:" > $sam_dir'uniquemapped.sam'
samtools view -@ $cpu -Sb $sam_dir'uniquemapped.sam' > $sam_dir'uniquemapped.bam'
samtools sort -@ $cpu -o $sam_dir'uniquemapped_sort.bam' $sam_dir'uniquemapped.bam'
samtools index -@ $cpu $sam_dir'uniquemapped_sort.bam'

