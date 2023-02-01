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

echo "mapping starts:"
## bam

if ! [ -x "$(command -v bowtie2)" ]; then
  echo 'Error: bowtie2 is not installed.' >&2
  exit 1
fi

# 1 比对
echo "  mapping--1 bowtie2 mapping:"
if [ $read2 ]; then
  bowtie2 -1 $read1 -2 $read2 -x $bowtie_index -S $sam_dir'original.sam' --threads $cpu &> $txt_dir'bowtie2_out.txt'
else
  bowtie2 -U $read1 -x $bowtie_index -S $sam_dir'original.sam' --threads $cpu &> $txt_dir'bowtie2_out.txt'
fi
cat $txt_dir'bowtie2_out.txt' | awk '{print"    "$0}'

# 2 split sam
echo "  mapping--2 split mapping reads:"