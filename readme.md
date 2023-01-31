## nASAP  
Quality Assessment and Comprehensive Analyses of Nascent RNA Profiling Datasets

## Conda Environments to Manage Dependencies for nASAP
```bash
cat > env.yml <<EOF
name: nasap_env1
dependencies:
  - libcurl
  - bioconda::deeptools
  - bioconda::pybigwig
  - bioconda::bioawk
  - bioconda::fastp=0.22.0
  - bioconda::flash=1.2.11
  - bioconda::bowtie2
  - bioconda::fastq-pair=1.0
  - bioconda::samtools=1.13
  - bioconda::bedtools=2.30.0
  - pip:
    - --index-url https://pypi.tuna.tsinghua.edu.cn/simple
    - nasap
EOF
conda env create -f ./env.yml
```

## Test the package installed 
```bash 
conda activate nasap_env
nasap --help 
```


## build a docker mirror and start a container 
```bash 
cd Docker 
docker build -t biodancer/nasap:latest . 
sudo docker run --rm -it -v `pwd`/data:/tmp -w /tmp biodancer/nasap:latest nasap --help 
```


## Usage 
```bash 
nasap all --read1 your_fastq.fq.gz --bowtie_index your_bowtie_index --gtf your_gtf --output_root output_dir
```
more usage information, refer to the document  
http://grobase.top/nasap_doc_en/1overview/  
