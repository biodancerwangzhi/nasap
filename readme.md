## nASAP  
Quality Assessment and Comprehensive Analyses of Nascent RNA Profiling Datasets

## Conda Environments to Manage Dependencies for nASAP
```bash
cat > env.yml <<EOF
name: nasap_env
dependencies:
  - libcurl=7.87.0
  - scipy=1.5.2
  - pandas=1.1.5
  - hvplot=0.7.3
  - pip=21.2.2
  - bioconda::deeptools=2.5.7
  - bioconda::pybigwig=0.3.17
  - bioconda::pysam=0.8.3
  - bioconda::bioawk=1.0
  - bioconda::bedtools=2.30.0
  - bioconda::bowtie2=2.3.5.1
  - bioconda::fastp=0.22.0
  - bioconda::flash=1.2.11
  - bioconda::fastq-pair=1.0
  - bioconda::samtools=1.13
  - pip:
    - nasap==0.2.10
EOF
conda env create -f ./env.yml
```

## Test the package installed 
```bash 
conda activate nasap_env
nasap --help 
```


## pull a docker mirror and start a container 
```bash 
docker pull biodancer/nasap:latest 
sudo docker run --rm -it -v `pwd`/data:/tmp -w /tmp biodancer/nasap:latest nasap --help 
```

## Usage 
```bash 
nasap all --read1 your_fastq.fq.gz --bowtie_index your_bowtie_index --gtf your_gtf --output_root output_dir
```
more usage information, refer to the document  
http://grobase.top/nasap_doc_en/1overview/  

## A example for using both docker container with test data. 
1 download test data here.
[onedrive](https://1drv.ms/u/s!AvDRT-KhJYcpgwuo35GDStZ5A5LI?e=bmsEZE)  

2 unzip the compressed files 
```bash
unzip test_nasap_data.zip
```
and the files directory 
```
test_nasap_data  
|--Homo_sapiens.GRCh38.93.gtf  
|--Homo_sapiens.GRCh38_tf_target.txt  
|--Homo_sapiens.GRCh38_enhancer_target.txt  
|--test_r1.fq.gz  
|--hg38_bowtie2_index  
|--|--hg38_bowtie2_index.1.bt2  
|--|--hg38_bowtie2_index.2.bt2  
|--|--hg38_bowtie2_index.3.bt2  
|--|--hg38_bowtie2_index.4.bt2  
|--|--hg38_bowtie2_index.rev.1.bt2  
|--|--hg38_bowtie2_index.rev.2.bt2  
```

3 init a docker container 
```bash
cd test_nasap_data && docker run --rm -it -v $(pwd):/home -w /home --name nasap_container biodancer/nasap /bin/bash
```

in docker container run the script  
```bash
nasap assessment --read1 ./test_r1.fq.gz  --adapter1 TGGAATTCTCGGGTGCCAAGG --bowtie_index ./hg38_bowtie2_index/hg38_bowtie2_index --gtf ./Homo_sapiens.GRCh38.93.gtf --output_root ./test_out
```