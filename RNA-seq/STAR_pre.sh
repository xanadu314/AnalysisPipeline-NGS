#!/bin/bash

module load sratoolkit/2.9.6 
module load STAR/2.7.0c 
modle load FastQC/0.11.5-Java-1.8.0_92
module load MultiQC/1.6

## Download genome file
wget ftp://ftp.ensemblgenomes.org/pub/bacteria/release-44/fasta/bacteria_3_collection/salmonella_enterica_subsp_enterica_serovar_typhimurium_str_st4_74/dna/Salmonella_enterica_subsp_enterica_serovar_typhimurium_str_st4_74.ASM18873v1.dna.chromosome.Chromosome.fa.gz
wget ftp://ftp.ensemblgenomes.org/pub/bacteria/release-44/gff3/bacteria_3_collection/salmonella_enterica_subsp_enterica_serovar_typhimurium_str_st4_74/Salmonella_enterica_subsp_enterica_serovar_typhimurium_str_st4_74.ASM18873v1.37.gff3.gz
gzip -d *.gz

## Download SRA file
for ((i=1027;i<=1035;i++));do echo SRR95$i;done > SRAID.txt
nohup prefetch --option-flie SRAID.txt &

## convert SRA to fastq
fastq-dump --split-3 --gzip *.sra

## Fastqc
for ((i=1027;i<=1035;i++));do fastqc SRR95$i.fastq.gz ;done
# fastqc -t 10 filename.fastq.gz可以加快速度
# 也可以用命令 ls *9510*.gz |xargs fastqc -t 10
#xargs是用来传递变量
multiqc ./ #合并fastqc文件

## Genome index generate
STAR --runMode genomeGenerate --runThreadN 1 --genomeFastaFiles st4_74.fa --sjdbGTFfile st4_74.gff3 --sjdbGTFtagExonParentTranscript Parent --sjdbOverhang 49 --genomeDir ~/genome/ --genomeSAindexNbases 1
# --runThreadN 线程数，--sjdbGTFtagExonParentTranscript Parent 针对GFF3文件，--sjdbOverhang 平均reads长度-1 --genomeSAindexNbases min(14, log2(GenomeLength)/2 - 1)

