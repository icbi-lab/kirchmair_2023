#!/bin/bash

in=$(readlink -f $1)
out=$(readlink -f $2)

wd=$(readlink -f ~/myScratch/tmp)
config=$(readlink -f env/nf.config)

fasta="/data/genomes/hg38/fasta/gencode/GRCh38.primary_assembly.genome.fa"
gtf="/data/genomes/hg38/annotation/gencode/gencode.v38.primary_assembly.annotation.gtf"
star_idx="/data/genomes/hg38/index/STAR/2.7.9a/gencode/gencode.v38.GRCh38.primary_assembly.genome/200"

nextflow run nf-core/rnaseq -r 3.4 \
  -profile singularity,cd8 \
  -w $wd \
  -c $config \
  --input $in \
  --outdir $out \
  --fasta $fasta \
  --gtf $gtf \
  --skip_trimming \
  --star_index $star_idx
