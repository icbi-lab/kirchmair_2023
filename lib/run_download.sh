#!/usr/bin/bash

#$ -S /usr/bin/bash
#$ -pe smp 1
#$ -cwd
#$ -V
#$ -N DOWNLOAD
#$ -o logs/$JOB_NAME-$JOB_ID.log
#$ -e logs/$JOB_NAME-$JOB_ID.err

acc=$1
name=$2
outdir=$(readlink -f $3)
tmpdir=$(readlink -f $4)

fasterq-dump $acc --threads 1 --temp $tmpdir --outdir $outdir
mv $outdir/$acc.fastq $outdir/$name.fastq
gzip $outdir/$name.fastq
