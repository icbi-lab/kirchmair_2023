#!/usr/bin/bash

#$ -S /usr/bin/bash
#$ -pe smp 1
#$ -cwd
#$ -V
#$ -N SNPCHECK
#$ -o logs/$JOB_NAME-$JOB_ID.log
#$ -e logs/$JOB_NAME-$JOB_ID.err

# NGSCheckMate
in=$(readlink -f $1)
out=$(readlink -f $2)
bed=$(readlink -f lib/NGSCheckMate/SNP/SNP_GRCh38_hg38_wChr.bed)
snp=$(readlink -f lib/NGSCheckMate/SNP/SNP.pt)

# BAM/VCF mode
NGSCheckMate=$(readlink -f lib/NGSCheckMate/ncm.py)
python $NGSCheckMate -B -l $in -bed $bed -O $out
