## Kirchmair et al., Frontiers in Immunology 2023

Data analyses for '<sup>13</sup>C tracer analysis reveals the landscape of metabolic checkpoints in human CD8<sup>+</sup> T cell differentiation and exhaustion' (Kirchmair et al., Frontiers in Immunology 2023)


 
## Environment setup
```bash
# source lib/make_env.sh # initial code to make conda envs
conda env create -n cd8 -f env/cd8.yml # recreate conda env
conda env create -n cd8_ngscm -f env/cd8_ngscm.yml # recreate conda env

conda activate cd8
Rscript -e 'devtools::install_github("AlexanderKirchmair/datamisc")'
Rscript -e 'devtools::install_github("AlexanderKirchmair/c13ms")'
Rscript -e 'install.packages("qualpalr", repos = "https://cran.wu.ac.at/")'
Rscript -e 'devtools::install_github("AlexanderKirchmair/DeLuciatoR")' # version forked from https://github.com/infotroph/DeLuciatoR

mkdir logs
```
 
Set up [NGSCheckMate-1.0.0](https://github.com/parklab/NGSCheckMate)
```bash
cd lib
git clone https://github.com/parklab/NGSCheckMate.git
echo 'SAMTOOLS=samtools' > lib/NGSCheckMate/ncm.conf
echo 'BCFTOOLS=bcftools' >> lib/NGSCheckMate/ncm.conf
echo 'REF=/data/genomes/hg38/fasta/gencode/GRCh38.primary_assembly.genome.fa' >> lib/NGSCheckMate/ncm.conf
cd ..
```


## RNA sequencing data

### Raw data download

```bash
conda activate cd8
```
 
Memory differentiation samples ([GSE234099](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE234099)):
```bash
mkdir -p data/rnaseq/MEM/00_RAW
accs=$(awk 'NR>1 {print $2 "-" $1}' "tables/GSE234099.txt")
for acc in $accs
do
  qsub lib/run_download.sh ${acc%-*} ${acc#*-} data/rnaseq/MEM/00_RAW ~/myScratch/tmp
  while [ $(qstat -s pr | grep -w -c "DOWNLOAD") -gt 3 ]; do sleep 3; done
done
```
 
Exhaustion samples ([GSE234100](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE234100)):
```bash
mkdir -p data/rnaseq/EXH/00_RAW
accs=$(awk 'NR>1 {print $2 "-" $1}' "tables/GSE234100.txt")
for acc in $accs
do
  qsub lib/run_download.sh ${acc%-*} ${acc#*-} data/rnaseq/EXH/00_RAW ~/myScratch/tmp
  while [ $(qstat -s pr | grep -w -c "DOWNLOAD") -gt 3 ]; do sleep 3; done
done
```
 
### Preprocessing

Trimming:
```bash
mkdir data/rnaseq/MEM/01_TRIMMED
for file in data/rnaseq/MEM/00_RAW/*fastq.gz
do
  qsub lib/run_trimming.sh $file data/rnaseq/MEM/01_TRIMMED
done

mkdir data/rnaseq/EXH/01_TRIMMED
for file in data/rnaseq/EXH/00_RAW/*fastq.gz
do
  qsub lib/run_trimming.sh $file data/rnaseq/EXH/01_TRIMMED
done
```

Read alignment and quantification using the [nf-core/rnaseq-3.4](https://nf-co.re/rnaseq/3.4) pipeline (set genome paths in `lib/run_rnaseq.sh`):
```bash
bash -i lib/run_rnaseq.sh 'tables/samplesheet_mem.csv' 'data/rnaseq/MEM/02_NF_results'
mv .nextflow.log logs/mem.nextflow.log

bash -i lib/run_rnaseq.sh 'tables/samplesheet_exh.csv' 'data/rnaseq/EXH/02_NF_results'
mv .nextflow.log logs/exh.nextflow.log
```


Check if the paired samples are matching with [NGSCheckMate-1.0.0](https://github.com/parklab/NGSCheckMate):
```bash
conda activate cd8_ngscm

mkdir data/rnaseq/MEM/samplecheck
ls -d data/rnaseq/MEM/02_NF_results/star_salmon/*bam > data/rnaseq/MEM/samplecheck/files.txt
qsub lib/run_NGSCheckMate.sh 'data/rnaseq/MEM/samplecheck/files.txt' 'data/rnaseq/MEM/samplecheck'
mv r_script.r.Rout data/rnaseq/MEM/samplecheck/

mkdir data/rnaseq/EXH/samplecheck
ls -d data/rnaseq/EXH/02_NF_results/star_salmon/*bam > data/rnaseq/EXH/samplecheck/files.txt
qsub lib/run_NGSCheckMate.sh 'data/rnaseq/EXH/samplecheck/files.txt' 'data/rnaseq/EXH/samplecheck'
mv r_script.r.Rout data/rnaseq/EXH/samplecheck/

Rscript lib/plot_NGSCheckMate.R
```

Gene sets were prepared by running `Rscript lib/prepare_genesets.R`.


## Metabolomics data
<sup>13</sup>C metabolomics data: `data/metabolomics`
 

## Seahorse data
Seahorse data: `data/seahorse`
 

## Analysis

The main analyses can be reproduced by rendering the the .Rmd files:

```bash
conda activate cd8
Rscript -e "rmarkdown::render('analyses/01-RNA-Differentiation.Rmd')"
Rscript -e "rmarkdown::render('analyses/02-13C-Differentiation.Rmd')"
Rscript -e "rmarkdown::render('analyses/03-RNA-Exhaustion.Rmd')"
Rscript -e "rmarkdown::render('analyses/04-13C-Exhaustion.Rmd')"
Rscript -e "rmarkdown::render('analyses/05-RNA-Exhaustion-Public.Rmd')"
Rscript -e "rmarkdown::render('analyses/06-RNA-Mitochondria.Rmd')"
Rscript -e "rmarkdown::render('analyses/07-Public-Dataset-Comparison.Rmd')"
```


## Results

To reproduce the final figures and tables, run:

```bash
conda activate cd8
Rscript -e "rmarkdown::render('analyses/08-Results.Rmd')"
```

