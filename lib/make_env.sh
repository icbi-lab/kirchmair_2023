
# Main conda env
envname="cd8"
installconda="nextflow=21.10.0 cutadapt=4.0 r-base=4.1.3 r-essentials \
              bioconductor-isocorrector bioconductor-org.hs.eg.db bioconductor-clusterprofiler bioconductor-deseq2 bioconductor-ihw \
              bioconductor-geoquery bioconductor-complexheatmap bioconductor-biomart bioconductor-limma bioconductor-tximport bioconductor-tcgabiolinks \
              r-nlme r-glmmtmb r-contrast r-multcomp r-betareg r-msigdbr r-r.utils r-dendsort r-pals \
              r-ggpubr r-dplyr r-tibble r-openxlsx r-ggplot2 r-ggrepel r-devtools r-crayon"

conda create -n $envname -y
conda activate $envname
conda config --env --add channels defaults
conda config --env --add channels r
conda config --env --add channels bioconda
conda config --env --add channels conda-forge

conda install $installconda -v -y
Rscript -e 'devtools::install_github("AlexanderKirchmair/datamisc")'
Rscript -e 'devtools::install_github("AlexanderKirchmair/c13ms")'

conda env export | grep -v "^prefix: " > env/$envname'.yml'
conda deactivate


# Conda env for NGSCheckMate
envname="cd8_ngscm"
conda create -n $envname -y
conda activate $envname
conda install python=2.7 samtools=1.3.1 -v -y
conda env export | grep -v "^prefix: " > env/$envname'.yml'
conda deactivate

