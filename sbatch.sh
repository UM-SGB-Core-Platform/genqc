#!/bin/bash
#SBATCH --account=def-mliu
#SBATCH --time=24:00:00 # 2 hours
#SBATCH --ntasks=1
#SBATCH --mem=32G

# Please don't make any changes to this section
##############################################
export NXF_SINGULARITY_CACHEDIR=$HOME/projects/def-mliu/NXF_CACHE_DIR

module load singularity
module load java/jdk13.0.1
module load admixture/1.3.0
module load gcc/11.2
# This version of python includes required package 'numpy' and 'pandas' on the UofM Grex cluster.
module load python/3.10.4

pip install openpyxl
##############################################


# From here, make changes to the pipline parameters to provide:
# --dataDir=The Location of your plink input files.
# --dataName=The name of your plink files.  Plink files share the same base name with .fam, .bim, .bed, OR .ped and .map
# --referenceDataDir=The location of your reference dataset files.
# --referenceDataName=The base name of the reference dataset files.  Same plink format as --refereceDataname
~/UofM_QC/nextflow ~/UofM_QC/QC.nf -profile singularity --dataDir=/home/umstua02/UKBioBank/ --dataName=ukbiobank_genotype --referenceDataDir=/home/umstua02/ --referenceDataName=referenceDataset -resume
