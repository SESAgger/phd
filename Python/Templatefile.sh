#!/bin/bash -l
# NOTE the -l flag!

#SBATCH -J J_Allcrc
#SBATCH -o O_Allcrc
#SBATCH -e E_Allcrc
#SBATCH -t 00:02:0
#SBATCH -A snic2017-7-384
#SBATCH -p core -n 1
#Mail updates when job starts and finishes
#SBATCH --mail-user sophie.agger@imbim.uu.se
#SBATCH --mail-type=ALL

module load bioinfo-tools  GATK/4.1.1.0

gatk CollectReadCounts -L canfam3.1.preprocessed.interval_list -R /proj/uppstore2017228/KLT.05.GRUS/DATA_merge/canFam3.1.fa -imr OVERLAPPING_ONLY -I  NAME1               -O NAME1
