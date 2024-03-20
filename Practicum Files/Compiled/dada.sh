#!/bin/bash
#Bash script to generate ASVs data
#sbatch -o $PROJECT/Practicum/16s_analysis/ASV/dada.out -t 6:00:00 -N 1 --ntasks-per-node=40 -J dada $PROJECT/Practicum/Scripts/dada.sh
export TZ='UTC'
cd $PROJECT/Practicum/16s_analysis/ASV
singularity exec -B "$HOME/qiime2:/home/qiime2" -B "$SCRATCH" -B "$PROJECT"  qiime2-practicum.sif\
  qiime dada2 denoise-paired\
  --i-demultiplexed-seqs paired-end-karelia.qza\
  --o-table asv-karelia-table\
  --o-representative-sequences asv-karelia-rep-seqs\
  --p-trim-left-f 0\
  --p-trim-left-r 0\
#trim length was decided by visual inspection
  --p-trunc-len-f 150\
  --p-trunc-len-r 150\
  --o-denoising-stats karelia-dada-stats\
  --verbose
exit