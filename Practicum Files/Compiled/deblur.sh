#!/bin/bash
#Bash script to generate OTU data
#sbatch -o $PROJECT/Practicum/16s_analysis/OTU/deblur.out -t 20:00:00 -N 1 --ntasks-per-node=40 -J deblur $PROJECT/Practicum/Scripts/deblur.sh
export TZ='UTC'
cd $PROJECT/Practicum/16s_analysis/OTU
singularity exec -B "$HOME/qiime2:/home/qiime2" -B "$SCRATCH" -B "$PROJECT"  qiime2-practicum.sif\
#Deblur requires the scores to be filtered first
  qiime quality-filter q-score \
  --i-demux paired-end-karelia.qza \
  --o-filtered-sequences karelia-filtered.qza \
  --o-filter-stats karelia-filter-stats.qza
singularity exec -B "$HOME/qiime2:/home/qiime2" -B "$SCRATCH" -B "$PROJECT"  qiime2-practicum.sif\
  qiime deblur denoise-16S \
  --i-demultiplexed-seqs karelia-filtered.qza \
#trim length was decided by visual inspection
  --p-trim-length 150 \
  --o-representative-sequences karelia-rep-seqs-deblur.qza \
  --o-table karelia-table-deblur.qza \
  --p-sample-stats \
  --o-stats karelia-deblur-stats.qza
exit