#!/bin/bash
#Bash script to generate taxonomy data from the ASV/OTUs using the silva classifier
#collapses based on genus level,
#and then generates differential abundance data from the ASV/OTU table
#Then generates relative frequencies 
#and then exports the table
#sbatch -o $PROJECT/Practicum/16s_analysis/ASV/taxdiff.out -t 10:00:00 -N 1 --ntasks-per-node=40 -J taxadiff $PROJECT/Practicum/Scripts/taxonomy-and-differential.sh
export TZ='UTC'
cd $PROJECT/Practicum/16s_analysis/ASV
singularity exec -B "$HOME/qiime2:/home/qiime2" -B "$SCRATCH" -B "$PROJECT"  qiime2-practicum.sif\
  qiime feature-classifier classify-sklearn \
  --i-classifier silva_132_99_v3v4_q2_2019-7.qza \
  --i-reads asv-karelia-rep-seqs.qza \
  --o-classification karelia-taxonomy-otu.qza
singularity exec -B "$HOME/qiime2:/home/qiime2" -B "$SCRATCH" -B "$PROJECT"  qiime2-practicum.sif\
  qiime taxa collapse \
  --i-table asv-karelia-table.qza \
  --i-taxonomy karelia-taxonomy-otu.qza \
  --p-level 6 \
  --o-collapsed-table asv-genus-table.qza
singularity exec -B "$HOME/qiime2:/home/qiime2" -B "$SCRATCH" -B "$PROJECT"  qiime2-practicum.sif\
  qiime composition add-pseudocount \
  --i-table asv-genus-table.qza \
  --o-composition-table comp-table-asv.qza
#Didnt work because "composition != freqeuncy"
#singularity exec -B "$HOME/qiime2:/home/qiime2" -B "$SCRATCH" -B "$PROJECT"  qiime2-practicum.sif\
 # qiime feature-table relative-frequency \
  #--i-table comp-table-otu.qza \
  #--o-relative-frequency-table otu-pseudo-rel-table.qza
singularity exec -B "$HOME/qiime2:/home/qiime2" -B "$SCRATCH" -B "$PROJECT"  qiime2-practicum.sif\
  qiime tools export \
  --input-path comp-table-asv.qza \
  --output-path exported-pseudo-asv-table
exit