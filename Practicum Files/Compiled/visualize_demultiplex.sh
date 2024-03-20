#!/bin/bash
#Bash script to generate demuxed data
#sbatch -o $PROJECT/Practicum/16s_analysis/ASV/demult.out -p debug -t 01:00:00 -N 1 --ntasks-per-node=40 -J vis_demultiplex $PROJECT/Practicum/Scripts/visualize_demultiplex.sh
export TZ='UTC'
cd $PROJECT/Practicum/16s_analysis/ASV
singularity exec -B "$HOME/qiime2:/home/qiime2" -B "$SCRATCH" -B "$PROJECT"  qiime2-practicum.sif\
  qiime demux summarize \
  --i-data paired-end-karelia.qza \
  --o-visualization paired-end-karelia.qzv
exit

 