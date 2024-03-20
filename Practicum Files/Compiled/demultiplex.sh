#!/bin/bash
#Bash script to generate demuxed data
#sbatch -o $PROJECT/Practicum/16s_analysis/ASV/demult.out -t 03:00:00 -N 1 --ntasks-per-node=40 -J demultiplex $PROJECT/Practicum/Scripts/demultiplex.sh
export TZ='UTC'
cd $PROJECT/Practicum/16s_analysis/ASV
singularity exec -B "$HOME/qiime2:/home/qiime2" -B "$SCRATCH" -B "$PROJECT"  qiime2-practicum.sif\
  qiime tools import\
  --type 'SampleData[PairedEndSequencesWithQuality]'\
  --input-path /project/o/oespinga/segurahi/Practicum/16s_analysis/ASV/Karelia_manifest.tsv\
  --output-path paired-end-karelia.qza\
  --input-format PairedEndFastqManifestPhred33V2
exit

 