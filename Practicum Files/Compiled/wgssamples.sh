#!/bin/bash
#sbatch -o $PROJECT/Practicum/Shotgun_analysis/wgs.out -p debug -t 1:00:00 -N 1 --ntasks-per-node=1 -J wgssubsample $PROJECT/Practicum/Scripts/wgssamples.sh
module load NiaEnv/2019b intelpython3  
source activate myPythonEnv
module load gcc/8.3.0
module load tbb/2019u9
module load bowtie2/2.4.4
module load gnu-parallel/20191122
cd $PROJECT/Practicum/DATA/KareliaMetagenomic/Complete/Metaphlan
#save the current directory
CWD="$(pwd)"
#make the directories
while read p; do  mkdir "$p" ; done <sampling_methods.txt                                                          
#move the required files to the correcti directory
while read p; do
	while read q; do
	 cp "$q"_profile.txt  "$p"; 
	done <"$p".txt 
done <sampling_methods.txt
#Go into each directory, and merge the files as needed
while read p; do
  cd "$p" ;
  merge_metaphlan_tables.py *_profile.txt > "$p"abundance_table.txt
  cd "$CWD" #Change back to old directory
  done <sampling_methods.txt
exit