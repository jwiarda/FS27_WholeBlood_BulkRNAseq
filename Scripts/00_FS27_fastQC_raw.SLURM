#!/bin/bash

#SBATCH --job-name=FS27_fastQC_raw                              # name of the job submitted
#SBATCH -p short                                        # name of the queue you are submitting to
#SBATCH -N 1                                            # number of nodes in this job
#SBATCH -n 40                                             # number of cores/tasks in this job
#SBATCH -t 48:00:00                                     # time allocated for this job hours:mins:seconds
#SBATCH -o "stdout.%j.%N"                               # standard out %j adds job number to outputfile name and %N adds the node name
#SBATCH -e "stderr.%j.%N"                               # optional but it prints our standard error
#SBATCH --mem=32G   					# memory e.g.: 100G ; 250G ; 100M etc... (need mem node for some reqs)
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jewiarda@gmail.com					# enter your email address here USDA emails sometimes dont work, I use gmail

# ENTER COMMANDS HERE:

# Run fastQC on raw .fastq files from FS27 whole blood; Jayne Wiarda 20221026

# Load fastQC module:
module load fastqc/0.11.9

# Make directory for output files:
mkdir /project/fsepru/jwiarda/FS27/fastQC_raw

# Run fastQC:
for file in /90daydata/fsepru/Kristen_FS23X_FS27/lane1/*-FS27-*
do
fastqc -f fastq -t 6 --noextract --nogroup --outdir /project/fsepru/jwiarda/FS27/fastQC_raw ${file} 
done

for file in /90daydata/fsepru/Kristen_FS23X_FS27/lane2/*-FS27-*
do
fastqc -f fastq -t 6 --noextract --nogroup --outdir /project/fsepru/jwiarda/FS27/fastQC_raw ${file}
done

#End of file
