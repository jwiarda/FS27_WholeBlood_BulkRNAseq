#!/bin/bash

#SBATCH --job-name=FS27_Trimmomatic                              # name of the job submitted
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

# Run Trimmomatic on raw .fastq files from FS27 whole blood; Jayne Wiarda 20221027

# Load Trimmomatic module:
module load java/1.8.0_121 # Trimmomatic requires java
module load trimmomatic/0.38

# Make directory for output files:
mkdir /project/fsepru/jwiarda/FS27/Trimmomatic

# Run Trimmomatic:
cd /90daydata/fsepru/Kristen_FS23X_FS27/lane1

for file in $(ls *-FS27*.fastq.gz | cut -d. -f1)
do
java -jar /software/7/apps/trimmomatic/0.38/trimmomatic-0.38.jar SE -threads 8 ${file}.fastq.gz /project/fsepru/jwiarda/FS27/Trimmomatic/${file}_trimmed.fastq.gz ILLUMINACLIP:/software/7/apps/trimmomatic/0.38/adapters/TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:25
done

cd /90daydata/fsepru/Kristen_FS23X_FS27/lane2

for file in $(ls *-FS27*.fastq.gz | cut -d. -f1)
do
java -jar /software/7/apps/trimmomatic/0.38/trimmomatic-0.38.jar SE -threads 8 ${file}.fastq.gz /project/fsepru/jwiarda/FS27/Trimmomatic/${file}_trimmed.fastq.gz ILLUMINACLIP:/software/7/apps/trimmomatic/0.38/adapters/TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:25
done

#End of file
