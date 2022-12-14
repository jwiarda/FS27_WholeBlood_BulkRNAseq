#!/bin/bash

#SBATCH --job-name=FS27_STAR_alignment                              # name of the job submitted
#SBATCH -p short                                        # name of the queue you are submitting to
#SBATCH -N 1                                            # number of nodes in this job
#SBATCH -n 40                                             # number of cores/tasks in this job
#SBATCH -t 48:00:00                                     # time allocated for this job hours:mins:seconds
#SBATCH -o "stdout.%j.%N"                               # standard out %j adds job number to outputfile name and %N adds the node name
#SBATCH -e "stderr.%j.%N"                               # optional but it prints our standard error
#SBATCH --mem=100G   					# memory e.g.: 100G ; 250G ; 100M etc... (need mem node for some reqs)
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jewiarda@gmail.com					# enter your email address here USDA emails sometimes dont work, I use gmail

# ENTER COMMANDS HERE:

# Run STAR alignment on trimmed .fastq files from FS27 whole blood; Jayne Wiarda 20221213

# Obtained genome file http://ftp.ensembl.org/pub/release-97/fasta/sus_scrofa/dna/ & annotation file modified from http://ftp.ensembl.org/pub/release-97/gtf/sus_scrofa/

# Load STAR module:
module load star/2.7.9a

# Generate genome index
## Make directory to store genome index in:
mkdir /project/fsepru/jwiarda/FS27/GenomeIndex_STAR

## Create index:
STAR --runThreadN 20 --runMode genomeGenerate --genomeDir /project/fsepru/jwiarda/FS27/GenomeIndex_STAR --genomeFastaFiles /project/fsepru/jwiarda/FS27/GeneAnnotationFiles/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa --sjdbGTFfile /project/fsepru/jwiarda/FS27/GeneAnnotationFiles/Sus_scrofa.Sscrofa11.1.97_modified06302021_JEW_SKS.gtf --sjdbOverhang 99

# Run STAR alignment
## Make directory for output files:
mkdir /project/fsepru/jwiarda/FS27/Alignments_STAR

## Run alignment:
cd /90daydata/fsepru/Kristen_FS23X_FS27/lane1

for file in $(ls *-FS27*.fastq.gz)
do
STAR --runThreadN 20 --genomeDir /project/fsepru/jwiarda/FS27/GenomeIndex_STAR --readFilesIn ${file} --readFilesCommand zcat --outFileNamePrefix "${file}"_ --outSAMtype BAM Unsorted --outReadsUnmapped Fastx --outFilterMultimapNmax 20 --outFilterType BySJout --outFilterMismatchNmax 30 
done

mv *{Aligned.out.bam,Aligned.out.sam,Log.final.out,Log.out,Log.progress.out,SJ.out.tab,Unmapped.out.mate1,Unmapped.out.mate2} /project/fsepru/jwiarda/FS27/Alignments_STAR

cd /90daydata/fsepru/Kristen_FS23X_FS27/lane2

for file in $(ls *-FS27*.fastq.gz)
do
STAR --runThreadN 20 --genomeDir /project/fsepru/jwiarda/FS27/GenomeIndex_STAR --readFilesIn ${file} --readFilesCommand zcat --outFileNamePrefix "${file}"_ --outSAMtype BAM Unsorted --outReadsUnmapped Fastx --outFilterMultimapNmax 20 --outFilterType BySJout --outFilterMismatchNmax 30
done

mv *{Aligned.out.bam,Log.final.out,Log.out,Log.progress.out,SJ.out.tab,Unmapped.out.mate1} /project/fsepru/jwiarda/FS27/Alignments_STAR

#End of file
