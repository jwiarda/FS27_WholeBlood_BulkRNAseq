#!/bin/bash

#SBATCH --job-name=FS27_FeatureCounts                              # name of the job submitted
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

# Run FeatureCounts on .bam files from FS27 whole blood; Jayne Wiarda 20221214

# Load subread module:
module load subread/2.0.2

# Make directory for output files:
mkdir /project/fsepru/jwiarda/FS27/FeatureCounts

# Run FeatureCounts:
cd /project/fsepru/jwiarda/FS27/Alignments_STAR
files=(`ls *.bam`)

featureCounts -a /project/fsepru/jwiarda/FS27/GeneAnnotationFiles/Sus_scrofa.Sscrofa11.1.97_modified06302021_JEW_SKS.gtf -G /project/fsepru/jwiarda/FS27/GeneAnnotationFiles/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa -F GTF -g gene_name -T 10 -s 1 -o FS27_GeneCountTable_LanesSeparate.txt ${files[@]}

mv FS27_GeneCountTable_LanesSeparate.txt /project/fsepru/jwiarda/FS27/FeatureCounts

#End of file
