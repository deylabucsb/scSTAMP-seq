#!/bin/bash
                             
#PBS -l nodes=1:ppn=1 
                           
#PBS -l walltime=48:00:00

#User Input
TXT_DIR="/home/achialastri/scMAT-seq/InTSSByGeneCell"
START_DIR="/home/achialastri/perlscripts/Genomic_Location_Search"
#GENOME_LOCI="/home/achialastri/perlscripts/Genomic_Location_Search/LocationFiles/mm10"
GENOME_LOCI="/home/achialastri/perlscripts/Genomic_Location_Search/LocationFiles/hg19"
TXT_FILE_USE=$1

#Arg 3 Upstream nt Genomic Loci.  Arg 4 Downstream nt Genomic Loci.  Arg 5 Name of Genomic Loci
perl $START_DIR/PullOut_Near_GenomicElements.pl $GENOME_LOCI/$2 $TXT_DIR/$TXT_FILE_USE 3000 3000 $3