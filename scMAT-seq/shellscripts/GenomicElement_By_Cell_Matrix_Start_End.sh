#!/bin/bash
                             
#PBS -l nodes=1:ppn=1 
                           
#PBS -l walltime=24:00:00

#User Input
TXT_DIR="/home/achialastri/scMAT-seq/InTSSByGeneCell"
START_DIR="/home/achialastri/perlscripts/Genomic_Location_Search"
#GENOME_LOCI="/home/achialastri/perlscripts/Genomic_Location_Search/LocationFiles/mm10"
GENOME_LOCI="/home/achialastri/perlscripts/Genomic_Location_Search/LocationFiles/hg19"
TXT_FILE_USE=$1


#Arg 3 nt Downstream TSE.  Arg 4 Upstream TSS.  Arg 5 Name of Genomic Loci
perl $START_DIR/Create_GenomicLocation_Counts_By_Cell_Matrix_Start_End.pl $GENOME_LOCI/$2 $TXT_DIR/$TXT_FILE_USE 3000 3000 $3
