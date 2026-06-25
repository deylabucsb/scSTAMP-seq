#!/bin/bash

#PBS -l nodes=1:ppn=1

#PBS -l walltime=24:00:00


#User Input
START_DIR="/home/piscopio/Run54/RobRun54-357662604/RPSD071422-CE6_P1S4-395424012"


#Standard Usage per person
CEL_BARCODES="/home/piscopio/perlscripts/mRNA_Mapping/cel-seq_barcodes.csv"
PERL_DIR="/home/piscopio/perlscripts/PC_Oligo"
PERL_DIR2="/home/piscopio/perlscripts/mRNA_Mapping"
HOME_DIR="/home/piscopio"
PC_Barcodes="PC_Barcodes_All.txt"

#Do not change
PAST_DIR=${START_DIR%/*}
OUT_NAME=${START_DIR##*/}
RUN_NAME=${PAST_DIR##*/}
R1="_L001-4_R1_001.fastq"
R2="_L001-4_R2_001.fastq"
FASTQ_R1=$OUT_NAME$R1
FASTQ_R2=$OUT_NAME$R2
LOG=${FASTQ_R1%??????}
LOG2=${FASTQ_R2%??????}

CelSeq_UMI_Length=6
PC_UMI_Length=5
MP_UMI_Length=5


# (BNL) Use these lines for when data comes back from Nextseq and is in a 4 lane format. Will delete these temp files after Triming
cat $START_DIR/*_L001_R1_001.fastq $START_DIR/*_L002_R1_001.fastq $START_DIR/*_L003_R1_001.fastq $START_DIR/*_L004_R1_001.fastq > $START_DIR/$OUT_NAME-_1_temp.fq
cat $START_DIR/*_L001_R2_001.fastq $START_DIR/*_L002_R2_001.fastq $START_DIR/*_L003_R2_001.fastq $START_DIR/*_L004_R2_001.fastq > $START_DIR/$OUT_NAME-_2_temp.fq

# (NOVOGENE) use these lines if in a 2 file format (1 for read 1, 1 for read 2) Aka from Novogene on Highseq
#cp $START_DIR/*_1.fq $START_DIR/$OUT_NAME-_1_temp.fq
#cp $START_DIR/*_2.fq $START_DIR/$OUT_NAME-_2_temp.fq

#Trim Files to desired length (25 Bp R1 is all thats needed. R2 76 mapping bases as historical)
perl /home/achialastri/perlscripts/MakeFastqShorter.pl $START_DIR/$OUT_NAME-_1_temp.fq $START_DIR/$FASTQ_R1.all 30
perl /home/achialastri/perlscripts/MakeFastqShorter.pl $START_DIR/$OUT_NAME-_2_temp.fq $START_DIR/$FASTQ_R2.all 46

#Run PC Oligo Fastq to final text files script
perl $PERL_DIR/PC_Barcode_Detection_WithMultiplexTags.pl --FASTQ_R1 $START_DIR/$FASTQ_R1.all --FASTQ_R2 $START_DIR/$FASTQ_R2.all --CELSEQ_BC $CEL_BARCODES --PC_BC $PERL_DIR/$PC_Barcodes --NumUMICelSeq $CelSeq_UMI_Length --NumUMIPC $PC_UMI_Length --OutName $START_DIR/$OUT_NAME --MultiPlex_BC $PERL_DIR/Multiplex_6p_Barcodes_3Hamming.txt --NumUMIMP $MP_UMI_Length

# Remove other unimportant files
rm $START_DIR/$FASTQ_R1.all
rm $START_DIR/$FASTQ_R2.all
rm $START_DIR/$OUT_NAME-_1_temp.fq
rm $START_DIR/$OUT_NAME-_2_temp.fq


#copy all important files to a single location for that run
TRANSCRIPTS="_PC_UMIs"
COUNTS="_PC_Raw"
MP_Raw="_MP_Raw"
MP_UMI="_MP_UMIs"

mkdir $PAST_DIR/$RUN_NAME$TRANSCRIPTS
cp $START_DIR/$OUT_NAME-PC_UMIs.txt $PAST_DIR/$RUN_NAME$TRANSCRIPTS

mkdir $PAST_DIR/$RUN_NAME$COUNTS
cp $START_DIR/$OUT_NAME-PC_Raw.txt $PAST_DIR/$RUN_NAME$COUNTS

mkdir $PAST_DIR/$RUN_NAME$MP_Raw
cp $START_DIR/$OUT_NAME-MP_Raw.txt $PAST_DIR/$RUN_NAME$MP_Raw

mkdir $PAST_DIR/$RUN_NAME$MP_UMI
cp $START_DIR/$OUT_NAME-MP_UMIs.txt $PAST_DIR/$RUN_NAME$MP_UMI
