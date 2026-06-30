#!/bin/bash                                                                                                                                                  
                                       
#PBS -l nodes=1:ppn=6                                                                                                                                        
                                     
#PBS -l walltime=2:00:00:00


#User Input    Make OUT_NAME the same as the fastq files before _L00#_R#_001.fastq
START_DIR="/home/achialastri/Plate11-H9-5mC-Adpt-Optimization/ACSD082718-P11-H9-5mC-1000nM"




#Standard Usage, no input required
BARCODES="/home/achialastri/perlscripts/MspJI/mspj1"
GENOME="/home/sdey/genomes/human_gene_models/hg19.fa"
PERL_DIR="/home/achialastri/perlscripts/MspJI"
ABASI_BC="/home/achialastri/perlscripts/5hmC/aba_barcodes.csv"
CEL_BARCODES="/home/achialastri/perlscripts/mRNA_Mapping/cel-seq_barcodes.csv"

#Do not change
PAST_DIR=${START_DIR%/*}
OUT_NAME=${START_DIR##*/}-MSPJI
RUN_NAME=${PAST_DIR##*/}
R1="_L001-4_R1_001.fastq"
R2="_L001-4_R2_001.fastq"
FASTQ_R1=$OUT_NAME$R1
FASTQ_R2=$OUT_NAME$R2

intermediateR1=${FASTQ_R1%??????}
intermediateR2=${FASTQ_R2%??????}
OUT_NAME_R1=$intermediateR1-MSPJI
OUT_NAME_R2=$intermediateR2-MSPJI

#Cat Fastq Files
cat $START_DIR/*L001_R1* $START_DIR/*L002_R1* $START_DIR/*L003_R1* $START_DIR/*L004_R1* > $START_DIR/$FASTQ_R1


#Pull out only MspJI Lines (Set --scMATseq == 1 if done with RNA, set --scTHseq == 1 if done with ABASI
perl $PERL_DIR/ExtractingMSPJIReads_UserInput_NoBarcodeCollisions.pl --FASTQ_R1 $START_DIR/$FASTQ_R1 --CELSEQ_BC $CEL_BARCODES --scMATseq 1 --MSPJI_BC $BARCODES --scTHseq 0 --ABASI_BC $ABASI_BC

#Mapping
/home/cwangsanuwat/bwa/bwa-0.7.15/bwa aln -q 0 -n 0.04 -k 2 -l 200 -t 6 -B 11 $GENOME $START_DIR/$OUT_NAME_R1.fastq > $START_DIR/$OUT_NAME_R1.sai


/home/cwangsanuwat/bwa/bwa-0.7.15/bwa samse -n 100 $GENOME $START_DIR/$OUT_NAME_R1.sai $START_DIR/$OUT_NAME_R1.fastq > $START_DIR/$OUT_NAME-se.sam

#Analysis prior to PCR can comment out if not interested in pre PCR information
perl $PERL_DIR/process_scmspji_MCviPIwithQC.pl $GENOME $START_DIR/$OUT_NAME-se.sam $BARCODES
/home/cwangsanuwat/src/samtools/samtools flagstat $START_DIR/$OUT_NAME-se.sam > $START_DIR/$OUT_NAME-se.flagstat


#Same as Bam File
/home/cwangsanuwat/src/samtools/samtools view -bS $START_DIR/$OUT_NAME-se.sam > $START_DIR/$OUT_NAME-se.bam

#Call Type of Methyaltion
perl $PERL_DIR/Make_Final_Version_MCvipI.pl $START_DIR/$OUT_NAME-se.faba



#Remove excess files created
rm $START_DIR/$OUT_NAME-se.sam

rm $START_DIR/$OUT_NAME_R1.sai

rm $START_DIR/$OUT_NAME_R1.fastq 
rm $START_DIR/$OUT_NAME_R2.fastq 
 
rm $START_DIR/$FASTQ_R1
rm $START_DIR/$FASTQ_R2

#copy all important files to a single location for that run
FABA="_FABA"
FLAGSTAT="_FLAGSTAT"
RABA="_RABA"
ACCESS="_FABA_ACCESS"
CG="_FABA_WCG"
CA="_FABA_CA"
CT="_FABA_CT"
CC="_FABA_CC"
AMBIG="_FABA_AMBIGUOUS"


mkdir $PAST_DIR/$RUN_NAME$FABA
cp $START_DIR/$OUT_NAME-se.faba $PAST_DIR/$RUN_NAME$FABA

mkdir $PAST_DIR/$RUN_NAME$FLAGSTAT
cp $START_DIR/$OUT_NAME-se.flagstat $PAST_DIR/$RUN_NAME$FLAGSTAT

mkdir $PAST_DIR/$RUN_NAME$RABA
cp $START_DIR/$OUT_NAME-se.raba $PAST_DIR/$RUN_NAME$RABA

mkdir $PAST_DIR/$RUN_NAME$ACCESS
cp $START_DIR/$OUT_NAME-se_Full5mC_Accessibilty_GC_Final_Rmdup.txt $PAST_DIR/$RUN_NAME$ACCESS

mkdir $PAST_DIR/$RUN_NAME$CG
cp $START_DIR/$OUT_NAME-se_Full5mC_Endogenous_WCG_Final_Rmdup.txt $PAST_DIR/$RUN_NAME$CG

mkdir $PAST_DIR/$RUN_NAME$CA
cp $START_DIR/$OUT_NAME-se_Full5mC_Endogenous_CA_Final_Rmdup.txt $PAST_DIR/$RUN_NAME$CA

mkdir $PAST_DIR/$RUN_NAME$CT
cp $START_DIR/$OUT_NAME-se_Full5mC_Endogenous_CT_Final_Rmdup.txt $PAST_DIR/$RUN_NAME$CT

mkdir $PAST_DIR/$RUN_NAME$CC
cp $START_DIR/$OUT_NAME-se_Full5mC_Endogenous_CC_Final_Rmdup.txt $PAST_DIR/$RUN_NAME$CC

mkdir $PAST_DIR/$RUN_NAME$AMBIG
cp $START_DIR/$OUT_NAME-se_Full5mC_Ambiguous_GCG_Final_Rmdup.txt $PAST_DIR/$RUN_NAME$AMBIG
cp $START_DIR/$OUT_NAME-se_Full5mC_Ambiguous_CCG_Final_Rmdup.txt $PAST_DIR/$RUN_NAME$AMBIG
