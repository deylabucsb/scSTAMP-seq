# This script is designed to take sequencing results from the photocleavable linker design and assign a cell and detection amount according to the UMIs seen on both the PC Oligo and on the Cel-seq2 like primers
# A secondary output file gives back all results seen regardless of UMI
# A similar output is created for Multiplexing Oligos. The barcode indicating a multiplex oligo is hardcoded in. - AJC, 10/12/21
# Created 6/4/2021 by Alex Chialastri. Modified 10/12/21 by AJC, to also account for Multiplexing Oligos containing PCR handle 1.

use warnings;
use Getopt::Long;
use List::MoreUtils qw(uniq);

GetOptions (   
	    "FASTQ_R1=s" => \$FASTQ_R1,
		"FASTQ_R2=s" => \$FASTQ_R2,
        "CELSEQ_BC=s" => \$BARCODES,
		"PC_BC=s" => \$PC_BARCODES,
		"NumUMICelSeq=i" => \$NumUMICelSeq,
		"NumUMIPC=i" => \$NumUMIPC,
		"OutName=s" => \$OutName,
		"MultiPlex_BC=s" => \$MP_BARCODES,
		"NumUMIMP=i" => \$NumUMIMP
    )   
 or die("Error in command line arguments\n usage: All Arguments in this order are needed: --FASTQ_R1 is the input name of your input R1 fastq file && --FASTQ_R2 is the input name of your input R2 fastq file && --CELSEQ_BC is a .csv file of Cel-seq cell specific barcodes RNA && --NumUMICelSeq is how many random bases are used in the Celseq2 primer && --NumUMIPC is how many random bases are used in the PC Oligo && --PC_BC is a text file with 3 columns the second being the Positional barcode and the third being the positional primer && --OutName is the naming convention for the output files && --MultiPlex_BC is a text file containing 6 bp barcodes for the multiplex oligo && --NumUMIMP is how many random bases are used in the Multiplex Oligo)\n");
 

# Read in CEL SEQ Barcodes
open($BC, "$BARCODES"); #opens celseq barcode file
	$i = 0;
	while ($BCline = <$BC>)
	{
		chomp $BCline;   
		@BarC = split("\t",$BCline);
		$BarCode[$i] = $BarC[1];
		$i++;     
	}
$CelSeqBC_Length = length($BarCode[0]);

	
# Read in PC Barcodes and Positional Primer
open($PC, "$PC_BARCODES"); #opens PC barcode file
	$i = 0;
	while ($PCline = <$PC>)
	{
		chomp $PCline;   
		@PC_BarC = split("\t",$PCline);
		$PC_BarCode_Label[$i] = $PC_BarC[0];
		$PC_BarCode[$i] = $PC_BarC[1];
		$PC_Primer[$i] = $PC_BarC[2];
		$i++;
	}
$PC_BC_Length = length($PC_BarCode[0]);
$PC_Primer_Length = length($PC_Primer[0]);

# Read in Multiplex Barcodes
open($MP_BC, "$MP_BARCODES"); #opens celseq barcode file
	$i = 0;
	while ($MP_BCline = <$MP_BC>)
	{
		chomp $MP_BCline;   
		@MP_BarC = split("\t",$MP_BCline);
		$MP_BarCode[$i] = $MP_BarC[0];
		$i++;     
	}
$MP_BC_Length = length($MP_BarCode[0]);


#Open the two fastq files
open($fastafileR1, "$FASTQ_R1");
open($fastafileR2, "$FASTQ_R2");

#Create our two output files for PC Oligo Detection
my $outputRawName = $OutName."-PC_Raw.txt";
open($outputRaw, '>', "$outputRawName");
my $outputUMIName = $OutName."-PC_UMIs.txt";
open($outputUMI, '>', "$outputUMIName");

#Create our two output files for Multiplex Oligo Detection
my $MP_outputRawName = $OutName."-MP_Raw.txt";
open($MP_outputRaw, '>', "$MP_outputRawName");
my $MP_outputUMIName = $OutName."-MP_UMIs.txt";
open($MP_outputUMI, '>', "$MP_outputUMIName");

#Loop through fastq files 4 lines at a time
my %PC_Hash_UMI;

#Setup the UMI PC Hash
for my $i (0 .. $#BarCode) {
	$Cell_BC_Temp = $BarCode[$i];
		for my $j (0 .. $#PC_BarCode) {	
			$PC_BC_Temp = $PC_BarCode[$j];
			$Holder = "a";
			push @{ $PC_Hash_UMI{$Cell_BC_Temp}{$PC_BC_Temp} }, $Holder; #Gives each set of Cells and PC Barcodes an inital key value of a, incase one is not detected, will subtract out a value of 1 at the end to account for the 1 that we added here.
		}
	}


#Setup the UMI MP Hash
my %MP_Hash_UMI;
$MP_MoleculeCode = "ACTG"; #This is the code that flags to us that this is not an inner PC barcode, but instead is from the Multiplex oligo. It was designed with a Hamming distance of >3 from the inner PC barcdoe
for my $i (0 .. $#BarCode) {
	$Cell_BC_Temp = $BarCode[$i];
		for my $j (0 .. $#MP_BarCode) {	
			$MP_BC_Temp = $MP_BarCode[$j];
			$Holder = "a";
			push @{ $MP_Hash_UMI{$Cell_BC_Temp}{$MP_BC_Temp} }, $Holder; #Gives each set of Cells and MP Barcodes an inital key value of a, incase one is not detected, will subtract out a value of 1 at the end to account for the 1 that we added here.
		}
	}
	

while ($header = <$fastafileR1>)
{
    $sequence = <$fastafileR1>;
    $plus     = <$fastafileR1>;
    $qual     = <$fastafileR1>;

	$headerR2   = <$fastafileR2>;
    $sequenceR2 = <$fastafileR2>;
    $plusR2     = <$fastafileR2>;
    $qualR2     = <$fastafileR2>;

	#Get sectional information from each peice of the reads
    $RNA_BC_Section = substr($sequence,0,$CelSeqBC_Length); #Get CELSEQ Barcode Position sequence
	$RNA_UMI_Section = substr($sequence,$CelSeqBC_Length,$NumUMICelSeq); #Get CELSEQ UMI Position sequence
    
	$PC_Primer_Section = substr($sequenceR2,0,$PC_Primer_Length); #Get PC Primer Position sequence
	$PC_BC_Section = substr($sequenceR2,$PC_Primer_Length,$PC_BC_Length); #Get PC Barcode Position sequence (Is also the MP Flag Position)
	$PC_UMI_Section =  substr($sequenceR2,$PC_Primer_Length+$PC_BC_Length,$NumUMIPC); #Get PC UMI Position sequence
	
	$MP_BC_Section = substr($sequenceR2,$PC_Primer_Length+$PC_BC_Length,$MP_BC_Length);  #Get MP Barcode Position sequence
	$MP_UMI_Section = substr($sequenceR2,$PC_Primer_Length+$PC_BC_Length+$MP_BC_Length,$NumUMIMP); #Get the MP UMI Position
	
	### Detection Section ###
	
	$Match_CelSeqBC_Flag = 0; #Check if CELSEQ Barcode is a exact match, Change Flag to True
	for my $i (0 .. $#BarCode) {
		if ($BarCode[$i] eq $RNA_BC_Section){
			$Match_CelSeqBC_Flag = 1;
			$Cell_BC_Temp = $BarCode[$i];
			}
	}
	
	$Match_PCBC_Flag = 0; #Check if PC Barcode is a exact match, and if the PCR Primer is only off by a max of 1 base. If so Change Flag to True
	for my $i (0 .. $#PC_BarCode) {
		$MatchingPCPrimerAmount = 0;
		if ($PC_BarCode[$i] eq $PC_BC_Section){
			for ($k=0;$k<$PC_Primer_Length;$k++){
			$Temp_PC_Primer_SequencedBase = substr($PC_Primer_Section,$k,1);
			$Temp_PC_Primer_RightBase = substr($PC_Primer[$i],$k,1);
		    if ($Temp_PC_Primer_SequencedBase eq $Temp_PC_Primer_RightBase){$MatchingPCPrimerAmount++;}
			}
		}
		
		if ($MatchingPCPrimerAmount >= $PC_Primer_Length-1){
			$Match_PCBC_Flag = 1;
			$PC_BC_Temp = $PC_BarCode[$i];
			}
	}
	
	
	$Match_MPBC_Flag = 0; #Check if MP Barcode is a exact match, and if the PCR Primer is only off by a max of 1 base. If so Change Flag to True
	if ($PC_BC_Section eq $MP_MoleculeCode){ #Only do this section if you see the Multiplex codex, "ACTG"
		for my $i (0 .. $#MP_BarCode) {
			$MatchingPCPrimerAmount = 0;
			if ($MP_BarCode[$i] eq $MP_BC_Section){
				for ($k=0;$k<$PC_Primer_Length;$k++){
				$Temp_PC_Primer_SequencedBase = substr($PC_Primer_Section,$k,1);
				$Temp_PC_Primer_RightBase = substr($PC_Primer[0],$k,1); #Has to match the inner PCR Handle so hardcoded here.
				if ($Temp_PC_Primer_SequencedBase eq $Temp_PC_Primer_RightBase){$MatchingPCPrimerAmount++;}
				}
			}
			
			if ($MatchingPCPrimerAmount >= $PC_Primer_Length-1){
				$Match_MPBC_Flag = 1;
				$MP_BC_Temp = $MP_BarCode[$i];
				}
		}
	}
	
	
	#If both Celseq Barcode and PC barcode section change their flags to True then store the UMI found into the correct cell and PC barcode section of the hash. PC and MP Flag should never be on at the same time but I will enforce that anyway.
	if ($Match_CelSeqBC_Flag == 1 && $Match_PCBC_Flag == 1 && $Match_MPBC_Flag == 0){
		$Total_UMI = $RNA_UMI_Section.$PC_UMI_Section;
		push @{ $PC_Hash_UMI{$Cell_BC_Temp}{$PC_BC_Temp} }, $Total_UMI; #Making an array for each cell PC combo that hold each total UMI found. This is the raw number found but if we make this list unique it is the number of UMIs found.
	}

	#If both Celseq Barcode and MP barcode section change their flags to True then store the UMI found into the correct cell and PC barcode section of the hash. PC and MP Flag should never be on at the same time but I will enforce that anyway.
	if ($Match_CelSeqBC_Flag == 1 && $Match_PCBC_Flag == 0 && $Match_MPBC_Flag == 1){
		$Total_UMI = $RNA_UMI_Section.$MP_UMI_Section;
		push @{ $MP_Hash_UMI{$Cell_BC_Temp}{$MP_BC_Temp} }, $Total_UMI; #Making an array for each cell PC combo that hold each total UMI found. This is the raw number found but if we make this list unique it is the number of UMIs found.
	}
}


#Print the Raw and UMI PC Hash
for my $i (0 .. $#BarCode) {
	$Cell_BC_Temp = $BarCode[$i];
	$Cell_Num_Temp = $i+1;
	print $outputUMI ("$Cell_Num_Temp ");
	print $outputRaw ("$Cell_Num_Temp ");
		for my $j (0 .. $#PC_BarCode) {	
			$PC_BC_Temp = $PC_BarCode[$j];
			$RawHere = @{ $PC_Hash_UMI{$Cell_BC_Temp}{$PC_BC_Temp} }; #Counts how many enteries are in the array of this Cell by PC_BC key
			$Raw_Print = $RawHere - 1; #Each PC_Hash_UMI key was initialized with one entry in an array sub subtract that one entry
			$RawToUnique = uniq(@{ $PC_Hash_UMI{$Cell_BC_Temp}{$PC_BC_Temp} }); #Make the array held in the key unique (aka remove douplicate UMIs), then count.
			$UMI_Print = $RawToUnique - 1; #Each PC_Hash_UMI key was initialized with one entry in an array sub subtract that one entry
			print $outputUMI ("$PC_BarCode_Label[$j] $UMI_Print ");
			print $outputRaw ("$PC_BarCode_Label[$j] $Raw_Print ");
		}
	print $outputUMI ("\n");
	print $outputRaw ("\n");
}


#Print the Raw and UMI MP Hash
for my $j (0 .. $#MP_BarCode) {
	$MP_Num = $j+1;
	$MP_Name = "MP_".$MP_Num;
	print $MP_outputUMI ("$MP_Name ");
	print $MP_outputRaw ("$MP_Name ");
	}
	print $MP_outputUMI ("\n");
	print $MP_outputRaw ("\n");
for my $i (0 .. $#BarCode) {
	$Cell_BC_Temp = $BarCode[$i];
	$Cell_Num_Temp = $i+1;
	#print $MP_outputUMI ("$Cell_Num_Temp ");
	#print $MP_outputRaw ("$Cell_Num_Temp ");
		for my $j (0 .. $#MP_BarCode) {	
			$MP_BC_Temp = $MP_BarCode[$j];
			$RawHere = @{ $MP_Hash_UMI{$Cell_BC_Temp}{$MP_BC_Temp} }; #Counts how many enteries are in the array of this Cell by MP_BC key
			$Raw_Print = $RawHere - 1; #Each MP_Hash_UMI key was initialized with one entry in an array sub subtract that one entry
			$RawToUnique = uniq(@{ $MP_Hash_UMI{$Cell_BC_Temp}{$MP_BC_Temp} }); #Make the array held in the key unique (aka remove douplicate UMIs), then count.
			$UMI_Print = $RawToUnique - 1; #Each MP_Hash_UMI key was initialized with one entry in an array sub subtract that one entry
			print $MP_outputUMI ("$UMI_Print ");
			print $MP_outputRaw ("$Raw_Print ");
		}
	print $MP_outputUMI ("\n");
	print $MP_outputRaw ("\n");
}
