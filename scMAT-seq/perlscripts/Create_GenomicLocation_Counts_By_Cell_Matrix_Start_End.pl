#!/usr/bin/perl
 
#use strict;
use warnings;

$DistancePriorToTSS = $ARGV[2]; #This max distance downstream the TES (3rd argument in command line)
$DistanceAfterToTSS = $ARGV[3]; #This max distance upstream the TSS (4th argument in command line)

$filename = $ARGV[0]; #This is the file giving gene TSS and TSE locations in a vector -> LongestGenesLocations.txt (Chr Start End Strand Name)
open($GeneLocations, '<', $filename);

$filename = $ARGV[1]; #This is the txt file for GCH or WCG
open($AllGCLoci, '<', $filename);

$GenomicElement = $ARGV[4]; #This argument is simply for naming what genomic element you are looking near

$outputfilename      = substr($filename,0,-4)."Matrix_ByCell_in_$GenomicElement-Start_to_End-with-$DistanceAfterToTSS-Upstream_$DistancePriorToTSS-Downstream.txt"; #ouput file name based on distances used and input file used.
open($output, '>', $outputfilename);

my $startTime = time; #To time the perl script
%HighValHash = (); #Will contain Lower value to searh for each gene with chromosomes as keys
%LowValHash = (); #Will contain Higher value to searh for each gene with chromosomes as keys

%StartHash = (); #Will contain start values with chromosomes as keys
%EndHash = (); #Will contain end values (highest number) with chromosomes as keys.  Not used but would need to be impleneted to create similar table but for within gene body, not just TSS.
%NameHash = (); #Will contain gene names with chromosomes as keys (corresponds with StartHash)
@NameArray = (); #List of all gene names

#### This section reads in each gene and TSS location and creates the boundries for the search along with an array of each gene used.
$GeneNumber = 0;
while ($line = <$GeneLocations>) { #read in all of the gene loctaions into a vector
	@x = split(" ",$line);
	$start = $x[1];
	$End = $x[2];
		if ($start > $End){ #this if statment assures the high number is greater than the low number accounting for if the TSS is on + or - strand
		$High = $start+$DistanceAfterToTSS;
		$Low = $End-$DistancePriorToTSS;
		}
		else {
		$High = $End+$DistancePriorToTSS;
		$Low = $start-$DistanceAfterToTSS;
		}
	$chr = 'chr'.$x[0];
	
	push @{ $StartHash{$chr} }, $start; #push the start position into the end of the hash with the right chromosome as the key
	push @{ $NameHash{$chr} }, $x[4]; #push name of gene into end of hash so it is lined up with its other start, high and low values
	push @{ NameArray }, $x[4]; #array of gene names read in order presented in txt file
	push @{ $HighValHash{$chr} }, $High; #Will always be the higher number, for + strand genes, this is downstream TSS.  Vice versa for - strand genes
	push @{ $LowValHash{$chr} }, $Low; #Will always be the lower number, for + strand genes, this is upstream TSS (true promoter).  Vice versa for - strand genes
	
	$GeneNumber++;	
	
}
#####


##### This section creates a hash within a hash (gene name and cell number) that contains an initialized value of 0 for that gene cell TSS
%GeneByCellHash = ();

for ($i=1;$i<=96;$i++) {
	$iKey = "Cell".$i;
		foreach (@NameArray){
		$GeneByCellHash{$_}{$iKey}=0;
		
		}
	
}
#####



##### This section reads through the file containing the found GC or WCG locations.  If a read is within the TSS it will add it to the GeneByCellHash we created above.

while ($line2 = <$AllGCLoci>) { #read in sequenced locations (cell chr loci strand)
	@x2 = split(" ",$line2);
	#$ChrVal = $x2[1]; $ChrKey = 'chr'.$ChrVal; # Use if does not have chr already written
	$ChrKey = $x2[1]; # Use if do have chr already written
	$CellKey = "Cell".$x2[0];
	$LociVal = $x2[2];
	#$StrandVal = $x2[3]; #strand value not used but could be if want to check strand info at TSS
	
		for($i=0; $i<@{$HighValHash{$ChrKey}}; $i++) { #instead of searching if the loci value is within the TSS of all genes and then checking if it is on the right chromosome.  Since we read in the gene TSS as a hash, we can search only those on the same chromosome as the read.
			if(@{$HighValHash{$ChrKey}}[$i]>$LociVal && $LociVal>@{$LowValHash{$ChrKey}}[$i]){ #checks if the read is within defined TSS
			$gname=@{$NameHash{$ChrKey}}[$i]; #if it is within TSS, grab name of TSS it is within
			$GeneByCellHash{$gname}{$CellKey}++; #add one to that gene for that cell
			}
		} 

}


#### This section prints out the Gene by cell matrix containing the number of times a site was found in the TSS

#creates header
print $output "GeneName\t";
for ($i=1;$i<=96;$i++){
#$iKey = "Cell".$i;
$iKey = $i;
print $output "$iKey\t"
}
print $output "\n";

#for each gene cycle through each cell and print out how many were in that TSS.
#Here I have also removed genes that contained 0 elements in TSS for all cells.
use List::Util 'sum'; #needed to sum a list
foreach (@NameArray) { #cycle through all genes
	$gname=$_;
	$value_count_gname = sum values %{$GeneByCellHash{$gname}}; #sums all counts for that gene
	
		if ($value_count_gname > 0 ) { #if a mark is present in at least 1 cell, print that gene for all cells.
		
			print $output "$gname\t";
				for ($i=1;$i<=96;$i++){
				$iKey = "Cell".$i;
				print $output "$GeneByCellHash{$gname}{$iKey}\t";
				
				}
			print $output "\n";
		}
}

#####

# To print out the script execution time
my $duration = time - $startTime;
print "Execution time: $duration s\n";


