# Make a txt file that is compatible to earlier versions from AvO's *.faba file
# AvO's script has a bug in the exact nucleotide positions. This is corrected in this script


$filename = $ARGV[0];
open($file, '<', $filename);

$outputfilenameCG      = substr($filename,0,-5)."_Full5mC_Endogenous_WCG_Final_Rmdup.txt";
open($outputCG, '>', $outputfilenameCG);

$outputfilenameGC      = substr($filename,0,-5)."_Full5mC_Accessibilty_GC_Final_Rmdup.txt";
open($outputGC, '>', $outputfilenameGC);

$outputfilenameGCG      = substr($filename,0,-5)."_Full5mC_Ambiguous_GCG_Final_Rmdup.txt";
open($outputGCG, '>', $outputfilenameGCG);

$outputfilenameCCG      = substr($filename,0,-5)."_Full5mC_Ambiguous_CCG_Final_Rmdup.txt";
open($outputCCG, '>', $outputfilenameCCG);

$outputfilenameCA      = substr($filename,0,-5)."_Full5mC_Endogenous_CA_Final_Rmdup.txt";
open($outputCA, '>', $outputfilenameCA);

$outputfilenameCC      = substr($filename,0,-5)."_Full5mC_Endogenous_CC_Final_Rmdup.txt";
open($outputCC, '>', $outputfilenameCC);

$outputfilenameCT      = substr($filename,0,-5)."_Full5mC_Endogenous_CT_Final_Rmdup.txt";
open($outputCT, '>', $outputfilenameCT);

#open($file, "/hpc/hub_oudenaarden/s.dey/5mC_RNA_Emb_Hyb_CAST_B6_sc_20161116/HybEmb5mCRNA-Lib1_5mCReads/AlleleSpecificData/FinalManuscript/HybEmb5mCRNA-Lib1_5mCReads_B6_se.allele.faba");
#open($output, ">/hpc/hub_oudenaarden/s.dey/5mC_RNA_Emb_Hyb_CAST_B6_sc_20161116/HybEmb5mCRNA-Lib1_5mCReads/AlleleSpecificData/FinalManuscript/HybEmb5mCRNA-Lib1_5mCReads_B6_se_Full5mC_CG_Final_Rmdup.txt");

while ($line = <$file>)
{
    @x = split("\t",$line);
    $x[0] = $x[0] + 1;
    $x[1] = $x[1] + 1;
    if ($x[3] == 1)
    {
	$x[2] = $x[2] + 1;
    }
    if ($x[3] == -1)
    {
        $x[2] =$x[2] +3;
		
    }
    
	$y = substr($x[8],0,3);
    if ($y eq "GCG")
    {
	#print $outputCG ("$x[0] $x[1] $x[2] $x[4] $x[3] $x[5]\n");
	print $outputGCG ("$x[0] chr$x[1] $x[2] $x[3]\n");
    }
    elsif ($y eq "CCG")
    {
	print $outputCCG ("$x[0] chr$x[1] $x[2] $x[3]\n");
    }
	elsif ($y =~ m/.CG/)
	{
		print $outputCG ("$x[0] chr$x[1] $x[2] $x[3]\n");

	}
	elsif ($y =~ m/GC./)
	{
		print $outputGC ("$x[0] chr$x[1] $x[2] $x[3]\n");

	}
	elsif ($y =~ m/.CA/)
	{
		print $outputCA ("$x[0] chr$x[1] $x[2] $x[3]\n");

	}
	elsif ($y =~ m/.CT/)
	{
		print $outputCT ("$x[0] chr$x[1] $x[2] $x[3]\n");

	}
	elsif ($y =~ m/.CC/)
	{
		print $outputCC ("$x[0] chr$x[1] $x[2] $x[3]\n");

	}
	}
	
	
	
