# Feature type	Chromosome/scaffold name	Start (bp)
$filename = $ARGV[0];
open($file, '<', $filename);

$outputfilenameCT      = substr($filename,0,-4)."_TopQuartile.txt";
open($outputCT, '>', $outputfilenameCT);

while ($line = <$file>)
{
@x = split(" ",$line);

#$x[5] >= 0.072 
#$x[5] < 0.072 & $x[5] >= 0.039
#$x[5] < 0.039 & $x[5] > 0.031
#$x[5] <= 0.031

	if ($x[5] >= 0.072) {
	print $outputCT ("$x[0] $x[1] $x[2] $x[3] $x[4] $x[5]\n");
	}

}
