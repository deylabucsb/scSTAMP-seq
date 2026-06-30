# Feature type	Chromosome/scaffold name	Start (bp)
$filename = $ARGV[0];
open($file, '<', $filename);

$outputfilenameCT      = substr($filename,0,-4)."_Reformated.txt";
open($outputCT, '>', $outputfilenameCT);

$count=0;
while ($line = <$file>)
{
@x = split("\t",$line);


$MidPoint = ($x[2]+$x[3])/2;
$chr=substr($x[1],3);
if ($x[1] eq chrX) {$chr=23;}
if ($x[1] eq chrY) {$chr=24;}
$SiteName="HeLaHypersenstivitySite".$count;

print $outputCT ("$chr $MidPoint $MidPoint 1 $SiteName $x[7]\n");
$count++;
}
