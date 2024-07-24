#A00404:158:HV2LNDSXX:2:2475:14660:4210/1        chr1    10098   10198
my %RE;
my %OUT;
open FP,"$ARGV[0]";
my @record=();
$record_total=-1;
while(<FP>)
{
	chomp();
	@F=split(/[:\-\t]/,$_);
	$pos=int(($F[1]+$F[2])/2);
	push(@F,$pos);
	push @record,[$F[0],$F[1],$F[2],$F[3],$F[4]];
	$record_total++;
}
close FP;
$i=0;
#A00404:158:HV2LNDSXX:3:2445:14714:20901 147     chr1    111142365       60      150M    =       111141832       -683    ACCCTCCCCTTTTTGAGTTTCTTAAAGGTATGGAATGTTATATTTCAGCCTTTCTGTTTCATCAGAGTCTGGAACCACTAAGTGTTTATTAAATACTAGACGAATAATATTTGTTAATTGTTCAGTAAGTACTGATTACCTACTAGCAGG        +???????????+???????????+???????????????????????????????????????????????????????????????????????5?????????????????????????????????????????????????????        MC:Z:150M       PG:Z:MarkDuplicates           MQ:i:60 AS:i:150        XS:i:0  MD:Z:150        NM:i:0  RG:Z:HG01529_TACTCATA-CCTGTGGC_HV2LNDSXX_L003
my %READ;
$last_pos=0;
$record_i=0;
$read_c=0;
print "$record[$record_i][0]:$record[$record_i][1]-$record[$record_i][2]\t$record[$record_i][3]\n";
$d_chr=$record[$record_i][0];$d_start=$record[$record_i][1];$d_end=$record[$record_i][2];

foreach $read_i (<STDIN>)
{
	$read_c++;
	chomp($read_i);
	@F=split(/\t/,$read_i);
	$cigar=$F[5];
	$map_len=0;
	while($cigar=~/[DM]/)
	{
		$cigar=~s/(\d+)[DM]//;
		$map_len=$map_len+$1;
	}
	$r_chr=$F[2];$r_start=$F[3];$r_end=$F[3]+$map_len;
	#	print abs($record[$record_i][4]-$pos)."\t".abs($record[$record_i+1][4]-$pos)."\n";
	if($F[3] lt $last_pos or $r_chr ne $d_chr or $r_start>$d_end)
	{
		$record_i++;
		$d_chr=$record[$record_i][0];$d_start=$record[$record_i][1];$d_end=$record[$record_i][2];
		print "\n$record[$record_i][0]:$record[$record_i][1]-$record[$record_i][2]\t$record[$record_i][3]\n";
		$read_c=0;
	}
	if($read_c<1000)
	{
		print "$read_i##123##456##";
	}
	$last_pos=$F[3];
}
print "\n";
