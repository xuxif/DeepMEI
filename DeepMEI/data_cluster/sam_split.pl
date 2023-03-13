#12      117526485       0       HG002   ALU
open FP,"$ARGV[0]";
@record_raw=<FP>;
close FP;
#HISEQ1:19:H8VDAADXX:2:1103:19682:15966  147     12      117526289       60      148M    =       117525989       -448    CAGCCTGGCCAACATGGCAAAACCCCGTCTCT
$r=0;
@record=split(/\t/,$record_raw[$r]);
open OUT,">HG002_$record[0]:".($record[1]-50)."-".($record[1]+50).".sam";
while(<STDIN>)
{
	@read=split(/\t/,$_);
	if("$read[2]" eq "$record[0]" and abs($read[3]-$record[1])<2000)
	{
		print OUT join("\t",@read);
	}
	else
	{
		$r++;
		@record=split(/\t/,$record_raw[$r]);
		open OUT,">HG002_$record[0]:".($record[1]-50)."-".($record[1]+50).".sam";
	}
}
close OUT;
