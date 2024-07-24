use List::Util qw/sum/;
use POSIX;
%OUT;
%OUTIN;
while(<STDIN>)
{
#1       149228  149328  ALU,ALU 1
#1       149467  149567  ALU,ALU	2

	chomp();
	@F=split(/\t/,$_);
	$OUT{"$F[0]\t$F[1]"}=$F[4];
	$OUTIN{"$F[0]\t$F[1]"}=$F[3];
}
close FP;

foreach $i (keys %OUT)
{
	#$count=$OUT{"$i"}+$CN{"$i"};
	$count=$OUT{"$i"};
	if($count>1)
	{
		my %te;
		@tes=split(/,/,$OUTIN{"$i"});
		$te{"ALU"}=0;
		$te{"LINE1"}=0;
		$te{"SVA"}=0;
		foreach $i (@tes) 
		{ 
			$te{$i}++;
		}
		if($te{"ALU"}>=$te{"LINE1"} and $te{"ALU"}>=$te{"SVA"} )
		{
			print "$i\tALU\n";
		}
		elsif($te{"ALU"}<=$te{"LINE1"} and $te{"LINE1"}>=$te{"SVA"} )
		{
			print "$i\tLINE1\n";
		}
		else 
		{
			print "$i\tSVA\n";
		}
	}
}
