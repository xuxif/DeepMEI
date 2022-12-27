use strict;
use warnings;
my @line=();
my @name=();
my $as="";
while(<STDIN>)
{
	@line=split;
	@name=split /\|/, $line[0];
	$as="NA";
#if ($tmp1[2]%4 >=2
	unless($line[2] eq "*" )
	{
		if( $name[2]%4 >=2)
		{
			for(my $i=11;$i<@line;$i++)
			{
				if($line[$i]=~s/AS:i://){$as=$line[$i];}
			}
			if($as eq "NA" || $as<20){next;}
			print "$name[3]\t$name[4]\t$name[6]\t$line[2]\t$line[3]\t$line[4]\t$line[5]\t$as\n";
		}
	}
}
