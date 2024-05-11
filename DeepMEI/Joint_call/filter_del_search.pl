
@seq_contex=<STDIN>;
chomp(@seq_contex);
$seq_line=join("#",@seq_contex);
chomp();
open FP,"$ARGV[0]";
$num=0;
while(<FP>)
{
	if($_=~/^@/) {next;};
	chomp();
	@F=split(/\t/,$_);
	if($F[5]=~/S$/)
	{
#		print "$F[5]\t";
		$F[5]=~s/.*?(\d+)S$/\1/;
		$seq=substr($F[9],-1*$F[5]);
#		print "$F[5]\t$seq\n";
		if($F[5]>10 and $seq_line=~/$seq/)
		{
			$num++;
		}
	}
	else
	{
		$F[5]=~s/^(\d+)S.*/\1/;
		$seq=substr($F[9],0,$F[5]);
		if($F[5]>10 and $seq_line=~/$seq/)
		{
			$num++;
		}
	}
}
print "$num\n";
