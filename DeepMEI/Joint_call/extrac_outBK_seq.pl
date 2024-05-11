$left=$ARGV[0];
$right=$ARGV[1];
$pad=10;
open FP,$ARGV[2];
while(<FP>)
{
	if($_=~/^@/)
	{next;}
	chomp();
	@F=split(/\t/,$_);
	if($F[5]=~/S/)
	{
		next;
	}
	if($F[3]+150<=$left-$pad)
	{
		print "left:$F[9]\n";
	}
	elsif($F[3]<=$left)
	{
		$seq=substr($F[9],0,$left-$F[3]+10);
		print "left:$seq\n";
	}

	if($right<$F[3])
	{
		print "right:$F[9]\n";
	} 
	elsif($right-$F[3]>=0)
	{
		if($right-$F[3]<10)
		{
			$seq=substr($F[9],0);
		}
		else
		{
			$seq=substr($F[9],$right-$F[3]-10);
		}

		print "right:$seq\n";
	} 
}
