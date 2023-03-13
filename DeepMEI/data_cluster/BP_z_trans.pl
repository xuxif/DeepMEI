while(<STDIN>)
{
	chomp();
	@F=split(/\t/,$_);
	$F[0]=~s/:.*//;
	if(($F[1] ==0 and $F[2] ==0) or (abs($F[1]-$F[2])>200 and $F[1]*$F[2]>0))
	{
		next;
	}
	if($F[1] == 0 ) 
	{
		$F[1]=$F[2]-10;
	} 
	if($F[2] == 0 )
	{ 
		$F[2]=$F[1]+10;
	} 
	if($F[1]>$F[2]) 
	{
		print "$F[0]\t$F[2]\t$F[1]\n";
	} 
	else 
	{
		print "$F[0]\t$F[1]\t$F[2]\n";
	}
}
