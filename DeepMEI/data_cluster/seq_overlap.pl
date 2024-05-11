chomp($seq_clip=<STDIN>);
chomp($seq_ref=<STDIN>);
$seq_clip=~s/(.)/\1\t/g;
$seq_clip=~s/\t$//g;
@seq_clip=split(/\t/,$seq_clip);

$seq_ref=~s/(.)/\1\t/g;
$seq_ref=~s/\t$//g;
@seq_ref=split(/\t/,$seq_ref);
$match=0;
$max_match=0;
$match_one_max=0;
for($i=0;$i<=$#seq_clip-5;$i++)
{
	for($j=$i,$k=0;$j<=$#seq_ref;$j++,$k++)
	{
		if($seq_ref[$j] eq $seq_clip[$k])
		{
			$match++;
		}
		elsif($match_one_max<$match)
		{
			$match_one_max=$match;
			$match=0;
		}
		else
		{
			$match=0;
		}
	}
	if($match_one_max>$max_match)
	{
		$max_match=$match_one_max;
	}
	$match_one_max=0;
	$match=0;
}
print "$max_match\n";
