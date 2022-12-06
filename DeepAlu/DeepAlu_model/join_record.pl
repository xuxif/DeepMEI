$last_chr=1;
$last_pos=1;
$distence=50;
$record='';
$num=1;
while(<STDIN>)
{
	chomp();
	@line=split(/\t/,$_);
	$new_chr=$line[0];	
	$new_pos=$line[1];	
	if($last_chr ne $new_chr or abs($new_pos-$last_pos) > $distence)
	{
		print "$_\n";
		$last_chr=$new_chr;
		$last_pos=$new_pos;
	}
	else
	{
		#		print "$last_pos\t$new_pos\n";		
	}
}
