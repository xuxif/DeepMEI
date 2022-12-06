$last_direction=1;
$last_pos=1;
$distence=10;
$record='';
$num=1;
$count=1;
while(<STDIN>)
{
	chomp();
	@line=split(/\t/,$_);
	$new_direction=$line[0];	
	$new_pos=$line[1];	
	if($last_direction ne $new_direction or abs($new_pos-$last_pos) > $distence)
	{
		print "$_\t$count\n";
		$count=1;
	}
	else
	{
		$count++;
	}
	$last_direction=$new_direction;
	$last_pos=$new_pos;
}
