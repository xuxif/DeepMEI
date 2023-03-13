%RE;
@name=`cat tmp.bed|cut -f3`;
chomp(@name);
foreach $i (@name)
{
	$RE{"$i"}=1;
}
open FP,"$ARGV[0]";
while(<FP>)
{
	chomp();
	$tmp_name=$_;
	$tmp_name=~s/\t.*//;
	if($RE{"$tmp_name"} eq 1)
	{
		print "$_\n";
#		$RE{"$tmp_name"}=0;
	}
}
