use List::Util qw/sum/;
use POSIX; 
$last_chr=100;
$last_pos=1;
$distence=10;
$hg_pos_mean=0;
$record='';
$num=0;
@cluster_pos='';
#1       320115 115M33S ALU     27      12      12S21M  21 132
#
while(<STDIN>)
{
	chomp();
	@line=split(/\t/,$_);
	$new_chr=$line[0];	
	$new_pos=$line[1];	
	push @cluster_pos,$new_pos;
	if($last_pos ==1)
	{
		$last_chr=$new_chr;
		$last_pos=$new_pos;
		next;
	}

#	if($last_chr eq $new_chr and abs($new_pos-$last_pos) < $distence )
	if($last_chr eq $new_chr and abs($new_pos-$last_pos) < $distence )
	{
		$num++;
	}

	if($last_chr ne $new_chr or abs($new_pos-$last_pos) > $distence or eof(STDIN))
	{
		if($last_chr ne $new_chr or abs($new_pos-$last_pos) > $distence )
		{
			pop @cluster_pos;
		}
		if($num>=0) 
		{
			$hg_pos_mean=ceil(sum(@cluster_pos)/@cluster_pos);
#			print "@cluster_pos\n";
			$chr=$last_chr;
			$start=$hg_pos_mean-100;
			$end=$hg_pos_mean+100;
#			print "$chr:$start-$end\t$num\n";
#			print "$chr\t$start\t$end\t$line[8]\n";
#			print "$chr\t$hg_pos_mean\t0\tHG002\t$num\n";
			print "$chr\t$hg_pos_mean\t0\tHG002\n";
#			print "$chr\t$hg_pos_mean\t$clip_part\n";
		}
		@cluster_pos=();
		push @cluster_pos,$new_pos;
		$num=0;
	}
	$last_chr=$new_chr;
	$last_pos=$new_pos;
}
