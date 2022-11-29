use List::Util qw/sum/;
use POSIX; 
#HISEQ1:18:H8VC6ADXX:1:1102:3956:71617   177     2       191250351       60      20S128M 1       49775580 
@left_clip_all;
@right_clip_all;
while(<STDIN>)
{
	chomp();
	@line=split(/\t/,$_);
	unless($line[5]=~/S/ or $line[4]<45)
	{
		next;
	}
	$chr=$line[2];
	$pos=$line[3];
	$line[5]=~s/(\d+)/\1\t/g;
	$line[5]=~s/([A-Z])/\1\t/g;
	$line[5]=~s/\t$//g;
	@hg_cigar=split(/\t/,$line[5]);
	$hg_max_length=0;
	for($i=1;$i<@hg_cigar;$i=$i+2)
	{
		if($hg_cigar[$i]=~/S/ and $hg_cigar[$i-1]>$clip_count)
		{
			$clip_pos=$i;
			$clip_count=$hg_cigar[$i-1];
		}
		elsif($hg_cigar[$i]=~/M/ and $hg_cigar[$i-1]>$map_count)
		{
			$map_pos=$i;
			$map_count=$hg_cigar[$i-1];
		}
	}
	if($map_pos>$clip_pos)
	{
		push @left_clip_all,$pos;
	}
	else
	{
		$clip_len=0;
		for($i=0;$i<$clip_pos-1;$i=$i+2) {$clip_len=$clip_len+$hg_cigar[$i];}
#		print "pos:$pos\t$clip_len\t$line[0]\n";
		$pos=$pos+$clip_len;
		push @right_clip_all,$pos;
	}
	$clip_len=0;
	$map_count=0;
	$map_pos=-1;
	$clip_pos=-1;
	$clip_count=0;
}
$last_i=0;
$distence=15;
$count=1;
$count_max_left=0;
$count_max_left_pos=0;
$j=0;
$count_left=1;
foreach $i (sort { $a <=> $b } @left_clip_all)
{
	if(abs($i-$last_i)<$distence)
	{
		$count++
	}
	else
	{
		$count=1;
		$last_i=$i;
		$count_left++;
	}
	if($count>$count_max_left)
	{
		$count_max_left=$count;
		$count_max_left_pos=$last_i;
	}
	if($j==$#left_clip_all and $count_max_left_pos==0)
	{
		$count_max_left_pos=$i;
		$count_max_left=$count;
	}	
	$j++;
#		print "$i\n";
}
$last_i=0;
$count=1;
$count_max_right=-1;
$j=0;
$count_right=1;
foreach $i (sort { $a <=> $b } @right_clip_all)
{
	if(abs($i-$last_i)<$distence)
	{
		$count++;
	}
	else
	{
		$count=1;
		$last_i=$i;
		$count_right++;
	}
	if($count>$count_max_right)
	{
		$count_max_right=$count;
		$count_max_right_pos=$last_i;
	}
	
	if($j==$#right_clip_all and $count_max_right_pos==0 )
	{
		$count_max_right_pos=$i;
		$count_max_right=$count;
	}	
	$j++;
#		print "$i\n";
}
#print "$count_max_left\t$count_max_right\n";
#print "$count_max_left_pos\t$count_max_right_pos\n";
if($count_max_left_pos-8>$count_max_right_pos and $count_max_left>=5 and $count_max_right>=5)
{
	print "wrong_2clip_position_1\n";
}	
elsif($count_right >6 or $count_right >6)
{
	print "wrong_2clip_position_2\n";
}	
elsif($count_max_left<=7 and  $count_max_right<=7)
{
	print "wrong_2clip_position_3\n";
}	
elsif($count_max_left/$count_max_right<0.13 or $count_max_right/$count_max_left<0.13)
{
	print "wrong_2clip_position_4\n";
}	
else
{
	print "correct_2clip_position\n";
}	
