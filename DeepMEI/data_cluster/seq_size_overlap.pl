$seq_ref=$ARGV[0];
$dir=$ARGV[1];
$match=0;
$max_clip_len=0;
while(<STDIN>)
{
	chomp();
	$seq_read=$_;
	if(length($seq_read)>$max_clip_len)
	{
		$max_clip_len=length($seq_read);
	}
#	if(length($seq_read) <20)
#	{
#		next;
#	}
	if($dir eq 'left')
	{
		$seq_key=substr($seq_ref,-5);
		if($seq_read=~/(.*$seq_key)(.*)/ and $seq_ref=~/$1/)
		{
			$seq_read=~s/.*$seq_key//;
			$match=length($seq_read);
			last;
		}
	}
	elsif($dir eq 'right')
	{
		$seq_key=substr($seq_ref,0,5);
		if($seq_read=~/(.*)($seq_key.*)/ and $seq_ref=~/$2/)
		{
			$seq_read=~s/$seq_key.*//;
			$match=length($seq_read);
			last;
		}
	}
}

print "$match\t$max_clip_len\n";


