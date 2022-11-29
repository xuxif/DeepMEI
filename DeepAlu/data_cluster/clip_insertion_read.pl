use List::Util qw/sum/,qw/max/;
sub clip_insertion
{
	(my @read)=@_;
	$cigar_t=$read[5];
	$cigar_t=~s/.*?(\d+)I.*/\1/g;
	if($read[5]=~/I/ and $cigar_t>5)
	{
		$max_cg=&maxCigar($read[5],'I');
		print "read5:$read[5]\tmax_cg:$max_cg\n";
		print "$read[5]\t$read[9]\n";
		$cigar_t=$read[5];
		$cigar_t=~s/${max_cg}I.*//;
		$left_len=&mapLen($cigar_t);

		$cigar_t=$read[5];
		$cigar_t=~s/.*?${max_cg}I//;
		$right_len=&mapLen($cigar_t);
		
		if($left_len>=$right_len)
		{
			$read[5]=~s/(${max_cg})I(.*)/\1S/;
			$cigar_t=$read[5];
			$cigar_t=~s/\d+D//g;
			$cigar_t=~s/[A-Z]/\t/g;
			@cigar_len=split(/\t/,$cigar_t);
			$clip_len_read=sum(@cigar_len);
			$read[9]=substr($read[9],0,$clip_len_read);
			$read[10]=substr($read[10],0,$clip_len_read);
		}	
		else
		{
			$read[5]=~s/(.*?)($max_cg)I/\2S/;
			$cigar_t=$1;
			$cigar_t=~s/(\d+[DM])/\t\1/g;
			$cigar_t=~s/\t(\d+)[DM].*?\t?/\1\t/g;
			$cigar_t=~s/\t$//g;
			@cigar_len=split(/\t/,$cigar_t);
			$clip_len_ref=sum(@cigar_len);
			$read[3]=$read[3]+$clip_len_ref;

			$cigar_t=$read[5];
			$cigar_t=~s/\d+D//g;
			$cigar_t=~s/[A-Z]/\t/g;
			@cigar_len=split(/\t/,$cigar_t);
			$clip_len_read=sum(@cigar_len)*(-1);
			$read[9]=substr($read[9],$clip_len_read,);
			$read[10]=substr($read[10],$clip_len_read,);
		}
		print "$read[5]\t$read[9]\n";
	}
	return @read
}

sub clip_del
{
	(my @read)=@_;
	$cigar_t=$read[5];
	$cigar_t=~s/.*?(\d+)D.*/\1/g;
	if($read[5]=~/D/ and $cigar_t>5)
	{
		$max_cg=&maxCigar($read[5],'D');
		print "0###$read[3]\t$read[5]\t$max_cg\n";
		$cigar_t=$read[5];
		$cigar_t=~s/${max_cg}D.*//;
		$left_len=&mapLen($cigar_t);

		$cigar_t=$read[5];
		$cigar_t=~s/.*?${max_cg}D//;
		$right_len=&mapLen($cigar_t);
		print "left_len:$left_len\tright_len:$right_len\n";	
		if($left_len>=$right_len)
		{
			$read[5]=~s/(${max_cg})D(.*)/\1S/;
			print "$read[5]\n";
			$cigar_t=$read[5];
#			$cigar_t=~s/\d+S$//;
			$cigar_t=~s/\d+D//g;
			$cigar_t=~s/[A-Z]/\t/g;
			@cigar_len=split(/\t/,$cigar_t);
			print "$cigar_t\n";
			$clip_len_read=sum(@cigar_len);
			if($clip_len_read>length($read[9]))
			{
				my $extra_len=$clip_len_read-length($read[9]);
                                $max_cg_t=$max_cg-$extra_len;
                                $read[5]=~s/${max_cg}S$/${max_cg_t}S/;
				$clip_len_read=length($read[9]);

			}
			$read[9]=substr($read[9],0,$clip_len_read);
			$read[10]=substr($read[10],0,$clip_len_read);
			print "1###$clip_len_read\n";
		}	
		else
		{
			$read[5]=~s/(.*?)(${max_cg})D/\2S/;
			$clip_len_ref=$2;
			$cigar_t=$1;
			$cigar_t_2=$1;
			$cigar_t=~s/(\d+[DM])/\t\1/g;
			$cigar_t=~s/\t(\d+)[DM].*?\t?/\1\t/g;
			$cigar_t=~s/\t$//g;
			@cigar_len=split(/\t/,$cigar_t);
			$clip_len_ref+=sum(@cigar_len);
			$read[3]=$read[3]+$clip_len_ref;
			
			$cigar_t=$read[5];
#			$cigar_t=~s/\d+S//;
			$cigar_t=~s/\d+D//g;
			$cigar_t=~s/[A-Z]/\t/g;
			@cigar_len=split(/\t/,$cigar_t);
			$clip_len_read=sum(@cigar_len)*(-1);
			if($clip_len_read*(-1)>length($read[9]))
			{
				my $extra_len=(-1)*$clip_len_read-length($read[9]);
                                $max_cg_t=$max_cg-$extra_len;
                                $read[5]=~s/^${max_cg}S/${max_cg_t}S/;
				print "max_cg:$max_cg\tmax_cg_t:$max_cg_t\tclip_len_read:$clip_len_read\tread[9]:".length($read[9])."\n";
				$clip_len_read=-1*length($read[9]);
			}
			print "1###$clip_len_read\n";
			$read[9]=substr($read[9],$clip_len_read,);
			$read[10]=substr($read[10],$clip_len_read,);
		}
		print "$read[3]\t$read[5]\n";
	}
	return @read
}
sub maxCigar
{
	($cigar,$find_c)=@_;
	print "$cigar\t$find_c\n";
	$cigar=~s/(\d+$find_c)/\1\t/g;
	$cigar=~s/\t$//;
	@cigars=split(/\t/,$cigar);
	for($i=0;$i<@cigars;$i++)
	{
		unless($cigars[$i]=~/$find_c/)
		{
			$cigars[$i]="";
		}
		$cigars[$i]=~s/.*?(\d+)$find_c/\1/;
	}

	return max(@cigars)
}
sub mapLen
{
	($cigar_t)=@_;
	$cigar_t=~s/([A-Z])/\1\t/g;
	@cigar_len=split(/\t/,$cigar_t);
	for($i=0;$i<@cigar_len;$i++)
	{
		unless($cigar_len[$i]=~/M/)
		{
			$cigar_len[$i]=0;	
		}
	}
	return sum(@cigar_len)
}	
while(<STDIN>)
{
	chomp();
	@line=split(/\t/,$_);
	@line=&clip_insertion(@line);
	@line=&clip_del(@line);
	print join("\t",@line)."\n";

}
