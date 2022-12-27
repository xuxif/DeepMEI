use List::Util qw/sum/;
$ori_dir="regions";
$dir="alu_sam";
$out_dir="split_softclipped";
if(not -e $out_dir)
{
`mkdir $out_dir`;
}
if(not $ARGV[0])
{
	$discord_file="HG00096_chr10:104057406-104057506_discord.sam";
	$ss_file="HG00096_chr10:104057406-104057506_alu.sam";
	$ori_file="HG00096_chr10:104057406-104057506.bam";
}
else
{
	$discord_file="$ARGV[0]_discord.sam";
	$ss_file="$ARGV[0]_alu.sam";
	$ori_file="$ARGV[0].bam";
}
%ORI;
%SS;
%DISCORD;
@direction_ori,@direction_alu;
@clip_len_ori,@clip_len_alu;
@reads_ori,@reads_alu;
@max_match_ori,@max_match_alu;
$cigar_filter_len=20;
#$out_dir="discord_read";
$bwa="/public/home/swgenetics_4/software/bwa-0.7.15/bwa";
`cat head.sam >${out_dir}/$ARGV[0]_mapRef.sam`;
`cat head.sam >${out_dir}/$ARGV[0]_mapClipR.sam`;
`cat head.sam >${out_dir}/$ARGV[0]_mapClipL.sam`;
open(MPR,">>${out_dir}/$ARGV[0]_mapRef.sam");
open(MPCL,">>${out_dir}/$ARGV[0]_mapClipL.sam");
open(MPCR,">>${out_dir}/$ARGV[0]_mapClipR.sam");
($direction,$clip_len,$reads,$max_match)=&read_clip_check("$ori_dir/$ori_file",0);
@direction_ori=@{$direction};
@clip_len_ori=@{$clip_len};
@reads_ori=@{$reads};
@max_match_ori=@{$max_match};

($direction,$clip_len,$reads,$max_match)=&read_clip_check("$dir/$ss_file",1);
@direction_alu=@{$direction};
@clip_len_alu=@{$clip_len};
@reads_alu=@{$reads};
@max_match_alu=@{$max_match};

for(my $i=0;$i<@clip_len_ori;$i++)
{
	my ($max_match_alu_pad,$max_match_alu_m)=split(/:/,$max_match_alu[$i]);
	my ($max_match_ori_pad,$max_match_ori_m)=split(/:/,$max_match_ori[$i]);

	if($clip_len_alu[$i] >10 and $clip_len_ori[$i]>10 and ($direction_alu[$i]+$direction_ori[$i])==3)
	{
		my @tmp_read_alu=split(/\t/,$reads_alu[$i]);
		my @tmp_read_ori=split(/\t/,$reads_ori[$i]);
		if($max_match_alu_pad<$max_match_ori_pad and ($max_match_alu_pad+$max_match_alu_m) < ($max_match_ori_pad+$max_match_ori_m))
		{
			$tmp_read_ori[6]='chr1';
			$tmp_read_ori[7]=$max_match_alu_pad+$max_match_alu_m;
			print MPCL join("\t",@tmp_read_ori)."\n";
#			print "MPCL\n";
		}
		elsif($max_match_alu_pad>$max_match_ori_pad and ($max_match_alu_pad+$max_match_alu_m) > ($max_match_ori_pad+$max_match_ori_m))
		{
			$tmp_read_ori[6]='chr2';
			$tmp_read_ori[7]=$max_match_alu_pad;
			print MPCR join("\t",@tmp_read_ori)."\n";
#			print "MPCR\n";
		}
		else
		{
#			print "MPR\n";
			print MPR "$reads_ori[$i]\n";
		}
	}
	else
	{
#		print "MPR\n";
		print MPR "$reads_ori[$i]\n";
	}
#	print "read_name:$tmp_read_ori[0],max_match_alu:$max_match_alu[$i],alu_cigar:$tmp_read_alu[5]\n$tmp_read_alu[9]\n$tmp_read_alu[2],$tmp_read_alu[5],\nmax_match_alu_pad:$max_match_alu_pad,max_match_alu_m:$max_match_alu_m,max_match_ori:$max_match_ori[$i],ori_cigar:$tmp_read_ori[5]\n$tmp_read_ori[9]\n$tmp_read_ori[2],$tmp_read_ori[5]\n,\nmax_match_ori_pad:$max_match_ori_pad,max_match_ori_m:$max_match_ori_m\n";
}

close MPCR;
close MPCL;
close MPR;
sub read_clip_check
{
	my $file=$_[0];
	my @FP_ori;
	my $type=$_[1];
	my @max_match,@max_clip;
	my $max_clip;
	if($type==0)
	{
		@FP_ori=`samtools view  $file`;
	}
	else
	{
		@FP_ori=`samtools view -F 256 $file`;
	}
	my @reads,@cigars,@reads,@reads_total;
	my $i=0,$cigar_filter_len,$direction,$clip_len;
	foreach my $each_read (@FP_ori)
	{
		chomp($each_read);
		@reads=split(/\t/,$each_read);
		$reads[5]=~s/(\d+)/#\1#/g;
		$reads[5]=~s/^#//g;
		@cigars=split(/#/,$reads[5]);		
		my $max_m_j=-1,$max_c_j=-1,$j=1,$max_m=0,$max_c=0;
		for(my $j=0;($j+1)<=@cigars;$j=$j+2)
		{
			if($cigars[($j+1)] eq 'M')
			{
				$max_m=$max_m>$cigars[$j]?$max_m:$cigars[$j];
				$max_m_j=$max_m>$cigars[$j]?$max_m_j:$j;
			}
			elsif($cigars[($j+1)] eq 'S')
			{
				$max_c=$max_c>$cigars[$j]?$max_c:$cigars[$j];
				$max_c_j=$max_c>$cigars[$j]?$max_c_j:$j;
			}
		}
		my $max_m_pad=0,$max_c_pad=0;
		for(my $k=0;$k<$max_m_j;$k=$k+2)
		{
			if($cigars[($k+1)] ne 'D')
			{
				$max_m_pad=$max_m_pad+$cigars[$k];
			}
			else
			{
				$max_m_pad=$max_m_pad-$cigars[$k];
			}
			
		}
		
		$max_match[$i]=$max_m_pad.":".$max_m;
		$max_clip[$i]=$max_c_pad.":".$max_c;
		if($max_c_j <$max_m_j and $max_c_j != -1 and $cigars[$max_c_j] >$cigar_filter_len)
		{
			$direction[$i]="1";
			$clip_len[$i]=$cigars[$max_c_j];
		}
		elsif($max_c_j > $max_m_j and $cigars[$max_c_j]>$cigar_filter_len)
		{
			$direction[$i]="2";
			$clip_len[$i]=$cigars[$max_c_j];
		}
		else
		{
			$direction[$i]="0";
			$clip_len[$i]=0;
		}
		$reads[5]=~s/#//g;
		$reads_total[$i]=$each_read;
		$i++;
	}
	
	return (\@direction,\@clip_len,\@reads_total,\@max_match)
}

