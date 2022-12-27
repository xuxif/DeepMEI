use List::Util qw/sum/;
$ori_dir="regions";
$dir="alu_sam";
$out_dir="split_softclipped";
$bwa="/public/home/swgenetics_4/software/bwa-0.7.15/bwa";
if(not -e $out_dir)
{
	`mkdir $out_dir`;
}
if(not $ARGV[0])
{
	$file_prefix="HG00096_chr10:26573519-26573619";
	$file_prefix="HG002_chr1:102454924-102455024";
	$discord_file="${file_prefix}_discord.sam";
	$ss_file="${file_prefix}_alu.sam";
	$ori_file="${file_prefix}.bam";
	$ARGV[0]=$file_prefix;
}
else
{
	$discord_file="$ARGV[0]_discord.sam";
	$ss_file="$ARGV[0]_alu.sam";
	$ori_file="$ARGV[0].bam";
}
open(MPR,">>${out_dir}/$ARGV[0]_mapRef.sam");
open(MPCL,">>${out_dir}/$ARGV[0]_mapClipL.sam");
open(MPCR,">>${out_dir}/$ARGV[0]_mapClipR.sam");
`cat head.sam >${out_dir}/$ARGV[0]_mapRef.sam`;
`cat head.sam >${out_dir}/$ARGV[0]_mapClipR.sam`;
`cat head.sam >${out_dir}/$ARGV[0]_mapClipL.sam`;
@ori=`samtools view ${ori_dir}/$ori_file`;
my @read_ori;
my @left_map_start ,@left_map_end,@left_clip_seq,@left_clip_len,@left_read_i;
my @right_map_start ,@right_map_end,@right_clip_seq,@right_clip_len,@right_read_i;
my $i=0,$match_ref=0;
foreach (@ori)
{
	chomp();
	@read=split(/\t/,$_);
	$read_ori[$i]=$_;	
	$match_ref=0;	
	if($read[5]=~/S/)
	{
		$map_start=$read[3];
		my $tmp_cigar=$read[5];
		my $alu_pos=$ori_file;
		$alu_pos=~s/.*://;
		$alu_pos=~s/\-.*//;
		$alu_pos+=50;
		$tmp_cigar=~s/(\d+[^MD\d])//g;
		$tmp_cigar=~s/[MD]/#/g;
		$map_len=sum(split(/#/,$tmp_cigar));
		$map_end=$map_start+$map_len;	

		my $has_polyAT=0;
		$tmp_cigar=$read[5];
		if(($map_start+$map_end)/2 < $alu_pos)
		{
			$tmp_cigar=~s/.*M//;
			my $clip_len=0,$clip_seq;
			if($tmp_cigar=~/S/)
			{
				$tmp_cigar=~s/[A-Z]/#/g;
				$clip_len=sum(split(/#/,$tmp_cigar));
				$clip_seq=substr($read[9],-$clip_len);

				$right_map_start[$right_i]=$map_start ;
				$right_map_end[$right_i]=$map_end;
				$right_clip_seq[$right_i]=$clip_seq;
				$right_clip_len[$right_i]=$clip_len;
				$right_read_i[$right_i]=$i;
				$right_i++;
			}
			else
			{
				$match_ref=1;
			}
		}
		else
		{
			$tmp_cigar=~s/\d+M.*//;
			my $clip_len=0,$clip_seq;
			if($tmp_cigar=~/S/)
			{
				$tmp_cigar=~s/[A-Z]/#/g;
				$clip_len=sum(split(/#/,$tmp_cigar));
				$clip_seq=substr($read[9],0,$clip_len);

				$left_map_start[$left_i]=$map_start ;
				$left_map_end[$left_i]=$map_end;
				$left_clip_seq[$left_i]=$clip_seq;
				$left_clip_len[$left_i]=$clip_len;
				$left_read_i[$left_i]=$i;
				$left_i++;
			}
			else
			{
				$match_ref=1;
			}
		}
	}
	else
	{
		$match_ref=1;
	}
	if($match_ref==1)
	{
		print MPR "$_\n";
	}
	$i++;
}
my @left_polyAT,$left_polyAT_total=0,$left_polyAT_ratio=0,$polyAT_direction=0;
my @right_polyAT,$right_polyAT_total=0,$right_polyAT_ratio=0;
my $left_map_mid=&mid(@left_map_start);
my $right_map_mid=&mid(@right_map_end);
for($i=0;$i<@left_clip_seq;$i++)
{
	if($left_clip_seq[$i]=~/[AT]{7,}/ or (@{[$left_clip_seq[$i]=~m/[^AT]/g]} <= 1 and abs($left_map_mid-$left_map_start[$i])<5))
	{
		$left_has_polyAT[$i]=1;
		$left_poly_AT_total+=@{[$left_clip_seq[$i]=~m/[AT]{7,}/g]};	
	}	
	else
	{
		$left_has_polyAT[$i]=0;
	}	
#	print "$left_clip_seq[$i]\t$left_has_polyAT[$i]\n";
}	
for($i=0;$i<@right_clip_seq;$i++)
{
	if($right_clip_seq[$i]=~/[AT]{7,}/ or (@{[$right_clip_seq[$i]=~m/[^AT]/g]} <= 1 and abs($right_map_mid-$right_map_end[$i])<5))
	{
		$right_has_polyAT[$i]=1;
		$right_poly_AT_total+=@{[$right_clip_seq[$i]=~m/[AT]{7,}/g]};	
	}	
	else
	{
		$right_has_polyAT[$i]=0;
	}	
#	print "$right_clip_seq[$i]\t$right_has_polyAT[$i]\t".abs($right_map_mid-$right_map_end)."\n";
}

if($left_poly_AT_total>$right_poly_AT_total and sum(@left_has_polyAT)/@left_clip_seq>0.6)
{
	$polyAT_direction=-1;
	for($i=0;$i<@left_has_polyAT;$i++)
	{
		my @tmp_read_ori=split(/\t/,$read_ori[$left_read_i[$i]]);
		$tmp_read_ori[6]='chr1';
		$tmp_read_ori[7]=$left_clip_len[$i];
		if($left_has_polyAT[$i]==1)
		{
			print MPCL join("\t",@tmp_read_ori)."\n";
		}
		else
		{
			print MPR join("\t",@tmp_read_ori)."\n";
		}
	}
}
elsif($left_poly_AT_total<$right_poly_AT_total and sum(@right_has_polyAT)/@right_clip_seq>0.6)
{
	$polyAT_direction=1;
	for($i=0;$i<@right_has_polyAT;$i++)
	{
		my @tmp_read_ori=split(/\t/,$read_ori[$right_read_i[$i]]);
		$tmp_read_ori[6]='chr2';
		$tmp_read_ori[7]=$right_clip_len[$i];
		if($right_has_polyAT[$i]==1)
		{
			print MPCR join("\t",@tmp_read_ori)."\n";
		}
		else
		{
			print MPR join("\t",@tmp_read_ori)."\n";
		}
	}
}

my @left_map_alu;
my @right_map_alu;
my @blastout;
if($polyAT_direction!=-1)
{
	open T_FASTA,">tmp_left.fa";
	for($i=0;$i<@left_clip_seq;$i++)
	{
		if($left_clip_len[$i]>15)
		{
			print T_FASTA ">$i\n$left_clip_seq[$i]\n";
		}
		$left_map_alu[$i]=0;
	}
	close T_FASTA;
	@blast_out=`blastn -db ME_ref.fa -query tmp_left.fa -evalue 0.0001 -outfmt "6 qseqid	sstrand"  -out /dev/stdout |sort |uniq `;
	print @blast_out;
	foreach $i (@blast_out)
	{
		$i=~s/\t.*//;
		$left_map_alu[$i]=1;
	}	
	
}
if($polyAT_direction!=1)
{
	open T_FASTA,">tmp_right.fa";
	for($i=0;$i<@right_clip_seq,$right_clip_len[$i]>15;$i++)
	{
		if($right_clip_len[$i]>15)
		{
			print T_FASTA ">seq_$right_read_i[$i]\n$right_clip_seq[$i]\n";
		}
	}
	close T_FASTA;
}



sub mid
{
    my @list = sort @_;
    my $count = @list;
    if( $count == 0 )
    {
        return 0;
        die "mid:reads is empty ";
    }   
    return $list[($count-1)/2];
}
