use List::Util qw/sum/;
use List::Util qw/max/;
#print "$ARGV[0]\t$ARGV[1]\t$ARGV[2]\t$ARGV[3]\n";
$ran_num=$ARGV[1];
$ori_dir="regions_$ran_num";
$dir="alu_sam";
$out_dir="split_softclipped_$ran_num";
$right_poly_AT_total=0;
$left_poly_AT_total=0;
if(not -e $out_dir)
{
	`mkdir $out_dir`;
}
$record="$ARGV[0]";
$record=~s/[ \t]+.*//;
$me_type=$ARGV[0];
$me_type=~s/.*[ \t]+//;
$bam_file="$ARGV[2]";
$REF="$ARGV[3]";
$step="$ARGV[4]";
my @group_ref,@group_left,@group_right;
#open(MPCHL,">>${out_dir}/$ARGV[0]_mapClipHL.sam");
#open(MPCHR,">>${out_dir}/$ARGV[0]_mapClipHR.sam");
if($step == 1)
{
#	@ori=`samtools view -T $REF $bam_file $record |head --lines=200 `;
	@ori=`cat regions_$ran_num/HG002_${record}.sam |head --lines=200 `;
	`rm regions_$ran_num/HG002_${record}.sam `;
}
elsif($step ==2)
{
#	@ori=`samtools view -T $REF $bam_file $record |head --lines=200 |perl modify_read_base.pl |perl get_second_alignment.pl  $bam_file $REF `;
	@ori=`cat regions_$ran_num/HG002_${record}.sam |head --lines=200 |perl modify_read_base.pl |perl get_second_alignment.pl  $bam_file $REF `;
	`rm regions_$ran_num/HG002_${record}.sam `;
}
#@ori=`samtools view -T $REF $bam_file $record |head --lines=200 `; #|perl modify_read_base.pl |perl get_second_alignment.pl  $bam_file $REF `;
my @read_ori;
my @left_map_start ,@left_map_end,@left_clip_seq,@left_clip_len,@left_read_i,@left_hard_read_i;
my @right_map_start ,@right_map_end,@right_clip_seq,@right_clip_len,@right_read_i,@right_hard_read_i;
my $i=0,$match_ref=0,$left_i=0,$right_i=0,$left_hard_i=0,$right_hard_i=0;
my $match_hard_right=0,$match_hard_right=0;

#HISEQ1:20:H9V1RADXX:2:1214:8847:71161   83      1       9996    0       51S97M  =       10063   -30     CCAAGCCCGAACCCGAACCTGAACAGTAACCGTAGTCCAAACCTGTACCCTTCCGATATCCCTTACCCTTACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTA    (&250&)&((2)&0(4+((++((:((0(()(((4+((8+(43(((++0(5(2&2(+3+((5,((55:5((;9?5((.>;.(..2?6(==<D=(==;C=8(89HHFD9HG@HF?GIIFD:IIIHFCIJIIIHJJJGHHHHHFFDFFCCC    NM:i:3  MD:Z:7A4A5A78   MC:Z:87M61S     AS:i:82 XS:i:79

my $alu_pos=$record;
$alu_pos=~s/.*:(\d+)\-(\d+).*/\1\t\2/;
my @alu_poss=split(/\t/,$alu_pos);
$alu_pos=sum(@alu_poss)/2;
foreach (@ori)
{
	chomp();
	@read=split(/\t/,$_);
	$read_ori[$i]=$_;
	$match_ref=0;	
	if($read[5]=~/(\d+)S/) # and $1<50)
	{
			
		$map_start=$read[3];
		my $tmp_cigar=$read[5];
		$tmp_cigar=~s/(\d+[^MD\d])//g;
		$tmp_cigar=~s/[MD]/#/g;
		$map_len=sum(split(/#/,$tmp_cigar));
		$map_end=$map_start+$map_len;	
		my $has_polyAT=0;
		$map_mid_pos=(($map_start+$map_end)/2);

		my $tmp_cigar=$read[5];
		
		if($tmp_cigar=~/\d+S$/ and $tmp_cigar=~/^\d+S/)	
		{
			$read_clip_part=1;
			if(abs($map_start-$alu_pos) < abs($map_end-$alu_pos))
			{
				$read_clip_part=0;	
			}
		}
		elsif($tmp_cigar=~/\d+S$/)	
		{
			$read_clip_part=1;
		}
		elsif($tmp_cigar=~/^\d+S/)	
		{
			$read_clip_part=0;
		}

		if($read_clip_part==1)
		{
			$tmp_cigar=~s/.*M//;
			my $clip_len=0,$clip_seq;
			$clip_len=$tmp_cigar;
			$clip_len=~s/.*?(\d+)S/\1/;
			if($tmp_cigar=~/S/)
			{
				$tmp_cigar=~s/[A-Z]/#/g;
				$clip_len=~s/.*?(\d+)S.*/\1/;
				$clip_seq=substr($read[9],-$clip_len);
				$clip_len_ran=int(rand(50));
				if($clip_len>50) 
				{
					$read[5]=~s/(\d+)S$/50S/;
					$read[9]=substr($read[9],0,length($read[9])-($clip_len-50));
					$read[10]=substr($read[10],0,length($read[10])-($clip_len-50));
					$clip_len=50;
					$clip_seq=substr($clip_seq,0,$clip_len);
					$read_ori[$i]=join("\t",@read);
				}
				$right_map_start[$right_i]=$map_start ;
				$right_map_end[$right_i]=$map_end;
				$right_clip_seq[$right_i]=$clip_seq;
				$right_clip_len[$right_i]=$clip_len;
				$right_read_i[$right_i]=$i;
				$right_i++;
			}
			else
			{
	#			$match_ref=1;
			}
		}
		else
		{
			$tmp_cigar=~s/\d+M.*//;
			my $clip_len=0,$clip_seq;
			$clip_len=$tmp_cigar;
			if($tmp_cigar=~/S/)
			{
				$tmp_cigar=~s/[A-Z]/#/g;
				#$clip_len=sum(split(/#/,$tmp_cigar));
				$clip_len=~s/.*?(\d+)S.*/\1/;
				$clip_seq=substr($read[9],0,$clip_len);
				$clip_len_ran=int(rand(50));
				if($clip_len>50) 
				{
					$read[5]=~s/^(\d+)S/50S/;
					$read[9]=substr($read[9],$clip_len-50);
					$read[10]=substr($read[10],$clip_len-50);
					$clip_len=50;
					$clip_seq=substr($clip_seq,-$clip_len);
					$read_ori[$i]=join("\t",@read);
				}
				
				$left_map_start[$left_i]=$map_start ;
				$left_map_end[$left_i]=$map_end;
				$left_clip_seq[$left_i]=$clip_seq;
				$left_clip_len[$left_i]=$clip_len;
				$left_read_i[$left_i]=$i;
				$left_i++;
			}
			else
			{
	#			$match_ref=1;
			}
		}
	}
	elsif($read[5]=~/H/)
	{
		if($read[5]=~/H$/)
		{
			if($read[5]=~/^\d+H/)
			{
	#			$match_ref=1;
			}
			else
			{
				$right_hard_read_i[$right_hard_i]=$i;
				$right_hard_i++;
			}
		}
		else
		{
			if($read[5]=~/H$/)
			{
	#			$match_ref=1;
			}
			else
			{
				$left_hard_read_i[$left_hard_i]=$i;
				$left_hard_i++;
			}
		}
	}
	else
	{
		$match_ref=1;
	}
	if($match_ref==1)
	{
		@line=split(/\t/,$_);
		@line=&clip_insertion(@line);
		@line=&clip_del(@line);
		#print MPR join("\t",@line)."\n";
		push @group_ref,join("\t",@line);

	}
	$i++;
}
my @left_polyAT,$left_polyAT_total=0,$left_polyAT_ratio=0,$polyAT_direction=0;
my @right_polyAT,$right_polyAT_total=0,$right_polyAT_ratio=0;
my $left_map_mid=&mid(@left_map_start);
my $right_map_mid=&mid(@right_map_end);
my $hard_part2_seq='';

my @right_clip_len_t=@right_clip_len;
for($t_i=0;$t_i <=$#right_clip_len_t;$t_i++)
{
    if(abs($right_map_end[$t_i]-$right_map_mid)>5)
    {
        $right_clip_len_t[$t_i]=0;
    }
}
for($hard_j=0;$hard_j<@right_clip_seq;$hard_j++)
{
	if($right_clip_len[$hard_j] == max(@right_clip_len_t) and abs($right_map_end[$hard_j]-$right_map_mid)<5)
	{
		$hard_part2_seq=$right_clip_seq[$hard_j];
		last;
	}
}
$fake_seq=0;
if(length($hard_part2_seq)==0)
{
	$hard_part2_seq='TCTCCAGGAATAGGAGCAGTGGGACCCCCTTCAAGAATC';
	$fake_seq=1;
}
$left_hard_supp_count=0;
$right_hard_supp_count=0;
foreach $hard_i (@right_hard_read_i)
{
	$hard_part2_len=length($hard_part2_seq);
	@tmp_hard_read_ori=split(/\t/,$read_ori[$hard_i]);
	$hard_len=&mapLen($tmp_hard_read_ori[5]);
	$tmp_hard_read_ori[9]=$tmp_hard_read_ori[9].$hard_part2_seq;
	$tmp_hard_read_ori[4]=50;
	$tmp_hard_read_ori[10]=$tmp_hard_read_ori[10].'B'x$hard_part2_len;
	$tmp_hard_read_ori[5]=~s/\d+H$/${hard_part2_len}S/;
	#	print "right:$hard_part2_len\n";
	if(abs($tmp_hard_read_ori[3]+$hard_len-$right_map_mid)<5 or $fake_seq==1)
	{
		#print MPCR join("\t",@tmp_hard_read_ori)."\n";
		push @group_right,join("\t",@tmp_hard_read_ori);
		$right_hard_supp_count=$right_hard_supp_count+1;
	}
	else
	{
		#print MPR $read_ori[$hard_i]."\n";
		push @group_ref,$read_ori[$hard_i];
	}

}
my @left_clip_len_t=@left_clip_len;
for($t_i=0;$t_i <=$#left_clip_len_t;$t_i++)
{
    if(abs($left_map_start[$t_i]-$left_map_mid)>5)
    {
        $left_clip_len_t[$t_i]=0;
    }
}
for($hard_j=0;$hard_j<@left_clip_seq;$hard_j++)
{
	if($left_clip_len[$hard_j] == max(@left_clip_len_t) and abs($left_map_start[$hard_j]-$left_map_mid)<5)
	{
		$hard_part2_seq=$left_clip_seq[$hard_j];	
		last;
	}
}
$fake_seq=0;
if(length($hard_part2_seq)==0)
{
	$hard_part2_seq='TCTCCAGGAATAGGAGCAGTGGGACCCCCTTCAAGAATC';
	$fake_seq=1;
}
foreach $hard_i (@left_hard_read_i)
{
	$hard_part2_len=length($hard_part2_seq);
	@tmp_hard_read_ori=split(/\t/,$read_ori[$hard_i]);
	$tmp_hard_read_ori[9]=$hard_part2_seq.$tmp_hard_read_ori[9];
	$tmp_hard_read_ori[10]='B'x$hard_part2_len.$tmp_hard_read_ori[10];
	$tmp_hard_read_ori[4]=50;
	$tmp_hard_read_ori[5]=~s/^\d+H/${hard_part2_len}S/;
	#	print "left:$hard_part2_len\n";
	if(abs($tmp_hard_read_ori[3]-$left_map_mid)<5 or $fake_seq==1)
	{
		#print MPCL join("\t",@tmp_hard_read_ori)."\n";
		push @group_left,join("\t",@tmp_hard_read_ori);
		$left_hard_supp_count=$left_hard_supp_count+1;
	}
	else
	{
		#print MPR $read_ori[$hard_i]."\n";
		#push @group_ref,$read_ori[$hard_i];
	}
}
for($i=0;$i<@left_clip_seq;$i++)
{
	#print "seq:$left_clip_seq[$i]\t".abs($left_map_mid-$left_map_start[$i])."\t$left_read_i[$i]\t".max(@left_clip_len)."\n";
	if($left_clip_seq[$i]=~/[AT]{7,}/ or (@{[$left_clip_seq[$i]=~m/[^AT]/g]} <= 3 and abs($left_map_mid-$left_map_start[$i])<5))
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
	if($right_clip_seq[$i]=~/[AT]{7,}/ or (@{[$right_clip_seq[$i]=~m/[^AT]/g]} <= 3 and abs($right_map_mid-$right_map_end[$i])<5))
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

$left_print=0;
$right_print=0;
@group_left_clip;
@group_right_clip;
if($left_poly_AT_total>$right_poly_AT_total and sum(@left_has_polyAT)/($#left_clip_seq+1)>0.6)
{
	$polyAT_direction=-1;
	for($i=0;$i<@left_has_polyAT;$i++)
	{
		my @tmp_read_ori=split(/\t/,$read_ori[$left_read_i[$i]]);
		#$tmp_read_ori[6]='1';
		#if($tmp_read_ori[2]=~/^chr/)
		#{
		#	$tmp_read_ori[6]='chr1';
		#}	
		#$tmp_read_ori[7]=$left_clip_len[$i];
		if($left_has_polyAT[$i]==1)
		{
			#print MPCL join("\t",@tmp_read_ori)."\n";
			push @group_left,join("\t",@tmp_read_ori);
			push (@group_left_clip,($left_clip_seq[$i]));
			$left_print++;
		}
		else
		{
			@line=&clip_insertion(@tmp_read_ori);
			@line=&clip_del(@line);
			#print MPR join("\t",@line)."\n";
			#push @group_ref,join("\t",@line);
		}
	}
}
elsif($left_poly_AT_total<$right_poly_AT_total and sum(@right_has_polyAT)/($#right_clip_seq+1)>0.6)
{
	$polyAT_direction=1;
	for($i=0;$i<@right_has_polyAT;$i++)
	{
		my @tmp_read_ori=split(/\t/,$read_ori[$right_read_i[$i]]);
		#$tmp_read_ori[6]='2';
		#if($tmp_read_ori[2]=~/^chr/)
		#{
		#	$tmp_read_ori[6]='chr2';
		#}	
		#$tmp_read_ori[7]=$right_clip_len[$i];
		if($right_has_polyAT[$i]==1)
		{
			#print MPCR join("\t",@tmp_read_ori)."\n";
			push @group_right,join("\t",@tmp_read_ori);
			push (@group_right_clip,($right_clip_seq[$i]));
			$right_print++;
		}
		else
		{
			@line=&clip_insertion(@tmp_read_ori);
			@line=&clip_del(@line);
			#print MPR join("\t",@line)."\n";
			#push @group_ref,join("\t",@line);
		}
	}
}

my @left_map_alu;
my @right_map_alu;
my @blastout;
if($polyAT_direction!=-1)
{
	for(my $i=0;$i<@left_clip_seq;$i++)
	{
		$left_map_alu[$i]=0;
	}
	@blast_out=`cat $out_dir/${record}.rp_left 2>/dev/null`;
	my $tmp_seq,@left_map_alu_seq;
	foreach $i (@blast_out)
	{
		push @left_map_alu_seq,$left_clip_seq[$i];
		$left_map_alu[$i]=1;
	}	
	chomp(@left_map_alu_seq);	
	for(my $i=0;$i<@left_clip_seq;$i++)
	{
		my @tmp_read_ori=split(/\t/,$read_ori[$left_read_i[$i]]);
		#$tmp_read_ori[6]='1';
		#if($tmp_read_ori[2]=~/^chr/)
		#{
		#	$tmp_read_ori[6]='chr1';
		#}	
		#$tmp_read_ori[7]=$left_clip_len[$i];
#		print "read:$tmp_read_ori[0]\tleft_start:$left_map_start[$i]\t$left_map_mid\n";
		if($left_map_alu[$i]==0)
		{
#				print abs($right_map_end[$i]-$right_map_mid)."\t".&seq_compare($t_seq,$right_clip_seq[$i],1)."\t".$right_clip_seq[$i]."\t$t_seq\n";
			if(abs($left_map_start[$i]-$left_map_mid)<5)
			{
				$left_map_alu[$i]=1;
			}
		}
		if($left_map_alu[$i]==0)
		{
			@line=&clip_insertion(@tmp_read_ori);
			@line=&clip_del(@line);
			#print MPR join("\t",@line)."\n";
			#push @group_ref,join("\t",@line);
		}
		else
		{
			#print MPCL join("\t",@tmp_read_ori)."\n";
			push @group_left,join("\t",@tmp_read_ori);
			push (@group_left_clip,($left_clip_seq[$i]));
			$left_print++;
		}
	}
	
}
if($polyAT_direction!=1)
{
	for(my $i=0;$i<@right_clip_seq;$i++)
	{
		$right_map_alu[$i]=0;
	}
	@blast_out=`cat $out_dir/${record}.rp_right 2>/dev/null`;
	my $tmp_seq,@right_map_alu_seq=();
	foreach $i (@blast_out)
	{
		push @right_map_alu_seq,$right_clip_seq[$i];
		$right_map_alu[$i]=1;
	}	
	chomp(@right_map_alu_seq);
	for(my $i=0;$i<@right_clip_seq;$i++)
	{
		my @tmp_read_ori=split(/\t/,$read_ori[$right_read_i[$i]]);
		#$tmp_read_ori[6]='2';
		#if($tmp_read_ori[2]=~/^chr/)
		#{
		#	$tmp_read_ori[6]='chr2';
		#}	
		#$tmp_read_ori[7]=$right_clip_len[$i];
		if($right_map_alu[$i]==0)
		{
#			print abs($right_map_end[$i]-$right_map_mid)."\n";
			if(abs($right_map_end[$i]-$right_map_mid)<5)
			{
				$right_map_alu[$i]=1;
			}
		}
		if($right_map_alu[$i]==0)
		{
			@line=&clip_insertion(@tmp_read_ori);
			@line=&clip_del(@line);
			#print MPR join("\t",@line)."\n";
			#push @group_ref,join("\t",@line);
		}
		else
		{
			#print MPCR join("\t",@tmp_read_ori)."\n";
			push @group_right,join("\t",@tmp_read_ori);
			push (@group_right_clip,($right_clip_seq[$i]));
			$right_print++;
		}
	}
}

my $clipR_chrom=0;
my $clipL_chrom=0;
my $clip_reads_count=0;
my $reads_count=0;
foreach $read_i (@group_right)
{
	@read_info=split(/\t/,$read_i);
	if($read_info[6] eq "=")
	{
		$clipR_chrom++;
	}
	if($read_info[5] =~/^\d+S/ and $read_info[5] =~/S$/)
	{
		$clip_both_count++;
	}
	$reads_count++;
}
foreach $read_i (@group_left)
{
	@read_info=split(/\t/,$read_i);
	if($read_info[6] eq "=")
	{
		$clipL_chrom++;
	}
	if($read_info[5] =~/^\d+S/ and $read_info[5] =~/S$/)
	{
		$clip_both_count++;
	}
	$reads_count++;
}
if($reads_count !=0)
{
	$clip_both=$clip_both_count/$reads_count;
}
else
{
	$clip_both=0
}
#my $clipL_chrom=`bash clip_same_chrom.sh ${out_dir}/${record}_mapClipL.sam`;
#my $clipR_chrom=`bash clip_same_chrom.sh ${out_dir}/${record}_mapClipR.sam`;
#my $clip_both=`bash clip_both.sh ${out_dir}/${record}_mapClipL.sam ${out_dir}/${record}_mapClipR.sam `;
#                        -b <(cat ../data_cluster/tmp_3_tsd.bed |perl -F'\t' -alne 'if($F[1]>=-10 and $F[1] <=50 and $F[2]>=2 and $F[3]>=2 and $F[4]<0.3) {print "$_";}' |cut -f1|perl -npe "s/HG002_//"|cut -d':' -f 2-|perl -npe "s/:/\t/;s/\-/\t/" ) \
my $tsd_len=$right_map_mid-$left_map_mid;
if($tsd_len>=-20 and  $tsd_len<=70 and $clipL_chrom>=2 and $clipR_chrom>=2 and $clip_both<0.3)
{
	open(BPINF,">${out_dir}/HG002_${record}_BPinfo.txt");
	print BPINF "$record\t$left_map_mid\t$right_map_mid\t";
	#	`cat head_$ran_num.sam >${out_dir}/HG002_${record}_mapRef.sam`;
	#`cat head_$ran_num.sam >${out_dir}/HG002_${record}_mapClipR.sam`;
	#`cat head_$ran_num.sam >${out_dir}/HG002_${record}_mapClipL.sam`;
	@record_split=split(/[:\-]/,$record);
	#1	527089	0	HG002	ALU
#	print join("\n",@group_left_clip)."\n";
	my $iden_sum_left=&clip_seq_iden(-1,@group_left_clip);
#	print join("\n",@group_right_clip)."\n";
	my $iden_sum_right=&clip_seq_iden(1,@group_right_clip);
	#print "$record_split[0]\t".(($record_split[1]+$record_split[2])/2)."\t0\tHG002\t$me_type\t$iden_sum_left\t$iden_sum_right\n";
	my $tsd_seq="";
	if($left_map_mid+1<$right_map_mid)
	{
		my $right_map_mid_t=$right_map_mid-1;
		$tsd_seq=`bash getfastq.sh $REF $record_split[0] $left_map_mid $right_map_mid_t`;
	}
	if($tsd_seq=~/^.?A{10,}.?$/ or $tsd_seq=~/^.?T{10,}.?$/ or $tsd_seq=~/^.?C{10,}.?$/ or $tsd_seq=~/^.?G{10,}.?$/  or $tsd_seq=~/^.?(AT){4,}.?$/ )
	{
	}	
	elsif($iden_sum_left>=0.8 and  $iden_sum_right>=0.8)
	{
        my $map_mid_all=int(($left_map_mid+$right_map_mid)/2);
    	open(MPR,">${out_dir}/HG002_${record}_mapRef.sam");
    	open(MPCL,">${out_dir}/HG002_${record}_mapClipL.sam");
    	open(MPCR,">${out_dir}/HG002_${record}_mapClipR.sam");
    	print "$record_split[0]\t".int(($record_split[1]+$record_split[2])/2)."\t0\tHG002\t$me_type\n";#$iden_sum_left\t$iden_sum_right\n";
    	foreach $read_i (@group_ref)
    	{
		@F=split(/\t/,$read_i);
		if($F[5] ne '*' and $F[4] == 0 )
		{
#			$F[4]=60;
		}
		print MPR join("\t",@F);
		print MPR "\n";
    	}
	$i=0;
    	foreach $read_i (@group_left)
    	{
		$i++;
		@F=split(/\t/,$read_i);
		if($F[5] ne '*' and $F[4] == 0 )
		{
#			$F[4]=60;
		}
		print MPCL join("\t",@F);
		print MPCL "\n";
		if($i>15)
		{
			last;
		}
    	}
	$i=0;
    	foreach $read_i (@group_right)
    	{
		$i++;
		@F=split(/\t/,$read_i);
		if($F[5] ne '*' and $F[4] == 0 )
		{
#			$F[4]=60;
		}
		print MPCR join("\t",@F);
		print MPCR "\n";
		if($i>15)
		{
			last;
		}
    	}
    	close(MPCR);
    	close(MPCL);
    	close(MPR);
	}
	print BPINF "$tsd_len\t$clipL_chrom\t$clipR_chrom\t$clip_both\t$left_hard_supp_count\t$right_hard_supp_count\t$polyAT_direction\n";
	close(BPINF);
}
#print "left:$left_print\tright:$right_print\n";
sub clip_seq_iden
{
	my ($dir,@seq_long)=@_;
	my %BC;
	my $base="";
	my @iden_sum;
	if($dir eq '-1')
	{
		for(my $i=0;$i<5;$i++)
		{
			for(my $j=0;$j<=$#seq_long;$j++)
			{
				$base=substr($seq_long[$j],-1*($i+1),1);	
				if(length($seq_long[$j])>=($i+1))
				{
					$BC{$base}++;
				}
			}
			if(sum(values %BC) eq 0 )
			{
				next;
			}
			push(@iden_sum,(max(values %BC)/sum(values %BC)));
			%BC="";
		}
	}
	else
	{
		for(my $i=0;$i<5;$i++)
		{
			for(my $j=0;$j<=$#seq_long;$j++)
			{
				$base=substr($seq_long[$j],$i,1);	
				if(length($seq_long[$j])>=($i+1))
				{
					$BC{$base}++;
				}
			}
			if(sum(values %BC) eq 0 )
			{
				next;
			}
			push(@iden_sum,(max(values %BC)/sum(values %BC)));
			%BC="";
		}
	}
	return sum(@iden_sum)/($#iden_sum+1);
}


sub seq_compare
{
	my $seq_long=$_[0];
	my $seq_short=$_[1];
	my $direction=$_[2];
	my $compare_len=length($seq_short)/3,$seq_long;
	if(length($seq_short)<=3)
	{
		$compare_len=1;
	}
#	print "start\n";
	for(my $i=0;$i<$compare_len;$i++)
	{
		if($direction==-1)
		{
			$seq_short_part_r=substr($seq_short,0,length($seq_short)-$i);
			$seq_short_part_l=substr($seq_short,$i);
			if($seq_long=~/$seq_short_part_r$/ or $seq_long=~/$seq_short_part_l$/)
			{
				#print "left match\n";
				return 1;
			}
		}
		else
		{
			$seq_short_part_r=substr($seq_short,0,length($seq_short)-$i);
			$seq_short_part_l=substr($seq_short,$i);
			#print "$seq_short_part_l\t$seq_short_part_r\t$seq_long\n";
			if($seq_long=~/^$seq_short_part_r/ or $seq_long=~/^$seq_short_part_l/)
			{
				print "right match\n";
				return 1;
			}
		}
	}
	return 0;
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

sub clip_insertion
{
	my $cigar_t,my $max_cg ,my $left_len,my $right_len ,my @cigar_len,my $clip_len_read,my $clip_len_ref;
	(@read)=@_;
	if($read[5]=~/[^MDSI\d]/)
	{
		return @read
	}
	$cigar_t=$read[5];
	$cigar_t=~s/.*?(\d+)I.*/\1/g;
	if($read[5]=~/I/ and $cigar_t>5)
	{
		$max_cg=&maxCigar($read[5],'I');
#		print "$read[5]\t$read[9]\n";
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
#		print "$read[5]\t$read[9]\n";
	}
	return @read
}

sub clip_del
{
	my $cigar_t,my $max_cg ,my $left_len,my $right_len ,my @cigar_len,my $clip_len_read,my $clip_len_ref;
	(@read)=@_;
	if($read[5]=~/[^MDSI\d]/)
	{
		return @read
	}
	$cigar_t=$read[5];
	$cigar_t=~s/.*?(\d+)D.*/\1/g;
	if($read[5]=~/D/ and $cigar_t>5)
	{
		$max_cg=&maxCigar($read[5],'D');
#		print "0###$read[3]\t$read[5]\t$max_cg\n";
		$cigar_t=$read[5];
		$cigar_t=~s/${max_cg}D.*//;
		$left_len=&mapLen($cigar_t);

		$cigar_t=$read[5];
		$cigar_t=~s/.*?${max_cg}D//;
		$right_len=&mapLen($cigar_t);
#		print "left_len:$left_len\tright_len:$right_len\n";	
		if($left_len>=$right_len)
		{
			$read[5]=~s/(${max_cg})D(.*)/\1S/;
			$cigar_t=$read[5];
#			$cigar_t=~s/\d+S$//;
			$cigar_t=~s/\d+D//g;
			$cigar_t=~s/[A-Z]/\t/g;
			@cigar_len=split(/\t/,$cigar_t);
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
			#print "1###$clip_len_read  max_cg:$max_cg \n";
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
                                $clip_len_read=-1*length($read[9]);
			}
#			print "1###$clip_len_read\n";
			$read[9]=substr($read[9],$clip_len_read,);
			$read[10]=substr($read[10],$clip_len_read,);
		}
#		print "$read[3]\t$read[5]\n";
	}
	return @read
}
sub maxCigar
{
	(my $cigar,my $find_c)=@_;
	my @cigars;
	$cigar=~s/(\d+$find_c)/\1\t/g;
	$cigar=~s/\t$//;
	@cigars=split(/\t/,$cigar);
	for(my $i=0;$i<@cigars;$i++)
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
	(my $cigar_t)=@_;
	my  @cigar_len;
	$cigar_t=~s/([A-Z])/\1\t/g;
	@cigar_len=split(/\t/,$cigar_t);
	for(my $i=0;$i<@cigar_len;$i++)
	{
		unless($cigar_len[$i]=~/[MD]/)
		{
			$cigar_len[$i]=0;	
		}
	}
	return sum(@cigar_len)
}	

sub std
{
	my @pos=@_;
	my $mean=&mid(@pos);
	my $sum_2=0;
	my $std_j=0;
	#	print "\nmean:$mean:".length(@pos)."\n";
	for(my $std_i=0;$std_i<=$#pos;$std_i++)
	{
		#	print "$std_i:$pos[$std_i]\t";
		if(abs($pos[$std_i]-$mean)<10)
		{
			$sum_2=$sum_2+($pos[$std_i]-$mean)*($pos[$std_i]-$mean);
			$std_j++;
		}
	}
	#if($std_j==0){	return 0 }
	if($std_j==0){	return 1 }
	my $std_value=sqrt($sum_2/$std_j);
	my $soft_val=$std_j/($#pos+1);
	return $soft_val
	#	return $std_value
}
