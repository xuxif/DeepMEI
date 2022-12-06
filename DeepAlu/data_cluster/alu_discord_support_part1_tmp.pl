use List::Util qw/sum/;
$ori_dir="regions_extra";
$dir="alu_sam";
$out_dir="split_softclipped";
$bwa="/public/home/swgenetics_4/software/bwa-0.7.15/bwa";
$bwa="bwa";
if(not -e $out_dir)
{
	`mkdir $out_dir`;
}
$ori_file="$ARGV[0]";
$ARGV[0]=~s/\..*?$//;
open(MPR,">${out_dir}/$ARGV[0]_mapRef.sam");
open(MPCL_ME_LONG,">${out_dir}/$ARGV[0]_mapClipL_long.fa");
open(MPCL_ME_SHORT,">${out_dir}/$ARGV[0]_mapClipL_short.fa");
open(MPCR_ME_LONG,">${out_dir}/$ARGV[0]_mapClipR_long.fa");
open(MPCR_ME_SHORT,">${out_dir}/$ARGV[0]_mapClipR_short.fa");
if($ori_file=~/bam$/)
{
	@ori=`samtools view ${ori_dir}/$ori_file`;
}
else
{
	@ori=`cat ${ori_dir}/$ori_file`;
}
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
				if($clip_len<15)
				{
					print MPCR_ME_SHORT ">$ARGV[0]:$i\n$clip_seq\n";
				}
				else
				{
					print MPCR_ME_LONG ">$ARGV[0]:$i\n$clip_seq\n";
				}
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
				if($clip_len<15)
				{
					print MPCL_ME_SHORT ">$ARGV[0]:$i\n$clip_seq\n";
				}
				else
				{
					print MPCL_ME_LONG ">$ARGV[0]:$i\n$clip_seq\n";
				}
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
close  MPCL_ME_SHORT;
close  MPCL_ME_LONG;
close  MPCR_ME_SHORT;
close  MPCR_ME_LONG;
