#HISEQ1:20:H9V1RADXX:1:2211:5155:41672   113     10      101340285       8       43S105M X       136864278       0       GAGAATTCTTTTTTTTTTTTTTGAGACGGAGTCTCGCTCTTTCGCCCAGGCCGGACTGCAGTGGCGCTATCTCAGCTCACTGCAAGCTCTGCCTCATCGTGGAGAATTCTATTGCCTCTCAAAGTGTGGTCCATAGACCAGCAGAATC    DDCCCDDBDDDDDDDDDDDDDEDDDDDDDDDDDDDDDDDDDDDBDBDDDDDDDDDDDDDDDDDDDDDEEEDEFFFFFFHHHHHJHHIJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJIJJJJJHHHHHFFFFFCCC    NM:i:3  MD:Z:9A2G4A87   MC:Z:20S128M    AS:i:90 XS:i:85 SA:Z:8,133030379,+,83S58M7S,8,2;
$REF=$ARGV[0];
$bam_file=$ARGV[1];
$file_pre=$ARGV[2];
%dis_name=();
%dis_strand=();
%dis_r=();
$dis_region="";
$max_clip_L=0;
$max_clip_L_seq='';
open CL,"${file_pre}_mapClipL.sam";
while(<CL>)
{
	chomp();
	if($_=~/^@/)
	{
		next;
	}
	@F=split(/\t/,$_);
	#read map reverse strand
	if($F[1]%32 >=16)
	{
		$dis_name{$F[0]}=1;
		$dis_region="$F[6]:".($F[7]-1)."-".($F[7]+1)." $dis_region";

		$dis_strand{$F[0]}='-';
		$dis_r{$F[0]}=$F[1];
			
	}
	$F[5]=~/^(\d+)S/;
	$max_clip_t=$1;
	if($max_clip_L<$max_clip_t)
	{
		$max_clip_L=$max_clip_t;
		$max_clip_L_seq=substr($F[9],0,$max_clip_L)
	}
}
close CL;
$max_clip_R=0;
$max_clip_R_seq='';
open CR,"${file_pre}_mapClipR.sam";
while(<CR>)
{
	chomp();
	if($_=~/^@/)
	{
		next;
	}
	@F=split(/\t/,$_);
	#read map reverse strand
	if($F[1]%32 <16)
	{
		$dis_name{$F[0]}=1;
		$dis_region="$F[6]:".($F[7]-1)."-".($F[7]+1)." $dis_region";

		$dis_strand{$F[0]}='+';
		$dis_r{$F[0]}=$F[1];
	}
	$F[5]=~/^(\d+)S/;
	$max_clip_t=$1;
	if($max_clip_R<$max_clip_t)
	{
		$max_clip_R=$max_clip_t;
		$max_clip_R_seq=substr($F[9],0,$max_clip_R)
	}
}

@dis_read_raw=`samtools view $bam_file $dis_region`;

@dis_seq_l=();
@dis_seq_r=();
foreach $read_i (@dis_read_raw)
{
	chomp($read_i);
	@F=split(/\t/,$read_i);
	if($dis_name{$F[0]}==1 and (($dis_r{$F[0]}%256 <128 and  $F[1]%256 >=128) or ($dis_r{$F[0]}%256 >= 128 and  $F[1]%256 < 128)))
	{
		$seq=$F[9];
		if(($F[1]%32 >=16 and $dis_r{$F[0]}%32 >=16) or ($F[1]%32 <16 and $dis_r{$F[0]}%32 <16) )
		{
			$seq=reverse($seq);$seq=~s/A/M/g;$seq=~s/T/A/g;$seq=~s/M/T/g;$seq=~s/C/M/g;$seq=~s/G/C/g;$seq=~s/M/G/g;
		}
		if($dis_strand{$F[0]} eq '+')
		{
			push(@dis_seq_r,$seq);
		}
		else
		{
			push(@dis_seq_l,$seq);
		}
		$dis_name{$F[0]}=0;
	}
}
foreach $seq_i (@dis_seq_l)
{
	if($max_clip_L_seq

