#HISEQ1:19:H8VDAADXX:1:2212:8084:89459	401	7	138791545	0	68M80H	=	138791246	-367	GTTTTAGAGATACTAAAGTACTTAGTAATCAAATGGTTGCTGAGATTTGCTTTAAAATAACACATTGG	CDDDDEDEEEDCDDEEEDDCDDEDDEDDEEEDDDDDDDDDDDDDCC>DDDDDDEEDEDC@CACDDCDD	NM:i:0	MD:Z:68	MC:Z:148M	AS:i:68	XS:i:22	SA:Z:2,124793723,+,83M65S,0,0;
$bam_file=$ARGV[0];
$ref=$ARGV[1];
while(<STDIN>)
{
	$print=0;
	chomp();
	@F=split(/\t/,$_);
	if($F[5]=~/H$/ and $F[5]=~/^\d+H/)
	{
		next;
	}
	if($F[5]=~/H/)
	{
		foreach $i (@F)
		{
			$sa_str=$i;
			if($sa_str=~s/SA:Z://)
			{
				@sa_all=split(/;/,$sa_str);
				foreach $j (@sa_all)
				{ 
					#print "$j\n";
					if($j=~/H/) {}
					else{
						@sa=split(/,/,$j);
						$sa[3]=~/(\d+)M/;
						$sa_mp_len=$1;
						$sa_start=$sa[1]-3+$sa_mp_len;
						$sa_end=$sa_start+2;
						$se_pos="$sa[0]:$sa_start-$sa_end";
						@se_ori=split(/\t/,`samtools view $bam_file -T $ref $se_pos |grep -m 1 $F[0]|cut -f10,11|head -n 1`);
						#print "samtools view $bam_file -T $ref $se_pos |grep -m 1 $F[0]|cut -f10,11|head -n 1";
						$se_seq=$se_ori[0];
						$se_qual=$se_ori[1];
						$se_qual=~s/\n$//;
						$se_seq_rev=reverse $se_seq;
						$se_seq_rev=~s/A/M/g;$se_seq_rev=~s/T/A/g;$se_seq_rev=~s/M/T/g;$se_seq_rev=~s/C/M/g;$se_seq_rev=~s/G/C/g; $se_seq_rev=~s/M/G/g;
						if($se_seq=~/$F[9]/){}
						elsif($se_seq_rev=~/$F[9]/)
						{
							$se_qual=reverse $se_qual;
							$se_seq=$se_seq_rev;
						}
						else {next;}
						$cigar_m=$F[5];
						$cigar_m=~s/(.*?)(\d+)H(.*)/\1\2S\3/;
						$clip_len=$2;
						if($F[5]=~/H$/)
						{
							$se_seq=~/$F[9](.*)/;
							$clip_seq=$1;
							$clip_seq_len=length($clip_seq);
							$clip_qual=substr($se_qual,-1*$clip_seq_len);
							$F[9]="$F[9]$clip_seq";
							$F[10]="$F[10]$clip_qual";
						}
						else
						{
							$se_seq=~/(.*)$F[9]/;
							$clip_seq=$1;
							$clip_seq_len=length($clip_seq);
							$clip_qual=substr($se_qual,0,$clip_seq_len);
							$F[9]="$clip_seq$F[9]";
							$F[10]="$clip_qual$F[10]";
						}
						if($clip_len != length($clip_seq)) { next;}
						$F[5]=$cigar_m;
						print join("\t",@F); print "\n";
						$print=1; last;
						#print "$F[0]\t$print\n";
					}
				}
				last;
			}
		}
	}
	if($print ==0)
	{
		print "$_\n";
	}
}
