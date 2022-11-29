#A00265:490:HGYGLDSXY:3:2228:22878:33223 353     1       223911529       0       46M104H 8       9526688 0       CCAGGAATAGGAGCAGTGGGACCCCCTTCAAGAATCCCCCATTTCC  @BDCDC?BAABCAABAAABBC?A@AABCABCACCBAAAAABBCCAA  SA:Z:3,46550021,-,111M39S,0,0;  MC:Z:150M       MD:Z:46 PG:Z:MarkDuplicates     RG:Z:4  NM:i:0  AS:i:46 XS:i:21
#A00265:490:HGYGLDSXY:3:2658:25753:17503 433     1       223911558       0       30H30M  8       42185462        0       AAGAATCCCCCATTTCCCATATACTTTCCT  CBA1BCA<<<ACCC..ABA/A?A?ED>CA  SA:Z:2,27662856,+,21S39M,0,1;   MC:Z:150M       MD:Z:30 PG:Z:MarkDuplicates     RG:Z:4  NM:i:0  AS:i:30 XS:i:0>>>>))))?/???
while(<STDIN>)
{
	chomp();
	@read=split(/\t/,$_);
	my $clip_len=0,@sa;
	if($read[5]=~/H$/)
	{
		$read[5]=~/(\d+)H$/;
		$part_1_len=$1;
		$read[5]=~s/(\d+)H$/\1S/;
		$clip_len=$1;
		for($sa_i=10;$sa_i<@read;$sa_i++)
		{
			if($read[$sa_i]=~/^SA:/)
			{
				$sa_seq=$read[$sa_i];
				$sa_seq=~s/.*://;
				@sa=split(/,/,$sa_seq);
				last;
			}
		}
		$sa[3]=~/(\d+)M/;
		$part_2_len=$1;
		$part_2_start=$sa[1]-1;
		$part_2_end=$sa[1]+$part_2_len-1;
		$part_2_seq=`echo  "$sa[0]\t$part_2_start\t$part_2_end"|bedtools getfasta -fi /DeepAlu/DeepAlu_model/reference/hs37d5.fa -bed /dev/stdin|grep -v "^>"`;
		$part_1_strand='+';
		$part_2_strand=$sa[2];
		if($read[1]%32>15)
		{
			$part_1_strand='-';
		}
		if($part_1_strand ne $part_2_strand)
		{
			$part_2_seq=reverse($part_2_seq);
			$part_2_seq=&comp_seq($part_2_seq);
		}
		$overlap_len=$part_2_len-$part_1_len+1;
		$part_2_seq=substr($part_2_seq,$overlap_len);
		$seq_read=`samtools view /DeepAlu/download_bam/Fudan_DNA_D5_1.recal.bam $sa[0]:$part_2_start\-$part_2_end 2>/dev/null|grep $read[0] |cut -f10`;
		print "start###############\n";
		print "$seq_read\n$part_1_strand\t$read[5]\t$sa[3]\n$read[9]  $part_2_seq\n";
		print "end###############\n";
		#print "$part_2_seq\n";
	}
}
	
sub comp_seq
{
        (my $seq)=@_;
	$seq=~s/A/M/g;
	$seq=~s/T/A/g;
	$seq=~s/M/T/g;
	$seq=~s/C/M/g;
	$seq=~s/G/C/g;
	$seq=~s/M/G/g;
	return $seq;
}

