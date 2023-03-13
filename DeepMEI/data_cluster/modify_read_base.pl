#HISEQ1:19:H8VDAADXX:2:2212:10501:78842  161     7       46282872        39      139M9S  hs37d5  14660836        0       TGCAACCACTTTAGAAAACAGGTTGGCCATTTCTGGCGTGAACCCGGGAGGCGGAGCTTGCAGTGAGCTGAGATTGCGCCACTGCACTCCAGCTTGGGCAACAGAGCGAGACTCCGTCTCAAAAAAAAAAAAAAAAAAAAAAAAAAAA    @@@DDDDDHFBBFEHIIIEHECA<A@EEGGFGIGGGGIEBGIEHGGHBHEEED@B=BBCCBCCCCCCCCCCBCCCCCBBBBBCCCCACCCCCCCCCBBBCBBBBCCCBBBBBCCCBBBBCDCCCBBBBBBBBBBBBBBBBBBBBBB>B    NM:i:3  MD:Z:26T66C21A23        MC:Z:148M
while(<STDIN>)
{
	chomp();
	@F=split(/\t/,$_);
	if($F[5]=~/[DI]/ or $_=~/^@/)
	{
		print "$_\n";
		next;
	}
	@md="";
	for($i=0;$i<=$#F;$i++)
	{
		if($F[$i]=~/MD/)
		{
			$md_tmp=$F[$i];
			$md_tmp=~s/MD:Z://;
			$md_tmp=~s/([A-Z])/\t\1\t/g;
			@md=split(/\t/,$md_tmp);
			last;
		}
	}
	if($#md eq 0)
	{
		print "$_\n";
		next;
	}
	$seq_t="";
	$seq=$F[9];
	$start=0;
	if($F[5]=~/^(\d+)S/)
	{
		$seq_t=substr($seq,0,$1);
		$start=$1;
	}
	for($i=0;$i<$#md;$i=$i+2)
	{
		$seq_t=$seq_t.substr($seq,$start,$md[$i]);
		$seq_t=$seq_t."$md[($i+1)]";
		$start=$start+$md[$i]+1;
	}
	$seq_t=$seq_t.substr($seq,$start);
	$F[9]=$seq_t;
	print join("\t",@F);
	print "\n";
}
