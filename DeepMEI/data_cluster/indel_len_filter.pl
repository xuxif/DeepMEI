#HISEQ1:18:H8VC6ADXX:1:1213:2189:87579   99      1       10001   9       105M1I42M       =       10177   333     TAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCCAACCCCAACCCTAACCCTAACCCCAACCCTAACCCCAACCCAAAC    ???DDBDDDDDDDIDEE:C3ACEEEEI3EEDDDDIID@D>BCADID@BADI<6=@AC9AD;?C6;?(;AAA>6>>(,(,5=92=(5958<(9<<<>(39<<808<<99&(2()91991<8((3+21)(2(&&(+9(((&+2(&0((+(    NM:i:5  MD:Z:102T22T11T5T3      MC:Z:54M1D16M8D78M      AS:i:121        XS:i:124
$bam_file=$ARGV[0];
$REF=$ARGV[1];
$record=$ARGV[2];

@record=split(/[:\-]/,$record);

$chr=$record[0];
$start=$record[1];
$end=$record[2];

$indel_read_count=0;

@reads=`samtools view -T $REF $bam_file $ARGV[2]`;
foreach $read (@reads)
{
	chomp($read);
	@read_detail=split(/\t/,$read);
	if($read_detail[5]=~/[DI]/)
	{
		$read_detail[5]=~s/([A-Z])/\t\1\t/g;
		$read_detail[5]=~s/^\t+//;
		$read_detail[5]=~s/\t+$//;
		#print "$read_detail[5]\n";
		@cigar=split(/\t/,$read_detail[5]);
		$match_len=0;
		for($i=0;$i<@cigar;$i=$i+2)
		{
			if(($cigar[$i+1] eq 'I') and $cigar[$i] >4 and ($read_detail[3]+$match_len)<$end+5 and  ($read_detail[3]+$match_len)>$start-5)
			{
				#				print (($read_detail[3]+$match_len)-($start+$end)/2);
				#				print "\n";
				$indel_read_count++
			}
			if(($cigar[$i+1] eq 'D') and $cigar[$i] >4 and ($read_detail[3]+$match_len)<$end+5 and  ($read_detail[3]+$match_len)>$start-5)
			{
				#print (($read_detail[3]+$match_len)-($start+$end)/2);
				#print "\n";
				$indel_read_count++
			}
			if(($cigar[$i+1] eq 'D') and $cigar[$i] >4 and ($read_detail[3]+$match_len+$cigar[$i])<$end+5 and  ($read_detail[3]+$match_len+$cigar[$i])>$start-5)
			{
				#print (($read_detail[3]+$match_len)-($start+$end)/2);
				#print "\n";
				$indel_read_count++
			}
			if($cigar[$i+1] eq 'M' or $cigar[$i+1] eq 'D')
			{
				$match_len=$match_len+$cigar[$i];
			}
		}
	}

}
if($indel_read_count==0)
{
	print "$chr\t$start\t$end\n";
}
#else
#{
#	print "1\n";
#}
