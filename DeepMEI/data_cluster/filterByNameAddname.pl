#10:101340206-101340306	HISEQ1:20:H9V1RADXX:2:1204:9388:79215	1	15	76120636	76120638	left	
#10:101340206-101340306	HISEQ1:18:H8VC6ADXX:2:1209:9515:34168	1	15	76120631	76120633	left

#HISEQ1:20:H9V1RADXX:2:1214:8847:71161   83      CCAAGCCCGAACCCGAACCTGAACAGTAACCGTAGTCCAAACCTGTACCCTTCCGATATCCCTTACCCTTACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAA
#HISEQ1:18:H8VC6ADXX:1:1105:1751:76902   387     GATAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAAC      HHHHHHGFBDFCED@CCBDDDDDB?<ADDDD<ACDDDDBDCDDB?CCDABDD<BBDBCC8AA?C

%RE;
%REG;
%REGF;
%REGD;
@name=`cat $ARGV[0] `;
chomp(@name);
foreach $i (@name)
{
#10:101340234-101340334	HISEQ1:20:H9V1RADXX:1:2211:5155:41672	1
	@dis=split(/\t/,$i);
	$RE{"$dis[1]"}=1;
	$REG{"$dis[1]"}=$dis[0];
	$REGF{"$dis[1]"}=$dis[2];
	$REGD{"$dis[1]"}=$dis[6];
}
while(<STDIN>)
{
	chomp();
	@F=split(/\t/,$_);
	$name=$F[0];
	$flag=$F[1];
	$seq=$F[9];
	$qual=$F[10];
	$read_r=1;
	if($flag%256 >=128)
	{
		$read_r=2;
	}
	if($RE{"$name"} ==1)
	{
		if($REGF{"$name"} != $read_r)
		{
			$seq=&extract_high_quality_bases($seq, $qual, $flag);
			if(($flag%32>=16 and $REGD{"$name"} eq "right") or ($flag%32<16 and $REGD{"$name"} eq "left"))
			{
				$seq=&get_reverse_complement($seq);
			}
			if(length($seq)>15)
			{
				print ">".$REG{"$name"}."_".$REGD{"$name"}."\n$seq\n";
			}
			$RE{"$name"}=0;
		}
	}
}


# 解析 SAM flags
sub parse_sam_flags {
    my ($flags) = @_;
    return (
        paired     => $flags & 0x1,
        proper     => $flags & 0x2,
        unmapped   => $flags & 0x4,
        mate_unmapped => $flags & 0x8,
        reverse    => $flags & 0x10,
        mate_reverse => $flags & 0x20,
        read1      => $flags & 0x40,
        read2      => $flags & 0x80,
        secondary  => $flags & 0x100,
        qc_fail    => $flags & 0x200,
        duplicate  => $flags & 0x400,
        supplementary => $flags & 0x800,
    );
}

# 质量检查函数
sub quality_check {
    my ($index, $qualities, $threshold) = @_;
    my $start = $index - 5 >= 0 ? $index - 5 : 0;
    my $end = $index + 5 < length($qualities) ? $index + 5 : length($qualities) - 1;

    my $low_quality_count = 0;
    for my $i ($start .. $end) {
        my $quality_score = ord(substr($qualities, $i, 1)) - 33;  # 转码
        $low_quality_count++ if $quality_score < $threshold;
    }
    return $low_quality_count < 5;
}

# 提取高质量碱基序列
sub extract_high_quality_bases {
    my ($sequence, $qualities, $flags, $threshold) = @_;
    $threshold ||= 30;

    my %parsed_flags = parse_sam_flags($flags);
    my @high_quality_bases;

    if (($parsed_flags{read1} && !$parsed_flags{reverse}) || ($parsed_flags{read2} && $parsed_flags{reverse})) {
        # Read 1 on forward strand or Read 2 on reverse strand: left-to-right
        for my $i (0 .. length($sequence) - 1) {
            last if !quality_check($i, $qualities, $threshold);
            push @high_quality_bases, substr($sequence, $i, 1);
        }
    } elsif (($parsed_flags{read2} && !$parsed_flags{reverse}) || ($parsed_flags{read1} && $parsed_flags{reverse})) {
        # Read 2 on forward strand or Read 1 on reverse strand: right-to-left
        for my $i (reverse 0 .. length($sequence) - 1) {
            last if !quality_check($i, $qualities, $threshold);
            unshift @high_quality_bases, substr($sequence, $i, 1);
        }
    } else {
        die "Invalid read orientation.";
    }

    return join('', @high_quality_bases);
}


# 函数定义：获取反向互补序列
sub get_reverse_complement {
    my ($sequence) = @_;

    # 创建碱基互补对的哈希表
    my %complement = (
        'A' => 'T',
        'T' => 'A',
        'C' => 'G',
        'G' => 'C',
        'N' => 'N',
        'X' => 'X',
        'a' => 't',  # 支持小写字母
        't' => 'a',
        'n' => 'n',
        'c' => 'g',
        'x' => 'x',
        'g' => 'c'
    );

    # 将序列中的每个碱基替换为其互补碱基
    my $complement_sequence = join('', map { $complement{$_} } split('', $sequence));

    # 反转互补序列
    my $reverse_complement_sequence = reverse($complement_sequence);

    return $reverse_complement_sequence;
}
