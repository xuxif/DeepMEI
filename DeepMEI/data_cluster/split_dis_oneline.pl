#>10:101340206-101340306_left    ACTGGGATTACAGGCATGCACCACCATGCCCAGTTAATTTTGTATTTTTAGTAGAGATGGGGTTTCACCATGTTGGTCAGGCTGGTCTTAAACTCCTGACCTCAAGTGATCCACCCACTTCAGCCTCCCAAAATGCTGGGATTTCAGG
%RE;
$dir=$ARGV[0];
while(<STDIN>)
{
	chomp();
	@F=split(/\t/,$_);
	$F[0]=~s/_left$//;
	$F[0]=~s/_right$//;
	$F[0]=~s/^>//;
	$RE{$F[0]}=$F[1]."\n$RE{$F[0]}";
}
foreach $key (keys %RE)
{
	open FP,">$dir/ass_seq_$key";
	print FP "$RE{$key}";
	close FP;
}


