while(<STDIN>)
{
	$name=$_;
	$seq=<STDIN>;
	$sig=<STDIN>;
	$qual=<STDIN>;
	chomp($name);
	chomp($seq);
	chomp($sig);
	chomp($qual);
	$seq=reverse($seq);
	$qual=reverse($qual);
	print "$name\n";
	print "$seq\n";
	print "$sig\n";
	print "$qual\n";
}
