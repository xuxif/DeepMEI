open RE,"$ARGV[0]";
while(<RE>)
{
	chomp();
	($record,$region_t,$name)=split(/\t/,$_);
	@region=split(/[\-:]/,$region_t);
	$RENAME{$name.$region[1]}=$record;
}
while(<STDIN>)
{
	chomp();
	($name,$pos,$seq)=split(/\t/,$_);
	if(exists $RENAME{$name.$pos})
	{ 
		print "$RENAME{$name.$pos}\t$seq\n";
		#delete $RENAME{$name};
	}
}
