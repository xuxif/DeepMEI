use warnings;
use List::Util qw/sum/;
use POSIX;
my @line=();
my @name=();
my $as="";
#soft|HISEQ1:13:H8G92ADXX:2:2116:10736:89620|147|1|72988040|60|119M29S|*|72988159|-119|*|*
#soft|HISEQ1:20:H9V1RADXX:2:2207:13925:91280|163|12|10732293|60|15S133M|*|10732293|358|mate_chr|mate_pos       0       L1M2_5end       2785    8       15M     *       0       0       TAAAACAACCAGAAA CCCFFFFFHHHHGJJ AS:i:-5 XN:i:0  XM:i:1  XO:i:0  XG:i:0  NM:i:1  MD:Z:0C14       YT:Z:U
$ran_num=rand(10000);
open DIS,">tmp_$ARGV[1]/dis_${ran_num}.tsv";
while(<STDIN>)
{
	if($_=~/^@/){ next;};
	@line=split;
	@name=split /\|/, $line[0];
	$as="NA";
	$chr=$name[3];
	$clip_pos=$name[4];

        $name[6]=~s/([A-Z])/\t$1\t/g;
	$name[6]=~s/\t$//;
	@hg_cigar=split(/\t/,$name[6]);
	$name[6]=~s/\t//g;
	#$alu_cigar_len=length($line[9]);
        $clip_part='none';
        if($hg_cigar[1]=~/^\d+S/) # and $hg_cigar[0]==$alu_cigar_len)
        {
                $clip_part='left';
        }
        else #($hg_cigar[$#hg_cigar]=~/S/) # and $hg_cigar[($#hg_cigar-1)]==$alu_cigar_len)
        {
                $clip_part='right';
                $clip_pos=$name[8];
        }
	
#if ($tmp1[2]%4 >=2
	unless($line[2] eq "*" )
	{
		#		if( $name[2]%4 >=2)
		#		{
			for(my $i=11;$i<@line;$i++)
			{
				if($line[$i]=~s/AS:i://){$as=$line[$i];}
			}
#			if($as eq "NA" ){next;}
			if($as eq "NA" || ($as<20 and $as >0) || $as<-15 ){next;}
			print "$chr\t$clip_pos\t$name[6]\t$line[2]\t$line[3]\t$line[4]\t$line[5]\t$as\t".abs($name[9])."\n";
			#		}
	}
	elsif($name[9] eq 0 and $name[10] ne "*")
	{
			print DIS "$name[10]\t$name[11]\t".($name[11]+50)."\t$chr\t$clip_pos\t$name[6]\t$line[2]\t$line[3]\t$line[4]\t$line[5]\t$as\t".abs($name[9])."\n";
	}
	#	elsif($line[9]=~/TTTTTTTTTT/ or $line[9]=~/AAAAAAAAAA/)
	#	{
	#		print "$chr\t$clip_pos\t$name[6]\t$line[2]\t$line[3]\t$line[4]\t$line[5]\t$as\tployAT\n";
	#	}
}
close(DIS);
print `bash ins_can.sh tmp_$ARGV[1]/dis_${ran_num}.tsv $ARGV[0] $ARGV[2] $ARGV[3]`;
#`rm dis_${ran_num}.tsv`
