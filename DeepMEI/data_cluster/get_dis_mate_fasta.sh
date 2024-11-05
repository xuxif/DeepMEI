bam_file=$1
REF=$2
file=$3
#samtools view  -T $REF $bam_file `cat $file|cut -f4-6|bedtools merge|perl -npe "s/\t/:/;s/\t/-/;s/\n/ /"|perl -npe "s/$/\n/"`|samtools fasta |perl ../filterByNameAddname.pl $file > fa_$file
samtools view  -T $REF $bam_file `cat $file|cut -f4-6|bedtools merge|perl -npe "s/\t/:/;s/\t/-/;s/\n/ /"|perl -npe "s/$/\n/"`|perl ../filterByNameAddname.pl $file > fa_$file

