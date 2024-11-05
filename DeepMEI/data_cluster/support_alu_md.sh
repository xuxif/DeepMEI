#makeblastdb -in ME_ref.fa -input_type fasta -dbtype nucl
#rm blastout/*
#rm support_alt_reads/*
#rm suppport_alu.tar.gz
dir=regions
#dir=NA12878_validate
#ls regions/|while read file;do
file=$1
	out_file=`echo $file|perl -npe "s/\.bam//"`
	samtools fasta $dir/$file  | \
	blastn -db ME_ref.fa -query /dev/stdin -evalue 0.0001 -outfmt "6 qseqid"  -out /dev/stdout |perl -npe "s/\/[12]$//"|sort |uniq >blastout/${out_file}.txt
	comm -23 <(samtools view $dir/$file|cut -f1,6|grep "S"|cut -f1|sort) <(cat blastout/${out_file}.txt)  |while read id
	do 
		echo -ne "$id\t";
		samtools view $dir/$file |grep "$id" |cut -f6,10|grep "S"|head --lines=1|perl -F'\t' -alne 'if($F[0]=~/^\d+S/) { $F[0]=~s/^(\d+)S.*/\1/;$seq=substr($F[1],0,$F[0]);} elsif($F[0]=~/\d+S$/) { $F[0]=~s/.*?(\d+)S$/\1/;$seq=substr($F[1],-$F[0]);} if($seq=~/AAAAAAAAAA/ or $seq =~/TTTTTTTTTT/) {print "T";} else { print "F";} '  
	
	done |grep "T$"|cut -f1 >support_alt_reads/${out_file}.txt
	cat blastout/${out_file}.txt >>support_alt_reads/${out_file}.txt
	rm blastout/${out_file}.txt
