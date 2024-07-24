REF=$1
out_dir=$2
record=$3

cat <(samtools view -H $out_dir/dis_read.bam) <(samtools view $out_dir/dis_read.bam  $record |perl -npe "s/^/$record\t/"|perl -F'\t' -alne '$F[0]=~/(.*):(\d+)-(\d+)/;$chrom=$1;$pos=$2+500;if(abs($pos-$F[4])>550 or  $chrom ne $F[7]){print "$_";}' |cut -f2-) \
	|samtools fastq /dev/stdin >$out_dir/${record}.fastq  2>/dev/null 
count=`cat $out_dir/${record}.fastq|wc -l`
if [[ $count -ge 12 ]]
then
	samtools faidx $REF $record > $out_dir/${record}.fasta
	bwa index $out_dir/${record}.fasta 2>/dev/null
	bwa mem $out_dir/${record}.fasta $out_dir/${record}.fastq 2>/dev/null \
	    | samtools view -F 4 |perl -F'\t' -alne '$F[2]=~/(.*):(\d+)-(\d+)/;$F[2]=$1;$F[3]=$2+$F[3];print join("\t",@F);'  \
	    >$out_dir/${record}_bwa.sam
fi
