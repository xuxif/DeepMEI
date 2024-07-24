# 提取所有的 read pairs
REF=$1
input_bam=$2
ran_num=$3
record=$4

echo $record |perl -npe 's/(.*):(\d+)-(\d+)/\1 \2 \3/' |while read chrom start end;do 

out_dir=regions_$ran_num
samtools view -h -F 16 -b $input_bam $chrom:$((start-500))-$((start-50))  > $out_dir/${record}_left_region_positive.bam

samtools view -h -f 16 -b $input_bam $chrom:$((end+50))-$((end+500))  > $out_dir/${record}_right_region_negative.bam

cat <(samtools view -H $input_bam) \
	<(cat \
		<(samtools view  -F 2 $out_dir/${record}_left_region_positive.bam ) \
		<(samtools view  -F 2 $out_dir/${record}_right_region_negative.bam ) \
		<(samtools view  $out_dir/${record}_left_region_positive.bam | awk '$9 > 1000 || $9 < -1000 || $9 == 0'  ) \
		<(samtools view  $out_dir/${record}_right_region_negative.bam | awk '$9 > 1000 || $9 < -1000 || $9 == 0' )  \
		| sort -k1,1 -k2,2n |uniq) \
 >$out_dir/${record}_discordant.sam

# 提取 discordant reads 的 read names 和 mate 的比对位置
read_name=`samtools view $out_dir/${record}_discordant.sam |perl -F'\t' -alne 'if($F[6] eq "=") {$F[6]=$F[2];} $read_name=$read_name." $F[0]";$region=$region." $F[6]:$F[7]\-".($F[7]+1);if(eof()){$read_name=~s/^ //;print "$read_name";}'`
region=`samtools view $out_dir/${record}_discordant.sam |perl -F'\t' -alne 'if($F[6] eq "=") {$F[6]=$F[2];} print "$F[6]\t$F[7]";' |pos2bed 100 |sort -k1,1 -k2,2n |bedtools merge|perl -npe "s/\t/:/;s/\t/-/;s/\n/ /"|perl -npe "s/$/\n/"`

 samtools view --threads 2 -b $input_bam $region \
	| samtools view -O BAM  -N <(echo "$read_name"|perl -npe "s/ /\n/g" )  \
	|samtools fastq  /dev/stdin  -o $out_dir/${record}_discordant_and_mate_reads_R1.fastq 2>/dev/null


samtools faidx $REF $chrom:$((pos-200))-$((pos+200)) > $out_dir/${record}_insertion_site_region.fasta
bwa index $out_dir/${record}_insertion_site_region.fasta 2>/dev/null 

cat \
	<(samtools view -H $input_bam) \
	<(bwa mem $out_dir/${record}_insertion_site_region.fasta $out_dir/${record}_discordant_and_mate_reads_R1.fastq 2>/dev/null \
	| samtools view -F 4 |perl -F'\t' -alne '$F[2]=~/(.*):(\d+)-(\d+)/;$F[2]=$1;$F[3]=$2+$F[3];print join("\t",@F);' |sort -k1,1 -k2,2n |uniq )  \
 >$out_dir/${record}_insertion_dis_align.sam

done
