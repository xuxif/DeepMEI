sample='HG002'
bam_file=$1
input_gt=$2
REF=$3
output=`echo $bam_file|perl -npe "s/.*\///;s/\.final.cram//;s/\.bam//"`
ran_num=$output
if [[  -d "regions_$ran_num" ]]
then
	rm -rf regions_$ran_num split_softclipped_$ran_num split_softclipped_sort_$ran_num head_$ran_num.sam
fi
echo "$ran_num"
mkdir  split_softclipped_sort_$ran_num
samtools view -H $bam_file >head_$ran_num.sam
sample='HG002';cat $input_gt|perl -F'\t' -alne '$start=$F[1]-50;$start2=$F[1]-100;;$end=$F[1]+50;$end2=$F[1]+100;print "$F[0]:$start-$end\t$F[4]";' |xargs -n 2 -I{} -P 40 perl alu_discord_support_part2_sample_refine.pl {} $ran_num $bam_file $REF >${output}_vsoft.bed
tar -czf split_softclipped_sort_${output}.tar split_softclipped_sort_${output}/
