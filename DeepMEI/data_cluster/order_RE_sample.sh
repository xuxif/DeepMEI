sample='HG002'
bam_file=$1
ran_num=$2
input_gt=$3
ME_REF=$4
REF=$5
step=$6
if [[  -d "regions_$ran_num" ]]
then
	rm -rf regions_$ran_num split_softclipped_$ran_num split_softclipped_sort_$ran_num head_$ran_num.sam
fi

mkdir regions_$ran_num split_softclipped_$ran_num split_softclipped_sort_$ran_num
#samtools view -T $REF -H $bam_file|grep -f chr_list.txt  >head_$ran_num.sam
samtools view -T $REF -H $bam_file  >head_$ran_num.sam
cd regions_$ran_num

sort -k1,1 -k2,2n $input_gt >input_sort.bed
num=`cat input_sort.bed |wc -l|perl -F'\t' -alne '$num=int($F[0]/2);print "$num";'`
paste -d'\n' <(head -n $num input_sort.bed ) <(tail -n $num input_sort.bed ) |split -l 10000 /dev/stdin shuff_split_candidate

#ls split_candidate*|while read file;do  num=`cat $file|wc -l|perl -F'\t' -alne '$num=int($F[0]/2);print "$num";'`;paste -d'\n' <(sort -k1,1 -k2,2n $file|head -n $num ) <(sort -k1,1 -k2,2n $file|tail -n $num ) |head -n `cat $file|wc -l` >shuff_$file;done
ls shuff_split_candidate*|xargs -n 1 -P 20 -I{} bash ../extrac_region_batch.sh $bam_file $REF {}
cd ..

#sample='HG002';cat $input_gt|perl -F'\t' -alne '$start=$F[1]-50;$start2=$F[1]-100;;$end=$F[1]+50;$end2=$F[1]+100;print "$F[0]:$start-$end\t$F[4]";' |xargs -n 2 -I{} -P 40 perl alu_discord_support_part2_sample_refine.pl {} $ran_num $bam_file $REF >../DeepAlu_script/vsoft_pos/${ran_num}_vsoft.bed
#sample='HG002';cat $input_gt|perl -F'\t' -alne '$start=$F[1]-50;$start2=$F[1]-100;;$end=$F[1]+50;$end2=$F[1]+100;print "$F[0]:$start-$end";' |xargs -n 1 -P 20 -I record bash extract_region_sample.sh record $bam_file $sample regions_$ran_num $REF
sample='HG002';cat $input_gt|perl -F'\t' -alne '$start=$F[1]-50;$start2=$F[1]-100;;$end=$F[1]+50;$end2=$F[1]+100;print "$F[0]:$start-$end\t$F[4]";' |xargs -n 2 -I{} -P 40 perl alu_discord_support_part2_sample_refine_vt.pl {} $ran_num $bam_file $REF $step >../DeepMEI_script/vsoft_pos/${ran_num}_vsoft.bed #${ran_num}_vsoft.bed
echo "generate record finished!"

#ls regions_$ran_num |xargs -n 1 -I {} -P 40 bash discover_bam_to_candidate_sample.sh {} $ran_num $ME_REF

##根据RM结果和polyA/T信号，对所有reads进行分组输出
#ls regions_$ran_num|xargs -n 1 -I{} -P 40 perl alu_discord_support_part2_sample.pl {} $ran_num >${ran_num}_filter.bed
#ls regions_$ran_num|xargs -n 1 -I{} -P 40 perl alu_discord_support_part2_sample_refine_vt.pl {} $ran_num >${ran_num}_filter.bed
#bedtools intersect -a <(cat $input_gt|perl -F'\t' -alne '$start=$F[1]-50;$end=$F[1]+50;print "$F[0]\t$start\t$end\t$_";') -b <(cat ${ran_num}_filter.bed|perl -F'\t' -alne '$start=$F[1]-50;$end=$F[1]+50;print "$F[0]\t$start\t$end\t$_";')  -wa|cut -f4- |sort |uniq >../DeepAlu_script/vsoft_pos/${ran_num}_vsoft.bed
#rm ${ran_num}_filter.bed
cd split_softclipped_$ran_num
#mkdir ../split_softclipped_sort
#对每种sam文件进行排序
#ls |sed -n '/mapClipL.sam$/p'|while read record;do cat <(cat ../head_$ran_num.sam) <(grep -v "^@" $record |perl -F'\t' -alne 'print "$F[5]\t$_";'|perl -npe "s/(\d+)S.*?\t/\1\t/"|sort -k5,5n -k1,1nr|cut -f2- ) > ../split_softclipped_sort_$ran_num/$record;done &
#ls |sed -n '/mapClipL.sam$/p'|while read record;do record_pre_1=`echo $record|perl -npe "s/mapClipL\.sam/BPinfo\.txt/"`;new_name_1=`cat $record_pre_1|perl -npe "s/:.*?\t/\t/"|perl -F'\t' -alne '$pos=$F[1];$start=$pos-50;$end=$pos+50;print "$F[0]:$start-$end";'` ;cat <(cat ../head_$ran_num.sam) <(grep -v "^@" $record |perl -F'\t' -alne 'print "$F[5]\t$_";'|perl -npe "s/(\d+)S.*?\t/\1\t/" |perl -F'\t' -alne '$start=$F[4]-$F[0];print "$start\t$_";'|sort -k1,1n |cut -f3- ) > ../split_softclipped_sort_$ran_num/HG002_${new_name_1}_mapClipL.sam;rm $record ;done &

#ls |sed -n '/mapClipR.sam$/p'|while read record ;do record_pre_2=`echo $record|perl -npe "s/mapClipR.sam/BPinfo\.txt/"`;new_name_2=`cat $record_pre_2|perl -npe "s/:.*?\t/\t/"|perl -F'\t' -alne '$pos=$F[1];$start=$pos-50;$end=$pos+50;print "$F[0]:$start-$end";'` ; cat <(cat ../head_$ran_num.sam) <(grep -v "^@" $record |sort -k4,4n ) > ../split_softclipped_sort_$ran_num/HG002_${new_name_2}_mapClipR.sam;rm $record ;done &

#ls |sed -n '/mapRef.sam$/p'|while read record ;do record_pre_3=`echo $record|perl -npe "s/mapRef.sam/BPinfo\.txt/"`;new_name_3=`cat $record_pre_3|perl -npe "s/:.*?\t/\t/"|perl -F'\t' -alne '$pos=$F[1];$start=$pos-50;$end=$pos+50;print "$F[0]:$start-$end";'` ; cat <(cat ../head_$ran_num.sam) <(grep -v "^@" $record |sort -k4,4n ) > ../split_softclipped_sort_$ran_num/HG002_${new_name_3}_mapRef.sam;rm $record ;done &

#ls |sed -n '/mapClipL.sam$/p'|while read record;do cat <(cat ../head_$ran_num.sam) <(grep -v "^@" $record |perl -F'\t' -alne 'print "$F[5]\t$_";'|perl -npe "s/(\d+)S.*?\t/\1\t/"|sort -k5,5n -k1,1nr|cut -f2- ) > ../split_softclipped_sort_$ran_num/$record;done &

ls |sed -n '/mapClipL.sam$/p'|while read record;do cat <(cat ../head_$ran_num.sam) <(grep -v "^@" $record |perl -F'\t' -alne 'print "$F[5]\t$_";'|perl -npe "s/(\d+)S.*?\t/\1\t/" |perl -F'\t' -alne '$start=$F[4]-$F[0];print "$start\t$_";'|sort -k1,1n |cut -f3-) > ../split_softclipped_sort_$ran_num/$record;rm $record ;done &

ls |sed -n '/mapClipR.sam$/p'|while read record ;do  cat <(cat ../head_$ran_num.sam) <(grep -v "^@" $record |sort -k4,4n ) > ../split_softclipped_sort_$ran_num/$record;rm $record ;done &

ls |sed -n '/mapRef.sam$/p'|while read record ;do cat <(cat ../head_$ran_num.sam) <(grep -v "^@" $record |sort -k4,4n ) > ../split_softclipped_sort_$ran_num/$record;rm $record ;done &
#ls |sed -n '/mapClipL.sam$/p'|while read record;do cat <(cat ../head_$ran_num.sam) <(grep -v "^@" $record |perl -F'\t' -alne 'print "$F[5]\t$_";'|perl -npe "s/(\d+)S.*?\t/\1\t/" |perl -F'\t' -alne '$start=$F[4]-$F[0];print "$start\t$_";'|sort -k1,1n |cut -f3-|perl -F'\t' -alne  '$md=$F[12];$md=~s/.*://;$md=~s/[^ATCG\n]//g;$md=length($md);if($md<=5){print "$_";}' ) > ../split_softclipped_sort_$ran_num/$record;rm $record ;done &

#ls |sed -n '/mapClipR.sam$/p'|while read record ;do  cat <(cat ../head_$ran_num.sam) <(grep -v "^@" $record |sort -k4,4n|perl -F'\t' -alne  '$md=$F[12];$md=~s/.*://;$md=~s/[^ATCG\n]//g;$md=length($md);if($md<=5){print "$_";}' ) > ../split_softclipped_sort_$ran_num/$record;rm $record ;done &

#ls |sed -n '/mapRef.sam$/p'|while read record ;do cat <(cat ../head_$ran_num.sam) <(grep -v "^@" $record |sort -k4,4n|perl -F'\t' -alne  '$md=$F[12];$md=~s/.*://;$md=~s/[^ATCG\n]//g;$md=length($md);if($md<=5){print "$_";}' ) > ../split_softclipped_sort_$ran_num/$record;rm $record ;done &
wait

#cd ../split_softclipped_sort_$ran_num
#ls |sed -n '/mapRef.sam$/p'|perl -npe "s/_mapRef.sam//"|xargs -n 1 -P 20 -I{} perl ../downsample_sam.pl {} $ran_num

