sample='HG002'
file=$1
rm -rf regions split_softclipped split_softclipped_sort
#rm -rf split_softclipped split_softclipped_sort

mkdir regions split_softclipped split_softclipped_sort
input_gt=input_gt.txt
#sample='HG002';cat $input_gt|perl -F'\t' -alne '$start=$F[1]-50;$start2=$F[1]-100;;$end=$F[1]+50;$end2=$F[1]+100;print "$F[0]:$start-$end\t";' |xargs -n 1 -P 20 -I {}  samtools view $file  {} -o regions/${sample}_${record}.bam
samtools view -H $file >head.sam
#sample='HG002';cat $input_gt|perl -F'\t' -alne '$start=$F[1]-50;$start2=$F[1]-100;;$end=$F[1]+50;$end2=$F[1]+100;print "$F[0]:$start-$end\t$F[0]:$start2-$end2";' |while read record region;do  echo $record;samtools view $file  $record -O SAM -o regions/${sample}_${record}.bam;done
sample='HG002';cat input_gt.txt|perl -F'\t' -alne '$start=$F[1]-50;$start2=$F[1]-100;;$end=$F[1]+50;$end2=$F[1]+100;print "$F[0]:$start-$end";' |xargs -n 1 -P 20 -I record bash extract_region.sh record $file $sample

##遍历所有bam文件，找出切断的reads，并保存到fasta文件
##ls regions|xargs -n 1 -P 40 perl alu_discord_support_part1.pl 
##rm clip_right_seq.fa clip_left_seq.fa
##cd split_softclipped/
##ls |sed -n '/mapClipR_long/p'|while read file;do cat $file >>../clip_right_seq.fa;done
##ls |sed -n '/mapClipL_long/p'|while read file;do cat $file >>../clip_left_seq.fa;done
##RepeatMasker对所有提取出fasta文件进行注释
##cd ..
##conda activate repeatmasker
##RepeatMasker -no_is -alu -species human -a  -e rmblast -excln clip_left_seq.fa 
##RepeatMasker -no_is -alu -species human -a  -e rmblast -excln clip_right_seq.fa
##拆分RM注释的文件到每个位点
##cat clip_left_seq.fa.out |perl -npe "s/ +/\t/g"|perl -npe "s/^\t//"|tail --lines=+4|cut -f5|perl -npe "s/:(\d+)$/\t\1/"|while read file id;do echo $id >>split_softclipped/${file}.rp_left;done
##cat clip_right_seq.fa.out |perl -npe "s/ +/\t/g"|perl -npe "s/^\t//"|tail --lines=+4|cut -f5|perl -npe "s/:(\d+)$/\t\1/"|while read file id;do echo $id >>split_softclipped/${file}.rp_right;done
ls regions |while read record;do bash discover_bam_to_candidate.sh $record test $file ;done 
##根据RM结果和polyA/T信号，对所有reads进行分组输出
ls regions|xargs -n 1 -P 40 perl alu_discord_support_part2.pl 
#ls regions|while read file;do perl alu_discord_support_part2.pl $file ;echo $file ;done 
cd split_softclipped
#mkdir ../split_softclipped_sort
#对每种sam文件进行排序
ls |sed -n '/mapClipL.sam$/p'|while read record;do cat <(cat ../head.sam) <(grep -v "^@" $record |perl -F'\t' -alne '$start=$F[3]-$F[7];print "$start\t$_"'|sort -k1,1n |cut -f2- ) > ../split_softclipped_sort/$record;done

ls |sed -n '/mapClipR.sam$/p'|while read record ;do cat <(cat ../head.sam) <(grep -v "^@" $record |sort -k4,4n ) > ../split_softclipped_sort/$record;done

ls |sed -n '/mapRef.sam$/p'|while read record ;do cat <(cat ../head.sam) <(grep -v "^@" $record |sort -k4,4n ) > ../split_softclipped_sort/$record;done

#ls |sed -n '/mapClipL/p'|sed -n '/HG002/p'|while read file ;do cat <(cat ../head.sam) <(grep -v "^@" $file |perl -F'\t' -alne '$start=$F[3]-$F[7];$F[2]="chr$F[2]";print "$start\t".join("\t",@F);'|sort -k1,1n |cut -f2- ) > ../split_softclipped_sort/$file;done
#ls |sed -n '/mapClipL/p'|sed -n '/HG002/p'|while read file ;do cat <(cat ../head.sam) <(grep -v "^@" $file |perl -F'\t' -alne '$start=$F[3]-$F[7];$F[2]="chr$F[2]";print "$start\t".join("\t",@F);'|sort -k1,1n |cut -f2- ) > ../split_softclipped_sort/$file;done
#
#ls |sed -n '/mapRef/p'|sed -n '/HG002/p'|while read file ;do cat <(cat ../head.sam) <(grep -v "^@" $file |perl -F'\t' -alne '$F[2]="chr$F[2]";print join("\t",@F);'|sort -k4,4n ) > ../split_softclipped_sort/$file;done
