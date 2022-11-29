sample='HG002'
file='/root/HG002.sort.bam'
file='/mnt/disk2/GIAB/AshkenazimTrio/HG002.hs37d5.60x.1.bam'
vcf=/home/xuxf/DeepAlu/final_vcf/true_new.vcf 
rm -rf regions regions_${sample}_combine
mkdir regions_${sample}_combine
#wget -c ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG002.hs37d5.60x.1.bam -O HG002.hs37d5.60x.1.bam
#wget ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG002.hs37d5.60x.1.bam.bai -O HG002.hs37d5.60x.1.bam.bai
#sample='HG002';cat input_gt_HG002.txt |cut -f2,3|perl -F'\t' -alne '$start=$F[1]-50;$end=$F[1]+50;print "$F[0]:$start-$end";' |while read record;do samtools view $file -O BAM $record|samtools reheader  head.sam  /dev/stdin >regions_${sample}/${sample}_${record}.bam;done 
samtools view -H $file >head.sam 
bp_var=-50
while [[ $bp_var -le 50 ]]
do
rm -rf regions
mkdir regions
cat ${vcf} |grep -v "^#"|cut -f1,2|perl -npe "s/$/\t${bp_var}/"|perl -F'\t' -alne '$start=$F[1]-50-$F[2];$end=$F[1]+50-$F[2];print "$F[0]:$start-$end";' |while read record;do samtools view $file -b $record|samtools reheader  head.sam  /dev/stdin >regions/${sample}_${record}.bam;echo ${sample}_${record};done

#遍历所有bam文件，找出切断的reads，并保存到fasta文件
rm -rf split_softclipped
mkdir split_softclipped

ls regions|cut -d'.' -f1|while read file;do perl alu_discord_support_part1.pl $file ;echo $file ;done 

rm clip_right_seq.fa clip_left_seq.fa
rm -rf clip_right_seq.fa.* -rf clip_left_seq.fa.*
cd split_softclipped/
ls |sed -n '/mapClipR_long/p'|while read file;do cat $file >>../clip_right_seq.fa;done
ls |sed -n '/mapClipL_long/p'|while read file;do cat $file >>../clip_left_seq.fa;done
#RepeatMasker对所有提取出fasta文件进行注释
cd ..
RepeatMasker -no_is -alu -species human -a  -e rmblast -excln clip_left_seq.fa 
RepeatMasker -no_is -alu -species human -a  -e rmblast -excln clip_right_seq.fa
#拆分RM注释的文件到每个位点
cat clip_left_seq.fa.out |perl -npe "s/ +/\t/g"|perl -npe "s/^\t//"|tail --lines=+4|cut -f5|perl -npe "s/:(\d)+$/\t\1/"|while read file id;do echo $id >>split_softclipped/${file}.rp_left;done
cat clip_right_seq.fa.out |perl -npe "s/ +/\t/g"|perl -npe "s/^\t//"|tail --lines=+4|cut -f5|perl -npe "s/:(\d)+$/\t\1/"|while read file id;do echo $id >>split_softclipped/${file}.rp_right;done
#根据RM结果和polyA/T信号，对所有reads进行分组输出
ls regions|cut -d'.' -f1|while read file;do perl alu_discord_support_part2.pl $file ;echo $file ;done 
cd split_softclipped
rm -rf split_softclipped_sort
mkdir ../split_softclipped_sort
#对每种sam文件进行排序
ls |sed -n '/mapClipL.sam$/p'|while read file ;do cat <(cat ../head.sam) <(grep -v "^@" $file |perl -F'\t' -alne '$start=$F[3]-$F[7];print "$start\t$_"'|sort -k1,1n |cut -f2- ) > ../split_softclipped_sort/$file;done

ls |sed -n '/mapClipR.sam$/p'|while read file ;do cat <(cat ../head.sam) <(grep -v "^@" $file |sort -k4,4n ) > ../split_softclipped_sort/$file;done

ls |sed -n '/mapRef.sam$/p'|while read file ;do cat <(cat ../head.sam) <(grep -v "^@" $file |sort -k4,4n ) > ../split_softclipped_sort/$file;done
mv split_softclipped_sort split_softclipped_sort_${bp_var}
tar -czf split_softclipped_sort_${bp_var}.tar.gz split_softclipped_sort_${bp_var}
((bp_var=bp_var+1));
done


#ls |sed -n '/mapClipL/p'|sed -n '/HG002/p'|while read file ;do cat <(cat ../head.sam) <(grep -v "^@" $file |perl -F'\t' -alne '$start=$F[3]-$F[7];$F[2]="chr$F[2]";print "$start\t".join("\t",@F);'|sort -k1,1n |cut -f2- ) > ../split_softclipped_sort/$file;done
#ls |sed -n '/mapClipL/p'|sed -n '/HG002/p'|while read file ;do cat <(cat ../head.sam) <(grep -v "^@" $file |perl -F'\t' -alne '$start=$F[3]-$F[7];$F[2]="chr$F[2]";print "$start\t".join("\t",@F);'|sort -k1,1n |cut -f2- ) > ../split_softclipped_sort/$file;done
#
#ls |sed -n '/mapRef/p'|sed -n '/HG002/p'|while read file ;do cat <(cat ../head.sam) <(grep -v "^@" $file |perl -F'\t' -alne '$F[2]="chr$F[2]";print join("\t",@F);'|sort -k4,4n ) > ../split_softclipped_sort/$file;done
