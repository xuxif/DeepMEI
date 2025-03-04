input_file=$1
output=`echo $1|perl -npe "s/.*\///"`
if [[ ! -n "$input_file" ]]
then
	input_file=/DeepAlu/download_bam/Fudan_DNA_D5_1.recal.bam
	output='FudanD5'
fi


BWA=bwa
ME_REF=$2
output=$3
REF=$4
ME_bed=$5
hd_bed=$6
parallel=$7
split_len=15
if [[ -n "$8" ]]
then
	split_len=$8
	output=$output"_$split_len"
	echo "split_len:$split_len"
	echo "output:$output"
fi
extractSC=SE-MEI/extractSoftclipped

rm -rf tmp_$output 
mkdir tmp_$output
samtools view -T $REF -H $input_file |grep "SN:" |cut -f2,3|perl -npe "s/SN://;s/LN://"|perl -npe "s/\t/\t1\t/"|grep -f chr_list.txt|bedtools makewindows -b /dev/stdin -w 50000000 |perl -npe "s/\t/:/;s/\t/\-/"|xargs  -I{} -P $parallel bash extract_candidate.sh $input_file $split_len $ME_REF $parallel {} $REF ${output} $ME_bed
cat tmp_$output/soft_candidate_${output}*.txt >candidate_pos/${output}_soft_candidate.txt
#cat tmp_$output/bk_*.tsv  >candidate_pos/${output}_bk.txt
rm -rf tmp_$output 
#$extractSC -l $split_len $input_file > /dev/stdout | \
#$BWA mem -t $parallel -k 19 -r 1.5 -c 100000 -m 50 -T 20 -h 10000 -a -Y -M $ME_REF /dev/stdin  > /dev/stdout | \
#perl soft_format_polyAT.pl > candidate_pos/${output}_soft_candidate.txt
#cat candidate_pos/${output}_soft_candidate.txt | \
#sort -k1,1 -k2,2n candidate_pos/${output}_soft_candidate.txt |cut -f1,2,4|grep -f chr_list.txt|perl -F'\t' -alne '$start=$F[1]-1;$end=$F[1];$F[2]=~s/_.*//;print "$F[0]\t$start\t$end\t$F[2]";' |bedtools merge -d 20 -c 4,4 -o collapse,count|perl add_sr_count.pl  >vsoft_pos/${output}_raw_vsoft.bed

#sort -k1,1 -k2,2n candidate_pos/${output}_soft_candidate.txt |cut -f1,2,4|grep -f chr_list.txt |perl -F'\t' -alne '$start=$F[1]-1;$F[2]=~s/_.*//;print "$F[0]\t$start\t$F[1]\t$F[2]";'|bedtools merge -d 2 -c 4,4 -o collapse,count| perl -F'\t' -alne 'if($F[4]>1) {my %te;@tes=split(/,/,$F[3]);foreach $i (@tes) { $te{$i}++;}; if($te{"ALU"}>=$te{"LINE1"} and $te{"ALU"}>=$te{"SVA"} ){print "$F[0]\t$F[1]\t$F[2]\tALU";} elsif($te{"ALU"}<=$te{"LINE1"} and $te{"LINE1"}>=$te{"SVA"} ){print "$F[0]\t$F[1]\t$F[2]\tLINE1";} else {print "$F[0]\t$F[1]\t$F[2]\tSVA";} }' |perl -F'\t' -alne 'if($F[1]>50) { print "$_";}' >vsoft_pos/${output}_raw_vsoft.bed

#bedtools intersect -a <(cat vsoft_pos/${output}_raw_vsoft.bed|cut -f1-3|perl -F'\t' -alne '$start=$F[1]-50;$end=$F[1]+50;print "$F[0]\t$start\t$end\t$F[2]";') -b $hd_bed -v |perl -F'\t' -alne '$pos=int(($F[1]+$F[2])/2) ;print "$F[0]\t$pos\t$F[3]";' |perl -F'\t' -alne 'if($F[1]>100) {print "$F[0]\t$F[1]\t0\tHG002\t$F[2]";}'  >vsoft_pos/${output}_extract_vsoft.bed

cat candidate_pos/${output}_soft_candidate.txt |
sort -k1,1 -k2,2n candidate_pos/${output}_soft_candidate.txt |cut -f1,2,4|grep -f chr_list.txt |perl -F'\t' -alne '$start=$F[1]-1;$F[2]=~s/_.*//;print "$F[0]\t$start\t$F[1]\t$F[2]";'|bedtools merge -d 10 -c 4,4 -o collapse,count| perl -F'\t' -alne 'if($F[4]>1) {my %te;@tes=split(/,/,$F[3]);foreach $i (@tes) { $te{$i}++;}; if($te{"ALU"}>=$te{"LINE1"} and $te{"ALU"}>=$te{"SVA"} ){print "$F[0]\t$F[1]\t$F[2]\tALU";} elsif($te{"ALU"}<=$te{"LINE1"} and $te{"LINE1"}>=$te{"SVA"} ){print "$F[0]\t$F[1]\t$F[2]\tLINE1";} else {print "$F[0]\t$F[1]\t$F[2]\tSVA";} }' |perl -F'\t' -alne 'if($F[1]>50) { print "$_";}' >vsoft_pos/${output}_raw_vsoft.bed

bedtools intersect -a <(cat vsoft_pos/${output}_raw_vsoft.bed|cut -f1-4) -b <(cat $hd_bed|perl -F'\t' -alne '$start=$F[1]-50;$end=$F[2]+50;if($start<0) {$start=0;};print "$F[0]\t$start\t$end";') -v  |perl -F'\t' -alne '$mean=int(($F[1]+$F[2])/2);if($mean >100) {print "$F[0]\t$mean\t0\tHG002\t$F[3]";}'  >vsoft_pos/${output}_extract_vsoft.bed

#cat vsoft_pos/${output}_raw_vsoft.bed|cut -f1-3| perl -F'\t' -alne 'if($F[1]>100) {print "$F[0]\t$F[1]\t0\tHG002\t$F[2]";}'  >vsoft_pos/${output}_extract_vsoft.bed

#paste <(cat vsoft_pos/${output}_raw_vsoft.bed) <(bedtools getfasta -tab -fi $REF -bed <(cat vsoft_pos/${output}_raw_vsoft.bed|perl -F'\t' -alne '$start=$F[1]-30;$end=$F[1]+30;print "$F[0]\t$start\t$end\t$F[3]";' )) |grep -Ev "T{10,}|A{10,}|C{10,}|G{10,}"|cut -f1-4|perl -F'\t' -alne '$mean=int(($F[1]+$F[2])/2);if($mean >100) {print "$F[0]\t$mean\t0\tHG002\t$F[3]";}'  >vsoft_pos/${output}_extract_vsoft.bed
