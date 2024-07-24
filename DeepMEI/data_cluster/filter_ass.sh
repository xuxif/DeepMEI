bam_file=$1
REF=$2
ran_num=$3
out_dir=$4
parallel=$5

ls |grep "mapRef"|perl -npe "s/_mapRef.*//"|xargs -P $parallel -I{} bash ../../Joint_call/filter_before_cnn_re.sh $ran_num {} $bam_file $REF 2 
cat `ls |grep "mapRef"|perl -npe "s/_mapRef.*//" |perl -npe "s/$/_discordant.txt/"|perl -npe "s/\n/ /"|perl -npe "s/$/\n/"` >tmp_discordant_all.txt

rm tmp_discordant_split_* 2>/dev/null

cat tmp_discordant_all.txt |sort -k4,4 -k5,5n >tmp_discordant_all_sort.txt

split  -n l/20 tmp_discordant_all_sort.txt tmp_discordant_split_

ls tmp_discordant_split_* |xargs -P $parallel -I{} bash ../get_dis_mate_fasta.sh $bam_file $REF {}

cat fa_tmp_discordant_split_* >tmp_clip.fa

cat tmp_clip.fa|perl -npe "s/^(>.*)\n/\1\t/" >tmp_clip_oneline.fa
cat tmp_clip_oneline.fa|grep -v ">soft"|perl ../split_dis_oneline.pl ../split_softclipped_$ran_num/

cd ../split_softclipped_$ran_num/

ls ass_seq_* |xargs -n 1 -P $parallel -I{}  python ../get_dis_mate_sequence_Int.py $bam_file {} $REF >tmp_ass_in.fa 

bash ../repeatmask_para.sh tmp_ass_in.fa $parallel 2>&1 >/dev/null

cat tmp_ass_in.fa.masked |perl -npe "s/\n/\t/"|perl -npe "s/>/\n>/g"|grep -v "^$"|perl -npe "s/\t/ /;s/\t//g"|perl -npe "s/ /\t/"|perl -F'\t' -alne '$count=0;$count++ while $F[1] =~ /N/g;$len=length($F[1]);print "$F[0]\t$count\t$len";'|perl -npe "s/>//;s/[:\-_]/\t/g" |sort -k1,1 -k2,2n|bedtools merge -c 4,5,6 -o distinct,sum,sum |perl -F'\t' -alne '$rate=$F[4]/$F[5]; print "$F[0]\t$F[1]\t$F[2]\t$rate\t$F[5]";'|perl -F'\t' -alne 'if($F[4]>200) { if($F[3]>0.9){print "$_\tdis_only"} elsif($F[3]<0.7){ print "$_\tdis_related";} else {print "$_\tdis_moderate";}} else {print "$_\tdis_short";}' >$out_dir/filter_ass_rate.bed
cat tmp_ass_in.fa.out |perl -npe "s/^ +//;s/ +/\t/g"|cut -f5-7,10,11|perl -npe "s/_left//;s/_right//"|tail --lines=+4|grep -E "Alu|LINE|SVA"|perl -F'\t' -alne 'print "$F[0]_$F[3]\t$F[1]\t$F[2]\t$F[4]";' |sort -k1,1 -k2,2n -k3,3n |bedtools merge -c 4 -o distinct |perl -F'\t' -alne '$len=$F[2]-$F[1];print "$F[0]\t1\t50\t$len\t$F[3]";'|bedtools merge -c 4,5 -o sum,distinct |perl -npe "s/_/\t/"|sort -k1,1 -k5,5nr|perl -F'\t' -alne 'if(not exists $RE{$F[0]}) {print "$F[0]\t$F[1]#$F[5]";$RE{$F[0]}=1;}'|perl -npe "s/-/\t/;s/:/\t/" >$out_dir/filter_ass_me.bed

bedtools intersect -a $out_dir/filter_ass_rate.bed -b $out_dir/filter_ass_me.bed -loj |cut -f1-6,10,11 > $out_dir/filter_ass.bed

#cat <(bedtools intersect -a <(bcftools view -H -i 'ME_int="ME_only"' ~/DeepMEI_output/HG002_bwa_sort.bam/HG002_bwa_sort.bam.vcf |pos2bed) -b <(cat filter_dis.bed|grep "related") -v ) <(bedtools intersect -a <(bcftools view -H -i 'ME_int="ME_related"' ~/DeepMEI_output/HG002_bwa_sort.bam/HG002_bwa_sort.bam.vcf |pos2bed) -b <(cat filter_dis.bed|grep  "dis_only") -wa -u )  >tmp_rs.bed

