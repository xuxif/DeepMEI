output=$1
ran_num=$2
base=/home/hjb_xxf/ssd_1/DeepMEI/
REF=${base}/DeepMEI/DeepMEI_model/reference/hs37d5.fa
ME_REF=${base}/DeepMEI/DeepMEI_model/reference/ME_add_ALU.fa
bam_file=${base}/HG002_bwa_sort.bam
cd ${base}/DeepMEI/data_cluster
#bash indel_len_filter.sh batch_cdgc/deepalu_${output}.bed $bam_file $REF $output $ran_num $ME_REF >batch_cdgc/deepalu_${output}.vcf.bed
#bedtools getfasta -tab -fi $REF -bed <(cat ${base}/DeepMEI/final_vcf/batch_cdgc/deepalu_${output}.bed |perl -F'\t' -alne 'print "$F[0]:$F[1]\-$F[2]";'|xargs -n 1 -I {} -P 30  perl indel_len_filter.pl {} $bam_file ) |grep -Ev "TTTTTTTTTTTTTTTT|AAAAAAAAAAAAAAAA" |cut -f1|perl -npe "s/[:\-]/\t/g" |while read line;do grep "$line" ${base}/DeepMEI/final_vcf/batch_cdgc/deepalu_${output}.bed ;done |perl -F'\t' -alne '$pos=int(($F[1]+$F[2])/2);print "$F[0]\t$pos\t$F[3]\t$F[4]\t$F[5]";' > ${base}/DeepMEI/final_vcf/batch_cdgc/deepalu_${output}.vcf.bed

cat ${base}/DeepMEI/final_vcf/batch_cdgc/deepalu_${output}.bed  |perl -F'\t' -alne '$pos=int(($F[1]+$F[2])/2);print "$F[0]\t$pos\t0\tHG002\tME";' >input_${ran_num}.bed
bash order_RE_sample.sh $bam_file $ran_num input_${ran_num}.bed $ME_REF $REF
#1       29321412        29321512        1       63.76242115     ALU
#bedtools intersect -a ${base}/DeepMEI/final_vcf/batch_cdgc/deepalu_${output}.bed   -b <(ls split_softclipped_sort_$ran_num/|perl -npe "s/_map.*//"|sort |uniq |while read file;do has_indel=`cat split_softclipped_sort_$ran_num/${file}* |grep -v "^@"|perl indel_len_filter.pl $file`;echo -e "$file\t$has_indel" ;done |grep "0$"|perl -npe "s/.*HG002_//;s/[:\-]/\t/g;s/\.bam//" ) -wa > ${base}/DeepMEI/final_vcf/batch_cdgc/deepalu_${output}.vcf.bed
#cat ${base}/DeepMEI/final_vcf/batch_cdgc/deepalu_${output}.bed  |while read info;do has_indel=`perl indel_len_filter.pl $file`;echo -e "$file\t$has_indel" ;done |grep "0$"|perl -npe "s/.*HG002_//;s/[:\-]/\t/g;s/\.bam//" ) -wa > ${base}/DeepMEI/final_vcf/batch_cdgc/deepalu_${output}.vcf.bed
#bedtools intersect -a ${base}/DeepMEI/final_vcf/batch_cdgc/deepalu_${output}.bed  -b <(cat input_${ran_num}.bed|perl -npe "s/\t/:/g"|while read info;do has_indel=`perl indel_len_filter.pl $info`;echo -e "$file\t$has_indel" ;done |grep "0$"|perl -npe "s/.*HG002_//;s/[:\-]/\t/g;s/\.bam//" ) -wa > ${base}/DeepMEI/final_vcf/batch_cdgc/deepalu_${output}.vcf.bed
bedtools intersect -a ${base}/DeepMEI/final_vcf/batch_cdgc/deepalu_${output}.bed -b <(cat ${base}/DeepMEI/final_vcf/batch_cdgc/deepalu_${output}.bed |cut -f1-3|perl -npe "s/\t/:/;s/\t/\-/;" |xargs -n 1 -P 20 -I {} perl indel_len_filter.pl $bam_file $REF {}) -wa > ${base}/DeepMEI/final_vcf/batch_cdgc/deepalu_${output}.vcf.bed

rm input_${ran_num}.bed

cd split_softclipped_sort_$ran_num


ls |grep "_map"|perl -npe "s/_map.*//"|sort |uniq |while read file;do  cat <(../../DeepMEI_script/SE-MEI/extractSoftclipped -l 15 <(cat <(samtools view -H $bam_file) <(cat ${file}_mapClip*|grep -v "^@")))|perl -npe "s/^(\+.*)\n/\1\t/"|perl -npe "s/^\+\t.*(\n)?//;s/\@soft\|.*/>${file}/"|sort |uniq |perl ../join_clip.pl|grep -v "^$" ;done  >tmp_clip.fa

ls *mapClipR*|while read file;do cat <(samtools view -H $bam_file) <(cat $file|grep -v "^@") |samtools view -f 32  |perl -F'\t' -alne 'if($F[6] eq "=") { $F[6]=$F[2];} if($F[8]==0 or $F[8]< -1000 or $F[8]>1000) { print "$F[6]:$F[7]-".($F[7]+1)."\t$F[0]\t$F[1]";}'|xargs -n 1 -I{} -P 20 bash ../extract_region_filter.sh $bam_file {} $file $REF ;done  >tmp_dis.fa
ls *mapClipL*|while read file;do cat <(samtools view -H $bam_file) <(cat $file|grep -v "^@")|samtools view -f 16 |perl -F'\t' -alne 'if($F[6] eq "=") { $F[6]=$F[2];} if($F[8]==0 or $F[8]< -1000 or $F[8]>1000) { print "$F[6]:$F[7]-".($F[7]+1)."\t$F[0]\t$F[1]";}'|xargs -n 1 -I{} -P 20 bash ../extract_region_filter.sh $bam_file {} $file $REF;done  >>tmp_dis.fa

ls |grep "_mapRef"|perl -npe "s/_map.*//"|while read file;do
record=`echo $file|perl -npe "s/HG002_//"`
left=`samtools view -h ${file}_mapClipL.sam -O BAM|bedtools bamtobed|cut -f2|perl -F'\t' -alne '$i=$F[0]+$i;$j=$j+1; END { print int($i/$j) }'`
right=`samtools view -h ${file}_mapClipR.sam -O BAM|bedtools bamtobed|cut -f3|perl -F'\t' -alne '$i=$F[0]+$i;$j=$j+1; END { print int($i/$j) }'`
samtools view ${file}_mapRef.sam|grep -E `samtools view -h ${file}_mapRef.sam|perl -F'\t' -alne 'if(/^@/) {print $_;next;}; if($F[6] eq "=") { $F[6]=$F[2];} if($F[8]==0 or $F[7]< -1000 or $F[8]>1000) { print "$_";}'|samtools view -O BAM |bedtools bamtobed|perl -npe "s/$/\t$left\t$right/"|perl -F'\t' -alne 'if($F[5] eq "+" and $F[2] <= $F[7]) { print $F[3];} ;if($F[5] eq "-" and $F[1] >= $F[6]) { print $F[3];} '|perl -npe "s/\/[12]//"|perl -npe "s/\n/\|/" |perl -npe "s/\|$/\n/"` 2>/dev/null |perl -F'\t' -alne 'if($F[6] eq "=") { $F[6]=$F[2];} if($F[8]==0 or $F[8]< -1000 or $F[8]>1000) { print "$F[6]:$F[7]-".($F[7]+1)."\t$F[0]\t$F[1]";}'|xargs -n 1 -I{} -P 20 bash ../extract_region_filter.sh $bam_file {} ${record} $REF;done >>tmp_dis.fa

cat tmp_dis.fa >>tmp_clip.fa


echo "repeatmasking"
bash  ../repeatmask.sh tmp_clip.fa 1>&1 >/dev/null

echo "repeatmask finished"

bedtools intersect -a ${base}/DeepMEI/final_vcf/batch_cdgc/deepalu_${output}.vcf.bed -b  <(cat tmp_clip.fa.out |grep -E "SINE\/Alu|LINE\/L1|Retroposon\/SVA"|perl -npe "s/^ +//;s/ +/\t/g"|cut -f5|cut -d'_' -f1-3|perl -npe "s/HG002_//"|cut -d'_' -f1|perl -npe "s/[:\-]/\t/g"|sort |uniq )  -wa |sort |uniq > ${base}/DeepMEI/final_vcf/batch_cdgc/deepalu_${output}.filtered.vcf.bed
#echo "1"
ls HG002_* |perl -npe "s/HG002_//;s/_map.*//"|sort|uniq |while read name;do bedtools getfasta -tab -fi $REF -bed <(samtools view HG002_${name}_mapClipL.sam  -O BAM |bedtools bamtobed|cut -f1,2,3|perl -F'\t' -alne 'print "$_\t$F[1]";'|sort -k1,1 -k2,2n|bedtools merge -c 4 -o median|cut -f1,4 |perl -F'\t' -alne '$F[1]=int($F[1]);$start=$F[1]+30;print "$F[0]\t$F[1]\t$start";' ) ;done |grep -E "A{10}|C{10}|G{10}|T{10}|C{5}[ATCG]C{5}"|perl -npe "s/:/\t/;s/\-/\t/"|sort -k1,1 -k2,2n |bedtools merge -c 4 -o collapse >simple_repeat.bed

#echo "2"
ls HG002_* |perl -npe "s/HG002_//;s/_map.*//"|sort|uniq |while read name;do bedtools getfasta -tab -fi $REF -bed <(samtools view HG002_${name}_mapClipR.sam  -O BAM |bedtools bamtobed|cut -f1,2,3|perl -F'\t' -alne 'print "$_\t$F[2]";'|sort -k1,1 -k2,2n|bedtools merge -c 4 -o median|cut -f1,4 |perl -F'\t' -alne '$F[1]=int($F[1]);$start=$F[1]-30;print "$F[0]\t$start\t$F[1]";' ) ;done |grep -E "A{10}|C{10}|G{10}|T{10}|C{5}[ATCG]C{5}"|perl -npe "s/:/\t/;s/\-/\t/"|sort -k1,1 -k2,2n |bedtools merge -c 4 -o collapse >>simple_repeat.bed
#echo "3"
n_simple_repeat=`cat simple_repeat.bed|wc -l `
if [[ $n_simple_repeat -eq 0 ]]
then
	cp ${base}/DeepMEI/final_vcf/batch_cdgc/deepalu_${output}.filtered.vcf.bed  ${base}/DeepMEI/final_vcf/batch_cdgc/deepalu_${output}.final.vcf.bed
else
	bedtools intersect -a ${base}/DeepMEI/final_vcf/batch_cdgc/deepalu_${output}.filtered.vcf.bed -b  simple_repeat.bed  -v |sort |uniq > ${base}/DeepMEI/final_vcf/batch_cdgc/deepalu_${output}.final.vcf.bed
fi
#echo "4"

