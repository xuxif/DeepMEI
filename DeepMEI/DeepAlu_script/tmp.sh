output="HG002_vultr"

BWA=bwa
ME_reference=/DeepAlu/DeepAlu_model/reference/ME.fa
threads=40
split_len=15
extractSC=/DeepAlu/DeepAlu_script/SE-MEI/extractSoftclipped
REF=/DeepAlu/DeepAlu_model/reference/Homo_sapiens_assembly38.fasta
$extractSC -l $split_len <(samtools view -h -T $REF /DeepAlu/HG002_bwa_sort.bam |perl -F'\t' -alne 'if($F[5]=~/[S]/ or $F[0]=~/^\@/) {print "$_";}') > /dev/stdout | \
$BWA mem -t $threads -k 19 -r 1.5 -c 100000 -m 50 -T 20 -h 10000 -a -Y -M $ME_reference /dev/stdin  2>/dev/null > /dev/stdout | \
        perl soft_format_polyAT.pl | \
        sort -k1,1 -k2,2n |cut -f1,2,4 |perl -F'\t' -alne '$start=$F[1]-1;if($start>500) { print "$F[0]\t$start\t$F[1]\t$F[2]";}'|bedtools merge -d 10 -c 4,4 -o collapse,count|perl -F'\t' -alne 'if($F[4]>1) {print "$_";}'|cut -f1-4 |while read chr start end mes;do echo -ne "$chr\t$start\t$end\t";echo $mes|perl -npe "s/,/\n/g"|sort |uniq -c |perl -npe "s/^ +//;s/ +/\t/"|sort -k1,1nr|head -n 1 |cut -f2;done |perl -F'\t' -alne 'if($F[0]=~/^(chr)?[XY\d]+$/) {print "$_";}' | \
        perl -F'\t' -alne '$start=$F[1]-30;$end=$F[1]+30;print "$F[0]\t$start\t$end\t$F[3]";'  |while read line;do echo -ne "$line\t";bedtools getfasta -tab -fi $REF -bed <(echo -e "$line");done |grep -Ev "T{10,}|A{10,}|C{10,}|G{10,}"|cut -f1-4|perl -F'\t' -alne '$mean=int(($F[1]+$F[2])/2);print "$F[0]\t$mean\t0\tHG002\t$F[3]";' > candidate_position/${output}_extract_vsoft.bed
samtools view -T $REF /DeepAlu/HG002_bwa_sort.bam -L candidate_position/${output}_extract_vsoft.bed -O CRAM,level=9 -o ${output}.cram -@ 20

