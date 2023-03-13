score_filter=$3
if [[ ! -n "$score_filter"  ]]
then
	score_filter=27.63
fi
score_pri=20
#1       150940113       0       HG002   ALU     2       2.785756765493639       4.272384121478002       ALU     0.37390524      0.15262619      0.47346854
bedtools merge -i <(cat $1 |sort -k1,1 -k2,2n |perl -npe "s/^/$score_pri\t/"|perl -F'\t' -alne 'if($F[8]>$F[0]) { print "$_";}'|cut -f2-|cut -f1,2,5,6,8 |uniq |perl -F'\t' -alne '$start=$F[1]-50;$end=$F[1]+50;print "$F[0]\t$start\t$end\t$F[2]\t$F[3]\t$F[4]";') -d 50  -c 4,5,6 -o collapse,collapse,max |while read chr start end me gts score;do gts=`echo $gts|perl -npe "s/,/\n/g"|sort|uniq -c |perl -npe "s/^ +//;s/ +/\t/g"|sort -k1,1n|tail --lines=1|cut -f2`;me=`echo $me |perl -npe "s/,/\n/g"|sort|uniq -c |perl -npe "s/^ +//;s/ +/\t/g"|sort -k1,1n|tail --lines=1|cut -f2`;echo -e "$chr\t$start\t$end\t$gts\t$score\t$me\t$score_filter";done |perl -F'\t' -alne 'if($F[4]>=$F[6]) {print "$_";}'|cut -f1-6 >../final_vcf/$2
cat $1 |cut -f1,2,5,6,8 |perl -F'\t' -alne '$start=$F[1]-50;$end=$F[1]+50;print "$F[0]\t$start\t$end\t$F[2]\t$F[3]\t$F[4]";' 
