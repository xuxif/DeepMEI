input=$1
start=$2
step=$3
end=$4
re=$5
if [[ $re -eq 1 ]]
then
	bedtools merge -i <( sort -k1,1 -k2,2n $input|perl -npe "s/^/20\t/"|perl -F'\t' -alne 'if($F[8]>$F[0]) { print "$_";}'|cut -f2-|cut -f1,2,5,6,8 |uniq |perl -F'\t' -alne '$start=$F[1]-50;$end=$F[1]+50;print "$F[0]\t$start\t$end\t$F[2]\t$F[3]\t$F[4]";') -d 50  -c 4,5,6 -o collapse,collapse,max |while read chr start end me gts score;do gts=`echo $gts|perl -npe "s/,/\n/g"|sort|uniq -c |perl -npe "s/^ +//;s/ +/\t/g"|sort -k1,1n|tail --lines=1|cut -f2`;me=`echo $me |perl -npe "s/,/\n/g"|sort|uniq -c |perl -npe "s/^ +//;s/ +/\t/g"|sort -k1,1n|tail --lines=1|cut -f2`;echo -e "$chr\t$start\t$end\t$gts\t$score\t$me";done >tmp.bed
fi
seq $start $step $end|while read score;do pp=`cat $input|perl -npe "s/$/\t$score/"|perl -F'\t' -alne 'if($F[4]>=$F[6]) { print "$_";}'|wc -l `;tp=`bedtools intersect -b /DeepMEI/final_vcf/true_giab.bed -a <(cat $input|perl -npe "s/$/\t$score/"|perl -F'\t' -alne 'if($F[4]>=$F[6]) { print "$_";}' )  -wa |wc -l `;precision=`echo "scale=5;$tp/$pp"|bc`;recall=`echo "scale=5;$tp/2000"|bc`;f1=`echo "scale=5;2*$precision*$recall/($precision+$recall)"|bc`;echo -e "$score\t$pp\t$tp\t$precision\t$recall\t$f1";done

