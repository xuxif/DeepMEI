
bam_input=$1
bedtools intersect -a ../Joint_call/joint_call_hg19.bed -b <( cat $bam_input |perl -F'\t' -alne '$start=$F[1]-10;$end=$F[1]+10;print "$F[0]\t$start\t$end";') -v |perl -F'\t' -alne '$pos=int(($F[1]+$F[2])/2);print "$F[0]\t$pos\t0\tHG002\tME";'
	
  

