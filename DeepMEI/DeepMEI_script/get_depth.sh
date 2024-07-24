#samtools depth -@ 50 ~/HG002_bwa_sort.bam -r $1|perl -F'\t' -alne 'if($F[2]>1000){print "$_";} ' |pos2bed 50 |bedtools merge >tmp_$1
samtools depth  $1 -r $2|perl -F'\t' -alne 'if($F[2]>500){print "$_";} ' |perl -F'\t' -alne '$start=$F[1]-50;$end=$F[1]+50;print "$F[0]\t$start\t$end";' |bedtools merge >$3/tmp_depth_$2
