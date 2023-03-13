bam_file=$1
REF=$2
file=$3
samtools view --threads 2 $bam_file -T $REF `cat $file |perl -F'\t' -alne '$start=$F[1]-50;$end=$F[1]+50;print "$F[0]:$start-$end";'|perl -npe "s/\n/ /"|perl -npe "s/$/\n/" ` |perl ../sam_split.pl $file
