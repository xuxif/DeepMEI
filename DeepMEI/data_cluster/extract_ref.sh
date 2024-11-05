bam_file=/DeepMEI/HG002_bwa_sort.bam
ls |grep "_mapRef"|perl -npe "s/_map.*//"|while read file;do
record=`echo $file|perl -npe "s/HG002_//"`
left=`samtools view -h ${file}_mapClipL.sam -O BAM|bedtools bamtobed|cut -f2|perl -F'\t' -alne '$i=$F[0]+$i;$j=$j+1; END { print int($i/$j) }'`
right=`samtools view -h ${file}_mapClipR.sam -O BAM|bedtools bamtobed|cut -f3|perl -F'\t' -alne '$i=$F[0]+$i;$j=$j+1; END { print int($i/$j) }'`
samtools view ${file}_mapRef.sam|grep -E `samtools view -h ${file}_mapRef.sam|perl -F'\t' -alne 'if(/^@/) {print $_;next;}; if($F[6] eq "=") { $F[6]=$F[2];} if($F[8]==0 or $F[7]< -1000 or $F[8]>1000) { print "$_";}'|samtools view -O BAM |bedtools bamtobed|perl -npe "s/$/\t$left\t$right/"|perl -F'\t' -alne 'if($F[5] eq "+" and $F[2] <= $F[7]) { print $F[3];} ;if($F[5] eq "-" and $F[1] >= $F[6]) { print $F[3];} '|perl -npe "s/\/[12]//"|perl -npe "s/\n/\|/" |perl -npe "s/\|$/\n/"`  |perl -F'\t' -alne 'if($F[6] eq "=") { $F[6]=$F[2];} if($F[8]==0 or $F[8]< -1000 or $F[8]>1000) { print "$F[6]:$F[7]-".($F[7]+1)."\t$F[0]\t$F[1]";}'|xargs  -I{} -P 20 bash ../extract_region_filter.sh $bam_file {} ${record}
done
