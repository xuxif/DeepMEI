bam_file=$2
REF=$3
ran_num=$1
echo "Processing raw filtering ...."


bash ../filter_polyN_mut.sh $ran_num $REF >filter_polyN.bed &

ls |grep "_mapRef"|perl -npe "s/_map.*//"|while read file
do
        cat ../split_softclipped_$ran_num/${file}_BPinfo.txt
done |perl ../BP_z_trans.pl | perl -npe "s/\t/:/;s/\t/\-/;" |xargs  -P 20 -I {} perl ../indel_len_filter.pl $bam_file $REF {} >filter_indel.bed &

ls |grep "mapRef"|cut -d'_' -f2|xargs  -I{} -P 20 bash ../indel_size_filter.sh $ran_num $REF $bam_file {}  | perl -F'\t' -alne 'if(($F[1]==0 and $F[3]==0) or ($F[1]==0 and ($F[3] <10 or $F[3]>40)) or ($F[3]==0 and ($F[1] <10 or $F[1]>40)) or ($F[1]==0 and $F[2]>$F[4]) or ($F[3]==0 and $F[2]<$F[4]) ){ print "$_"}'|perl -npe "s/[:\-]/\t/g"  >filter_largeIndel.bed &

wait
bedtools intersect \
	-a <(bedtools intersect -a filter_largeIndel.bed -b filter_indel.bed -wa) \
	-b filter_polyN.bed -v  >filter_brm.bed
#chr10   103352681       103352781       0       53      0       37
bedtools intersect -a <(ls |grep "mapRef"|perl -npe "s/HG002_//;s/_mapRef.sam//;s/:/\t/;s/-/\t/") -b filter_brm.bed -v|perl -npe "s/\t/:/;s/\t/-/"|cut -f1|xargs  -P 10 -I{} rm HG002_{}_mapRef.sam HG002_{}_mapClipR.sam HG002_{}_mapClipL.sam

ls |grep "mapRef"|perl -npe "s/_mapRef.*//"|xargs  -P 20 -I{} bash ../../Joint_call/filter_before_cnn_re.sh $ran_num {} $bam_file $REF &

../../DeepMEI_script/SE-MEI/extractSoftclipped -l 15 <(cat <(samtools view -H $bam_file ) <(cat *_mapClip*.sam |grep -v "^@") ) |samtools import  /dev/stdin 2>/dev/null |samtools fasta |cut -d'|' -f1,4,9 |perl -npe "s/(>.*)\n/\1\t/"|perl -F'\t' -alne '$len=length($F[1]);print "$_\t$len";'|sort -k1,1 -k3,3nr |uniq |perl -F'\t' -alne 'if($last=~/$F[1]/) {} else { print "$F[0]\n$F[1]";} $last=$F[1];'  >tmp_clip.fa

wait
cat *_discordant.fa >>tmp_clip.fa
split -n 20 tmp_clip.fa tmp_clip_split
ls tmp_clip_split*|xargs  -P 20 -I{} bash  ../repeatmask.sh {} 

cat tmp_clip_split*.out |perl -npe "s/HG002_//"|grep -E "SINE\/Alu|LINE\/L1|Retroposon\/SVA"|perl -npe "s/^ +//;s/ +/\t/g" |cut -f5,10,11|perl -npe "s/_.*?\t/\t/;s/:/\t/;s/\-/\t/"|perl -F'\t' -alne 'if($F[0]=~/soft/) { @region=split(/\|/,$F[0]);$start=$region[2]-1;$end=$region[2]+1;print "$region[1]\t$start\t$end\t$F[1]\t$F[2]";} else {print $_;}' >filter_rm.bed

bedtools intersect \
	-a <(bedtools intersect \
		-a <(bedtools intersect -a filter_rm.bed -b filter_largeIndel.bed -wa ) \
		-b filter_indel.bed -wa) \
	-b filter_polyN.bed -v  >filter.bed
