bam_file=$2
REF=$3
ran_num=$1
parallel=$4
echo "Processing raw filtering ...."


#bash ../filter_polyN_mut.sh $ran_num $REF 2>/dev/null >filter_polyN.bed &
ls |grep "mapRef"|cut -d'_' -f2|xargs -P $parallel -I{} perl ../check_polyAT.pl {} |uniq  >filter_polyN.bed

ls |grep "_mapRef"|perl -npe "s/_map.*//"|while read file
do
        cat ../split_softclipped_$ran_num/${file}_BPinfo.txt
done |perl ../BP_z_trans.pl | perl -npe "s/\t/:/;s/\t/\-/;" |xargs -P $parallel -I {} perl ../indel_len_filter.pl $bam_file $REF {} >filter_indel.bed &

ls |grep "mapRef"|cut -d'_' -f2|xargs -I{} -P $parallel bash ../indel_size_filter.sh $ran_num $REF $bam_file {} 2>/dev/null | perl -F'\t' -alne 'if(($F[1]==0 and $F[3]==0) or ($F[1]==0 and ($F[3] <10 or $F[3]>40)) or ($F[3]==0 and ($F[1] <10 or $F[1]>40)) or ($F[1]==0 and $F[2]>$F[4]) or ($F[3]==0 and $F[2]<$F[4]) ){ print "$_"}'|perl -npe "s/[:\-]/\t/g"  >filter_largeIndel.bed &

ls |grep "mapRef"|perl -npe "s/_mapRef.*//"|xargs -P $parallel -I{} bash ../../Joint_call/filter_before_cnn_re.sh $ran_num {} $bam_file $REF 1 &

wait
bedtools intersect \
	-a <(bedtools intersect -a filter_largeIndel.bed -b filter_indel.bed -wa) \
	-b filter_polyN.bed -v  >filter_brm.bed
#chr10   103352681       103352781       0       53      0       37
bedtools intersect -a <(ls |grep "mapRef"|perl -npe "s/HG002_//;s/_mapRef.sam//;s/:/\t/;s/-/\t/") -b filter_brm.bed -v|perl -npe "s/\t/:/;s/\t/-/"|cut -f1|xargs -P $parallel -I{} rm HG002_{}_mapRef.sam HG002_{}_mapClipR.sam HG002_{}_mapClipL.sam

#date

cat `ls |grep "mapRef"|perl -npe "s/_mapRef.*//" |perl -npe "s/$/_discordant.txt/"|perl -npe "s/\n/ /"|perl -npe "s/$/\n/"` >tmp_discordant_all.txt

rm tmp_discordant_split_* 2>/dev/null

cat tmp_discordant_all.txt |sort -k4,4 -k5,5n >tmp_discordant_all_sort.txt

split  -n l/20 tmp_discordant_all_sort.txt tmp_discordant_split_

ls tmp_discordant_split_* |xargs -P $parallel -I{} bash ../get_dis_mate_fasta.sh $bam_file $REF {}

#../../DeepMEI_script/SE-MEI/extractSoftclipped -l 15 <(cat <(samtools view -H $bam_file ) <(cat *_mapClip*.sam |grep -v "^@") ) |samtools import  /dev/stdin |samtools fasta |cut -d'|' -f1,4,9 |perl -npe "s/(>.*)\n/\1\t/"|perl -F'\t' -alne '$len=length($F[1]);print "$_\t$len";'|sort -k1,1 -k3,3nr |uniq |perl -F'\t' -alne 'if($last=~/$F[1]/) {} else { print "$F[0]\n$F[1]";} $last=$F[1];'  >tmp_clip.fa
cat HG002*_mapClip*.sam |grep -v "^@" |cut -f1,3,4|sort -k2,2 -k3,3n|uniq >tmp_map_clip.txt
../../DeepMEI_script/SE-MEI/extractSoftclipped -l 15 <(samtools view -h $bam_file  -T $REF `cat tmp_map_clip.txt |perl -F'\t' -alne '$start=$F[2]-50;$end=$F[2]+50;print "$F[1]\t$start\t$end";'|bedtools merge |perl -npe "s/\t/:/;s/\t/-/;s/\n/ /"|perl -npe "s/$/\n/"`|perl ../select_ori_read.pl ) |samtools import  /dev/stdin |samtools fasta 2>/dev/null |cut -d'|' -f1,4,7,9  |perl -F'\|' -alne 'if($_=~/^>/){ if($F[2]=~/S$/){$F[2]="right";}else {$F[2]="left";} $name="$F[0]|$F[1]|$F[3]|$F[2]"; }else {print $seq="$name\t$_\t".length($_);}'|sort -k1,1 -k3,3nr |uniq |perl -F'\t' -alne 'if($last=~/$F[1]/) {} else { print "$F[0]\n$F[1]";} $last=$F[1];' >tmp_clip.fa
#../../DeepMEI_script/SE-MEI/extractSoftclipped -l 15 <(cat <(samtools view -H $bam_file ) <(cat *_mapClip*.sam |grep -v "^@") ) |samtools import  /dev/stdin |samtools fasta |cut -d'|' -f1,4,7,9  |perl -F'\|' -alne 'if($_=~/^>/){ if($F[2]=~/S$/){$F[2]="right";}else {$F[2]="left";} $name="$F[0]|$F[1]|$F[3]|$F[2]"; }else {print $seq="$name\t$_\t".length($_);}'|sort -k1,1 -k3,3nr |uniq |perl -F'\t' -alne 'if($last=~/$F[1]/) {} else { print "$F[0]\n$F[1]";} $last=$F[1];' >tmp_clip.fa

#date
wait
cat fa_tmp_discordant_split_* >>tmp_clip.fa

cat tmp_clip.fa|perl -npe "s/^(>.*)\n/\1\t/" >tmp_clip_oneline.fa

rm tmp_clip_split_* 2>/dev/null
split  -n l/20 tmp_clip_oneline.fa tmp_clip_split_

ls tmp_clip_split_*|while read file;do cat $file |perl -npe "s/\t/\n/" >rm_${file};done

#date
ls rm_tmp_clip_split*|xargs -P $parallel -I{} bash  ../repeatmask.sh {} 2>&1 >/dev/null

cat rm_tmp_clip_split*.out |perl -npe "s/HG002_//"|grep -E "SINE\/Alu|LINE\/L1|Retroposon\/SVA"|perl -npe "s/^ +//;s/ +/\t/g" |cut -f5,10,11|perl -npe "s/_.*?\t/\t/;s/:/\t/;s/\-/\t/"|perl -F'\t' -alne 'if($F[0]=~/soft/) { @region=split(/\|/,$F[0]);$start=$region[2]-1;$end=$region[2]+1;print "$region[1]\t$start\t$end\t$F[1]\t$F[2]";} else {print $_;}' >filter_rm.bed

cat rm_tmp_clip_split*.out |grep "soft" |perl -npe "s/HG002_//"|grep -E "SINE\/Alu|LINE\/L1|Retroposon\/SVA"|perl -npe "s/^ +//;s/ +/\t/g" |cut -f5,10-14|perl -npe "s/\t\(\d+\)//;s/left_\d+/left/;s/right_\d+/right/"|python ../clip_ME_clue.py >clip_ME_clue.txt
bedtools intersect \
	-a <(bedtools intersect \
		-a <(bedtools intersect -a filter_rm.bed -b filter_largeIndel.bed -wa ) \
		-b filter_indel.bed -wa) \
	-b filter_polyN.bed -v  >filter.bed
echo "Raw filtering finished!"
