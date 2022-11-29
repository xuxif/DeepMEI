ran_num=$1;
REF=$2;
bam_file=$3;
record=$4;
samtools view -T $REF $bam_file $record >tmp_${record}.sam
seq_left=`bedtools getfasta -fi $REF -tab -bed <(cat ../split_softclipped_${ran_num}/HG002_${record}_BPinfo.txt|cut -f1,2|perl -F'\t' -alne '$F[0]=~s/:.*//;$start=$F[1]-50;$end=$F[1]-1;print "$F[0]\t$start\t$end";')|cut -f2`
seq_left_match=`../../DeepAlu_script/SE-MEI/extractSoftclipped -l 10 <(cat <(samtools view -H HG002_${record}_mapClipL.sam ) <( samtools view HG002_${record}_mapClipL.sam |cut -f1,6|perl -npe "s/\d+S/S/"|while read name cigar ;do cat tmp_${record}.sam |grep "$name"|grep "$cigar	";done |perl -F'\t' -alne '$F[1]=0;print join("\t",@F);') )  2>/dev/null|samtools import /dev/stdin|cut -f10 |perl ../seq_size_overlap.pl $seq_left left `

seq_right=`bedtools getfasta -fi $REF -tab -bed <(cat ../split_softclipped_${ran_num}/HG002_${record}_BPinfo.txt|cut -f1,3|perl -F'\t' -alne '$F[0]=~s/:.*//;$start=$F[1]+1;$end=$F[1]+50;print "$F[0]\t$start\t$end";')|cut -f2`
seq_right_match=`../../DeepAlu_script/SE-MEI/extractSoftclipped -l 10 <(cat <(samtools view -H HG002_${record}_mapClipR.sam ) <( samtools view HG002_${record}_mapClipR.sam |cut -f1,6|perl -npe "s/\d+S$//"|while read name cigar ;do cat tmp_${record}.sam |grep "$name"|grep "	$cigar";done |perl -F'\t' -alne '$F[1]=0;print join("\t",@F);') ) 2>/dev/null |samtools import /dev/stdin|cut -f10|perl ../seq_size_overlap.pl $seq_right right`
echo -e "$record\t$seq_left_match\t$seq_right_match";

#bedtools getfasta -fi $REF -tab -bed <(cat ../split_softclipped_${ran_num}/HG002_${record}_BPinfo.txt|cut -f1,2|perl -F'\t' -alne '$F[0]=~s/:.*//;$start=$F[1]-50;$end=$F[1]-1;print "$F[0]\t$start\t$end";')|cut -f2
#../../DeepAlu_script/SE-MEI/extractSoftclipped -l 10 <(cat <(samtools view -H HG002_${record}_mapClipL.sam ) <( samtools view HG002_${record}_mapClipL.sam |cut -f1,6|perl -npe "s/\d+S/S/"|while read name cigar ;do cat tmp_${record}.sam |grep "$name"|grep "$cigar	";done |perl -F'\t' -alne '$F[1]=0;print join("\t",@F);') )  2>/dev/null|samtools import /dev/stdin|cut -f10 |perl ../seq_size_overlap.pl $seq_left left 

#bedtools getfasta -fi $REF -tab -bed <(cat ../split_softclipped_${ran_num}/HG002_${record}_BPinfo.txt|cut -f1,3|perl -F'\t' -alne '$F[0]=~s/:.*//;$start=$F[1]+1;$end=$F[1]+50;print "$F[0]\t$start\t$end";')|cut -f2
#../../DeepAlu_script/SE-MEI/extractSoftclipped -l 10 <(cat <(samtools view -H HG002_${record}_mapClipR.sam ) <( samtools view HG002_${record}_mapClipR.sam |cut -f1,6|perl -npe "s/\d+S$//"|while read name cigar ;do cat tmp_${record}.sam |grep "$name"|grep "	$cigar";done |perl -F'\t' -alne '$F[1]=0;print join("\t",@F);') ) 2>/dev/null |samtools import /dev/stdin|cut -f10|perl ../seq_size_overlap.pl $seq_right right
#echo -e "$record\t$seq_left_match\t$seq_right_match";
