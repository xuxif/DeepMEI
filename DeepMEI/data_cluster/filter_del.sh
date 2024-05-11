#cat ~/DeepMEI_output/HG002_test/deepmei_HG002_test.bed |cut -f1-3|perl -npe "s/\t/:/;s/\t/\-/"|while read record;do ../../DeepMEI_script/SE-MEI/extractSoftclipped -l 15 <(cat HG002_${record}_mapClipL.sam |perl -F'\t' -alne 'if($_=~/^@/) {print $_;next;} $F[1]=2;print join("\t",@F);' ) |samtools import -i /dev/stdin |cut -f10 |grep -v "^$"|while read seq;do grep "$seq" HG002_${record}_mapRef.sam |cut -f10|perl -npe "s/^/$seq\t/"|cut -f1 ;done |sort |uniq |wc -l |perl -npe "s/^/$record\t/";done  >tmp_left.txt  &

#cat ~/DeepMEI_output/HG002_test/deepmei_HG002_test.bed |cut -f1-3|perl -npe "s/\t/:/;s/\t/\-/"|while read record;do ../../DeepMEI_script/SE-MEI/extractSoftclipped -l 15 <(cat HG002_${record}_mapClipR.sam |perl -F'\t' -alne 'if($_=~/^@/) {print $_;next;} $F[1]=2;print join("\t",@F);' ) |samtools import -i /dev/stdin |cut -f10 |grep -v "^$"|while read seq;do grep "$seq" HG002_${record}_mapRef.sam |cut -f10|perl -npe "s/^/$seq\t/"|cut -f1 ;done |sort |uniq |wc -l |perl -npe "s/^/$record\t/";done  >tmp_right.txt  &

#cat tmp_left.txt tmp_right.txt  |perl -npe "s/[:\-]/\t/g"|sort -k1,1 -k2,2n|bedtools merge -c 4 -o collapse |grep -v $'\t'"0," |grep -v ",0"|wc -l 

bedtools intersect -a ~/DeepMEI_output/HG002_test/deepmei_HG002_test_BP.bed -b ~/DeepMEI_output/HG002_test/deepmei_HG002_test.bed -wa |sort |uniq |cut -f 1-5|while read chr start end left right
do
	ex_start=$(($start-150))
	ex_end=$(($end+150))
	perl ../../Joint_call/extrac_outBK_seq.pl $left $right HG002_$chr:${start}-${end}_mapRef.sam >tmp_$chr:${ex_start}-${ex_end}.tsv
	left_num=`grep "left" tmp_$chr:${ex_start}-${ex_end}.tsv |perl ../../Joint_call/filter_del_search.pl  HG002_${chr}:$start-${end}_mapClipL.sam` 
	right_num=`grep "right" tmp_$chr:${ex_start}-${ex_end}.tsv |perl ../../Joint_call/filter_del_search.pl  HG002_${chr}:$start-${end}_mapClipR.sam`
	wait
	if [[ $((${left_num}*${right_num})) -gt 0 ]]
	then
		echo -e "$chr\t$start\t$end\t$left_num\t$right_num"
	fi
	wait
done

