cat ~/DeepMEI_output/HG002_test/deepmei_HG002_test.bed |cut -f1-3|perl -npe "s/\t/:/;s/\t/\-/"|while read record;do ../../DeepMEI_script/SE-MEI/extractSoftclipped -l 10 HG002_${record}_mapClipL.sam|samtools import -i /dev/stdin 2> /dev/null |grep -v "^@"|cut -f10|perl -npe "s/^/>$record\n/"  ;done >tmp_left.fa &
cat ~/DeepMEI_output/HG002_test/deepmei_HG002_test.bed |cut -f1-3|perl -npe "s/\t/:/;s/\t/\-/"|while read record;do ../../DeepMEI_script/SE-MEI/extractSoftclipped -l 10 HG002_${record}_mapClipR.sam|samtools import -i /dev/stdin 2> /dev/null |grep -v "^@"|cut -f10|perl -npe "s/^/>$record\n/"  ;done >tmp_right.fa &
bedtools getfasta -fi ~/DeepMEI/DeepMEI_model/reference/hs37d5.fa -bed  <( cat ~/DeepMEI_output/HG002_test/deepmei_HG002_test_BP.bed |cut -f1,4,5 |perl -F'\t' -alne '$start=$F[1];$end=$F[2];if($start>$end) {$start=$F[2];$end=$F[1];};$start=$start-10;$end=$end+10;print "$F[0]\t$start\t$end";' ) >tmp_bkregion.fa &
wait
bash ../repeatmask.sh tmp_left.fa -noint 2>&1 >/dev/null &
bash ../repeatmask.sh tmp_left.fa  -noint 2>&1 >/dev/null &
bash ../repeatmask.sh tmp_bkregion.fa -noint 2>&1 >/dev/null  &
wait

cat tmp_right.fa.out |perl -npe "s/^ +//;s/ +/\t/g"|cut -f5,10,11|grep -E "Low_complexity|Simple_repeat"|perl -npe "s/_\d+\t/\t/;s/:/\t/;s/\-/\t/" >tmp_right.bed
cat tmp_left.fa.out |perl -npe "s/^ +//;s/ +/\t/g"|cut -f5,10,11|grep -E "Low_complexity|Simple_repeat"|perl -npe "s/_\d+\t/\t/;s/:/\t/;s/\-/\t/" >tmp_left.bed

bedtools intersect -a tmp_left.bed  -b tmp_right.bed -wa -wb |cut -f1-4,9|grep -Ev "\(A\)n|\(T\)n"   |perl -npe "s/\(//g;s/\)n//g"|perl -F'\t' -alne 'if($F[3]=~/$F[4]/ or $F[4]=~/$F[3]/) {print "$_";}' |cut -f1-3|sort|uniq  >tmp_soft_clip_low.bed 
bedtools intersect  -a  <(cat tmp_bkregion.fa.out |perl -npe "s/ +//;s/ +/\t/g"|cut -f5,10,11|grep -E "Low|Simple"|grep -Ev "\(A\)|\(T\)"|perl -npe "s/[:\-]/\t/g"|cut -f1-4) -b <(bedtools intersect -a tmp_right.bed -b tmp_left.bed -wa -u ) -wa -u  >tmp_bk_region_low.bed
