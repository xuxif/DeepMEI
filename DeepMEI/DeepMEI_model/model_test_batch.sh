
usage() {
    echo "Usage:"
    echo "  bash model_test_batch.sh [-i bamfile] [-r reference] [-b bedfile]  "
    exit -1
}

while getopts i:r:m:b:q:v:w:o:c:j:d:h option
do
   case "${option}"  in  
                i) bam_file=${OPTARG};;
                r) REF=${OPTARG};;
                m) ME_REF=${OPTARG};;
                b) input_gt=${OPTARG};;
                q) quick_model=${OPTARG};;
                o) output_prefix=${OPTARG};;
                v) docker=${OPTARG};;
                w) output_dir=${OPTARG};;
                d) depth=${OPTARG};;
                c) clean=${OPTARG};;
                j) joint=${OPTARG};;
                h) usage;;
                ?) usage;;
   esac
    
done
ran_num=$RANDOM
echo "Random seed is $ran_num"
echo "Joint analysis is $joint"

if [[ ! -n "$quick_model" ]]
then
	quick_model=0
fi
echo "quick_model:$quick_model";
input_bai=`echo $bam_file|perl -npe "s/\.bam/\.bai/;s/\.cram/\.crai/"`
if [[ ! -f "${bam_file}.bai" ]]  && [[ ! -f $input_bai ]] && [[ ! -f "${bam_file}.crai" ]]
then
	echo "Indexed file is missed. Generate index ..."
	samtools index $bam_file -@ 20
fi
if [[ ! -n "$docker" ]]
then
	base=`echo $PWD/$0 |perl -npe "s/\/DeepMEI\/DeepMEI_model\/model_test.*//"`

	if [[ ! -n "$REF" ]]
	then
		echo "-r reference.fa is required"
		exit
	fi
	if [[ ! -f $bam_file ]] && [[ ! -f "$PWD/$bam_file" ]] 
	then
		echo "-r $bam_file is not existed!"
		exit
	fi
	if [[  -f "$PWD/$bam_file" ]] 
	then
		bam_file=$PWD/$bam_file
	fi
	if [[ ! -f $REF ]] && [[ ! -f "$PWD/$REF " ]] 
	then
		echo "-r Reference:$REF is not existed!"
		exit
	fi
	if [[  -f "$PWD/$REF" ]] 
	then
		REF=$PWD/$REF
	fi

	if [ ! -n "$output_prefix" ]; then
		output=`echo $bam_file|perl -npe "s/.*\///"`
	else
		output=$output_prefix
	fi
	
	if [[ ! -n "$output_dir" ]]
	then
		output_dir=${base}
	fi
	ref_chr=`bash $base/DeepMEI/DeepMEI_model/reference_check.sh $bam_file $REF |cut -f1`
	ref_version=`bash $base/DeepMEI/DeepMEI_model/reference_check.sh $bam_file $REF|cut -f2`
	output_dir=`echo $output_dir|perl -F'\t' -alne 'use Cwd(getcwd,cwd);if($_=~/^\//) {print $_;} else { print getcwd()."\/".$_;}'`

	mkdir ${output_dir} 2>/dev/null
	mkdir ${output_dir}/DeepMEI_output 2>/dev/null
	mkdir ${output_dir}/DeepMEI_output/$output 2>/dev/null
	echo "Output folder is :$output_dir/DeepMEI_output/$output"


	if [[ ! -n "$depth" ]]
	then
		echo "Estimated sequencing depth in chr1: 1M-2M"
		reference_type=`samtools view --reference $REF -H $bam_file|grep "SN:chr"|wc -l`
		if [[ $reference_type -gt 0 ]]
		then
			region_depth='chr1:1M-2M'
		else
			region_depth='1:1M-2M'
		fi
		depth=`samtools coverage --reference $REF -r $region_depth $bam_file |cut -f7|tail -n 1|perl -npe "s/\.\d+//"`
		if [[ $depth -eq 0 ]]
		then
			depth=25
			echo "Use defalult sequence depth:25"
		fi
		if [[ $depth -ge 100 ]]
		then
			echo "In chr1:1M-2M the depth is more than 100X, the default setting is 100X"
			depth=100
		fi
	fi
else
	base=/
	ln -s ${base}/root/DeepMEI ${base}/DeepMEI
	export PATH=${base}/DeepMEI/anaconda3/bin:$PATH
	if [[ ! -n "$REF" ]]; then
		REF=${base}/DeepMEI/DeepMEI_model/reference/hs37d5.fa
		echo "Use reference:hs37d5 (hg19 or GRCh37)"
		depth=`samtools coverage -r 1:1M-2M --reference $REF $bam_file |cut -f7|tail -n 1|perl -npe "s/\.\d+//"`
	elif [[ $REF -eq '37' || $REF -eq '19' || $REF -eq 'GRCh37' || $REF -eq 'hg19' || $REF -eq 'HG19' ]]
	then
		REF=${base}/DeepMEI/DeepMEI_model/reference/hs37d5.fa
		echo "Use reference:hs37d5 (hg19 or GRCh37)"
		depth=`samtools coverage --reference $REF -r 1:1M-2M $bam_file |cut -f7|tail -n 1|perl -npe "s/\.\d+//"`
	else
		REF=${base}/DeepMEI/DeepMEI_model/reference/Homo_sapiens_assembly38.fasta
		echo "Use reference:hg38 (GRCh38)"
		depth=`samtools coverage  --reference $REF -r chr1:1M-2M $bam_file |cut -f7|tail -n 1|perl -npe "s/\.\d+//"`
	fi
fi



if [ ! -n "$ME_REF" ]; then
	ME_REF=${base}/DeepMEI/DeepMEI_model/reference/ME_add_ALU.fa
	echo "Use default mobile element reference:$ME_REF"
else
	echo "Load user provided Reference sequence :$ME_REF"
fi
#echo "$bam_file $REF $input_gt";
#exit
input_addr=`echo $bam_file|perl -npe "s/.*\.//"`
#if [[ $input_addr == 'cram' ]]
#then
#	echo "Convert cram file to sam ......"
#	REF=${base}/DeepMEI/DeepMEI_model/reference/Homo_sapiens_assembly38.fasta
#	new_bam=`echo $bam_file|perl -npe "s/cram$/bam/"`
#	/DeepMEI/DeepMEI_model/samtools-1.10/samtools view -T $REF $bam_file  -o $new_bam -@ 20
#	bam_file=$new_bam
#	samtools index $bam_file -@ 20
#fi
if [[ $ref_version -eq 38 ]]
then
        input_db=$base/DeepMEI/Joint_call/joint_call_hg38.bed
        ME_bed=$base/DeepMEI/DeepMEI_model/reference/ME_38.bed
        head_vcf=$base/DeepMEI/DeepMEI_model/head_38.vcf
fi
if [[ $ref_version -eq 19 ]] && [[ $ref_chr == "nonchr" ]]
then
        input_db=$base/DeepMEI/Joint_call/joint_call_hg19_nonchr.bed
        ME_bed=$base/DeepMEI/DeepMEI_model/reference/ME_nonchr.bed
        head_vcf=$base/DeepMEI/DeepMEI_model/head_19_nochr.vcf
fi
if [[ $ref_version -eq 19 ]] && [[ $ref_chr == "chr" ]]
then
        input_db=$base/DeepMEI/Joint_call/joint_call_hg19_chr.bed
        ME_bed=$base/DeepMEI/DeepMEI_model/reference/ME_chr.bed
        head_vcf=$base/DeepMEI/DeepMEI_model/head_19_chr.vcf
fi


if [ ! -n "$input_gt" ]; then
	echo "Start candidate insertion searching....."
	cd ${base}/DeepMEI/DeepMEI_script/
	bash discover_bam_to_candidate.sh $bam_file $ME_REF $output $REF $ME_bed
	input_gt=${base}/DeepMEI/DeepMEI_script/vsoft_pos/${output}_extract_vsoft.bed
else
	no_chr_record=`cat $input_gt|cut -f1|perl -npe "s/^(chr)?\d+//;s/^(chr)?[XY]//"|grep -v "^$"|wc -l`
	if [[ $no_chr_record -gt 0 ]];then
		echo "ERROR:Input record include non chr!"
		exit
	fi
	echo "Get candidate site from user provided:$input_gt!"
	joint=0
fi

input_gt_count=`cat $input_gt|wc -l`
if [[ $input_gt_count -eq 0 ]];
then
	exit 0
fi

if [[ "$joint" == "1" ]]; then
	echo "Candidates from joint calling ...."
	bedtools intersect -a $input_db -b <( cat $input_gt |perl -F'\t' -alne '$start=$F[1]-10;$end=$F[1]+10;print "$F[0]\t$start\t$end";') -v |perl -F'\t' -alne '$pos=int(($F[1]+$F[2])/2);print "$F[0]\t$pos\t0\tHG002\tME";' >${output_dir}/DeepMEI_output/$output/deepmei_${output}_rawinput.bed
	cat $input_gt >>${output_dir}/DeepMEI_output/$output/deepmei_${output}_rawinput.bed
	input_gt=${output_dir}/DeepMEI_output/$output/deepmei_${output}_rawinput.bed
fi

echo "Starting generate region file......"
cd ${base}/DeepMEI/data_cluster/
bash  order_RE_sample.sh $bam_file $ran_num $input_gt $ME_REF $REF 1
cd ${base}/DeepMEI/data_cluster/split_softclipped_sort_$ran_num

bash ../../Joint_call/filter_before_cnn.sh $ran_num $bam_file $REF

bedtools intersect -a <(cat ${base}/DeepMEI/DeepMEI_script/vsoft_pos/${ran_num}_vsoft.bed |perl -F'\t' -alne '$start=$F[1]-50;$end=$F[1]+50;print "$F[0]\t$start\t$end\t$_";') -b filter.bed -wa  -wb |cut -f1-3,12,13|sort -k1,1 -k2,2n|uniq -c |perl -npe "s/ +(.*) (.*)/\2\t\1/"|sort -k1,1 -k2,2n -k6,6nr -k4,4r|perl -F'\t' -alne 'print "$F[3]\t$F[4]\t$F[5]\t$F[0]\t$F[1]\t$F[2]";' |uniq -f3 |perl -F'\t' -alne 'print "$F[3]\t$F[4]\t$F[5]\t0\tHG002\t$F[0]:$F[1]";' >${base}/DeepMEI/DeepMEI_script/vsoft_pos/${ran_num}_vsoft_filter.bed

bedtools intersect \
	-a <(ls ${base}/DeepMEI/data_cluster/split_softclipped_${ran_num}/|grep "_BPinfo.txt"|while read tmp_bp;do cat ${base}/DeepMEI/data_cluster/split_softclipped_${ran_num}/$tmp_bp;done |perl -npe "s/:/\t/;s/\-/\t/" |perl -F'\t' -alne '$dis=abs((($F[1]+$F[2])-($F[3]+$F[4]))/2);print "$dis\t$F[5]\t$F[6]\t$F[7]\t$F[8]\t$F[9]\t$F[10]\t$F[11]\t$F[0]\t$F[1]\t$F[2]\t$F[3]\t$F[4]";'  |sort -k9,9 -k12,12n -k1,1n |uniq -f 11 |cut -f2-|perl -F'\t' -alne 'print "$F[7]\t$F[8]\t$F[9]\t$F[10]\t$F[11]\t$F[0]\t$F[1]\t$F[2]\t$F[3]\t$F[4]\t$F[5]\t$F[6]";') \
	-b filter.bed -wa -u  \
	>${output_dir}/DeepMEI_output/$output/deepmei_${output}_BP_raw.bed


#cat ${base}/DeepMEI/DeepMEI_script/vsoft_pos/${ran_num}_vsoft_filter.bed |while read chr start end info;do num=`grep -m 1 "^$chr"$'\t'"$start"$'\t'"$end"$'\t' ${output_dir}/DeepMEI_output/$output/deepmei_${output}_BP_raw.bed |wc -l`;pos=`grep -m 1 "^$chr"$'\t'"$start"$'\t'"$end"$'\t' ${output_dir}/DeepMEI_output/$output/deepmei_${output}_BP_raw.bed |perl -F'\t' -alne '$pos=int(($F[3]+$F[4])/2);print "$pos";'`;if [[ $num -gt 0 ]];then echo -e "$chr\t$pos\t$info";fi done   >${output_dir}/DeepMEI_output/$output/deepmei_${output}_input.bed
cat ${base}/DeepMEI/DeepMEI_script/vsoft_pos/${ran_num}_vsoft_filter.bed |while read chr start end info;do num=`grep -m 1 "^$chr"$'\t'"$start"$'\t'"$end"$'\t' ${output_dir}/DeepMEI_output/$output/deepmei_${output}_BP_raw.bed |wc -l`;pos=`grep -m 1 "^$chr"$'\t'"$start"$'\t'"$end"$'\t' ${output_dir}/DeepMEI_output/$output/deepmei_${output}_BP_raw.bed |perl -F'\t' -alne '$pos=int(($F[3]+$F[4])/2);print "$pos";'`;if [[ $num -gt 0 ]];then echo -e "$chr\t$pos\t$info";fi done   >${output_dir}/DeepMEI_output/$output/deepmei_${output}_input.bed

cat ${output_dir}/DeepMEI_output/$output/deepmei_${output}_BP_raw.bed|perl -F'\t' -alne '$pos=int(($F[3]+$F[4])/2);$F[1]=$pos-50;$F[2]=$pos+50;print join("\t",@F);' >${output_dir}/DeepMEI_output/$output/deepmei_${output}_BP.bed
rm ${output_dir}/DeepMEI_output/$output/deepmei_${output}_BP_raw.bed

input_gt=${output_dir}/DeepMEI_output/$output/deepmei_${output}_input.bed
cd ..
bash  order_RE_sample.sh $bam_file $ran_num $input_gt $ME_REF $REF 2
input_gt=${base}/DeepMEI/DeepMEI_script/vsoft_pos/${ran_num}_vsoft.bed

echo "Starting DeepMEI prediction ......"
cd ${base}/DeepMEI/DeepMEI_model

threshold_model=`cat depth_threshold.txt|grep "^$depth	"|cut -f2`
#cp  $input_gt input_gt.txt
#python model_test.py       -i $input_gt -o ${base}/DeepMEI/DeepMEI_model/batch_cdgc/deepmei_${output}_predict.txt -s ${base}/DeepMEI/data_cluster/split_softclipped_sort_$ran_num -r $REF  -q $quick_model
multi_process=`ls multi_process/|wc -l`
while [[ $multi_process -ge 5 ]]
do
	multi_process=`ls multi_process/|wc -l`
	sleep 10
done
echo "$output" >multi_process/${output}.txt

if [[ ! -f "weights/val_best_model/variables/variables.data-00000-of-00001" ]]
then
	cat weights/val_best_model/variables/variables.data-00000-of-00001a* >weights/val_best_model/variables/variables.data-00000-of-00001
fi
python model_test_refine.py -i $input_gt -o ${base}/DeepMEI/DeepMEI_model/batch_cdgc/deepmei_${output}_predict.txt -s ${base}/DeepMEI/data_cluster/split_softclipped_sort_$ran_num -r $REF  -q $quick_model -t $threshold_model

rm multi_process/${output}.txt

#chr10:100003075-100003175	100003144	100003131	-13	2	1	0	0	0	0

#bash transPredictTobed.sh batch_cdgc/deepmei_${output}_predict.txt batch_cdgc/deepmei_${output}.bed $threshold_model

cat batch_cdgc/deepmei_${output}_predict.txt |cut -f1,2,5,6,8 |perl -F'\t' -alne '$start=$F[1]-50;$end=$F[1]+50;print "$F[0]\t$start\t$end\t$F[2]\t$F[3]\t$F[4]";' |sort -k1,1 -k2,2n |uniq >../final_vcf/batch_cdgc/deepmei_${output}_cnn.bed

cd ${base}/DeepMEI/data_cluster/split_softclipped_sort_$ran_num
bedtools intersect -a ${output_dir}/DeepMEI_output/$output/deepmei_${output}_BP.bed -b ../../final_vcf/batch_cdgc/deepmei_${output}_cnn.bed -wa |sort |uniq |cut -f 1-5|while read chr start end left right
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
done >filter_del.bed

cat ../../final_vcf/batch_cdgc/deepmei_${output}_cnn.bed |cut -f1-3|perl -npe "s/\t/:/;s/\t/\-/"|while read record;do ../../DeepMEI_script/SE-MEI/extractSoftclipped -l 10 HG002_${record}_mapClipL.sam|samtools import -i /dev/stdin 2> /dev/null |grep -v "^@"|cut -f10|perl -npe "s/^/>$record\n/"  ;done >tmp_left.fa &
cat ../../final_vcf/batch_cdgc/deepmei_${output}_cnn.bed |cut -f1-3|perl -npe "s/\t/:/;s/\t/\-/"|while read record;do ../../DeepMEI_script/SE-MEI/extractSoftclipped -l 10 HG002_${record}_mapClipR.sam|samtools import -i /dev/stdin 2> /dev/null |grep -v "^@"|cut -f10|perl -npe "s/^/>$record\n/"  ;done >tmp_right.fa &
bedtools getfasta -fi $REF -bed  <( cat ${output_dir}/DeepMEI_output/$output/deepmei_${output}_BP.bed |cut -f1,4,5 |perl -F'\t' -alne '$start=$F[1];$end=$F[2];if($start>$end) {$start=$F[2];$end=$F[1];};$start=$start-10;$end=$end+10;print "$F[0]\t$start\t$end";' ) >tmp_bkregion.fa &
cat ../../final_vcf/batch_cdgc/deepmei_${output}_cnn.bed |cut -f1-3|perl -npe "s/\t/:/;s/\t/\-/"|while read record;do poly_C=`samtools view HG002_${record}_mapClipR.sam |cut -f6,10|perl -F'\t' -alne '$F[0]=~s/.*?(\d+)M.*/\1/;$seq=substr($F[1],$F[0]-10,20) ;print "$seq";'|grep "C\{10\}"|wc -l `;echo -e "$record\t$poly_C";done |grep -v "0$" >tmp_polyC_right.txt  &
cat ../../final_vcf/batch_cdgc/deepmei_${output}_cnn.bed |cut -f1-3|perl -npe "s/\t/:/;s/\t/\-/"|while read record;do poly_C=`samtools view HG002_${record}_mapClipR.sam |cut -f6,10|perl -F'\t' -alne '$F[0]=~s/^(\d+)S.*/\1/;$seq=substr($F[1],$F[0]-10,20) ;print "$seq";'|grep "C\{10\}"|wc -l `;echo -e "$record\t$poly_C";done |grep -v "	0$" >tmp_polyC_left.txt  &
wait
bash ../repeatmask.sh tmp_left.fa -noint 2>&1 >/dev/null &
bash ../repeatmask.sh tmp_right.fa  -noint 2>&1 >/dev/null &
bash ../repeatmask.sh tmp_bkregion.fa -noint 2>&1 >/dev/null  &
wait

cat tmp_right.fa.out |perl -npe "s/^ +//;s/ +/\t/g"|cut -f5,10,11|grep -E "Low_complexity|Simple_repeat"|perl -npe "s/_\d+\t/\t/;s/:/\t/;s/\-/\t/" >tmp_right.bed
cat tmp_left.fa.out |perl -npe "s/^ +//;s/ +/\t/g"|cut -f5,10,11|grep -E "Low_complexity|Simple_repeat"|perl -npe "s/_\d+\t/\t/;s/:/\t/;s/\-/\t/" >tmp_left.bed

bedtools intersect -a tmp_left.bed  -b tmp_right.bed -wa -wb 2>/dev/null |cut -f1-4,9|grep -Ev "\(A\)n|\(T\)n"   |perl -npe "s/\(//g;s/\)n//g"|perl -F'\t' -alne 'if($F[3]=~/$F[4]/ or $F[4]=~/$F[3]/) {print "$_";}' |cut -f1-3|sort|uniq  >tmp_soft_clip_low.bed
bedtools intersect  -a  <(cat tmp_bkregion.fa.out |perl -npe "s/ +//;s/ +/\t/g"|cut -f5,10,11|grep -E "Low|Simple"|grep -Ev "\(A\)|\(T\)"|perl -npe "s/[:\-]/\t/g"|cut -f1-4) -b <(bedtools intersect -a tmp_right.bed -b tmp_left.bed -wa -u ) -wa -u  |cut -f1-3 >tmp_bk_region_low.bed

cat tmp_polyC_left.txt tmp_polyC_right.txt |perl -npe "s/[:\-]/\t/g"|sort -k1,1 -k2,2n |bedtools merge -c 4 -o sum |perl -F'\t' -alne 'if($F[3]>2){print $_;}' >filter_polyC.bed

filter_num=`cat filter_del.bed tmp_soft_clip_low.bed tmp_bk_region_low.bed filter_polyC.bed|wc -l`
if [[ $filter_num -gt 0 ]]
then
bedtools intersect -a ../../final_vcf/batch_cdgc/deepmei_${output}_cnn.bed -b <(cat filter_del.bed tmp_soft_clip_low.bed tmp_bk_region_low.bed filter_polyC.bed|cut -f1-3) -v >../../final_vcf/batch_cdgc/deepmei_${output}.bed
else
	cat ../../final_vcf/batch_cdgc/deepmei_${output}_cnn.bed > ../../final_vcf/batch_cdgc/deepmei_${output}.bed
fi
cd ${base}/DeepMEI/DeepMEI_model

#1       FRAM:SINE/Alu   2       65.93868112158563       11031133        11031153        20      14      25      0.0136986301369863      0       9       0

cat <(cat $head_vcf |perl -npe "s/SAMPLE$/${output}/") >${base}/DeepMEI/final_vcf/batch_cdgc/${output}.vcf
cat ../final_vcf/batch_cdgc/deepmei_${output}.bed |while read chr start end info;do echo -ne "$chr\t$info\t"; grep "^$chr"$'\t'"$start"$'\t'"$end"$'\t'  ${output_dir}/DeepMEI_output/$output/deepmei_${output}_BP.bed |cut -f4-  ;done |perl -F'\t' -alne '$F[2]=~s/1/0\/1/;$F[2]=~s/2/1\/1/;$pos=$F[4];$me=$F[1];$me=~s/.*://;if($F[4]>$F[5]){$pos=$F[5];};print "$F[0]\t$pos\t.\tN\t<INS:ME:$me>\t.\tPASS\tSVTYPE=INS;clipPosLeft=$F[4];clipPosRight=$F[5];PS=$F[3];ME=$F[1];TSD_len=$F[6];clipLeftNum=$F[7];clipRightNum=$F[8];clipRate=$F[9];HclipLeftNum=$F[10];HclipRightNum=$F[11];polyN_direction=$F[12]\tGT\t$F[2]";'|sort -k1,1 -k2,2n >>${base}/DeepMEI/final_vcf/batch_cdgc/${output}.vcf


#cat <(cat ../head.vcf|perl -npe "s/SAMPLE$/${output}/") \
#<(bedtools intersect -a ${base}/DeepMEI/final_vcf/batch_cdgc/deepmei_${output}.final.vcf.bed -b <(cat ${base}/DeepMEI/data_cluster/split_softclipped_${ran_num}/*_BPinfo.txt |perl -npe "s/:.*?\t/\t/"|perl -F'\t' -alne 'if($F[1]>$F[2]) {print "$F[0]\t$F[2]\t$F[1]\t$_";} else { print "$F[0]\t$F[1]\t$F[2]\t$_";}') -wa -wb |cut -f4-9,13-|perl -F'\t' -alne '$F[0]=~s/1/0\/1/;$F[0]=~s/2/1\/1/;print "$F[3]\t$F[4]\t.\tN\t<INS:ME:$F[2]>\t.\tPASS\tSVTYPE=INS;clipPosLeft=$F[4];clipPosRight=$F[5];PS=$F[1];ME=$F[2];TSD_len=$F[6];clipLeftNum=$F[7];clipRightNum=$F[8];clipRate=$F[9];HclipLeftNum=$F[10];HclipRightNum=$F[11];polyN_direction=$F[12]\tGT\t$F[0]";'|sort -k1,1 -k2,2n) >${base}/DeepMEI/final_vcf/batch_cdgc/${output}.vcf

cd ${base}/DeepMEI/final_vcf/
num_output=`cat batch_cdgc/deepmei_${output}.bed|wc -l `
if [[ $num_output -gt 10 ]]
then
	samtools view -T $REF $bam_file -O BAM --threads 20  `cat batch_cdgc/deepmei_${output}.bed|cut -f1-3|sort -k1,1 -k2,2n |perl -npe "s/\t/:/;s/\t/\-/;s/\n/ /"|perl -npe "s/$/\n/"`  |samtools sort --threads 20 >batch_cdgc/deepmei_${output}.bam 
	samtools index batch_cdgc/deepmei_${output}.bam
fi

mv ${base}/DeepMEI/final_vcf/batch_cdgc/*${output}*.bed ${output_dir}/DeepMEI_output/$output/
mv ${base}/DeepMEI/final_vcf/batch_cdgc/${output}.vcf ${output_dir}/DeepMEI_output/$output/
mv ${base}/DeepMEI/final_vcf/batch_cdgc/*${output}*.bam ${output_dir}/DeepMEI_output/$output/
mv ${base}/DeepMEI/final_vcf/batch_cdgc/*${output}*.bam.bai ${output_dir}/DeepMEI_output/$output/
if [ ! -n "$clean" ]; then
	cd ${base}/DeepMEI/data_cluster/ 
	echo "#$clean#"
	echo "cleaning tmp file ......"
	rm -rf split_softclipped_sort_$ran_num regions_$ran_num split_softclipped_$ran_num #head_$ran_num.sam
#	if [[ $input_addr == 'cram' ]]
#	then
#		rm $bam_file
#	fi	
else
	echo "Random number is $ran_num Temp file not be remove....."
	echo "${base}/DeepMEI/data_cluster/split_softclipped_sort_$ran_num"
	echo "${base}/DeepMEI/DeepMEI_model/batch_cdgc/deepmei_${output}_predict.txt"
	echo "$input_gt"
fi
