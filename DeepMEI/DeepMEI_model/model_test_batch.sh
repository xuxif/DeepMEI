
usage() {
    echo "Usage:"
    echo "  bash model_test_batch.sh [-i bamfile] [-r reference] [-b bedfile]  "
    exit -1
}

while getopts i:r:m:b:q:v:w:o:c:d:h option
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
                h) usage;;
                ?) usage;;
   esac
    
done
ran_num=$RANDOM
echo "Random seed is $ran_num"

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
	else
		REF=$REF
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
output=`echo $bam_file|perl -npe "s/.*\///"`
if [ ! -n "$output_prefix" ]; then
	output=`echo $bam_file|perl -npe "s/.*\///"`
else
	output=$output_prefix
fi

if [ ! -n "$output_dir" ]; then
	output_dir=${base}
fi
echo "Output folder is :$output_dir/DeepMEI_output/$output"

if [ ! -n "$input_gt" ]; then
	echo "Start candidate insertion searching....."
	cd ${base}/DeepMEI/DeepMEI_script/
	bash discover_bam_to_candidate.sh $bam_file $ME_REF $output $REF
	input_gt=${base}/DeepMEI/DeepMEI_script/vsoft_pos/${output}_extract_vsoft.bed
else
	no_chr_record=`cat $input_gt|cut -f1|perl -npe "s/^(chr)?\d+//;s/^(chr)?[XY]//"|grep -v "^$"|wc -l`
	if [[ $no_chr_record -gt 0 ]];then
		echo "ERROR:Input record include non chr!"
		exit
	fi
	echo "Get candidate site from user provided:$input_gt!"
fi

input_gt_count=`cat $input_gt|wc -l`
if [[ $input_gt_count -eq 0 ]];
then
	exit 0
fi
echo "Starting generate region file......"
cd ${base}/DeepMEI/data_cluster/
bash  order_RE_sample.sh $bam_file $ran_num $input_gt $ME_REF $REF
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

bash transPredictTobed.sh batch_cdgc/deepmei_${output}_predict.txt batch_cdgc/deepmei_${output}.bed $threshold_model


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
cd ${base}/DeepMEI/data_cluster
#bash indel_len_filter.sh batch_cdgc/deepmei_${output}.bed $bam_file $REF $output $ran_num $ME_REF >batch_cdgc/deepmei_${output}.vcf.bed
#bedtools getfasta -tab -fi $REF -bed <(cat ${base}/DeepMEI/final_vcf/batch_cdgc/deepmei_${output}.bed |perl -F'\t' -alne 'print "$F[0]:$F[1]\-$F[2]";'|xargs -n 1 -I {} -P 30  perl indel_len_filter.pl {} $bam_file ) |grep -Ev "TTTTTTTTTTTTTTTT|AAAAAAAAAAAAAAAA" |cut -f1|perl -npe "s/[:\-]/\t/g" |while read line;do grep "$line" ${base}/DeepMEI/final_vcf/batch_cdgc/deepmei_${output}.bed ;done |perl -F'\t' -alne '$pos=int(($F[1]+$F[2])/2);print "$F[0]\t$pos\t$F[3]\t$F[4]\t$F[5]";' > ${base}/DeepMEI/final_vcf/batch_cdgc/deepmei_${output}.vcf.bed

cat ${base}/DeepMEI/final_vcf/batch_cdgc/deepmei_${output}.bed  |perl -F'\t' -alne '$pos=int(($F[1]+$F[2])/2);print "$F[0]\t$pos\t0\tHG002\tME";' >input_${ran_num}.bed
bash order_RE_sample.sh $bam_file $ran_num input_${ran_num}.bed $ME_REF $REF
#1       29321412        29321512        1       63.76242115     ALU
#bedtools intersect -a ${base}/DeepMEI/final_vcf/batch_cdgc/deepmei_${output}.bed   -b <(ls split_softclipped_sort_$ran_num/|perl -npe "s/_map.*//"|sort |uniq |while read file;do has_indel=`cat split_softclipped_sort_$ran_num/${file}* |grep -v "^@"|perl indel_len_filter.pl $file`;echo -e "$file\t$has_indel" ;done |grep "0$"|perl -npe "s/.*HG002_//;s/[:\-]/\t/g;s/\.bam//" ) -wa > ${base}/DeepMEI/final_vcf/batch_cdgc/deepmei_${output}.vcf.bed
#cat ${base}/DeepMEI/final_vcf/batch_cdgc/deepmei_${output}.bed  |while read info;do has_indel=`perl indel_len_filter.pl $file`;echo -e "$file\t$has_indel" ;done |grep "0$"|perl -npe "s/.*HG002_//;s/[:\-]/\t/g;s/\.bam//" ) -wa > ${base}/DeepMEI/final_vcf/batch_cdgc/deepmei_${output}.vcf.bed
#bedtools intersect -a ${base}/DeepMEI/final_vcf/batch_cdgc/deepmei_${output}.bed  -b <(cat input_${ran_num}.bed|perl -npe "s/\t/:/g"|while read info;do has_indel=`perl indel_len_filter.pl $info`;echo -e "$file\t$has_indel" ;done |grep "0$"|perl -npe "s/.*HG002_//;s/[:\-]/\t/g;s/\.bam//" ) -wa > ${base}/DeepMEI/final_vcf/batch_cdgc/deepmei_${output}.vcf.bed
bedtools intersect -a ${base}/DeepMEI/final_vcf/batch_cdgc/deepmei_${output}.bed -b <(cat ${base}/DeepMEI/final_vcf/batch_cdgc/deepmei_${output}.bed |cut -f1-3|perl -npe "s/\t/:/;s/\t/\-/;" |xargs -n 1 -P 20 -I {} perl indel_len_filter.pl $bam_file $REF {}) -wa > ${base}/DeepMEI/final_vcf/batch_cdgc/deepmei_${output}.vcf.bed

rm input_${ran_num}.bed

cd split_softclipped_sort_$ran_num

ls |grep "_map"|perl -npe "s/_map.*//"|sort |uniq |while read file;do  cat <(../../DeepMEI_script/SE-MEI/extractSoftclipped -l 15 <(cat <(samtools view -H -T $REF $bam_file) <(cat ${file}_mapClip*|grep -v "^@")))|perl -npe "s/^(\+.*)\n/\1\t/"|perl -npe "s/^\+\t.*(\n)?//;s/\@soft\|.*/>${file}/"|sort |uniq |perl ../join_clip.pl|grep -v "^$" ;done  >tmp_clip.fa &

ls *mapClipR*|while read file;do cat <(samtools view -H -T $REF $bam_file) <(cat $file|grep -v "^@") |samtools view -f 32  |perl -F'\t' -alne 'if($F[6] eq "=") { $F[6]=$F[2];} if($F[8]==0 or $F[8]< -1000 or $F[8]>1000) { print "$F[6]:$F[7]-".($F[7]+1)."\t$F[0]\t$F[1]";}'|xargs -n 1 -I{} -P 20 bash ../extract_region_filter.sh $bam_file {} $file $REF ;done  >tmp_dis_1.fa &
ls *mapClipL*|while read file;do cat <(samtools view -H -T $REF $bam_file) <(cat $file|grep -v "^@")|samtools view -f 16 |perl -F'\t' -alne 'if($F[6] eq "=") { $F[6]=$F[2];} if($F[8]==0 or $F[8]< -1000 or $F[8]>1000) { print "$F[6]:$F[7]-".($F[7]+1)."\t$F[0]\t$F[1]";}'|xargs -n 1 -I{} -P 20 bash ../extract_region_filter.sh $bam_file {} $file $REF;done  >>tmp_dis_2.fa &

ls |grep "_mapRef"|perl -npe "s/_map.*//"|while read file;do
	record=`echo $file|perl -npe "s/HG002_//"`
	left=`samtools view -h ${file}_mapClipL.sam -O BAM|bedtools bamtobed|cut -f2|perl -F'\t' -alne '$i=$F[0]+$i;$j=$j+1; END { print int($i/$j) }'`
	right=`samtools view -h ${file}_mapClipR.sam -O BAM|bedtools bamtobed|cut -f3|perl -F'\t' -alne '$i=$F[0]+$i;$j=$j+1; END { print int($i/$j) }'`
	samtools view ${file}_mapRef.sam|grep -E `samtools view -h ${file}_mapRef.sam|perl -F'\t' -alne 'if(/^@/) {print $_;next;}; if($F[6] eq "=") { $F[6]=$F[2];} if($F[8]==0 or $F[7]< -1000 or $F[8]>1000) { print "$_";}'|samtools view -O BAM |bedtools bamtobed|perl -npe "s/$/\t$left\t$right/"|perl -F'\t' -alne 'if($F[5] eq "+" and $F[2] <= $F[7]) { print $F[3];} ;if($F[5] eq "-" and $F[1] >= $F[6]) { print $F[3];} '|perl -npe "s/\/[12]//"|perl -npe "s/\n/\|/" |perl -npe "s/\|$/\n/"` 2>/dev/null |perl -F'\t' -alne 'if($F[6] eq "=") { $F[6]=$F[2];} if($F[8]==0 or $F[8]< -1000 or $F[8]>1000) { print "$F[6]:$F[7]-".($F[7]+1)."\t$F[0]\t$F[1]";}'|xargs -n 1 -I{} -P 20 bash ../extract_region_filter.sh $bam_file {} ${record} $REF;done >tmp_dis_3.fa &

wait
cat tmp_dis_*.fa >>tmp_clip.fa


echo "repeatmasking"
bash  ../repeatmask.sh tmp_clip.fa 2>&1 >/dev/null

echo "repeatmask finished"

bedtools intersect -a ${base}/DeepMEI/final_vcf/batch_cdgc/deepmei_${output}.vcf.bed -b  <(cat tmp_clip.fa.out |grep -E "SINE\/Alu|LINE\/L1|Retroposon\/SVA"|perl -npe "s/^ +//;s/ +/\t/g"|cut -f5|cut -d'_' -f1-3|perl -npe "s/HG002_//"|cut -d'_' -f1|perl -npe "s/[:\-]/\t/g"|sort |uniq )  -wa |sort |uniq > ${base}/DeepMEI/final_vcf/batch_cdgc/deepmei_${output}.filtered.vcf.bed
#echo "1"
#ls HG002_* |perl -npe "s/HG002_//;s/_map.*//"|sort|uniq |while read name;do bedtools getfasta -tab -fi $REF -bed <(samtools view HG002_${name}_mapClipL.sam  -O BAM |bedtools bamtobed|cut -f1,2,3|perl -F'\t' -alne 'print "$_\t$F[1]";'|sort -k1,1 -k2,2n|bedtools merge -c 4 -o median|cut -f1,4 |perl -F'\t' -alne '$F[1]=int($F[1]);$start=$F[1]+30;print "$F[0]\t$F[1]\t$start";' ) ;done |grep -E "A{10}|C{10}|G{10}|T{10}|C{5}[ATCG]C{5}"|perl -npe "s/:/\t/;s/\-/\t/"|sort -k1,1 -k2,2n |bedtools merge -c 4 -o collapse >simple_repeat.bed

#echo "2"
#ls HG002_* |perl -npe "s/HG002_//;s/_map.*//"|sort|uniq |while read name;do bedtools getfasta -tab -fi $REF -bed <(samtools view HG002_${name}_mapClipR.sam  -O BAM |bedtools bamtobed|cut -f1,2,3|perl -F'\t' -alne 'print "$_\t$F[2]";'|sort -k1,1 -k2,2n|bedtools merge -c 4 -o median|cut -f1,4 |perl -F'\t' -alne '$F[1]=int($F[1]);$start=$F[1]-30;print "$F[0]\t$start\t$F[1]";' ) ;done |grep -E "A{10}|C{10}|G{10}|T{10}|C{5}[ATCG]C{5}"|perl -npe "s/:/\t/;s/\-/\t/"|sort -k1,1 -k2,2n |bedtools merge -c 4 -o collapse >>simple_repeat.bed
#echo "3"
bash ../filter_polyN_mut.sh $ran_num $REF 2>/dev/null >tmp_polyN.bed
bedtools intersect -a ${base}/DeepMEI/final_vcf/batch_cdgc/deepmei_${output}.filtered.vcf.bed -b  tmp_polyN.bed  -v |sort |uniq > ${base}/DeepMEI/final_vcf/batch_cdgc/deepmei_${output}.polyN.vcf.bed

ls |grep "mapRef"|cut -d'_' -f2|xargs -n 1 -I{} -P 20 bash ../indel_size_filter.sh $ran_num $REF $bam_file {} 2>/dev/null >tmp_size.txt
cat tmp_size.txt |perl -F'\t' -alne 'if(($F[1]==0 and $F[3]==0) or ($F[1]==0 and ($F[3] <10 or $F[3]>40)) or ($F[3]==0 and ($F[1] <10 or $F[1]>40)) or ($F[1]==0 and $F[2]>$F[4]) or ($F[3]==0 and $F[2]<$F[4]) ){ print "$_"}'|perl -npe "s/[:\-]/\t/g"  >tmp_largeIndel.bed

bedtools intersect -a ${base}/DeepMEI/final_vcf/batch_cdgc/deepmei_${output}.polyN.vcf.bed -b  tmp_largeIndel.bed  -wa |sort |uniq > ${base}/DeepMEI/final_vcf/batch_cdgc/deepmei_${output}.final.vcf.bed
#echo "4"
cat <(cat ../head.vcf|perl -npe "s/SAMPLE$/${output}/") \
<(bedtools intersect -a ${base}/DeepMEI/final_vcf/batch_cdgc/deepmei_${output}.final.vcf.bed -b <(cat ${base}/DeepMEI/data_cluster/split_softclipped_${ran_num}/*_BPinfo.txt |perl -npe "s/:.*?\t/\t/"|perl -F'\t' -alne 'if($F[1]>$F[2]) {print "$F[0]\t$F[2]\t$F[1]\t$_";} else { print "$F[0]\t$F[1]\t$F[2]\t$_";}') -wa -wb |cut -f4-9,13-|perl -F'\t' -alne '$F[0]=~s/1/0\/1/;$F[0]=~s/2/1\/1/;print "$F[3]\t$F[4]\t.\tN\t<INS:ME:$F[2]>\t.\tPASS\tSVTYPE=INS;clipPosLeft=$F[4];clipPosRight=$F[5];PS=$F[1];ME=$F[2];TSD_len=$F[6];clipLeftNum=$F[7];clipRightNum=$F[8];clipRate=$F[9];HclipLeftNum=$F[10];HclipRightNum=$F[11];polyN_direction=$F[12]\tGT\t$F[0]";'|sort -k1,1 -k2,2n) >${base}/DeepMEI/final_vcf/batch_cdgc/${output}.vcf

cd ${base}/DeepMEI/final_vcf/
samtools view -T $REF $bam_file -O BAM --threads 20  `cat batch_cdgc/deepmei_${output}.bed|cut -f1-3|sort -k1,1 -k2,2n |perl -npe "s/\t/:/;s/\t/\-/;s/\n/ /"|perl -npe "s/$/\n/"`  |samtools sort --threads 20 >batch_cdgc/deepmei_${output}.bam 
samtools index batch_cdgc/deepmei_${output}.bam

mkdir ${base}/DeepMEI_output 2>/dev/null
mkdir ${base}/DeepMEI_output/$output 2>/dev/null
mv ${base}/DeepMEI/final_vcf/batch_cdgc/*${output}*.bed ${output_dir}/DeepMEI_output/$output/
mv ${base}/DeepMEI/final_vcf/batch_cdgc/${output}.vcf ${output_dir}/DeepMEI_output/$output/
mv ${base}/DeepMEI/final_vcf/batch_cdgc/*${output}*.bam ${output_dir}/DeepMEI_output/$output/
mv ${base}/DeepMEI/final_vcf/batch_cdgc/*${output}*.bam.bai ${output_dir}/DeepMEI_output/$output/
if [ ! -n "$clean" ]; then
	cd ${base}/DeepMEI/data_cluster/ 
	echo "#$clean#"
	echo "cleaning tmp file ......"
	rm -rf split_softclipped_sort_$ran_num regions_$ran_num split_softclipped_$ran_num head_$ran_num.sam
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
