bam_file=/DeepMEI/HG002_bwa_sort.bam
#bam_file=/DeepMEI/test.bam
#input_gt=/DeepMEI/DeepMEI_script/HG002_vsoft_discord.bed
input_gt=/DeepMEI/DeepMEI_model/test.bed
cd /DeepMEI/data_cluster/
input_gt_new=/DeepMEI/data_cluster/input_gt.txt
cat  $input_gt|perl -F'\t' -alne 'for($i=-20;$i<20;$i++) {print "$F[0]\t".($i+$F[1])."\t$F[2]\t$F[3]\t$F[4]";}'> $input_gt_new
bash  order_RE_sample.sh $bam_file 'test' $input_gt_new /DeepMEI/DeepMEI_model/reference/ALU.fa
cd -
#cat $input_gt|while read record;
#do
rm -rf not_match_img_one/*
python3.6 model_test_image.py  not_match_img_one $input_gt_new #
#done
