bam_file=/DeepAlu/HG002_bwa_sort.bam
input_gt=/DeepAlu/final_vcf/all_tool_igv_chekck.bed
output=allTool
cd /DeepAlu/data_cluster/
cp  $input_gt input_gt.txt
bash  order_RE.sh $bam_file
cd -
cp  $input_gt input_gt.txt
python3.6 model_test.py $input_gt batch_cdgc/deepalu_${output}_predict.txt
bash transPredictTobed.sh batch_cdgc/deepalu_${output}_predict.txt batch_cdgc/deepalu_${output}.bed
