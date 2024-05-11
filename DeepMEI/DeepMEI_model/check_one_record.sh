i=$PWD/$1
o=/home/hjb_xxf/DeepMEI/DeepMEI_model/batch_cdgc/deepmei_test_predict.txt
s=/home/hjb_xxf/DeepMEI/data_cluster/split_softclipped_sort_test
r=/home/hjb_xxf/DeepMEI/DeepMEI_model/reference/hs37d5.fa
q=0
t=40.3
base=/home/hjb_xxf
bam_file=/home/hjb_xxf/ssd_2/HG002_bwa_sort.bam
ran_num=test
r_me=~/DeepMEI/DeepMEI_model/reference/ME.bed
cd /home/hjb_xxf/DeepMEI/data_cluster/
bash  order_RE_sample.sh $bam_file $ran_num $i $r_me $r
cd /home/hjb_xxf/DeepMEI/DeepMEI_model
i=~/DeepMEI/DeepMEI_script/vsoft_pos/test_vsoft.bed
python model_test_refine.py -i $i -o $o -s $s -r $r -q $q -t $t
cat /home/hjb_xxf/DeepMEI/DeepMEI_model/batch_cdgc/deepmei_test_predict.txt

