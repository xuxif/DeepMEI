
ran_num=13518
base=/home/hjb_xxf/ssd_1/DeepMEI/DeepAlu/DeepAlu_model
cd $base
python model_test_refine.py -i test2.bed -o batch_cdgc/deepalu_test.txt -s ../data_cluster/split_softclipped_sort_$ran_num/ -r reference/hs37d5.fa  -q 0 -t 40.12
cat batch_cdgc/deepalu_test.txt 

#cd /home/hjb_xxf/ssd_1/DeepMEI/DeepAlu/data_cluster/split_softclipped_25501
#record_check=HG002_13:41252374-41252474
#ls |grep "$record_check"|sed -n '/mapClipL.sam$/p'|while read record;do cat <(cat ../head_$ran_num.sam) <(grep -v "^@" $record |perl -F'\t' -alne 'print "$F[5]\t$_";'|perl -npe "s/(\d+)S.*?\t/\1\t/" |perl -F'\t' -alne '$start=$F[4]-$F[0];print "$start\t$_";'|sort -k1,1n |cut -f3- ) > ../split_softclipped_sort_$ran_num/$record;done &
#
#ls |grep "$record_check"|sed -n '/mapClipR.sam$/p'|while read record ;do  cat <(cat ../head_$ran_num.sam) <(grep -v "^@" $record |sort -k4,4n ) > ../split_softclipped_sort_$ran_num/$record;done &
#
#ls |grep "$record_check"|sed -n '/mapRef.sam$/p'|while read record ;do cat <(cat ../head_$ran_num.sam) <(grep -v "^@" $record |sort -k4,4n ) > ../split_softclipped_sort_$ran_num/$record;done &

#wait
