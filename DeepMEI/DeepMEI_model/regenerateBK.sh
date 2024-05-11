ran_num=$RANDOM
bam_file=$1
ME_REF=~/DeepMEI/DeepMEI_model/reference/ME_add_ALU.fa
REF=~/DeepMEI/DeepMEI_model/reference/Homo_sapiens_assembly38.fasta
output=`echo $bam_file|perl -npe "s/.*\///;s/deepalu_//;s/\.bam//"`
base=~
cd ~/DeepMEI/data_cluster/
cat /NAS/DeepMEI/DeepMEI_output/$output/*.final.vcf.bed |perl -F'\t' -alne '$pos=int(($F[1]+$F[2])/2);print "$F[0]\t$pos\t0\tHG002\tME";' >input_${ran_num}.bed
bash order_RE_sample.sh $bam_file $ran_num input_${ran_num}.bed $ME_REF $REF
cat ${base}/DeepMEI/data_cluster/split_softclipped_${ran_num}/* >/NAS/DeepMEI/DeepMEI_output/$output/${output}_BP.txt
rm -rf split_softclipped_sort_$ran_num split_softclipped_$ran_num input_${ran_num}.bed

