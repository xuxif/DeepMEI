file_dis=$1
n_REF=`echo $2|grep "38"|wc -l`
ME_bed=../DeepAlu_model/reference/ME.bed
if [[ $n_REF -eq 1 ]]
then
	ME_bed=../DeepAlu_model/reference/ME_38.bed
fi
bedtools intersect -a <(sort -k1,1 -k2,2n $1|grep -f chr_list.txt ) -b $ME_bed  -wa -wb |cut -f4-12,19 |perl -F'\t' -alne 'print "$F[0]\t$F[1]\t$F[2]\t$F[9]\t$F[4]\t$F[5]\t$F[6]\t$F[7]\t$F[8]";' 

