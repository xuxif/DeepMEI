file_dis=$1
REF=$2
bam_file=$3

ref_chr=`bash ../DeepMEI_model/reference_check.sh $bam_file $REF |cut -f1`
ref_version=`bash ../DeepMEI_model/reference_check.sh $bam_file $REF|cut -f2`
if [[ $ref_version -eq 38 ]]
then
        ME_bed=../DeepMEI_model/reference/ME_38.bed
fi
if [[ $ref_version -eq 19 ]] && [[ $ref_chr == "nonchr" ]]
then
        ME_bed=../DeepMEI_model/reference/ME_nonchr.bed
fi
if [[ $ref_version -eq 19 ]] && [[ $ref_chr == "chr" ]]
then
        ME_bed=../DeepMEI_model/reference/ME_chr.bed
fi

bedtools intersect -a <(sort -k1,1 -k2,2n $1|grep -f chr_list.txt ) -b $ME_bed  -wa -wb |cut -f4-12,19 |perl -F'\t' -alne 'print "$F[0]\t$F[1]\t$F[2]\t$F[9]\t$F[4]\t$F[5]\t$F[6]\t$F[7]\t$F[8]";' 

