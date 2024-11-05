input_file=$1
BWA=bwa
ME_REF=$3
threads=20
split_len=15
extractSC=SE-MEI/extractSoftclipped
input_record=`echo $input_file|perl -npe "s/\..*$//";`
ran_num=${RANDOM}_$1
while [[ -f test_num_${ran_num}.sam ]];
do
        ran_num=${RANDOM}_$1
done

# extract @soft|HISEQ1:20:H9V1RADXX:2:1214:8847:71161|83|1|9996|0|51S97M|*|9996|-30|1|10063
#                                                 .*0  1 2 3    4  5     6  7    8  9  10
cat head_${2}.sam <(cat regions_${2}/$input_file|grep -v "^@"|grep -n "" |perl -npe "s/^(\d+):(.*?)\t/\2:\1\t/")  >test_num_${ran_num}.sam
$extractSC -l $split_len test_num_${ran_num}.sam > /dev/stdout | \
bowtie2 -p 10 -x $ME_REF -  |samtools view -q 10 > /dev/stdout  | \
perl -F'\t' -alne 'print "$F[0]\t".length($F[9]); '|perl -npe "s/.*://;s/\|/\t/g"|cut -f1,6,12|perl -F'\t' -alne 'if($F[1]=~/^$F[2]S/) { print "$F[0]\tleft";} else { print "$F[0]\t\tright";}'|while read pos direction;do echo $pos >>split_softclipped_${2}/${input_record}.rp_${direction};done

rm test_num_${ran_num}.sam
