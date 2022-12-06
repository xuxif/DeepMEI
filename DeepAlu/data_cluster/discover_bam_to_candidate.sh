input_file=$1
output=$2
BWA=bwa
ME_reference=/DeepAlu/DeepAlu_model/reference/ALU.fa
threads=20
split_len=15
extractSC=SE-MEI/extractSoftclipped
sample_bam=$3
input_record=`echo $input_file|perl -npe "s/\..*$//";`
cat <(samtools view -H $sample_bam ) <(cat regions/$input_file|grep -v "^@"|grep -n "" |perl -npe "s/^(\d+):(.*?)\t/\2:\1\t/") |samtools view -b >test_num.bam
$extractSC -l $split_len test_num.bam > /dev/stdout | \
$BWA mem -t $threads -k 19 -r 1.5 -c 100000 -m 50 -T 20 -h 10000 -a -Y -M $ME_reference /dev/stdin 2>/dev/null |samtools view -q 10 > /dev/stdout 2>/dev/null | \
perl -F'\t' -alne 'print "$F[0]\t".length($F[9]); '|perl -npe "s/.*://;s/\|/\t/g"|cut -f1,6,12|perl -F'\t' -alne 'if($F[1]=~/^$F[2]S/) { print "$F[0]\tleft";} else { print "$F[0]\t\tright";}'|while read pos direction;do echo $pos >>split_softclipped/${input_record}.rp_${direction};done
