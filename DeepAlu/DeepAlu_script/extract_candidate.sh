input_file=$1
file_name=`echo $input_file|perl -npe 's/.*\///'`
split_len=$2
ME_REF=$3
threads=$4
chr=$5
REF=$6
output=$7
extractSC=SE-MEI/extractSoftclipped
BWA=bwa
$extractSC -l $split_len <(samtools view -T $REF -h $input_file ${chr} ) > /dev/stdout | \
$BWA mem -t $threads -k 19 -r 1.5 -c 100000 -m 50 -T 20 -h 10000 -a -Y -M $ME_REF /dev/stdin  > /dev/stdout 2>/dev/null | \
perl soft_format_polyAT.pl $REF $output > tmp_$output/soft_candidate_${output}_${chr}.txt

