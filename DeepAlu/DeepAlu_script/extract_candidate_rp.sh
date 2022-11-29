input_file=$1
file_name=`echo $input_file|perl -npe 's/.*\///'`
split_len=$2
ME_REF=$3
threads=$4
chr=$5
REF=$6
extractSC=SE-MEI/extractSoftclipped
BWA=bwa
$extractSC -l $split_len <(samtools view -T $REF -h $input_file ${chr} ) > /dev/stdout | perl -F'\t' -alne 'if($num%4 le 1) {print "$_";} $num++;'|perl -npe "s/\@soft/>soft/" > soft_candidate_${file_name}_${chr}.txt

