REF=$1
ME_REF=$2
region=$3
ran_num=$4
input_file=$5
threads=20
split_len=15
extractSC=SE-MEI/extractSoftclipped
region_s=`echo $region|perl -npe "s/:/_/;s/\-/_/"`
$extractSC -l $split_len <(samtools view -h $input_file $region -@ 10) > /dev/stdout | \
bowtie2 -p 20 -x $ME_REF -  > /dev/stdout  | \
perl soft_format_polyAT.pl > candidate_pos/${output}_${ran_num}_${region_s}_soft_candidate.txt
cat candidate_pos/${output}_${ran_num}_${region_s}_soft_candidate.txt|sort  -k1,1 -k2,2n|grep -v '*'|perl join_record.pl|perl -F'\t' -alne 'if($F[0]=~/^(chr)?[XY\d]+$/) {print "$_";}' >vsoft_pos/${output}_${ran_num}_${region_s}_vsoft.bed
rm candidate_pos/${output}_${ran_num}_${region_s}_soft_candidate.txt
