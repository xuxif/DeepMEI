ran_num=$1
REF=$2

ls |grep "mapRef"|cut -d'_' -f2|while read record
do
num=`cat <(
	paste <(../../DeepAlu_script/SE-MEI/extractSoftclipped -l 10 ../split_softclipped_sort_$ran_num/HG002_${record}_mapClipL.sam|samtools import /dev/stdin|cut -f10 |grep -E "A{5}$|T{5}$" |perl -npe "s/^.*(.{5})/\1/"|sort |uniq -c |perl -npe "s/^ +//;s/ +/\t/"|sort -k1,1nr|head -n 1|cut -f2) \
		<( ../../DeepAlu_script/SE-MEI/extractSoftclipped -l 10 ../split_softclipped_sort_$ran_num/HG002_${record}_mapClipR.sam|samtools import /dev/stdin|cut -f10|grep -E "^A{5}$|^T{5}" |perl -npe "s/(.{5}).*$/\1/"|sort |uniq -c |perl -npe "s/^ +//;s/ +/\t/"|sort -k1,1nr|head -n 1|cut -f2) ) \
	<(bedtools getfasta -fi $REF -bed <(cat ../split_softclipped_$ran_num/HG002_${record}_BPinfo.txt|cut -f1-3|perl -npe "s/:.*?\t/\t/") -tab|cut -f2|perl -npe "s/.$//" |grep -E "A{10}|T{10}") |wc -l `
	if [[ $num -eq 2 ]]
	then
		echo $record|perl -npe "s/[:\-]/\t/g"
	fi
done 
