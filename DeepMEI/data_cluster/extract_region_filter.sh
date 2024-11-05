file=$1
record=$3
REF=$4
region=`echo ${2}|perl -npe "s/ +/\t/g"|cut -f1`
name=`echo ${2}|perl -npe "s/ +/\t/g"|cut -f2`
flag=`echo ${2}|perl -npe "s/ +/\t/g"|cut -f3|perl -npe "s/^/scale=0;/;s/$/%128/"|bc`
#echo ${2}|perl -npe "s/ +/\t/g"|cut -f3|perl -npe "s/^/scale=0;/;s/$/%128/"
if [[ $flag -ge 64 ]]
then
	flag=128
else
	flag=64
fi
#samtools view -f $flag -T $REF $file $region -O SAM |grep $name  |cut -f10|perl -npe "s/^/>$record\n/"|perl -npe "s/_mapClip.*//"
samtools view -h  -f $flag -T $REF $file $region |samtools fasta -n /dev/stdin |grep -m 1 -A1 "$name"  |tail -n 1 |perl -npe "s/^/>$record\n/"
