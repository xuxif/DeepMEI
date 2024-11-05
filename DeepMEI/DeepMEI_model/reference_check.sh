bam_file=$1
REF=$2
chr=`samtools view -H $bam_file $REF |grep '^@SQ' |grep -E "SN:chr1"$'\t'"|SN:1"$'\t'|cut -f2|perl -npe "s/SN://"|perl -F'\t' -alne 'if($_=~/chr/) {print "chr";} else {print "nonchr";}'`
version=`samtools view -H $bam_file $REF |grep '^@SQ' |grep -E "SN:chr1"$'\t'"|SN:1"$'\t'|cut -f3|perl -npe "s/LN://"|perl -F'\t' -alne 'if($F[0]==248956422) {print "38";} elsif($F[0]==249250621) {print "19"} elsif($F[0] ==248387328) {print "13";} else {print "0";}'`
len_bam=`samtools view -H $bam_file $REF |grep '^@SQ' |grep -E "SN:chr1"$'\t'"|SN:1"$'\t'|cut -f3|perl -npe "s/LN://"`
dict_file=`echo $REF|perl -npe "s/\.fa$//;s/\.fasta$//"`
if [[ ! -f "${dict_file}.dict" ]] && [[ ! -f "${REF}.dict" ]]
	then
		echo "The dictionary file for the reference genome(${REF}.dict or ${dict_file}.dict) was not found."
	exit 2
fi
len_dict=`cat ${REF}.dict ${dict_file}.dict 2>/dev/null |grep '^@SQ' |grep -E "SN:chr1"$'\t'"|SN:1"$'\t'|cut -f3|perl -npe "s/LN://"|uniq`
dict_chr=`cat ${REF}.dict ${dict_file}.dict 2>/dev/null |grep '^@SQ' |grep -E "SN:chr1"$'\t'"|SN:1"$'\t'|cut -f2|perl -npe "s/SN://"|uniq|perl -F'\t' -alne 'if($_=~/chr/) {print "chr";} else {print "nonchr";}'`
if [[ "$version" -eq 0 ]]
then
	echo "2"
	exit
elif [[ "$len_bam" != "$len_dict" ]] 
then
	echo "3"
	exit

elif [[ "$chr" != "$dict_chr" ]] 
then
	echo "3"
	exit
else
	echo -e "$chr\t$version"
fi
