bam_file=$1
REF=$2
chr=`samtools view -H $bam_file $REF |grep '^@SQ' |grep -E "SN:chr1"$'\t'"|SN:1"$'\t'|cut -f2|perl -npe "s/SN://"|perl -F'\t' -alne 'if($_=~/chr/) {print "chr";} else {print "nonchr";}'`
version=`samtools view -H $bam_file $REF |grep '^@SQ' |grep -E "SN:chr1"$'\t'"|SN:1"$'\t'|cut -f3|perl -npe "s/LN://"|perl -F'\t' -alne 'if($F[0]==248956422) {print "38";} elsif($F[0]==249250621) {print "19"} else {print "0";}'`
len_bam=`samtools view -H $bam_file $REF |grep '^@SQ' |grep -E "SN:chr1"$'\t'"|SN:1"$'\t'|cut -f3|perl -npe "s/LN://"`
dict_file=`echo $REF|perl -npe "s/\.fa//;s/\.fasta//"`
if [[ ! -f "${dict_file}.dict" ]] && [[ ! -f "${REF}.dict" ]]
	then
		echo "The dictionary file for the reference genome(${REF}.dict or ${dict_file}.dict) was not found."
	exit
fi
len_dict=`cat ${REF}.dict ${dict_file}.dict 2>/dev/null|grep '^@SQ' |grep -E "SN:chr1"$'\t'"|SN:1"$'\t'|cut -f3|perl -npe "s/LN://"`
dict_chr=`cat ${REF}.dict ${dict_file}.dict 2>/dev/null|grep '^@SQ' |grep -E "SN:chr1"$'\t'"|SN:1"$'\t'|cut -f2|perl -npe "s/SN://"|perl -F'\t' -alne 'if($_=~/chr/) {print "chr";} else {print "nonchr";}'`
if [[ $version -eq 0 ]]
then
	echo "Reference sequence length is not match GRCh38, GRCh37 or hg19."
	echo "If you need to use DeepMEI for other species or non-GRCh38, GRCh37 and hg19 reference genomes, please contact us."
	exit
fi
if [[ $len_bam -ne $len_dict ]] 
then
	echo "The provided reference genome is inconsistent with the reference genome used for bam file mapping"
	exit
fi

if [[ "$chr" != "$dict_chr" ]] 
then
	echo "The provided reference genome is inconsistent with the reference genome used for bam file mapping"
	exit
fi
echo -e "$chr\t$version"
