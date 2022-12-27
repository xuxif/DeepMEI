bam_file=$1
strs=$2
id=`echo -e "$strs"|cut -f2`
region=`echo -e "$strs"|cut -f1`
ran=$RANDOM
samtools view $bam_file $region |grep -E "$id" >> ${bam_file}_dis_${ran}.sam

