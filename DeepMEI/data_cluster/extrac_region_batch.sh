bam_file=$1
REF=$2
file=$3
ran_num=$4
step=$5
insize_in=$6
#samtools view --threads 10 $bam_file -T $REF `cat $file|cut -f1 |perl -npe "s/\n/ /"|perl -npe "s/$/\n/" `  | perl read2region.pl $file |perl alu_discord_support_part2_pipe.pl $ran_num $bam_file $REF $step  >${file}_candidate.txt
#samtools view --threads 10 $bam_file -T $REF `cat $file|cut -f1 |perl -npe "s/\n/ /"|perl -npe "s/$/\n/" `  | perl read2region.pl $file >${file}_candidate.sam
record_count=`cat $file|wc -l`
if [[ $record_count -ge 1 ]]
then
	if [[ $step == 1 ]]
	then
		samtools view --threads 10 $bam_file -T $REF `cat $file|cut -f1 |perl -npe "s/\n/ /"|perl -npe "s/$/\n/" `  | perl read2region.pl $file |perl alu_discord_support_part2_pipe.pl $ran_num $bam_file $REF $step $insize_in
	fi

	if [[ $step == 2 ]]
	then
		 cat $file |perl alu_discord_support_part2_pipe.pl $ran_num $bam_file $REF $step $insize_in
	fi
fi
