ran_num=$1
record=$2
bam_file=$3
REF=$4
left=`cat ../split_softclipped_$ran_num/${record}_BPinfo.txt |cut -f2`
right=`cat ../split_softclipped_$ran_num/${record}_BPinfo.txt |cut -f3`



region=`echo $record|perl -npe "s/HG002_//;s/_map.*//"`
#cat <(samtools view -H $bam_file) <(cat ${record}_mapClipL.sam ${record}_mapClipR.sam  ${record}_mapRef.sam |grep -v "^@")  |perl -F'\t' -alne 'if($F[8]>1000 or $F[8]< -1000 or $F[8]==0 or $_=~/^@/) {print "$_";} '|bedtools bamtobed |perl -npe "s/$/\t$left\t$right/"|perl -F'\t' -alne 'if(($F[5] eq "+" and $F[2]<$F[$#F]+5) or ($F[5] eq "-" and $F[1]>$F[($#F-1)]-5)) {print $_;}'|cut -f4|perl -npe "s/\/(\d)/\t\1/;s/^/$region\t/" >${record}_discordant.txt

#cat ${record}_map*.sam |grep -f <(cut -f2 ${record}_discordant.txt) |cut -f3,7,8|perl -F'\t' -alne 'if("$F[1]" eq "="){$F[1]=$F[0];} print "$F[1]\t$F[2]\t$F[2]";'|sort -k1,1 -k2,2n |uniq|perl -npe "s/\t/:/;s/\t/-/" >${record}_discordant_region.txt
#num=`cat ${record}_discordant_region.txt|wc -l`
#if [[ $num -gt 0 ]]
#then
#samtools view -h  -T $REF $bam_file `cat ${record}_discordant_region.txt|perl -npe "s/\n/ /"|perl -npe "s/$/\n/"`|samtools fasta 2>/dev/null|perl ../filterByNameAddname.pl ${record}_discordant.txt >${record}_discordant.fa
#fi


python ../get_dis_mate_sequence.py $bam_file $region $left $right $REF >${record}_discordant.txt
#num=`cat ${record}_discordant_region.txt|wc -l`
#if [[ $num -gt 0 ]]
#then
#samtools view -h  -T $REF $bam_file `cat ${record}_discordant.txt|cut -f4|perl -npe "s/\n/ /"|perl -npe "s/$/\n/"`|samtools fasta 2>/dev/null|perl ../filterByNameAddname.pl ${record}_discordant.txt >${record}_discordant.fa
#fi


