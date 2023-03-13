ran_num=$1
record=$2
bam_file=$3
REF=$4
left=`cat ../split_softclipped_$ran_num/${record}_BPinfo.txt |cut -f2`
right=`cat ../split_softclipped_$ran_num/${record}_BPinfo.txt |cut -f3`



region=`echo $record|perl -npe "s/HG002_//;s/_map.*//"`
cat  ${record}_mapClipL.sam  <(cat ${record}_mapClipR.sam  ${record}_mapRef.sam |grep -v "^@")  |perl -F'\t' -alne 'if($F[8]>1000 or $F[8]< -1000 or $F[8]==0 or $_=~/^@/) {print "$_";} '|bedtools bamtobed |perl -npe "s/$/\t$left\t$right/"|perl -F'\t' -alne 'if(($F[5] eq "+" and $F[2]<$F[$#F]+5) or ($F[5] eq "-" and $F[1]>$F[($#F-1)]-5)) {print $_;}'|cut -f4|perl -npe "s/\/(\d)/\t\1/;s/^/$region\t/" >${record}_discordant.txt

samtools view -h -P $bam_file $region|samtools fasta 2>/dev/null|perl ../filterByNameAddname.pl ${record}_discordant.txt >${record}_discordant.fa


