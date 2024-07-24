# 提取所有的 read pairs
REF=$1
input_bam=$2
bam_file=$2
ran_num=$3
input_gt=$4

if [[  -d "regions_dis_$ran_num" ]]
then
        rm -rf regions_dis_$ran_num
fi

mkdir regions_dis_$ran_num
out_dir=regions_dis_$ran_num

sort -k1,1 -k2,2n $input_gt |perl -F'\t' -alne '$start=$F[1]-50;$end=$F[1]+50;print "$F[0]:$start-$end";' >$out_dir/input_sort.bed
split -n l/100  --additional-suffix=_input $out_dir/input_sort.bed  $out_dir/split_

ls $out_dir/|grep "^split"|grep "input$" | xargs -n 1 -I{} -P 20 bash remap_discordant_extract_dis.sh $input_bam {} $out_dir 
ls $out_dir/|grep "^split"|grep "input$" | while read file;do cat $out_dir/tmp_${file}_dis_list.txt;rm $out_dir/tmp_${file}_dis_list.txt $out_dir/$file ;done >$out_dir/read_dis_list.txt


samtools view -T $REF -H $bam_file |grep "SN:" |cut -f2,3|perl -npe "s/SN://;s/LN://"|perl -npe "s/\t/\t1\t/"|grep -f chr_list.txt|bedtools makewindows -b /dev/stdin -w 50000000 |perl -npe "s/\t/:/;s/\t/\-/"|xargs -n 1 -P 20 -I{} samtools view $bam_file -N $out_dir/read_dis_list.txt -o $out_dir/{}.sam {}


 cat <(samtools view -H $bam_file) <(cat $out_dir/*.sam |perl -F'\t' -alne '$chrom=$F[2];$pos=$F[3];if($F[6] ne "=") {$F[2]=$F[6];} $F[3]=$F[7];$F[6]=$chrom;$F[7]=$pos;print join("\t",@F);')|samtools sort -@ 20 --write-index -o $out_dir/dis_read.bam 2>/dev/null

rm $out_dir/*.sam

sort -k1,1 -k2,2n $input_gt |cut -f1,2 |sort -k1,1 -k2n |perl -F'\t' -alne '$start=$F[1]-500;$end=$F[1]+500;print "$F[0]\t$start\t$end";' |bedtools merge |perl -npe "s/\t/:/;s/\t/-/" >$out_dir/dis_read.txt
#record	read_name

cat $out_dir/dis_read.txt | grep -f chr_list.txt| xargs -n 1 -I{} -P 20 bash remap_discordant_bwa.sh $REF $out_dir {}

cat  <(samtools view -H $bam_file) <( ls $out_dir/|grep "_bwa.sam$"|while read file;do cat $out_dir/$file;done) |samtools sort -@ 20 --write-index -o $out_dir/dis_read_remap.bam

