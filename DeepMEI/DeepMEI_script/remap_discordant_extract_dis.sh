input_bam=$1
split_file=$2
out_dir=$3
samtools view $input_bam `cat $out_dir/$split_file |perl -npe "s/\n/ /"|perl -npe "s/$/\n/" ` |perl -F'\t' -alne 'if(($F[1]%4<2 and $F[1]%16 <8) or abs($F[8])>1000) { print "$F[0]";}' >$out_dir/tmp_${split_file}_dis_list.txt
