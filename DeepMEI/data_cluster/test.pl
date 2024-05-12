$line=`cat regions_2203/split_ac_input_candidate.sam|tail -n 1`;
$line=`grep -A1 "chr13:21242625-21242725" regions_2203/split_ac_input_candidate.sam|tail -n 1`;
@read=`echo "$line"`;
print "$line\t@read\n";
