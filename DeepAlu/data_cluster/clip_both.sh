clip_b=`cat $1 $2|grep -v "^@"|cut -f6|perl -npe "s/^\d+//"|grep "^S"|grep "S$"|wc -l`
total=`cat $1 $2|grep -v "^@"|wc -l`
echo -e "$clip_b\t$total"|perl -F'\t' -alne 'if($F[1] !=0) { print $F[0]/$F[1];} else { print "0";}' |perl -npe "s/\n$//"
