grep -v "^@" $1|perl -F'\t' -alne 'if($F[6] eq "=") { print "$_";}' |wc -l|perl -npe "s/\n//"
