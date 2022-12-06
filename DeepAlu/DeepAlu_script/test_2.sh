samtools view $1 $2 |perl -F'\t' -alne 'if($F[5]=~/S/){$i++;} if(eof()){print "$i";}'
#samtools view $1 $2 |perl -F'\t' -alne '$i++; if(eof()){print "$i";}'
