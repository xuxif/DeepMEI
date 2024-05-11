
samtools view -H $1 |grep "SN:" |cut -f2,3|perl -npe "s/SN://;s/LN://"|perl -npe "s/\t/\t1\t/"|grep "^1	"|bedtools makewindows -b /dev/stdin -w 2000000 |perl -npe "s/\t/:/;s/\t/\-/"|xargs -n 1 -I{} -P 20 bash test_2.sh $1 {}
