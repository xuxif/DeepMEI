bedtools getfasta -fi $1 -bed <(echo -e "$2\t$3\t$4") -tab |cut -f2
