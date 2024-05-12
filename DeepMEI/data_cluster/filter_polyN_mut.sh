ran_num=$1
REF=$2

ls |grep "mapRef"|cut -d'_' -f2|xargs -n 1 -P 20 -I{} perl ../check_polyAT.pl {}
