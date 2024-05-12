num_gt=`cat input_sort.bed|wc -l `
if [[ $num_gt -gt 10000 ]]
then
        split -l 1000  --additional-suffix=_input input_sort.bed  split_
else
        split -n l/100  --additional-suffix=_input input_sort.bed  split_
fi



