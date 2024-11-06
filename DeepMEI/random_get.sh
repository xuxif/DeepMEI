ran=$RANDOM
dir=$1
while [[ -f "$dir/random_num/ran_$ran" ]]
do
	sleep 1
	ran=${RANDOM}
done
touch $dir/random_num/ran_$ran
echo $ran
