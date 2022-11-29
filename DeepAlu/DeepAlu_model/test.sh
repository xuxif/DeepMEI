multi_process=`ls multi_process/|wc -l`
while [[ $multi_process -ge 5 ]]
do
        multi_process=`ls multi_process/|wc -l`
        sleep 10
done
a=2
if [[ ! -f "test.sha"  ]]
then
	echo "yes"
fi
