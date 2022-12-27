input_bai="test.sh"
if [[ ! -f "${input_bai}.sh" ]]  && [[ ! -f $input_bai ]]
then
	echo "not existed!"
fi

