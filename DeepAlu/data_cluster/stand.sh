#PBS -N simulate
#PBS -o log_out_simulate
#PBS -e log_err_simulate
#PBS -q SWH
#PBS -l nodes=1:ppn=1,mem=20GB
#PBS -l walltime=500:00:00
#PBS -S /usr/bin/bash
cd /public/home/swgenetics_1/mobileElement/TypeTE/test_data/simulate
bash order.sh
