samtools view -T ${5} ${2} ${1} -O SAM |head --lines=200 > ${4}/${3}_${1}.bam 
