picard FilterSamReads \
I=input.bam \ 
O=output.bam \ 
READ_LIST_FILE=read_names.txt \ 
FILTER=includeReadList

Filter by interval

java -jar picard.jar FilterSamReads \ 
I=input.bam \ 
O=output.bam \ 
INTERVAL_LIST=regions.interval_list \ 
FILTER=includePairedIntervals

