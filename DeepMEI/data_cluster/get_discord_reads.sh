apt-get install samtools -y

reference=reference/Homo_sapiens_assembly38.fasta
cp /content/drive/Shareddrives/ibicqupt/deepalu_data/reference/Homo_sapiens_assembly38* .
cp drive/Shareddrives/ibicqupt/deepalu_data/software/head.sam .
cp drive/Shareddrives/ibicqupt/deepalu_data/software/sample_list.txt .
cp drive/Shareddrives/ibicqupt/deepalu_data/input_data/regions.tar.gz  .
tar -xzf regions.tar.gz
cat sample_list.txt|while read sample
do
	cp drive/Shareddrives/ibicqupt/deepalu_data/download_bam/${sample}.final.cram* .
	ls regions |grep "$sample" |while read file
	do
		samtools view -F 3854  regions/$file |cut -f1,3,7,8|perl -npe "s/(.*?)\t(.*?)\t(=)\t(.*)/\1\t\2\t\2\t\4/"|cut -f1,3,4|perl -F'\t' -alne '$start=$F[2]-1;$end=$F[2]+1;print "$F[0]\t$F[1]:$start\-$end";' |while read read_name region
		do
			cat <(cat head.sam) <(samtools view -h ${sample}.final.cram  -T $reference $region |grep "$read_name" )|samtools view -O BAM -o discord_read/$file 
		done
	done
	rm ${sample}.final.cram*
done


sample='HG002'
	cp drive/Shareddrives/ibicqupt/deepalu_data/download_bam/HG002.GRCh38.60x.1.ba* .
	ls regions |grep "$sample" |while read file
	do
		samtools view -F 3854  regions/$file |cut -f1,3,7,8|perl -npe "s/(.*?)\t(.*?)\t(=)\t(.*)/\1\t\2\t\2\t\4/"|cut -f1,3,4|perl -F'\t' -alne '$start=$F[2]-1;$end=$F[2]+1;print "$F[0]\t$F[1]:$start\-$end";' |while read read_name region
		do
			cat <(cat head.sam) <(samtools view -h HG002.GRCh38.60x.1.bam  $region |grep "$read_name" )|samtools view -O BAM -o discord_read/$file 
		done
	done
	rm HG002.GRCh38.60x.1.ba*

tar -czf discord_read.tar.gz discord_read/
