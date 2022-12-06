record=$1
((record=$record+1))
chr=`sed -n "${record}p" input_gt_HG002.txt |cut -f1`
sample=`sed -n "${record}p" input_gt_HG002.txt |cut -f4`
suffix=`sed -n "${record}p" input_gt_HG002.txt |cut -f2|perl -F'\t' -alne '$start=$F[0]-50;$end=$F[0]+50;print "$start-$end"'`
ls -lh split_softclipped/${sample}_${chr}:${suffix}*
samtools view split_softclipped/${sample}_${chr}:${suffix}_mapClipL.sam|wc -l 
samtools view split_softclipped/${sample}_${chr}:${suffix}_mapClipR.sam|wc -l 
samtools view split_softclipped/${sample}_${chr}:${suffix}_mapRef.sam|wc -l 
perl alu_discord_support.pl ${sample}_${chr}:${suffix}
samtools view split_softclipped/${sample}_${chr}:${suffix}_mapClipL.sam|wc -l 
samtools view split_softclipped/${sample}_${chr}:${suffix}_mapClipR.sam|wc -l 
samtools view split_softclipped/${sample}_${chr}:${suffix}_mapRef.sam|wc -l 


