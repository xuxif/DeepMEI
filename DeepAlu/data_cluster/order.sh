#python define-landscape_random-insertions-freq-range.py --chassis chr1.seq --te-seqs ME_ref.fa --N 4 --insert-count 10000 --output output.pgd --min-distance 300 --min-freq 0.1 --max-freq 1 
#python build-population-genome.py --chassis chr1.seq --te-seqs ME_ref.fa --pgd output.pgd --output output.pg
python read_individual_illumina-PE.py --reads 45000000 --pg output.pg --read-length 150 --inner-distance 20 --std-dev 10 --error-rate 0.0001 --fraction-chimera 0.02  --fastq-prefix sample_
