import pysam
import subprocess
import os
import shutil


def create_directory(directory_path):
    if os.path.exists(directory_path):
        shutil.rmtree(directory_path)  # 删除已存在的目录及其内容
    os.makedirs(directory_path)  # 创建新的目录

def extract_reads_with_mates(bam_file, region, max_insert_size):
    soft_clipped_reads = []
    discordant_reads = []
    mate_reads = []

    # 打开 BAM 文件并遍历 region 中的 reads
    bamfile = pysam.AlignmentFile(bam_file, "rb")

    for read in bamfile.fetch(region=region):
#        if read.is_unmapped or read.is_secondary or read.is_supplementary:
#            continue
        
        # Soft-clipped reads
        if read.cigartuples and any(op == 4 for op, length in read.cigartuples):
            try:
                mate_read = bamfile.mate(read)
                soft_clipped_reads.append(read)
                mate_reads.append(mate_read)
                if read.is_reverse and not mate_read.is_reverse :
                    print
            except ValueError:
                print("Do not find mate")
                continue
        
        # Discordant reads and their mates
        elif (read.is_paired and not read.is_proper_pair) or abs(read.template_length) > max_insert_size:
            try:
                mate_read = bamfile.mate(read)
                discordant_reads.append(read)
                mate_reads.append(mate_read)
            except ValueError:
                print("Do not find mate")
                continue

    bamfile.close()

    return soft_clipped_reads, discordant_reads, mate_reads

def filter_reads(reads, min_quality=20):
    return [read for read in reads if read.mapping_quality >= min_quality]

def classify_reads_by_flag(reads):
    r1_reads = []
    r2_reads = []

    for read in reads:
        if read.is_read1:
            r1_reads.append(read)
        elif read.is_read2:
            r2_reads.append(read)

    return r1_reads, r2_reads

def save_paired_reads_to_fastq(forward_reads, reverse_reads, forward_output, reverse_output):
    with open(forward_output, "w") as forward_handle, open(reverse_output, "w") as reverse_handle:
        for f_read in forward_reads:
            forward_handle.write(f"@{f_read.query_name}\n")
            forward_handle.write(f"{f_read.query_sequence}\n")
            forward_handle.write("+\n")
            forward_handle.write(f"{''.join([chr(q + 33) for q in f_read.query_qualities])}\n")
            
        for r_read in reverse_reads:
            reverse_handle.write(f"@{r_read.query_name}\n")
            reverse_handle.write(f"{r_read.query_sequence}\n")
            reverse_handle.write("+\n")
            reverse_handle.write(f"{''.join([chr(q + 33) for q in r_read.query_qualities])}\n")

bam_file = "/home/xuxf/HG002_bwa_sort.bam"
region = "10:108991227-108991327"
max_insert_size = 1000
soft_clipped_reads, discordant_reads, mate_reads = extract_reads_with_mates(bam_file, region, max_insert_size)


# 获取所有的 mate reads

# 分类 reads 为 R1 和 R2
r1_soft_clipped_reads, r2_soft_clipped_reads = classify_reads_by_flag(soft_clipped_reads)
r1_discordant_reads, r2_discordant_reads = classify_reads_by_flag(discordant_reads)
r1_mate_reads, r2_mate_reads = classify_reads_by_flag(mate_reads)

# 合并所有 R1 和 R2 reads
r1_reads = r1_soft_clipped_reads + r1_discordant_reads + r1_mate_reads
r2_reads = r2_soft_clipped_reads + r2_discordant_reads + r2_mate_reads

directory_path = "spades_output_"+region

# 创建目录
create_directory(directory_path)
# 将 reads 分开保存
save_paired_reads_to_fastq(r1_reads, r2_reads, 'spades_output_'+region+'/reads_1.fastq', 'spades_output_'+region+'/reads_2.fastq')



spades_command = [
    "spades.py",
    "-1", "spades_output_"+region+"/reads_1.fastq",
    "-2", "spades_output_"+region+"/reads_2.fastq",
    "-o", "spades_output_"+region
]
subprocess.run(spades_command)

