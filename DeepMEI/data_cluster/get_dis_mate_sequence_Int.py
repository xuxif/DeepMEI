import pysam
import sys

def extract_high_quality_bases(read, threshold=30):
    sequence = read.query_sequence
    qualities = read.query_qualities
    
    def quality_check(index, qualities, threshold):
        start = max(0, index - 5)
        end = min(len(qualities), index + 6)  # +6 because range end is exclusive
        low_quality_count = sum(1 for i in range(start, end) if qualities[i] < threshold)
        return low_quality_count < 5    

    high_quality_bases = []
    
    # Determine the direction based on read type and strand orientation
    if (read.is_read1 and not read.is_reverse) or (read.is_read2 and read.is_reverse):
        # Read 1 on forward strand or Read 2 on reverse strand: left-to-right
        for i in range(len(sequence)):
            if not quality_check(i, qualities, threshold):
                break
            high_quality_bases.append(sequence[i])
    elif (read.is_read2 and not read.is_reverse) or (read.is_read1 and read.is_reverse):
        # Read 2 on forward strand or Read 1 on reverse strand: right-to-left
        for i in range(len(sequence) - 1, -1, -1):
            if not quality_check(i, qualities, threshold):
                break
            high_quality_bases.append(sequence[i])
        # Since we collected bases from right-to-left, reverse them
        high_quality_bases = high_quality_bases[::-1]
    else:
        raise ValueError("Invalid read orientation.")
    
    return ''.join(high_quality_bases)


def find_and_extend(short_seq, long_seqs,direction='left', min_extend=15):
    def get_extension(seq, long_seqs, direction,rm_index):
        for i, long_seq in enumerate(long_seqs):
            pos = long_seq.find(seq)
            if pos != -1 and i not in rm_index:
                if direction == 'left' and pos > 0:
                    return long_seq[:pos], i  # 返回从起始位置到找到位置之前的序列以及索引
                elif direction == 'right' and pos + len(seq) < len(long_seq):
                    return long_seq[pos + len(seq):], i  # 返回从找到位置之后到结尾的序列以及索引
        return None, None
    ori_short_seq=short_seq
    extended_seq = short_seq
    rm_index=[]
    while True:
        if direction == 'left':
            short_seq = extended_seq[:min_extend]  # 提取新的短序列
        elif direction == 'right':
            short_seq = extended_seq[-min_extend:]  # 提取新的短序列

        extension, index = get_extension(short_seq, long_seqs,direction,rm_index)

        if not extension :
            if len(extended_seq)<100:
                if direction == 'left':
                    extended_seq=extended_seq[5:]  # 提取新的短序列
                    short_seq = extended_seq[:min_extend]  # 提取新的短序列
                elif direction == 'right':
                    extended_seq= extended_seq[:-5]  # 提取新的短序列
                    short_seq = extended_seq[-min_extend:]  # 提取新的短序列
                if len(short_seq) <10:
                    break
            else:
                break
        else:
            if direction == 'left':
                extended_seq = extension + extended_seq
            elif direction == 'right':
                extended_seq = extended_seq + extension
        # 从列表中移除已匹配的长序列
        rm_index.append(index)
    if len(ori_short_seq)>len(extended_seq):
        extended_seq=ori_short_seq
    return extended_seq

def reverse_complement(dna_sequence):
    # 定义碱基互补对
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C','N':'N','X':'X'}
    
    # 生成互补序列
    complement_sequence = ''.join(complement[base] for base in dna_sequence)
    
    # 反转互补序列
    reverse_complement_sequence = complement_sequence[::-1]
    
    return reverse_complement_sequence
def get_mate_sequence(bam_file, read):
    """ Get the mate read's sequence """
    with pysam.AlignmentFile(bam_file, "rb") as bam_t:
        try:
            mate = bam_t.mate(read)
            return mate.query_sequence
        except ValueError:
            # Mate not found or read is unmapped
            return None

def extract_reads(bam_file, target_region, left, right,REF):
    """ Extract reads from a BAM file within a target region """
    file_type="rb"
    if bam_file.lower().endswith('.cram') :
        file_type="rc"
    mate_reads = []
    seqs=[]
    with pysam.AlignmentFile(bam_file, file_type,reference_filename=REF) as bam:
        chrom, region = target_region.split(':')
        start, end = map(int, region.split('-'))
        read_reverse=0
        read_forward=0
        for read in bam.fetch(chrom, start-1000, end+1000):
            if abs(read.template_length) > 1000 or (read.next_reference_id != -1 and read.reference_id != read.next_reference_id):
                end_position = read.reference_end
                start_position=read.reference_start
                if (not read.is_reverse and end_position <= right+5) or (read.is_reverse and start_position >= left-5):
                    mate_reads.append(read)
        for read in mate_reads:
            try:
                mate_read = bam.mate(read)
                sequence=extract_high_quality_bases(mate_read)
                if (read.is_reverse and mate_read.is_reverse) or ((not read.is_reverse) and (not mate_read.is_reverse)):
                    sequence=reverse_complement(sequence)

                seqs.append(sequence)
#                print(sequence)
            except ValueError:
                continue
    return seqs

#            if read_reverse>=10 and read_forward>=10:
#                break

# Example usage
if len(sys.argv) != 4:
    print("Usage: python script.py <bam_file> <bp file> <reference file>")
    sys.exit(1)

bam_file = sys.argv[1]
dis_file= sys.argv[2]
REF= sys.argv[3]
parts=[]
#6:158638748-158638848	158638840	158638856	16	21	18	0	4	7	-1	ACCCAAGAATTATCAATAAAAAAATAAATTAAAAAAAAAAAAAAAAAAA	AGTCTTTGCCGCCGCGCCGGCGAGCGCCGCCCGGGAGGCAGCGGCTGGAGGAGCGGACGGGC

#HG002_X:35027996-35028096_BPinfo.txt

with open('HG002_'+dis_file[8:]+'_BPinfo.txt', "r") as file:
    for line in file:
        parts = line.strip().split()

target_region = parts[0]
left_position = int(parts[1])
right_position = int(parts[2])
short_seq_left=parts[10]
short_seq_right=parts[11]

long_seqs=[]

try: 
    with open(dis_file, "r") as file:
        for line in file:
            long_seqs.append(line.strip())
except FileNotFoundError:
    sys.exit(0)
except Exception as e:
    sys.exit(0)

#long_seqs=extract_reads(bam_file, target_region, left_position, right_position,REF)
longest_extended_seq_left = find_and_extend(short_seq_left, long_seqs, direction='left')

right_comp=0
left_comp=0

pos = longest_extended_seq_left.find(short_seq_right[:15])
if pos != -1:
    longest_extended_seq_left=longest_extended_seq_left[pos:]
    left_comp=1

longest_extended_seq_right= find_and_extend(short_seq_right, long_seqs, direction='right')

pos = longest_extended_seq_right.find(short_seq_left[-15:])
if pos != -1:
    longest_extended_seq_right=longest_extended_seq_right[:pos+len(short_seq_left)]
    right_comp=1



print('>'+target_region+'_left\t'+str(left_comp)+'\n'+longest_extended_seq_left)
print('>'+target_region+'_right\t'+str(right_comp)+'\n'+longest_extended_seq_right)
