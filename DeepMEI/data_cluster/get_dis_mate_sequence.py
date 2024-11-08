import pysam
import sys

def get_mate_sequence(bam_file, read):
    """ Get the mate read's sequence """
    with pysam.AlignmentFile(bam_file, "rb") as bam_t:
        try:
            mate = bam_t.mate(read)
            return mate.query_sequence
        except ValueError:
            # Mate not found or read is unmapped
            return None

def extract_reads(bam_file, target_region, left, right,REF,step=1):
    """ Extract reads from a BAM file within a target region """
    file_type="rb"
    dis_count=10
    dis_distance=0
    if bam_file.lower().endswith('.cram') :
        file_type="rc"
    with pysam.AlignmentFile(bam_file, file_type,reference_filename=REF) as bam:
        chrom, region = target_region.split(':')
        if step == 2:
            dis_count=40
            dis_distance=1000
            if left < right:
                start=int(left)
                end=int(right)
            else:
                start=int(right)
                end=int(left)
        else:
            start, end = map(int, region.split('-'))
        
        read_reverse=0
        read_forward=0
        for read in bam.fetch(chrom, start-dis_distance, end+dis_distance):
            if abs(read.template_length) > 1000 or (read.next_reference_id != -1 and read.reference_id != read.next_reference_id):
                end_position = read.reference_end
                start_position=read.reference_start
                if (not read.is_reverse and end_position <= right+5) or (read.is_reverse and start_position >= left-5):
                    read_order=2
                    if read.is_read1:
                        read_order=1
                    if read.is_reverse: 
                        read_reverse+=1
                    else:
                        read_forward+=1
                    mate_pos_start=read.next_reference_start-1
                    mate_pos_end=read.next_reference_start+1
                    if mate_pos_start<=0:
                        mate_pos_start=1
                    if mate_pos_end<=0:
                        mate_pos_end=1
                    read_strand='left'
                    if read.is_reverse :
                       read_strand='right'
                    print(f"{target_region}\t{read.query_name}\t{read_order}\t{read.next_reference_name}\t{mate_pos_start}\t{mate_pos_end}\t{read_strand}")
            if read_reverse>=dis_count and read_forward>=dis_count:
                break

# Example usage
if len(sys.argv) != 7:
    print("Usage: python script.py <bam_file> <target_region> <left_position> <right_position>")
    sys.exit(1)

bam_file = sys.argv[1]
target_region = sys.argv[2]
left_position = int(sys.argv[3])
right_position = int(sys.argv[4])
REF= sys.argv[5]
step=int(sys.argv[6])

extract_reads(bam_file, target_region, left_position, right_position,REF,step)

