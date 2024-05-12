import sys

#soft|10|101340285|left	AluYk4	SINE/Alu	261	298 30
#soft|10|101340326|right	Alu	SINE/Alu	176	219 50
#soft|10|10135357|left	Alu	SINE/Alu	1	50 20
def read_saved_file(filename):
    saved_records = {}
    with open(filename, 'r') as file:
        for line in file:
            parts = line.strip().split()
            saved_records[parts[0]] = line.strip()
    return saved_records

def process_input_and_find_matches(saved_records):
    for line in sys.stdin:
        parts = line.strip().split()
        key_left = f"soft|{parts[0]}|{parts[3]}|left"
        key_right = f"soft|{parts[0]}|{parts[4]}|right"

        match_left = saved_records.get(key_left)
        match_right = saved_records.get(key_right)
        direct=[]
        left_cord=[]
        right_cord=[]
        dis_left=0
        dis_right=0
        me_type="none"
        cord=[]
        if match_left or match_right:
            if match_left:
               direct.append('left')
               match_left=match_left.split()
               cord.append(int(match_left[3]))
               cord.append(int(match_left[4]))
               dis_left=int(match_left[5])
               me_type_left=f"{match_left[1]}:{match_left[2]}"
            if match_right:
               direct.append('right')
               match_right=match_right.split()
               cord.append(int(match_right[3]))
               cord.append(int(match_right[4]))
               dis_right=int(match_right[5])
               me_type_right=f"{match_right[1]}:{match_right[2]}"

            if dis_left>dis_right:
                me_type=me_type_left
            else:
                me_type=me_type_right
            cord_max=max(cord)
            cord_min=min(cord)
            direct='-'.join(direct)
            print(f"{line.strip()}\t{me_type}:{cord_min}-{cord_max}:{direct}")
        else:
            print(f"{line.strip()}\tnone")


# Example usage
saved_file = "clip_ME_clue.txt"  # Replace with the path to your saved file
#saved_file = "clip_ME_clue_t.txt"  # Replace with the path to your saved file
saved_records = read_saved_file(saved_file)
process_input_and_find_matches(saved_records)

