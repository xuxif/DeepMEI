import sys

def process_data():
    records = {}

    for line in sys.stdin:
        parts = line.strip().split()
        key = parts[0]  # First column as the key
        fourth, fifth = map(int, parts[3:5])  # Convert the fourth and fifth columns to integers
        distance = abs(fifth - fourth)

        if key not in records:
            # Store the entire record and the min/max values
            records[key] = [line.strip(), min(fourth, fifth), max(fourth, fifth), distance]
        else:
            # Update the record if a longer distance is found
            if distance > records[key][3]:
                records[key][0] = line.strip()
                records[key][3] = distance
            # Always update min and max values
                records[key][1] = min(fourth, fifth)
                records[key][2] = max( fourth, fifth)

    for key, value in records.items():
        parts = value[0].split()
        # Replace the fourth and fifth columns with the min and max values
        parts[3] = str(value[1])
        parts[4] = str(value[2])
        parts.append(str(value[3]))
        print('\t'.join(parts))

# Call the function
process_data()

