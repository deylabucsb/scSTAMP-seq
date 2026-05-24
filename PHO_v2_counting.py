import sys
import matplotlib.pyplot as plt

def count_barcode_and_sequence_occurrences(trimmed_file, barcode_file, inner_sequence, outer_sequence):
    barcode_counts = {}
    sequence_counts = {'inner_sequence': 0, 'outer_sequence': 0, 'none': 0}
    inner_correct_next_4_bp = 0
    outer_correct_next_4_bp = 0
    inner_incorrect_next_4_bp = 0
    outer_incorrect_next_4_bp = 0

    # Read the barcodes from the text file
    barcode_to_number = {}
    with open(barcode_file, 'r') as f:
        for line in f:
            number, barcode = line.strip().split('\t')
            barcode_to_number[barcode] = number
            barcode_counts[barcode] = {
                'count': 0, 'IN_PHO1': 0,'OUT_PHO1': 0,'IN_PHO2': 0, 'OUT_PHO2': 0,
                'IN_PHO9': 0, 'OUT_PHO9': 0, 'IN_PHO10': 0,
                'OUT_PHO10': 0, 'incorrect_count': 0,
                'UMIs_IN': set(), 'UMIs_OUT': set()  # Set to store unique UMIs for each barcode
            }

    # Read the trimmed file and process each line
    with open(trimmed_file, 'r') as f:
        lines = f.readlines()

    for idx, line in enumerate(lines):
        first_8bp = line[:8]
        
        if first_8bp in barcode_counts:
            # First UMI: 6 bp after the initial 8 bp
            first_umi = line[8:14]
            
            # Initialize unique UMI tuple
            unique_umi_tuple_IN = None
            unique_umi_tuple_OUT = None
            
            # Check if the inner sequence is present on the same line
            if inner_sequence in line:
                sequence_counts['inner_sequence'] += 1

                # Get next 4 bp and second UMI after inner sequence
                next_4_bp_start = line.index(inner_sequence) + len(inner_sequence)
                next_4_bp = line[next_4_bp_start:next_4_bp_start + 4]
                second_umi = line[next_4_bp_start + 4:next_4_bp_start + 10]

                unique_umi_tuple_IN = (first_umi, second_umi)

                # Only count if the UMI pair is unique
                if unique_umi_tuple_IN not in barcode_counts[first_8bp]['UMIs_IN']:
                    barcode_counts[first_8bp]['count'] += 1
                    barcode_counts[first_8bp]['UMIs_IN'].add(unique_umi_tuple_IN)

                    # Increment specific inner sequence counts based on next 4 bp
                    if next_4_bp == 'TGAC':
                        barcode_counts[first_8bp]['IN_PHO1'] += 1
                        inner_correct_next_4_bp += 1
                    elif next_4_bp == 'AGGT':
                        barcode_counts[first_8bp]['IN_PHO2'] += 1
                        inner_correct_next_4_bp += 1
                    elif next_4_bp == 'GACA':
                        barcode_counts[first_8bp]['IN_PHO9'] += 1
                        inner_correct_next_4_bp += 1
                    elif next_4_bp == 'CTGA':
                        barcode_counts[first_8bp]['IN_PHO10'] += 1
                        inner_correct_next_4_bp += 1
                    else:
                        inner_incorrect_next_4_bp += 1
                        barcode_counts[first_8bp]['incorrect_count'] += 1

            # Check if the outer sequence is present on the same line
            elif outer_sequence in line:
                sequence_counts['outer_sequence'] += 1

                # Get next 4 bp and second UMI after outer sequence
                next_4_bp_start = line.index(outer_sequence) + len(outer_sequence)
                next_4_bp = line[next_4_bp_start:next_4_bp_start + 4]
                second_umi = line[next_4_bp_start + 4:next_4_bp_start + 10]

                unique_umi_tuple_OUT = (first_umi, second_umi)

                # Only count if the UMI pair is unique
                if unique_umi_tuple_OUT not in barcode_counts[first_8bp]['UMIs_OUT']:
                    barcode_counts[first_8bp]['count'] += 1
                    barcode_counts[first_8bp]['UMIs_OUT'].add(unique_umi_tuple_OUT)

                    # Increment specific outer sequence counts based on next 4 bp
                    if next_4_bp == 'GTGA':
                        barcode_counts[first_8bp]['OUT_PHO1'] += 1
                        inner_correct_next_4_bp += 1
                    elif next_4_bp == 'CAGA':
                        barcode_counts[first_8bp]['OUT_PHO2'] += 1
                        outer_correct_next_4_bp += 1
                    elif next_4_bp == 'ATCG':
                        barcode_counts[first_8bp]['OUT_PHO9'] += 1
                        outer_correct_next_4_bp += 1
                    elif next_4_bp == 'TAGT':
                        barcode_counts[first_8bp]['OUT_PHO10'] += 1
                        outer_correct_next_4_bp += 1
                    else:
                        outer_incorrect_next_4_bp += 1
                        barcode_counts[first_8bp]['incorrect_count'] += 1

            else:
                sequence_counts['none'] += 1
                barcode_counts[first_8bp]['incorrect_count'] += 1

    # Remove barcode counts with counts less than 10000
    barcode_counts_filtered = {barcode: counts for barcode, counts in barcode_counts.items() if counts['count'] >= 0}

    return barcode_counts_filtered, sequence_counts, inner_correct_next_4_bp, outer_correct_next_4_bp, inner_incorrect_next_4_bp, outer_incorrect_next_4_bp, barcode_to_number

# File paths
trimmed_file = ''
barcode_file = ''
inner_sequence = 'ACGTAACCGCGCCGGA'
outer_sequence = 'TTGAGGTAGTGTGGAG'

# Count occurrences of each barcode and sequences
barcode_counts, sequence_counts, inner_correct_next_4_bp, outer_correct_next_4_bp, inner_incorrect_next_4_bp, outer_incorrect_next_4_bp, barcode_to_number = count_barcode_and_sequence_occurrences(trimmed_file, barcode_file, inner_sequence, outer_sequence)

# Save output to a file
with open('GradP8L4_counts_UMI.txt', 'w') as f:
    for barcode, counts in barcode_counts.items():
        barcode_number = barcode_to_number.get(barcode, "Unknown")
        f.write(f"{barcode_number}\t{counts['count']}\t{counts['IN_PHO1']}\t{counts['OUT_PHO1']}\t{counts['IN_PHO2']}\t{counts['OUT_PHO2']}\t{counts['IN_PHO9']}\t{counts['OUT_PHO9']}\t{counts['IN_PHO10']}\t{counts['OUT_PHO10']}\t{counts['incorrect_count']}\n")
  