
import pandas as pd
import numpy as np

# Define the wild type sequence
wild_type_sequence = 'gtgagcgtcaacgacggggtcgatcagatgggcgccgagcccgacatcatggaattcgtcgaacagatgggcggctatttcgagtccaggagtttgactcggttggcgggtcgattgttgggctggctgctggtgtgtgatcccgagcggcagtcctcggaggaactggcgacggcgctggcggccagcagcggggggatcagcaccaatgcccggatgctgatccaatttgggttcattgagcggctcgcggtcgccggggatcggcgcacctatttccggttgcggcccaacgctttcgcggctggcgagcgtgaacgcatccgggcaatggccgaactgcaggacctggctgacgtggggctgagggcgctgggcgacgccccgccgcagcgaagccgacggctgcgggagatgcgggatctgttggcatatatggagaacgtcgtctccgacgccctggggcgatacagccagcgaaccggagaggacgactga'

# Function to apply insertions or deletions to the wild type sequence
def apply_mutations(wild_type, mutation_type, position, inserted_deleted):
    try:
        position = int(position)  # Convert position to integer
    except ValueError:
        print(f"Error converting position {position} to integer. Skipping mutation.")
        return np.nan  # Return NaN for invalid position

    mutated_sequence = list(wild_type)

    if mutation_type == 'ins':
        if not pd.isna(position):  # Check for NaN before inserting
            mutated_sequence.insert(position - 1, inserted_deleted)
    elif mutation_type == 'del':
        if not pd.isna(position) and not pd.isna(inserted_deleted):  # Check for NaN before deleting
            del mutated_sequence[position - 1:position - 1 + len(inserted_deleted)]

    return ''.join(str(x) for x in mutated_sequence)

# Read the Excel file
file_path = '/Users/jihadsnobre1/Library/CloudStorage/OneDrive-ITG/PhD/computationalproteomics/databases/MGIT_database.xlsx'
df = pd.read_excel(file_path)

# Apply mutations and create a new column for mutated sequences
df['Mutated_Sequence'] = df.apply(lambda row: apply_mutations(wild_type_sequence, row['Type_indel'], row['Nucleotide'], row['inserted_deleted']), axis=1)

# Save the mutated data to a new Excel file
output_file_path = '/Users/jihadsnobre1/Library/CloudStorage/OneDrive-ITG/PhD/computationalproteomics/databases/MGIT_databaset_sequence.xlsx'
df.to_excel(output_file_path, index=False)

print(f"Mutated data saved to {output_file_path}")
