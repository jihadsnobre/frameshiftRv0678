import pandas as pd

# Function to substitute nucleotides in the sequence
def substitute_nucleotides(wild_type_seq, position, ref, alt):
    """
    Substitute the nucleotide in the wild type sequence at the given position
    with the alternative nucleotide. The position is assumed to be 1-based indexing.
    """
    # Adjust for 0-based indexing in Python
    position -= 1

    # Check if the reference nucleotide matches the wild type sequence at the given position
    if wild_type_seq[position:position+len(ref)] == ref:
        # Perform the substitution
        return wild_type_seq[:position] + alt + wild_type_seq[position+len(ref):]
    else:
        # If the reference nucleotide does not match, return the original sequence
        # This might indicate an error in the data
        return wild_type_seq

# Load the Excel file
file_path = '/Users/jihadsnobre1/Library/CloudStorage/OneDrive-ITG/PhD/computationalproteomics/databases/WHO_nt_info_fs.xlsx'
df = pd.read_excel(file_path)

# Wild type sequence provided
wild_type_sequence = "gtgagcgtcaacgacggggtcgatcagatgggcgccgagcccgacatcatggaattcgtcgaacagatgggcggctatttcgagtccaggagtttgactcggttggcgggtcgattgttgggctggctgctggtgtgtgatcccgagcggcagtcctcggaggaactggcgacggcgctggcggccagcagcggggggatcagcaccaatgcccggatgctgatccaatttgggttcattgagcggctcgcggtcgccggggatcggcgcacctatttccggttgcggcccaacgctttcgcggctggcgagcgtgaacgcatccgggcaatggccgaactgcaggacctggctgacgtggggctgagggcgctgggcgacgccccgccgcagcgaagccgacggctgcgggagatgcgggatctgttggcatatatggagaacgtcgtctccgacgccctggggcgatacagccagcgaaccggagaggacgactga"

# Apply the substitution function to each row in the dataframe
df['MUTATED_SEQUENCE'] = df.apply(lambda row: substitute_nucleotides(wild_type_sequence, row['Nucleotide'], row['REF'], row['ALT']), axis=1)

# Optionally, save the modified dataframe to a new Excel file
output_file_path = '/Users/jihadsnobre1/Library/CloudStorage/OneDrive-ITG/PhD/computationalproteomics/databases/WHO_nt_info_fs_sequence.xlsx'
df.to_excel(output_file_path, index=False)