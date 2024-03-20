import pandas as pd
from Bio import pairwise2
from Bio.Seq import Seq

# Load the dataset from Excel
excel_file = "/Users/jihadsnobre1/Library/CloudStorage/OneDrive-ITG/PhD/computationalproteomics/databases/WHO_nt_info_fs_sequence_sixframes_stop.xlsx"
df = pd.read_excel(excel_file)

# Define the wild-type sequence
wild_type_sequence = "VSVNDGVDQMGAEPDIMEFVEQMGGYFESRSLTRLAGRLLGWLLVCDPERQSSEELATALAASSGGISTNARMLIQFGFIERLAVAGDRRTYFRLRPNAFAAGERERIRAMAELQDLADVGLRALGDAPPQRSRRLREMRDLLAYMENVVSDALGRYSQRTGEDD"
wild_type_sequence = Seq(wild_type_sequence)

# Function to calculate conservation percentage
def calculate_conservation(sequence):
    alignment = pairwise2.align.globalxx(wild_type_sequence, sequence, one_alignment_only=True)[0]
    conservation_percentage = (alignment.score / len(wild_type_sequence)) * 100
    return conservation_percentage

# List of frame columns to process
frame_columns = ["frame 1", "frame 2", "frame 3", "frame 4", "frame 5", "frame 6"]

# Loop over each frame column and apply the alignment function
for frame_column in frame_columns:
    # Apply the function to the current frame column
    df[f"{frame_column}_alignment"] = df[frame_column].apply(calculate_conservation)

# Write the updated dataframe to a new Excel file
output_excel_file = "/Users/jihadsnobre1/Library/CloudStorage/OneDrive-ITG/PhD/computationalproteomics/databases/WHO_nt_info_fs_sequence_sixframes_stop_align.xlsx"
df.to_excel(output_excel_file, index=False)

print(f"Results written to {output_excel_file}")
