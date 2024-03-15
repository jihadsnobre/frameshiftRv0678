#!/usr/bin/python3
import pandas as pd
import pandas as pd
from Bio.Seq import Seq

# Define a function to count amino acids before stop codon
def count_aa_before_stop(seq):
    if pd.isna(seq):
        return float("nan")
    else:
        seq_obj = Seq(seq)
        stop_index = seq_obj.find('*')
        if stop_index == -1:
            return len(seq_obj)
        else:
            return stop_index

# Read in the Excel file
df = pd.read_excel('/Users/jihadsnobre1/Library/CloudStorage/OneDrive-ITG/PhD/computationalproteomics/databases/Claudio_fs_sequence_sixframes.xlsx')

# Create new columns for the amino acid counts
df['frame 1 aa count'] = df['frame 1'].apply(count_aa_before_stop)
df['frame 2 aa count'] = df['frame 2'].apply(count_aa_before_stop)
df['frame 3 aa count'] = df['frame 3'].apply(count_aa_before_stop)
df['frame 4 aa count'] = df['frame 4'].apply(count_aa_before_stop)
df['frame 5 aa count'] = df['frame 5'].apply(count_aa_before_stop)
df['frame 6 aa count'] = df['frame 6'].apply(count_aa_before_stop)

# Write the output to a new Excel file
df.to_excel('/Users/jihadsnobre1/Library/CloudStorage/OneDrive-ITG/PhD/computationalproteomics/databases/Claudio_fs_sequence_sixframes_stop.xlsx', index=False)

