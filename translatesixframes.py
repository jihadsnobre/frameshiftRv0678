#!/usr/bin/python3
from Bio.Seq import Seq
import pandas as pd
import warnings
from Bio import BiopythonWarning
warnings.simplefilter('ignore', BiopythonWarning)


# Load the original excel file
df = pd.read_excel("/Users/jihadsnobre1/Library/CloudStorage/OneDrive-ITG/PhD/computationalproteomics/databases/WHO_nt_info_fs_sequence.xlsx")
df["Mutated_Sequence"] = df["Mutated_Sequence"].astype(str)


# Define a function that translates a nucleotide sequence to a protein sequence
from Bio.Seq import Seq

def translate_seq(seq):
    # Convert input sequence to string if not already a string
    if not isinstance(seq, str):
        seq = str(seq)
    # Use the BioPython library to translate the sequence
    protein_seq = Seq(seq).translate(table="Bacterial")
    return str(protein_seq)


# Iterate through each row in the dataframe and translate the sequence in all six frames
for i, row in df.iterrows():
    seq = row["Mutated_Sequence"]
    frame_1 = translate_seq(seq)
    frame_2 = translate_seq(seq[1:])
    frame_3 = translate_seq(seq[2:])
    seq = Seq(seq)
    frame_4 = translate_seq(seq.reverse_complement())
    frame_5 = translate_seq(seq.reverse_complement()[1:])
    frame_6 = translate_seq(seq.reverse_complement()[2:])
    
    # Add the translated sequences as new columns in the dataframe
    df.at[i, "frame 1"] = frame_1
    df.at[i, "frame 2"] = frame_2
    df.at[i, "frame 3"] = frame_3
    df.at[i, "frame 4"] = frame_4
    df.at[i, "frame 5"] = frame_5
    df.at[i, "frame 6"] = frame_6

# Save the modified dataframe to a new excel file
df.to_excel("/Users/jihadsnobre1/Library/CloudStorage/OneDrive-ITG/PhD/computationalproteomics/databases/WHO_nt_info_fs_sequence_sixframes.xlsx", index=False)
