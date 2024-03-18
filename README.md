Overview

This repository contains the scripts and datasets utilized in the study "Frameshift Mutations in Rv0678 Preserve Bedaquiline Susceptibility in Mycobacterium tuberculosis by Maintaining Protein Integrity." The research investigates the role of frameshift mutations in the Rv0678 gene and their impact on the efficacy of bedaquiline treatment for tuberculosis. This README outlines the analytical processes and datasets involved.

R Script

Frameshift_MS_submission.R: This script encompasses all analyses presented in the manuscript, segmented according to the paper's results sections. It includes data preprocessing, statistical analyses, and visualizations pertinent to the study's findings.

Datasets
The analysis uses several datasets, including:

WHO_mut: Genotypic data featuring sample IDs and detailed mutation characteristics, obtained from WHO_git_mut.txt on the WHO catalogue GitHub.
WHO_pheno: Phenotypic data corresponding to the samples, sourced from WHO_pheno.txt on the WHO catalogue github.
WHO-UCN-TB-2023.5-eng.xlsx: excel file extracted from the WHO catalogue, detailing protein-level mutations related to bedaquiline and their interpretations based on the WHO catalogue algorithm.
WHO_nt_info.xlsx: Can be found in the WHO catalogue 2023 update, this file provides  nucleotide information necessary for analyzing late stop codons and calculating alignment scores.
MGIT_database (second sheet of the supplementary tables) contains isolates extracted from the literature with available MIC in MGIT

Python Scripts for analysis on stop codons and alignement scores

Additional analyses on stop codon counts and alignment scores are performed using the following python scripts, in the specified order:

mutate_sequence.py for the MGIT database: Generates mutated nucleotide sequences for each mutation.
mutate_sequence_REFALTNucleotide.py for the WHO database: Similar to the above, but tailored for the WHO database.
translatesixframes.py: Translates mutated nucleotide sequences into six amino acid reading frames.
countaabeforestop.py: Counts the number of codons before the stop codon in each reading frame.
align.py: Calculates alignment scores compared to the wild type for each reading frame.
