######WHO catalogue- Impact of the position of frameshift mutations in the Rv0678 gene on BDQ phenotype
library(dplyr)
library(ggplot2)
library(readxl)
library(tidyverse)

WHO_mut <- read.csv("WHO_git_mut.txt", header = TRUE, sep = ",")
WHO_pheno <- read.csv("WHO_pheno.txt", header = TRUE, sep = ",")

WHO = merge(WHO_mut, WHO_pheno, by = "sample_id")
WHO <- WHO %>%
  filter(phenotypic_category == "ALL" | phenotypic_category == "WHO") 
WHO <- WHO %>% 
  filter(resolved_symbol == "Rv0678")

WHO_solo <- WHO %>%
  group_by(sample_id) %>%
  filter(n_distinct(variant_category) == 1)

WHO_solo <- WHO_solo %>% 
  filter(max.af. >= 0.75)

table(WHO_solo$predicted_effect)
table(WHO_solo$phenotypic_category)

ggplot(WHO_solo, aes(x = predicted_effect)) +
  geom_bar() +
  labs(title = "WHO catalogue mutations in Rv0678",
       x = "Predicted effect",
       y = "Isolates count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

WHO_solo <- WHO_solo %>% 
  filter(predicted_effect %in% c("frameshift"))

WHO_solo <- WHO_solo %>%
  mutate(aa_pos = as.numeric(gsub("[^0-9]", "", variant_category)))

ggplot(WHO_solo, aes(x = aa_pos)) +
  geom_bar(stat = "count") +
  labs(title = "Position frameshift mutations in Rv0678 gene reported in WHO catalogue",
       x = "Rv0678 codons",
       y = "Count isolates n")

WHO_solo <- WHO_solo %>%
  mutate(domain = case_when(
    aa_pos < 16 ~ "before_dimerization",
    between(aa_pos, 16, 33) ~ "first_dimerization",
    between(aa_pos, 34, 99) ~ "DNA_binding",
    between(aa_pos, 100, 160) ~ "second_dimerization",
    aa_pos >= 161 ~ "after_dimerization",
    TRUE ~ "unknown"  # If none of the conditions are met
  ))

table(WHO_solo$domain)
table(WHO_solo$phenotype)

###statistical analysis domains

WHO_solo <- WHO_solo %>%
  mutate(DNA_binding_helix = ifelse(between(aa_pos, 36, 62), "yes", "no"))

WHO_solo <- WHO_solo %>%
  mutate(DNA_binding = ifelse(between(aa_pos, 34, 99), "yes", "no"))

WHO_solo <- WHO_solo %>%
  mutate(functional = ifelse(between(aa_pos, 16, 160), "yes", "no"))

WHO_solo <- WHO_solo %>%
  mutate(dimerization = ifelse(between(aa_pos, 16, 34) | between(aa_pos, 100, 160), "yes", "no"))

WHO_solo <- WHO_solo %>%
  mutate(first_dimerization = ifelse(between(aa_pos, 16, 34), "yes", "no"))

WHO_solo <- WHO_solo %>%
  mutate(second_dimerization = ifelse(between(aa_pos, 100, 160), "yes", "no"))

WHO_solo$functional <- factor(WHO_solo$functional, levels = c("yes", "no"))
contingency_table <- table(WHO_solo$functional, WHO_solo$phenotype)
fisher_result <- fisher.test(contingency_table)

WHO_solo$DNA_binding <- factor(WHO_solo$DNA_binding, levels = c("yes", "no"))
contingency_table <- table(WHO_solo$DNA_binding, WHO_solo$phenotype)
fisher_result <- fisher.test(contingency_table)

WHO_solo$dimerization <- factor(WHO_solo$dimerization, levels = c("yes", "no"))
contingency_table <- table(WHO_solo$dimerization, WHO_solo$phenotype)
fisher_result <- fisher.test(contingency_table)

WHO_solo$first_dimerization <- factor(WHO_solo$first_dimerization, levels = c("yes", "no"))
contingency_table <- table(WHO_solo$first_dimerization, WHO_solo$phenotype)
fisher_result <- fisher.test(contingency_table)

WHO_solo$second_dimerization <- factor(WHO_solo$second_dimerization, levels = c("yes", "no"))
contingency_table <- table(WHO_solo$second_dimerization, WHO_solo$phenotype)
fisher_result <- fisher.test(contingency_table)

WHO_solo$DNA_binding_helix <- factor(WHO_solo$DNA_binding_helix, levels = c("yes", "no"))
contingency_table <- table(WHO_solo$DNA_binding_helix, WHO_solo$phenotype)
fisher_result <- fisher.test(contingency_table)

##association analysis unique frameshift mutations with WHO interpretation (ass w resistance/uncertain significance)

WHO_unique <- read_excel("WHO-UCN-TB-2023.5-eng.xlsx", skip = 2)
WHO_unique <- WHO_unique %>%
  filter(drug == "Bedaquiline" & gene == "Rv0678")

WHO_unique <- WHO_unique %>% 
  filter(effect %in% c("frameshift"))

table(WHO_unique$`INITIAL CONFIDENCE GRADING`)

ggplot(WHO_unique, aes(x = `FINAL CONFIDENCE GRADING`)) +
  geom_bar() +
  labs(title = "Bar Plot of FINAL CONFIDENCE GRADING",
       x = "FINAL CONFIDENCE GRADING WHO",
       y = "Count")

ggplot(WHO_unique, aes(x = `FINAL CONFIDENCE GRADING`)) +
  geom_bar() +
  labs(title = "Bar Plot of FINAL CONFIDENCE GRADING",
       x = "FINAL CONFIDENCE GRADING",
       y = "Count") +
  theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1))

WHO_unique <- WHO_unique %>%
  mutate(aa_pos = as.numeric(gsub("[^0-9]", "", mutation)))

table(WHO_unique$domain)

WHO_unique$domain_category <- ifelse(WHO_unique$aa_pos < 16 | WHO_unique$aa_pos > 160, "no", "yes")
WHO_unique$domain_category <- factor(WHO_unique$domain_category, levels = c("yes", "no"))
contingency_table <- table(WHO_unique$domain_category, WHO_unique$`INITIAL CONFIDENCE GRADING`)
fisher_test_result <- fisher.test(contingency_table)
print(fisher_test_result)

WHO_unique$DNA_binding <- ifelse(WHO_unique$aa_pos >= 34 & WHO_unique$aa_pos <= 99, 
                                 "yes", 
                                 "no")
WHO_unique$DNA_binding <- factor(WHO_unique$DNA_binding, levels = c("yes", "no"))
contingency_table <- table(WHO_unique$DNA_binding, WHO_unique$`INITIAL CONFIDENCE GRADING`)
fisher_test_result <- fisher.test(contingency_table)


WHO_unique$helix <- ifelse(WHO_unique$aa_pos >= 36 & WHO_unique$aa_pos <= 62, 
                           "yes", 
                           "no")
WHO_unique$helix <- factor(WHO_unique$helix, levels = c("yes", "no"))
contingency_table <- table(WHO_unique$helix, WHO_unique$`INITIAL CONFIDENCE GRADING`)
fisher_test_result <- fisher.test(contingency_table)

WHO_unique$dimerization <- ifelse((WHO_unique$aa_pos >= 16 & WHO_unique$aa_pos <= 32) | 
                                    (WHO_unique$aa_pos >= 100 & WHO_unique$aa_pos <= 160), 
                                  "yes", 
                                  "no")
WHO_unique$dimerization <- factor(WHO_unique$dimerization, levels = c("yes", "no"))
contingency_table <- table(WHO_unique$dimerization, WHO_unique$`INITIAL CONFIDENCE GRADING`)
fisher_test_result <- fisher.test(contingency_table)
print(fisher_test_result)

####WHO catalogue- Late stop codons and conserved reading frames 

library(readxl)
library(dplyr)
library(ggplot2)
library(tidyr)

#adapt input from WHO raw file to be able to run python scripts to get codon counts before stop and aligment scores
WHO_frames <- read_excel("WHO_nt_info.xlsx")
WHO_frames <- WHO_frames %>%
  filter(startsWith(variant, "Rv0678"), endsWith(variant, "fs"))

WHO_frames <- WHO_frames %>%
  mutate(Nucleotide = position - 778989)

WHO_frames <- WHO_frames %>%
  mutate(REF = tolower(reference_nucleotide), 
         ALT = tolower(alternative_nucleotide))

writexl::write_xlsx(WHO_frames, "WHO_nt_info_fs.xlsx")

###now run python scripts to get codon counts before stop and aligment scores -> output is called WHO_nt_info_fs_sequence_sixframes_stop_align.xlsx

WHO_frames <- read_excel("WHO_nt_info_fs_sequence_sixframes_stop_align.xlsx")

WHO_frames$highest_count <- pmax(WHO_frames$`frame 1 aa count`, 
                                 WHO_frames$`frame 2 aa count`, 
                                 WHO_frames$`frame 3 aa count`, 
                                 WHO_frames$`frame 4 aa count`, 
                                 WHO_frames$`frame 5 aa count`, 
                                 WHO_frames$`frame 6 aa count`, 
                                 na.rm = TRUE)

WHO_frames$highest_alignment <- pmax(WHO_frames$`frame 1_alignment`, 
                                     WHO_frames$`frame 2_alignment`, 
                                     WHO_frames$`frame 3_alignment`, 
                                     WHO_frames$`frame 4_alignment`, 
                                     WHO_frames$`frame 5_alignment`, 
                                     WHO_frames$`frame 6_alignment`, 
                                     na.rm = TRUE)


stop_long <- WHO_frames %>%
  pivot_longer(cols = ends_with("count"), names_to = "frame", values_to = "stop") %>%
  mutate(frame = sub("aa count$", "", frame))

stop_long$stop
ggplot(stop_long, aes(x = stop)) +
  geom_histogram(bins = 30, fill = "yellowgreen", color = "black") +
  geom_vline(xintercept = 160, linetype = "dotted", color = "red", size = 1) +
  theme_minimal() +
  labs(title = "Distribution of count before stop codon",
       x = "Count codons before stop codon",
       y = "Number of unique frameshift mutations") +
  facet_wrap(~frame, ncol = 3, labeller = as_labeller(c('frame 1 ' = "Frame 1 (ORF)",
                                                        'frame 2 ' = "Frame 2 (ARF forward)",
                                                        'frame 3 ' = "Frame 3 (ARF forward)",
                                                        'frame 4 ' = "Frame 4 (ARF reverse)",
                                                        'frame 5 ' = "Frame 5 (ARF reverse)",
                                                        'frame 6 '= "Frame 6 (ARF reverse)",
                                                        highest_count = "Highest count before stop codon"))) +
  theme(strip.text.x = element_text(size = 10, face = "bold"))

count_after_160 <- stop_long %>%
  group_by(frame) %>%
  summarise(CountAfter160 = sum(stop > 160))
print(count_after_160)

# Get unique values from the 'frame' column
unique(stop_long$frame)

alignment_long <- WHO_frames %>%
  pivot_longer(cols = ends_with("alignment"), names_to = "frame", values_to = "alignment_score") %>%
  mutate(frame = sub("_alignment$", "", frame))

count_above_85 <- alignment_long %>%
  group_by(frame) %>%
  summarise(CountAbove85 = sum(alignment_score > 85))

num_greater_than_85 <- WHO_frames %>%
  filter(highest_alignment > 85) %>%
  nrow()
num_between_50_and_85 <- WHO_frames %>%
  filter(highest_alignment > 50 & highest_alignment <= 85) %>%
  nrow()

num_less_than_50 <- WHO_frames %>%
  filter(highest_alignment < 50) %>%
  nrow()

variants_with_high_alignment_in_multiple_frames <- WHO_frames %>%
  mutate(HighScoreFrames = (`frame 1_alignment` > 85) +
           (`frame 2_alignment` > 85) +
           (`frame 3_alignment` > 85) +
           (`frame 4_alignment` > 85) +
           (`frame 5_alignment` > 85) +
           (`frame 6_alignment` > 85)) %>%
  filter(HighScoreFrames >= 2) 

table(alignment_long$frame, alignment_long$alignment_score)

ggplot(alignment_long, aes(x = alignment_score)) +
  geom_histogram(bins = 30, fill = "darkgreen", color = "black") +
  geom_vline(xintercept = 90, linetype = "dotted", color = "red", size = 1) +
  theme_minimal() +
  labs(title = "Distribution of Alignment Scores",
       x = "Alignment Score",
       y = "Number of unique frameshift mutations") +
  facet_wrap(~frame, ncol = 3, labeller = as_labeller(c('frame 1' = "Frame 1 (ORF)",
                                                        'frame 2' = "Frame 2 (ARF forward)",
                                                        'frame 3' = "Frame 3 (ARF forward)",
                                                        'frame 4' = "Frame 4 (ARF reverse)",
                                                        'frame 5' = "Frame 5 (ARF reverse)",
                                                        'frame 6'= "Frame 6 (ARF reverse)",
                                                        'highest' = "Highest aligment score"))) +
  theme(strip.text.x = element_text(size = 10, face = "bold"))

#####MIC MGIT database analysis - Impact of late stop codons and conserved reading frames on BDQ MIC 
library(readxl)
library(dplyr)
library(ggplot2)

MGIT <- read_excel("MGIT_databaset_sequence_sixframes_stop_align.xlsx")

MGIT_filtered_not_truncated <- MGIT %>%
  filter(Truncated == "no")

MGIT_filtered_not_truncated <- MGIT_filtered_not_truncated %>%
  mutate(aa_pos = Nucleotide / 3)

ggplot(MGIT_filtered_not_truncated, aes(x = aa_pos)) +
  geom_bar(stat = "count") +
  labs(title = "Position frameshift mutations reported in MGIT dataset in the Rv0678 gene",
       x = "Codon Rv0678",
       y = "Count isolates containing a frameshift mutation")

# Create the 'highest_alignment' column
MGIT_filtered_not_truncated <- MGIT_filtered_not_truncated %>%
  mutate(
    highest_alignment = pmax(!!!select(., contains("alignment"))),
    highest_count = pmax(!!!select(., contains("count")))
  )

MGIT_filtered_not_truncated$MGIT <- as.numeric(MGIT_filtered_not_truncated$MGIT)
MGIT_filtered_not_truncated$aa_pos <- as.numeric(MGIT_filtered_not_truncated$aa_pos)

ggplot(MGIT_filtered_not_truncated, aes(x = aa_pos)) +
  geom_bar(stat = "count") +
  labs(title = "Position frameshift mutations reported in MGIT dataset in the Rv0678 gene",
       x = "Codon Rv0678",
       y = "Count isolates containing a frameshift mutation")

categories <- c("<16", "16-33", "34-100", "101-160", ">161")

MGIT_filtered_not_truncated$stop_codon_category <- cut(
  MGIT_filtered_not_truncated$`frame 1 aa count`,
  breaks = c(0, 16, 33, 100, 160, 166),
  labels = c("beginning of the protein (aa 0-16)", "1st dimerization domain (aa 17-33)", "DNA binding domain (aa 34-99)", "2nd dimerization domain (aa 100-160)", "end of the protein (aa 161-165)"),
  include.lowest = TRUE
)

MGIT_filtered_not_truncated$position_category <- cut(
  MGIT_filtered_not_truncated$aa_pos,
  breaks = c(0, 16, 33, 100, 160, 166),
  labels = c("beginning of the protein (aa 0-16)", "1st dimerization domain (aa 17-33)", "DNA binding domain (aa 34-99)", "2nd dimerization domain (aa 100-160)", "end of the protein (aa 161-165)"),
  include.lowest = TRUE
)

MGIT_filtered_not_truncated$MGIT = as.numeric(MGIT_filtered_not_truncated$MGIT)

MGIT_filtered_not_truncated$MIC_interpretation <- ifelse(MGIT_filtered_not_truncated$MGIT < 1.1, "S", "R")

table(MGIT_filtered_not_truncated$stop_codon_category, MGIT_filtered_not_truncated$MIC_interpretation)

ggplot(MGIT_filtered_not_truncated, aes(x = stop_codon_category, y = MGIT)) +
  geom_boxplot(fill = "pink") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red", size = 0.5) +
  labs(title = "Boxplot position stop codon vs MIC",
       x = "Position stop codon induced by frameshift mutation",
       y = "BDQ MIC (µg/mL)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

ggplot(MGIT_filtered_not_truncated, aes(x = position_category, y = MGIT)) +
  geom_boxplot(fill = "pink") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red", size = 0.5) +
  labs(title = "Boxplot position frameshift mutation vs MIC",
       x = "Position  frameshift mutation",
       y = "MIC (µg/mL)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

#statistical analysis
MGIT_filtered_not_truncated$end <- ifelse(MGIT_filtered_not_truncated$stop_codon_category == "end of the protein (aa 161-165)", "yes", "no")
MGIT_filtered_not_truncated$end <- factor(MGIT_filtered_not_truncated$end, levels = c("yes", "no"))
MGIT_filtered_not_truncated$MIC_interpretation <- ifelse(MGIT_filtered_not_truncated$MGIT < 1.1, "S", "R")
MGIT_filtered_not_truncated$MIC_interpretation <- factor(MGIT_filtered_not_truncated$MIC_interpretation, levels = c("S", "R"))
contingency_table <- table(MGIT_filtered_not_truncated$end, MGIT_filtered_not_truncated$MIC_interpretation)
fisher.test(contingency_table)

MGIT_filtered_not_truncated$first_dimerization <- ifelse(MGIT_filtered_not_truncated$stop_codon_category == "1st dimerization domain (aa 17-33)", "yes", "no")
MGIT_filtered_not_truncated$first_dimerization <- factor(MGIT_filtered_not_truncated$first_dimerization, levels = c("yes", "no"))
contingency_table <- table(MGIT_filtered_not_truncated$first_dimerization, MGIT_filtered_not_truncated$MIC_interpretation)
fisher.test(contingency_table)

MGIT_filtered_not_truncated$second_dimerization <- ifelse(MGIT_filtered_not_truncated$stop_codon_category == "2nd dimerization domain (aa 100-160)", "yes", "no")
MGIT_filtered_not_truncated$second_dimerization <- factor(MGIT_filtered_not_truncated$second_dimerization, levels = c("yes", "no"))
contingency_table <- table(MGIT_filtered_not_truncated$second_dimerization, MGIT_filtered_not_truncated$MIC_interpretation)
fisher.test(contingency_table)

MGIT_filtered_not_truncated$DNA_binding <- ifelse(MGIT_filtered_not_truncated$stop_codon_category == "DNA binding domain (aa 34-99)", "yes", "no")
MGIT_filtered_not_truncated$sDNA_binding <- factor(MGIT_filtered_not_truncated$DNA_binding, levels = c("yes", "no"))
MGIT_filtered_not_truncated$MIC_interpretation <- factor(MGIT_filtered_not_truncated$MIC_interpretation, levels = c("R", "S"))
contingency_table <- table(DNA_binding = factor(MGIT_filtered_not_truncated$DNA_binding, levels = c("yes", "no")), 
                           MIC_interpretation = MGIT_filtered_not_truncated$MIC_interpretation)
fisher.test(contingency_table)

#plot alignment scores, to repeat for each frame, here example for frame 3
ggplot(MGIT_filtered_not_truncated, aes(x = `frame 3_alignment`, y = MGIT, color = position_category)) +
  geom_point(size = 2) +  # Increase point size
  geom_hline(yintercept = 1, color = "red", linetype = "dashed") +
  geom_vline(xintercept = 85, color = "blue", linetype = "dotted") +
  labs(title = "Frame 3 (ARF)",
       x = "Alignment score", y = "BDQ MIC (µg/mL)", color = "Position of frameshift mutation") +
  theme(legend.position = "right",
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 14),  
        axis.text.y = element_text(size = 14),  
        strip.text.x = element_text(size = 16),
        plot.title = element_text(hjust = 0.5, size = 24)) +
  scale_x_continuous(limits = c(0, 100))

table(MGIT_filtered_not_truncated$position_category)

n_isolates <- MGIT_filtered_not_truncated %>%
  filter(highest_alignment < 85, MGIT <= 1) %>%
  summarise(count = n())
print(n_isolates)

MGIT_filtered_not_truncated$MIC_interpretation <- ifelse(MGIT_filtered_not_truncated$MGIT < 1.1, "S", "R")

MGIT_filtered_not_truncated$aligment_conserved <- ifelse(MGIT_filtered_not_truncated$highest_alignment > 85, "yes", "no")

MGIT_filtered_not_truncated$aligment_conserved <- factor(MGIT_filtered_not_truncated$aligment_conserved, levels = c("yes", "no"))
MGIT_filtered_not_truncated$MIC_interpretation <- factor(MGIT_filtered_not_truncated$MIC_interpretation, levels = c("S", "R"))
contingency_table <- table(MGIT_filtered_not_truncated$aligment_conserved, MGIT_filtered_not_truncated$MIC_interpretation)
fisher.test(contingency_table)

#table of unique mutations for the manuscript
MGIT_table <- MGIT_filtered_not_truncated %>%
  group_by(rv0678) %>%
  summarise(
    aa_pos = first(aa_pos),
    position_category = first(position_category),
    stop_codon_category = first(stop_codon_category),
    frame_1_aa_count = first(`frame 1 aa count`),
    frame_1_alignment = first(`frame 1_alignment`),
    frame_2_alignment = first(`frame 2_alignment`),
    frame_3_alignment = first(`frame 3_alignment`),
    frame_4_alignment = first(`frame 4_alignment`),
    frame_5_alignment = first(`frame 5_alignment`),
    frame_6_alignment = first(`frame 6_alignment`),
    highest_alignment = first(highest_alignment),
    number_R_isolates = sum(MIC_interpretation == "R"),
    number_S_isolates = sum(MIC_interpretation == "S")
  ) %>%
  ungroup()

library(openxlsx)
write.xlsx(MGIT_table, "MGIT_table_MS.xlsx")

