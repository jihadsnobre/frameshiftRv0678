######WHO catalogue- Impact of the position of frameshift mutations in the Rv0678 gene on BDQ phenotype
library(dplyr)
library(ggplot2)

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

WHO_mutations_36_62 <- WHO_solo %>%
  filter(aa_pos >= 36 & aa_pos <= 62)

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

####analysis at a unique frameshift level 
library(readxl)
library(tidyverse)
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

####22012023 nt info 
library(readxl)
library(dplyr)
library(ggplot2)
library(tidyr)
WHO_frames <- read_excel("WHO_nt_info.xlsx")
WHO_frames <- WHO_frames %>%
  filter(startsWith(variant, "Rv0678"), endsWith(variant, "fs"))

WHO_frames <- WHO_frames %>%
  mutate(Nucleotide = position - 778989)

WHO_frames <- WHO_frames %>%
  mutate(reference_nucleotide = tolower(reference_nucleotide), 
         alternative_nucleotide = tolower(alternative_nucleotide))

writexl::write_xlsx(WHO_frames, "WHO_nt_info_fs.xlsx")
###now run python scripts 

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

# Print the counts
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

#####MIC MGIT database analysis 
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

####now analysis at an isolate level 
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

####now same analysis at a unique mutation level 
library(ggplot2)

MGIT_unique <- MGIT_filtered_not_truncated %>%
  group_by(rv0678) %>%
  summarise(mean_MGIT = mean(MGIT, na.rm = TRUE),
            mean_aa_pos = mean(aa_pos, na.rm = TRUE),
            mean_frame1_aa_count = mean(`frame 1 aa count`, na.rm = TRUE),
            mean_frame2_aa_count = mean(`frame 2 aa count`, na.rm = TRUE),
            mean_frame3_aa_count = mean(`frame 3 aa count`, na.rm = TRUE),
            mean_frame4_aa_count = mean(`frame 4 aa count`, na.rm = TRUE),
            mean_frame5_aa_count = mean(`frame 5 aa count`, na.rm = TRUE),
            mean_frame6_aa_count = mean(`frame 6 aa count`, na.rm = TRUE),
            mean_frame1_alignment = mean(`frame 1_alignment`, na.rm = TRUE),
            mean_frame2_alignment = mean(`frame 2_alignment`, na.rm = TRUE),
            mean_frame3_alignment = mean(`frame 3_alignment`, na.rm = TRUE),
            mean_frame4_alignment = mean(`frame 4_alignment`, na.rm = TRUE),
            mean_frame5_alignment = mean(`frame 5_alignment`, na.rm = TRUE),
            mean_frame6_alignment = mean(`frame 6_alignment`, na.rm = TRUE),
            mean_highest_alignment = mean(highest_alignment, na.rm = TRUE),
            mean_highest_count = mean(highest_count, na.rm = TRUE))

ggplot(MGIT_unique, aes(x = mean_aa_pos)) +
  geom_bar(stat = "count") +
  labs(title = "Position frameshift mutations reported in MGIT dataset in the Rv0678 gene",
       x = "Codon Rv0678",
       y = "Count isolates containing a frameshift mutation")

categories <- c("<16", "16-33", "34-100", "101-160", ">161")

MGIT_unique$stop_codon_category <- cut(
  MGIT_unique$mean_frame1_aa_count,
  breaks = c(0, 16, 33, 100, 160, 166),
  labels = c("begin of the protein (aa 0-16)", "1st dimerization domain (aa 17-33)", "DNA binding domain (aa 34-99)", "2nd dimerization domain (aa 100-160)", "end of the protein (aa 161-165)"),
  include.lowest = TRUE
)

MGIT_unique$position_category <- cut(
  MGIT_unique$mean_aa_pos,
  breaks = c(0, 16, 33, 100, 160, 166),
  labels = c("begin of the protein (aa 0-16)", "1st dimerization domain (aa 17-33)", "DNA binding domain (aa 34-99)", "2nd dimerization domain (aa 100-160)", "end of the protein (aa 161-165)"),
  include.lowest = TRUE
)

ggplot(MGIT_unique, aes(x = stop_codon_category, y = mean_MGIT)) +
  geom_boxplot(fill = "pink") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red", size = 0.5) +
  labs(title = "Boxplot position stop codon vs MIC",
       x = "Position stop codon induced by frameshift mutation",
       y = "Mean BDQ MIC (µg/mL)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

ggplot(MGIT_unique, aes(x = position_category, y = mean_MGIT)) +
  geom_boxplot(fill = "pink") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red", size = 0.5) +
  labs(title = "Boxplot position frameshift mutation vs MIC",
       x = "Position  frameshift mutation",
       y = "MIC (µg/mL)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

table(MGIT_unique$position_category)
table(MGIT_unique$stop_codon_category)

library(openxlsx)
write.xlsx(MGIT_unique, "MGIT_unique.xlsx", sheetName = "Sheet1", rowNames = FALSE)

MGIT_unique$end <- ifelse(MGIT_unique$stop_codon_category == "end of the protein (aa 161-165)", "yes", "no")
MGIT_unique$end <- factor(MGIT_unique$end, levels = c("yes", "no"))
MGIT_unique$MIC_interpretation <- ifelse(MGIT_unique$mean_MGIT < 1.1, "S", "R")
contingency_table <- table(MGIT_unique$end, MGIT_unique$MIC_interpretation)
fisher.test(contingency_table)
wilcox.test(mean_MGIT ~ end, data = MGIT_unique)

MGIT_unique$first_dimerization <- ifelse(MGIT_unique$stop_codon_category == "1st dimerization domain (aa 17-33)", "yes", "no")
MGIT_unique$first_dimerization <- factor(MGIT_unique$first_dimerization, levels = c("yes", "no"))
contingency_table <- table(MGIT_unique$first_dimerization, MGIT_unique$MIC_interpretation)
fisher.test(contingency_table)
wilcox.test(mean_MGIT ~ first_dimerization, data = MGIT_unique)

library(ggplot2)

library(dplyr)
count_columns <- paste("frame", 1:6, "aa count", sep = " ")
alignment_columns <- paste0("frame ", 1:6, "_alignment")

# Create the 'highest_alignment' column
MGIT_filtered_not_truncated <- MGIT_filtered_not_truncated %>%
  mutate(highest_alignment = pmax(!!!select(., all_of(alignment_columns))))

unique_data <- MGIT_filtered_not_truncated %>%
  group_by(rv0678) %>%
  summarise(mean_MGIT = mean(MGIT, na.rm = TRUE),
            mean_aa_pos = mean(aa_pos, na.rm = TRUE),
            mean_stop = mean(highest_count, na.rm = TRUE),
            mean_alignment = mean(highest_alignment, na.rm = TRUE))

library(ggplot2)
ggplot(unique_data, aes(x = mean_alignment, y = mean_MGIT)) +
  geom_point() +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +  
  labs(x = "Alignment Score", y = "MGIT (ug/mL)") +
  ggtitle("Scatterplot MGIT vs highest alignment score ORF and ARFs")

contingency_table <- table(unique_data$mean_MGIT < 1, unique_data$mean_alignment > 85)
fisher_test_result <- fisher.test(contingency_table)
print(fisher_test_result)

frame_labels <- setNames(c("Frame 1", "Frame 2", "Frame 3", "Highest Alignment"),
                         c("mean_frame1_alignment", "mean_frame2_alignment", "mean_frame3_alignment", "mean_highest_alignment"))

ggplot(MGIT_unique_long, aes(x = alignment_score, y = mean_MGIT, color = position_category)) +
  geom_point() +
  geom_hline(yintercept = 1, color = "red", linetype = "dashed") +
  geom_vline(xintercept = 85, color = "blue", linetype = "dotted") +
  labs(x = "Alignment score", y = "mean BDQ MIC (µg/mL)", color = "Position of frameshift mutation") +
  theme(legend.position = "right",
        axis.title.x = element_text(size = 20),  # Increase x-axis label text size
        axis.title.y = element_text(size = 20),  # Increase y-axis label text size
        strip.text.x = element_text(size = 16)) + # Increase facet label text size
  facet_wrap(~ frame, scales = 'free_x', labeller = labeller(frame = frame_labels)) +
  scale_x_continuous(limits = c(0, 100))

colnames(MGIT_unique)

MGIT_unique$max_al <- ifelse(MGIT_unique$mean_frame1_alignment > 85, "Yes", "No")
contingency_table = table(MGIT_unique$max_al, MGIT_unique$MIC_interpretation )
fisher.test(contingency_table)

alignment_long <- MGIT_unique %>%
  pivot_longer(cols = ends_with("alignment"), names_to = "frame", values_to = "alignment_score") %>%
  mutate(frame = sub("_alignment$", "", frame))

ggplot(alignment_long, aes(x = alignment_score)) +
  geom_histogram(bins = 30, fill = "pink", color = "black") +
  geom_vline(xintercept = 85, linetype = "dotted", color = "red", size = 1) +
  theme_minimal() +
  labs(title = "Distribution of Alignment Scores",
       x = "Alignment Score",
       y = "Number of Isolates") +
  facet_wrap(~frame, ncol = 3, labeller = as_labeller(c(mean_frame1 = "Frame 1 (ORF)",
                                                        mean_frame2 = "Frame 2 (ARF forward)",
                                                        mean_frame3 = "Frame 3 (ARF forward)",
                                                        mean_frame4 = "Frame 4 (ARF reverse)",
                                                        mean_frame5 = "Frame 5 (ARF reverse)",
                                                        mean_frame6 = "Frame 6 (ARF reverse)",
                                                        mean_highest = "Highest aligment score"))) +
  theme(strip.text.x = element_text(size = 10, face = "bold"))

#####Cryptic analysis 
library(tidyverse)
library(readxl)
ukmyc_phenotypes <- read.csv("UKMYC_PHENOTYPES.csv")
ukmyc_phenotypes <- ukmyc_phenotypes %>% filter(DRUG == "BDQ")
mutations_cryptic <-  read_excel("rv0678_cryptic_variants.xlsx")
mutations_cryptic <- mutations_cryptic %>% filter(MUTATION_TYPE == "INDEL")
mutations_cryptic <- mutations_cryptic %>% filter(IS_INDEL == "TRUE")
cryptic <- mutations_cryptic %>%
  left_join(ukmyc_phenotypes, by = "UNIQUEID")
cryptic_corr <- cryptic[!is.na(cryptic$MIC), ]
cryptic_corr <- cryptic[!is.na(cryptic$LOG2MIC), ]

###description
ggplot(cryptic_corr, aes(x = AMINO_ACID_NUMBER)) +
  geom_bar(stat = "count") +
  labs(title = "Position frameshift mutations reported in Cryptic dataset in the Rv0678 gene",
       x = "Codon Rv0678",
       y = "Count isolates containing a frameshift mutation")

ggplot(cryptic_corr, aes(x = MIC)) +
  geom_bar(stat = "count") +
  labs(title = "MIC distribution Cryptic dataset",
       x = "MIC BMD",
       y = "Count isolates")

table(cryptic_corr$BINARY_PHENOTYPE)

###make unique table

categories <- c("<16", "16-33", "34-100", "101-160", ">161")

# Create a new column category based on aa_pos
cryptic_corr$category <- cut(cryptic_corr$POSITION,
                             breaks = c(-Inf, 16, 33, 100, 160, Inf),
                             labels = categories,
                             include.lowest = TRUE)

# Create a violin plot
ggplot(cryptic_corr, aes(x = category, y = MIC)) +
  geom_violin(trim = FALSE, fill = "lightblue") +
  geom_boxplot(width = 0.2, fill = "lightblue", color = "black") +
  labs(x = "Frameshift mutation position Rv0678 codons", y = "MGIT (ug/mL)") +
  ggtitle("Violin Plot of MGIT across aa_pos Categories") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 14),   # Adjust x-axis label size
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),  # Adjust x-axis title size
        axis.title.y = element_text(size = 14))  

library(openxlsx)
write.xlsx(cryptic_corr, "cryptic_corr.xlsx")

####frames analysis description 
cryptic_frames <-  read_excel("cryptic_sequence_sixframes_stop_align.xlsx")

####position how many isolates in which position how many >160 

library(readxl)
library(ggplot2)
library(dplyr)

# 1. Count isolates with AMINO_ACID_NUMBER > 160 and plot AMINO_ACID_NUMBER vs MIC
cryptic_frames <- cryptic_frames %>% 
  filter(PHENOTYPE_QUALITY != "LOW")

cryptic_frames$MIC <- as.character(cryptic_frames$MIC)  # Ensure MIC is a character for replacement
cryptic_frames$MIC <- gsub("<=0.015", "0.015", cryptic_frames$MIC)  # Replace '<=0.015' with '0.015'
cryptic_frames$MIC <- gsub("<=0.008", "0.008", cryptic_frames$MIC)  # Replace '<=0.008' with '0.008'
cryptic_frames$MIC <- as.numeric(cryptic_frames$MIC)  # Convert the MIC column to numeric

isolates_with_high_amino_acid_number <- sum(cryptic_frames$AMINO_ACID_NUMBER > 160)
print(paste("Number of isolates with AMINO_ACID_NUMBER > 160:", isolates_with_high_amino_acid_number))

# Plotting AMINO_ACID_NUMBER vs MIC
cryptic_frames$MIC = as.numeric(cryptic_frames$MIC)
ggplot(cryptic_frames, aes(x = AMINO_ACID_NUMBER, y = MIC)) +
  geom_point() +
  geom_hline(yintercept = 0.25, color = "red") +
  labs(title = "position frameshift mutation vs MIC (ug/mL)", x = "position frameshift mutation", y = "MIC (ug/mL)")

# 2. Count isolates with stop codon > 160 in frame 1 and plot correlation with MIC
isolates_with_high_stop_codon_frame1 <- sum(cryptic_frames$`frame 1 aa count` > 160)
print(paste("Number of isolates with stop codon > 160 in frame 1:", isolates_with_high_stop_codon_frame1))

# Plotting stop codon count in frame 1 vs MIC
ggplot(cryptic_frames, aes(x = `frame 1 aa count`, y = MIC)) +
  geom_point() +
  geom_hline(yintercept = 0.25, color = "red") +
  labs(title = "stop codon position in Frame 1 (ORF) vs MIC", x = "Stop codon position in Frame 1 (ORF)", y = "MIC (ug/mL)")

cryptic_frames$max_stop_codon_count <- apply(cryptic_frames[, c('frame 1 aa count', 'frame 2 aa count', 'frame 3 aa count')], 
                                             1, max)


# Alignment score
ggplot(cryptic_frames, aes(x = `frame 1_alignment`, y = MIC)) +
  geom_point() +
  geom_hline(yintercept = 0.25, color = "red") +
  labs(title = "Alignment score ORF vs MIC", 
       x = "Alignment score ORF", 
       y = "MIC (ug/mL)")

cryptic_frames$max_alignment_score <- pmax(cryptic_frames$`frame 1_alignment`, 
                                           cryptic_frames$`frame 2_alignment`, 
                                           cryptic_frames$`frame 3_alignment`, 
                                           cryptic_frames$`frame 4_alignment`, 
                                           cryptic_frames$`frame 5_alignment`, 
                                           cryptic_frames$`frame 6_alignment`, 
                                           na.rm = TRUE)

ggplot(cryptic_frames, aes(x = max_alignment_score, y = MIC)) +
  geom_point() +
  geom_hline(yintercept = 0.25, color = "red") +
  labs(title = "Max alignment score vs MIC", 
       x = "Max alignment score ORF and ARFs", 
       y = "MIC (ug/mL)")


cryptic_unique <- cryptic_frames %>%
  group_by(Nucleotide, REF, ALT, AMINO_ACID_NUMBER) %>%
  summarise(
    mean_MIC = mean(MIC, na.rm = TRUE),
    mean_frame1_aa_count = mean(`frame 1 aa count`, na.rm = TRUE),
    mean_frame2_aa_count = mean(`frame 2 aa count`, na.rm = TRUE),
    mean_frame3_aa_count = mean(`frame 3 aa count`, na.rm = TRUE),
    mean_frame4_aa_count = mean(`frame 4 aa count`, na.rm = TRUE),
    mean_frame5_aa_count = mean(`frame 5 aa count`, na.rm = TRUE),
    mean_frame6_aa_count = mean(`frame 6 aa count`, na.rm = TRUE),
    mean_frame1_alignment = mean(`frame 1_alignment`, na.rm = TRUE),
    mean_frame2_alignment = mean(`frame 2_alignment`, na.rm = TRUE),
    mean_frame3_alignment = mean(`frame 3_alignment`, na.rm = TRUE),
    mean_frame4_alignment = mean(`frame 4_alignment`, na.rm = TRUE),
    mean_frame5_alignment = mean(`frame 5_alignment`, na.rm = TRUE),
    mean_frame6_alignment = mean(`frame 6_alignment`, na.rm = TRUE),
    mean_max_alignment = mean (`max_alignment_score`, na.rm = TRUE))

##Boxplot stop codons and statistical testing 
cryptic_unique$stop_codon_category <- cut(
  cryptic_unique$mean_frame1_aa_count,
  breaks = c(0, 16, 33, 100, 160, 166),
  labels = c("begin of the protein (aa 0-16)", "1st dimerization domain (aa 17-33)", "DNA binding domain (aa 34-99)", "2nd dimerization domain (aa 100-160)", "end of the protein (aa 161-165)"),
  include.lowest = TRUE
)

cryptic_unique$position_category <- cut(
  cryptic_unique$AMINO_ACID_NUMBER,
  breaks = c(0, 16, 33, 100, 160, 166),
  labels = c("begin of the protein (aa 0-16)", "1st dimerization domain (aa 17-33)", "DNA binding domain (aa 34-99)", "2nd dimerization domain (aa 100-160)", "end of the protein (aa 161-165)"),
  include.lowest = TRUE
)

ggplot(cryptic_unique, aes(x = stop_codon_category, y = mean_MIC)) +
  geom_boxplot(fill = "lightblue") +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "red", size = 0.5) +
  labs(title = "Boxplot position stop codon vs MIC",
       x = "Position stop codon induced by frameshift mutation",
       y = "MIC (µg/mL)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

ggplot(cryptic_unique, aes(x = position_category, y = mean_MIC)) +
  geom_boxplot(fill = "lightblue") +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "red", size = 0.5) +
  labs(title = "Boxplot position frameshift mutation vs MIC",
       x = "Position  frameshift mutation",
       y = "MIC (µg/mL)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

cryptic_unique$end <- ifelse(cryptic_unique$stop_codon_category == "end of the protein (aa 160-165)", "Yes", "No")
cryptic_unique$MIC_interpretation <- ifelse(cryptic_unique$mean_MIC < 0.25, "S", "R")
cryptic_unique$MIC_interpretation <- ifelse(cryptic_unique$mean_MIC < 0.25, "S", "R")
contingency_table <- table(cryptic_unique$end, cryptic_unique$MIC_interpretation)
fisher.test(contingency_table)
wilcox.test(mean_MIC ~ end, data = cryptic_unique)

cryptic_unique$DNA_binding <- ifelse(cryptic_unique$stop_codon_category == "DNA binding domain (aa 34-99)", "Yes", "No")
contingency_table <- table(cryptic_unique$DNA_binding, cryptic_unique$MIC_interpretation)
fisher.test(contingency_table)
wilcox.test(mean_MIC ~ first_dimerization, data = cryptic_unique)

ggplot(cryptic_unique, aes(x = mean_frame1_alignment, y = mean_MIC, color = stop_codon_category)) +
  geom_point() +
  geom_hline(yintercept = 0.25, color = "red", linetype = "dashed") +
  labs(title = "Alignment score ORF vs MIC", 
       x = "Alignment score ORF", 
       y = "MIC (ug/mL)",
       color = "Position stop codon") +
  theme(legend.position = "right")

cryptic_unique$max_al <- ifelse(cryptic_unique$mean_max_alignment > 90, "Yes", "No")
contingency_table = table(cryptic_unique$max_al, cryptic_unique$MIC_interpretation )
fisher.test(contingency_table)

###ORF and ARFs aligment score distribution 

alignment_long <- cryptic_unique %>%
  pivot_longer(cols = ends_with("alignment"), names_to = "frame", values_to = "alignment_score") %>%
  mutate(frame = sub("_alignment$", "", frame))

ggplot(alignment_long, aes(x = alignment_score)) +
  geom_histogram(bins = 30, fill = "skyblue", color = "black") +
  geom_vline(xintercept = 90, linetype = "dotted", color = "red", size = 1) +
  theme_minimal() +
  labs(title = "Distribution of Alignment Scores",
       x = "Alignment Score",
       y = "Number of Isolates") +
  facet_wrap(~frame, ncol = 3, labeller = as_labeller(c(mean_frame1 = "Frame 1 (ORF)",
                                                        mean_frame2 = "Frame 2 (ARF forward)",
                                                        mean_frame3 = "Frame 3 (ARF forward)",
                                                        mean_frame4 = "Frame 4 (ARF reverse)",
                                                        mean_frame5 = "Frame 5 (ARF reverse)",
                                                        mean_frame6 = "Frame 6 (ARF reverse)",
                                                        mean_max = "Highest aligment score across reading frames"))) +
  theme(strip.text.x = element_text(size = 10, face = "bold"))
