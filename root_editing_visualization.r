library(ggplot2)
library(dplyr)
library(openxlsx)
library(tidyr)

# Load the data
df <- read.xlsx("CRISPR_Genotyping_Results.xlsx", sheet = "roots_and_plants")

data_edits <- df %>%
  drop_na(Variant_Percentage) %>%
  filter(!Variant_Percentage == "-")

data_edits$Variant_Percentage <- as.numeric(data_edits$Variant_Percentage)
data_edits$Knockout_Score <- as.numeric(data_edits$Knockout_Score)

# Combine data into a single column for 'Value' and add a 'Type' column
data_long <- data_edits %>%
  pivot_longer(cols = c(Variant_Percentage, Knockout_Score), 
               names_to = "Type", 
               values_to = "Value")

# Reorder guides by ascending median Variant % or Knockout-Score
guide_order <- data_long %>%
  group_by(Guide) %>%
  summarize(median_value = median(Value, na.rm = TRUE)) %>%
  arrange(median_value) %>%
  pull(Guide)

data_long$Guide <- factor(data_long$Guide, levels = guide_order)

# Define dodge position for consistent alignment
dodge <- position_dodge(width = 0.8)

# if Conserved_Identifier contains plant, filter it out
data_long <- data_long %>%
  filter(!str_detect(Conserved_Identifier, "plant")) %>%
  filter(!str_detect(Conserved_Identifier, "Plant"))

# set order of x-axis 
guide_order <- c("G12", "G13", "G14", "G48", "G75", "G127", "G132")

data_long <- data_long %>%
  mutate(Guide = factor(Guide, levels = guide_order))

# Create the boxplot with points correctly aligned
data_long %>%
  filter(Type == "Variant_Percentage") %>%
  filter(!Filename == "Amp1R_RR2_G14_PREMIX_F01.ab1") %>% # failed 
  unique() %>%
  ggplot(., aes(x = Guide, y = Value)) +
  geom_boxplot(position = dodge, alpha = 0.7) + 
  geom_jitter(size = 2, width = 0.2, alpha = 0.7) + 
  #geom_point(position = dodge, size = 1, alpha = 0.9) + # Points aligned with dodge
  labs(
    title = "Variant Percentages by Guide in Marker Positive Roots",
    x = "Guide",
    y = "Percentage",
    fill = "Type",
    color = "Type"
  ) +
  theme_classic()


# save variant percentage graph
ggsave("Variant_Percentage_by_Guide.pdf", dpi=300, width = 5, height = 5)

# provide summary statistics for Variant_Percentage 
data_edits %>%
  group_by(Guide) %>%
  filter(!str_detect(Conserved_Identifier, "plant")) %>%
  filter(!str_detect(Conserved_Identifier, "Plant")) %>%
  summarize(
    mean_variant_percentage = mean(Variant_Percentage, na.rm = TRUE),
    median_variant_percentage = median(Variant_Percentage, na.rm = TRUE),
    sd_variant_percentage = sd(Variant_Percentage, na.rm = TRUE),
    n = n()
  )
