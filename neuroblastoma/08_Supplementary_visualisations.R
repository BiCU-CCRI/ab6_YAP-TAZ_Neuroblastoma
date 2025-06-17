library(stringr)
library(dplyr)
library(tidyr)
library(ggplot2)

url <- "https://plus.figshare.com/ndownloader/files/51065297"
file_path <- "~/workspace/neuroblastoma/data/temp_data/Model.csv"

download.file(url = url, destfile = file_path)
Model_dataframe <- read.csv(file = file_path, header = TRUE, sep = ",")
Model_dataframe <- Model_dataframe %>% dplyr::select(ModelID, StrippedCellLineName, OncotreeSubtype)
Model_dataframe$MetaGroup <- ifelse(Model_dataframe$OncotreeSubtype == "Neuroblastoma", "Neuroblastoma", "Other")
Model_dataframe$MetaGroup <- ifelse(Model_dataframe$StrippedCellLineName %in% c("GIMEN", "NB69", "KPNSI9S"), 
                                    Model_dataframe$StrippedCellLineName,
                                    Model_dataframe$MetaGroup)


url <- "https://plus.figshare.com/ndownloader/files/51064667"
file_path <- "~/workspace/neuroblastoma/data/temp_data/ScreenGeneEffect.csv"
download.file(url = url, destfile = file_path)
Screen_Gene_Effect <- read.csv(file = file_path, header = TRUE, sep = ",")
colnames(Screen_Gene_Effect) <- str_replace_all(colnames(Screen_Gene_Effect), 
                                                          pattern = "\\.\\..*$", replacement = "")
Screen_Gene_Effect_compact <- Screen_Gene_Effect[, c("X", "WWTR1", "YAP1", "TEAD1", "TEAD2", "TEAD3", "TEAD4", "JUN", "FOS", "FOSL1", "FOSL2")]
colnames(Screen_Gene_Effect_compact)[1] <- "ModelID"
expression_long <- Screen_Gene_Effect_compact %>%
  pivot_longer(cols = -ModelID, names_to = "Gene", values_to = "GeneEffect") %>%
  left_join(Model_dataframe, by = "ModelID") %>%
  mutate(order_flag = ifelse(MetaGroup != "Other", 1, 0)) %>%
  arrange(order_flag) %>% 
  mutate(order_flag = ifelse(MetaGroup %in% c("GIMEN", "NB69", "KPNSI9S"), 1, 0)) %>%
  arrange(order_flag)

ggplot(expression_long, aes(x = Gene, y = GeneEffect, color = MetaGroup)) +
  geom_point() +
  scale_color_manual(
    values = c(
      "GIMEN" = "firebrick",
      "NB69" = "green",
      "KPNSI9S" = "yellow",
      "Neuroblastoma" = "navy",
      "Other" = "lightgrey"
    )
  ) +
  scale_alpha_manual(values = c(
      "GIMEN" = 1,
      "NB69" = 1,
      "KPNSI9S" = 1,
      "Neuroblastoma" = 1,
      "Other" = 0.2
  )) +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Gene Effect per Sample",
       x = "Sample (X)", y = "Gene Effect")


table((expression_long %>% filter(Gene == "FOS"))$MetaGroup)
