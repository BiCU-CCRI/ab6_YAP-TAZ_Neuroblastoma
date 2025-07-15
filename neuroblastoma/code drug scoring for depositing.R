library(tidyverse)

setwd("C:/Users/soeren.strohmenger/OneDrive - CCRI/Bioinfo/Drug scoring revision/test")

#calculation of AUC scores SK-N-SH
viabilities <- read.csv("SK-N-SH_viabilities.csv")

viabilities <- viabilities %>% rename("ADR" = NumberOfCells_ADR, "MES" = NumberOfCells_MES)

controls <- viabilities %>% filter(Drug == "DMSO")

ref_ADR = mean(controls$ADR)
ref_MES = mean(controls$MES)

viabilities <- viabilities %>% mutate(ADR_rel = ADR/ref_ADR, MES_rel = MES/ref_MES)

druglist <- unique(viabilities$Drug)
druglist <- druglist[druglist != "DMSO"]


scores <- data.frame(drug = "", score_ADR =NA , score_MES =NA )

for (drug in druglist){
  temp2 <- filter(viabilities, Drug == drug)
  concentrations <- sort(unique(temp2$Concentration))
  
  adr_means <- c()
  mes_means <- c()
  for (i in concentrations){
    #gain mean values for viabilities at concentrations
    adr_means <- append(adr_means, c( mean(filter(temp2, Concentration == i)$ADR_rel)))
    mes_means <- append(mes_means, c( mean(filter(temp2, Concentration == i)$MES_rel)))
    
  }
  auc_adr <- 0
  auc_mes <- 0
  l <- length(concentrations)
  for (i in c(2:(l))){
    auc_adr <- auc_adr + ((1/(l-1)) * 0.5 * (adr_means[i] +adr_means[i-1]))
    auc_mes <- auc_mes + ((1/(l-1)) * 0.5 * (mes_means[i] +mes_means[i-1]))
  }
  #revert to area above curve
  score_adr <- (1-auc_adr)
  score_mes <- (1-auc_mes)
  
  
  #add to scores file
  scores <- add_row(scores, drug = drug, score_ADR = score_adr, score_MES = score_mes)
}

write.csv(scores, file ="SK-N-SH_reanalysis_v3.csv" )

#calculation of AUC scores STA-NB-8
viabilities <- read.csv("STA-NB-8_viabilities.csv")

viabilities <- viabilities %>% rename("ADR" = NumberOfCells_ADR, "MES" = NumberOfCells_MES)

controls <- viabilities %>% filter(Drug == "DMSO")

ref_ADR = mean(controls$ADR)
ref_MES = mean(controls$MES)

viabilities <- viabilities %>% mutate(ADR_rel = ADR/ref_ADR, MES_rel = MES/ref_MES)

druglist <- unique(viabilities$Drug)
druglist <- druglist[druglist != "DMSO"]


scores <- data.frame(drug = "", score_ADR =NA , score_MES =NA )

for (drug in druglist){
  temp2 <- filter(viabilities, Drug == drug)
  concentrations <- sort(unique(temp2$Concentration))
  
  adr_means <- c()
  mes_means <- c()
  for (i in concentrations){
    #gain mean values for viabilities at concentrations
    adr_means <- append(adr_means, c( mean(filter(temp2, Concentration == i)$ADR_rel)))
    mes_means <- append(mes_means, c( mean(filter(temp2, Concentration == i)$MES_rel)))
    
  }
  auc_adr <- 0
  auc_mes <- 0
  l <- length(concentrations)
  for (i in c(2:(l))){
    auc_adr <- auc_adr + ((1/(l-1)) * 0.5 * (adr_means[i] +adr_means[i-1]))
    auc_mes <- auc_mes + ((1/(l-1)) * 0.5 * (mes_means[i] +mes_means[i-1]))
  }
  #revert to area above curve
  score_adr <- (1-auc_adr)
  score_mes <- (1-auc_mes)
  
  
  #add to scores file
  scores <- add_row(scores, drug = drug, score_ADR = score_adr, score_MES = score_mes)
}

write.csv(scores, file ="STA-NB-8_reanalysis_v3.csv" )


#calculation of AUC scores STA-NB-10
viabilities <- read.csv("STA-NB-10_viabilities.csv")

viabilities <- viabilities %>% rename("ADR" = NumberOfCells_ADR, "MES" = NumberOfCells_MES)

controls <- viabilities %>% filter(Drug == "DMSO")

ref_ADR = mean(controls$ADR)
ref_MES = mean(controls$MES)

viabilities <- viabilities %>% mutate(ADR_rel = ADR/ref_ADR, MES_rel = MES/ref_MES)

druglist <- unique(viabilities$Drug)
druglist <- druglist[druglist != "DMSO"]


scores <- data.frame(drug = "", score_ADR =NA , score_MES =NA )

for (drug in druglist){
  temp2 <- filter(viabilities, Drug == drug)
  concentrations <- sort(unique(temp2$Concentration))
  
  adr_means <- c()
  mes_means <- c()
  for (i in concentrations){
    #gain mean values for viabilities at concentrations
    adr_means <- append(adr_means, c( mean(filter(temp2, Concentration == i)$ADR_rel)))
    mes_means <- append(mes_means, c( mean(filter(temp2, Concentration == i)$MES_rel)))
    
  }
  auc_adr <- 0
  auc_mes <- 0
  l <- length(concentrations)
  for (i in c(2:(l))){
    auc_adr <- auc_adr + ((1/(l-1)) * 0.5 * (adr_means[i] +adr_means[i-1]))
    auc_mes <- auc_mes + ((1/(l-1)) * 0.5 * (mes_means[i] +mes_means[i-1]))
  }
  #revert to area above curve
  score_adr <- (1-auc_adr)
  score_mes <- (1-auc_mes)
  
  
  #add to scores file
  scores <- add_row(scores, drug = drug, score_ADR = score_adr, score_MES = score_mes)
}

write.csv(scores, file ="STA-NB-10_reanalysis_v3.csv" )



#calculation of AUC scores CLB-Ma
viabilities <- read.csv("CLB-Ma_viabilities.csv")

viabilities <- viabilities %>% rename("ADR" = NumberOfCells_ADR, "MES" = NumberOfCells_MES)

controls <- viabilities %>% filter(Drug == "DMSO")

ref_ADR = mean(controls$ADR)
ref_MES = mean(controls$MES)

viabilities <- viabilities %>% mutate(ADR_rel = ADR/ref_ADR, MES_rel = MES/ref_MES)

druglist <- unique(viabilities$Drug)
druglist <- druglist[druglist != "DMSO"]


scores <- data.frame(drug = "", score_ADR =NA , score_MES =NA )

for (drug in druglist){
  temp2 <- filter(viabilities, Drug == drug)
  concentrations <- sort(unique(temp2$Concentration))
  
  adr_means <- c()
  mes_means <- c()
  for (i in concentrations){
    #gain mean values for viabilities at concentrations
    adr_means <- append(adr_means, c( mean(filter(temp2, Concentration == i)$ADR_rel)))
    mes_means <- append(mes_means, c( mean(filter(temp2, Concentration == i)$MES_rel)))
    
  }
  auc_adr <- 0
  auc_mes <- 0
  l <- length(concentrations)
  for (i in c(2:(l))){
    auc_adr <- auc_adr + ((1/(l-1)) * 0.5 * (adr_means[i] +adr_means[i-1]))
    auc_mes <- auc_mes + ((1/(l-1)) * 0.5 * (mes_means[i] +mes_means[i-1]))
  }
  #revert to area above curve
  score_adr <- (1-auc_adr)
  score_mes <- (1-auc_mes)
  
  
  #add to scores file
  scores <- add_row(scores, drug = drug, score_ADR = score_adr, score_MES = score_mes)
}

write.csv(scores, file ="CLB-Ma_reanalysis_v3.csv" )



#read all scores lists
SST055 <- read.csv("SK-N-SH_reanalysis_v3.csv")
SST062 <- read.csv("STA-NB-10_reanalysis_v3.csv")
SST064 <- read.csv("STA-NB-8_reanalysis_v3.csv")
SST070 <- read.csv("CLB-Ma_reanalysis_v3.csv")

#List of experiments with names
experiment_names <- c("SST055", "SST062", "SST064", "SST070")
experiments <- list(SST055, SST062, SST064, SST070)
names(experiments) <- experiment_names

#loop through all experiments for positive control correction
for (i in seq_along(experiments)) {
  # Access the dataframe
  df <- experiments[[i]]
  
  # Perform the calculations and update the dataframe
  df <- df %>%
    arrange(drug) %>%
    mutate(adj_score_ADR = if_else(score_ADR < 0, 0, score_ADR)) %>%
    mutate(adj_score_MES = if_else(score_MES < 0, 0, score_MES)) %>%
    mutate(adj_anti_MES = adj_score_MES - adj_score_ADR)
  
  
  # Save the modified dataframe back to the original variable
  assign(experiment_names[i], df)
}

write.csv(SST055, file ="SK-N-SH_scores.csv" )
write.csv(SST062, file ="STA-NB-10_scores.csv" )
write.csv(SST064, file ="STA-NB-8_scores.csv" )
write.csv(SST070, file ="CLB-Ma_scores.csv" )


#remove DMSO and empty row, outlier obatoclax (only results from 2 datapoints), and removal of retinoids: isotretinoin, tretinoin, bexarotene, alitretinoin (affects VIM staining), and Verteporfin (wrong concentration printed)
SST055 <-  SST055 %>% filter(drug != "DMSO") %>% filter(drug != "") %>% filter (drug != "Obatoclax") %>% filter(drug != "Isotretinoin") %>% filter(drug != "Tretinoin") %>% filter(drug != "Alitretinoin") %>% filter (drug != "Bexarotene") %>% filter (drug != "Verteporfin")
SST062 <-  SST062 %>% filter(drug != "DMSO") %>% filter(drug != "") %>% filter (drug != "Obatoclax") %>% filter(drug != "Isotretinoin") %>% filter(drug != "Tretinoin") %>% filter(drug != "Alitretinoin") %>% filter (drug != "Bexarotene") %>% filter (drug != "Verteporfin")
SST064 <-  SST064 %>% filter(drug != "DMSO") %>% filter(drug != "") %>% filter (drug != "Obatoclax") %>% filter(drug != "Isotretinoin") %>% filter(drug != "Tretinoin") %>% filter(drug != "Alitretinoin") %>% filter (drug != "Bexarotene") %>% filter (drug != "Verteporfin")
SST070 <-  SST070 %>% filter(drug != "DMSO") %>% filter(drug != "") %>% filter (drug != "Obatoclax") %>% filter(drug != "Isotretinoin") %>% filter(drug != "Tretinoin") %>% filter(drug != "Alitretinoin") %>% filter (drug != "Bexarotene") %>% filter (drug != "Verteporfin")

#SK-N-SH top 10 top/bottom
top_hits <- top_n(SST055, 10, adj_anti_MES)
bot_hits <- top_n(SST055, -10, adj_anti_MES)
tb_hits <- rbind(top_hits, bot_hits)
ggplot(tb_hits, aes(x=adj_anti_MES, y= reorder(drug, adj_anti_MES), fill=adj_anti_MES)) + geom_bar(stat="summary", fun = "mean") + scale_fill_gradient2(name = "anti-MES", high="#EB363C", low="blue", mid="grey95") +
  theme_minimal()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.y = element_text(size=8, color = "black", face = "bold"), aspect.ratio = 1/1) + 
  labs(x = "anti-MES score", y = "")
ggsave("SK-N-SH_individual_drugs_tb10_v3.pdf")

#STA-NB-10 top hits
top_hits <- SST062 %>% filter(adj_anti_MES > 0.00)
ggplot(top_hits, aes(x=adj_anti_MES, y= reorder(drug, adj_anti_MES), fill=adj_anti_MES)) + geom_bar(stat="summary", fun = "mean") + scale_fill_gradient2(name = "anti-MES", high="#EB363C", low="blue", mid="grey95") +
  theme_minimal()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.y = element_text(size=8, color = "black", face = "bold")) + 
  labs(x = "anti-MES score", y = "")
ggsave("NB-10_individual_drugs_tophits_v3.pdf")

#STA-NB-8 top hits
top_hits <- SST064 %>% filter(adj_anti_MES > 0.19)
ggplot(top_hits, aes(x=adj_anti_MES, y= reorder(drug, adj_anti_MES), fill=adj_anti_MES)) + geom_bar(stat="summary", fun = "mean") + scale_fill_gradient2(name = "anti-MES", high="#EB363C", low="blue", mid="grey95") +
  theme_minimal()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.y = element_text(size=6, color = "black", face = "bold")) + 
  labs(x = "anti-MES score", y = "")
ggsave("NB-8_individual_drugs_tophits_v3.pdf")

#CLB-Ma top hits
top_hits <- SST070 %>% filter(adj_anti_MES > 0.19)
ggplot(top_hits, aes(x=adj_anti_MES, y= reorder(drug, adj_anti_MES), fill=adj_anti_MES)) + geom_bar(stat="summary", fun = "mean") + scale_fill_gradient2(name = "anti-MES", high="#EB363C", low="blue", mid="grey95") +
  theme_minimal()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.y = element_text(size=3.3, color = "black", face = "bold")) + 
  labs(x = "anti-MES score", y = "")
ggsave("CLB-Ma_individual_drugs_tophits_v3.pdf")


# combine anti-MES values in one data frame and save
all_lines_anti_MES <- data.frame(drug = SST055$drug, SK_N_SH = SST055$adj_anti_MES, CLB_Ma = SST070$adj_anti_MES, STA_NB_8 = SST064$adj_anti_MES, STA_NB_10 = SST062$adj_anti_MES)

write.csv(all_lines_anti_MES, "all_lines_anti-MES_v3.csv" )

# combine AUC_ADR values in one data frame and save
all_lines_AUC_ADR <- data.frame(drug = SST055$drug, SK_N_SH = SST055$adj_score_ADR, CLB_Ma = SST070$adj_score_ADR, STA_NB_8 = SST064$adj_score_ADR, STA_NB_10 = SST062$adj_score_ADR)
write.csv(all_lines_AUC_ADR, "all_lines_AUC_ADR_v3.csv" )

# combine AUC_MES values in one data frame and save
all_lines_AUC_MES <- data.frame(drug = SST055$drug, SK_N_SH = SST055$adj_score_MES, CLB_Ma = SST070$adj_score_MES, STA_NB_8 = SST064$adj_score_MES, STA_NB_10 = SST062$adj_score_MES)
write.csv(all_lines_AUC_MES, "all_lines_AUC_MES_v3.csv" )

#combined plot anti-MES all cell lines
all_lines <- read.csv("all_lines_anti-MES_v3.csv")

numeric_cols <- sapply (all_lines, is.numeric)
numeric_cols["X"] <- FALSE  # Exclude column "X" from the condition
to_keep <- apply(all_lines[, numeric_cols], 1, function(row){
  all(row > 0, na.rm = TRUE) || all(row < 0, na.rm = TRUE)
})

only_consistent <- all_lines[to_keep,]

individual_long <- only_consistent %>%
  select(drug, SK_N_SH, CLB_Ma, STA_NB_8, STA_NB_10) %>%
  gather(key = "cell_line", value = "value", -drug)

summary_stats <- individual_long %>%
  group_by(drug) %>%
  summarise(
    avg_anti_MES = mean(value),
    sem = sd(value)/sqrt(4),
    .groups = "drop"
  )

top_hits <- summary_stats %>% top_n(-10, avg_anti_MES)
bot_hits <- summary_stats %>% top_n(10, avg_anti_MES)
tb_summary <- bind_rows(top_hits, bot_hits)

tb_drugs <- tb_summary$drug
individual_long_tb <- individual_long %>%
  filter(drug %in% tb_drugs)

tb_summary$drug <- factor(tb_summary$drug, levels = unique(tb_summary$drug))
individual_long_tb$drug <- factor(individual_long_tb$drug, levels = levels(tb_summary$drug))

tb_summary$drug <- factor(tb_summary$drug, levels = tb_summary %>% arrange(avg_anti_MES) %>% pull(drug))
individual_long_tb$drug <- factor(individual_long_tb$drug, levels = levels(tb_summary$drug))


# Definition custom shapes for geom_jitter
shape_values <- c(
  "SK_N_SH" = 16,    
  "CLB_Ma" = 17,     
  "STA_NB_8" = 15,   
  "STA_NB_10" = 18   
)


ggplot() +
  geom_bar(data = tb_summary, aes(y = drug, x = avg_anti_MES, fill = avg_anti_MES), stat = "identity") +
  geom_errorbarh(data = tb_summary, aes(y = drug, xmin = avg_anti_MES - sem, xmax = avg_anti_MES + sem), height = 0.3, linewidth = 0.2, color = "black") +
  scale_shape_manual(values = shape_values) +
  geom_jitter(data = individual_long_tb, aes(y = drug, x = value, shape = cell_line), width = 0.01, height = 0.01, color = "black", size = 0.75) +
  scale_fill_gradient2(name = "anti_MES", high = "#EB363C", low = "blue", mid = "grey95") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 8, color = "black", face = "bold"),
    aspect.ratio = 1 / 1
  ) +
  labs(y = "Drug", x = "anti-MES score")

ggsave("combined_drugs_tb_10_datapoints_symbols_consistent_error_bars_v3.pdf")
