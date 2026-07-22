library(ggplot2)
library(gdata)
library(readxl)
library(ggsignif)
library(plyr)
library(dplyr)
library(RColorBrewer)
library(ggrepel)
library(reshape2)
library(tidyr)
library(ggpubr)

#Figure 4A
bm_local = read_excel( "./fig2_data.xlsx", sheet = "Brain vs. Local")
nbm_local = read_excel( "./fig2_data.xlsx", sheet = "Brain vs. NBM")

bm_local_tmb <- bm_local %>%
  filter(`Gene/Biomarker` == 'TMB-High')

nbm_local_tmb <- nbm_local %>%
  filter(`Gene/Biomarker` == 'TMB-High')

bm_nbm_local_tmb <- bm_local_tmb %>%
  bind_cols(select(nbm_local_tmb, -any_of(names(bm_local_tmb)))) %>%
  mutate(`Brain Mets` = `% in Brain`) %>%
  mutate(`Local Biopsies` = `% in Local`) %>%
  mutate(`Non-Brain Mets` = `% in NBM`) 

bm_nbm_local_tmb <- bm_nbm_local_tmb %>%
  mutate(BM_total = `BM|Gene/Biomarker(+)` + `BM|Gene/Biomarker(-)`) %>%
  mutate(Local_total = `Local|Gene/Biomarker(+)` + `Local|Gene/Biomarker(-)`) %>%
  mutate(NBM_total =  `NBM|Gene/Biomarker(+)`+ `NBM|Gene/Biomarker(-)`) %>%
  mutate(BM_TMBH = `BM|Gene/Biomarker(+)`) %>%
  mutate(Local_TMBH = `Local|Gene/Biomarker(+)`) %>%
  mutate(NBM_TMBH = `NBM|Gene/Biomarker(+)`)

bm_nbm_local_tmb_percent <- bm_nbm_local_tmb %>%
  select(Disease, `Brain Mets`, `Local Biopsies`, `Non-Brain Mets`) %>%
  melt()

bm_nbm_local_tmb_percent$Disease <- factor(bm_nbm_local_tmb_percent$Disease, levels = c("Melanoma","NSCLC","Neuroendocrine","Esophagus","Breast","CRC"))

bm_nbm_local_tmb_numbers <- bm_nbm_local_tmb %>%
  select(Disease, BM_TMBH, BM_total, Local_TMBH, Local_total, NBM_TMBH, NBM_total)

bm_nbm_local_tmb_numbers <- bm_nbm_local_tmb_numbers %>%
  pivot_longer(cols = -Disease,names_to = c("variable", ".value"),names_sep = "_") %>%
  mutate(variable = recode(variable, BM = "Brain Mets", Local = "Local Biopsies", NBM = "Non-Brain Mets"), 
         variable = factor(variable, levels = c("Brain Mets", "Local Biopsies", "Non-Brain Mets")))

bm_nbm_local_tmb_numbers$Disease <- factor(bm_nbm_local_tmb_numbers$Disease, levels = c("Melanoma","NSCLC","Neuroendocrine","Esophagus","Breast","CRC"))

bm_nbm_local_tmb_numbers <- bm_nbm_local_tmb_numbers %>%
  left_join(bm_nbm_local_tmb_percent %>% select(Disease, variable, value), by = c("Disease", "variable"))

pdf("Figure4A.pdf", height = 7, width=9)
p <- ggplot(bm_nbm_local_tmb_percent, aes(x = Disease, y = value*100, fill = variable)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(data = bm_nbm_local_tmb_numbers, aes(x= Disease, y=value*100, label = sprintf("frac(%d,%d)", TMBH, total), group = variable),
    position = position_dodge(width = 0.8),
    vjust = -0.2,
    size = 3,
    lineheight = 0.9,
    inherit.aes = FALSE,
    parse = TRUE) +
  scale_fill_manual(values = c("Brain Mets" = "#FF4C00", "Non-Brain Mets" = "darkgreen", "Local Biopsies" = "#64CCC9")) +
  coord_cartesian(ylim = c(0, 75)) +
  labs(y = "Prevalence (%)", x = NULL, fill = NULL) +
  theme_classic() +
  theme(legend.position = "bottom")

print(p)
dev.off()

#Figure 4B
data_4B = read.table("./data_4B.csv", header=TRUE, sep=",")

data_4B <- data_4B %>%
  mutate(tTMB = as.numeric(tTMB))

data_4B$biopsy <- factor(data_4B$biopsy, levels = c('Local', "Brain met", "NBM"))

diseases = c('breast', 'CRC', 'esophagus', 'NSCLC', 'melanoma', 'neuroendocrine')
disease <- 'breast'
for (disease in diseases) {
  
  data_subset <- data_4B %>%
    filter(newDG == disease)
  
  summary <- data_subset %>%
    dplyr::group_by(biopsy) %>%
    dplyr::summarise(n=n())
  
  pdf(file = paste0(disease,"Figure_4B.pdf"), height=4, width=4)
  p <- ggplot(data_subset, aes(y=tTMB, x=biopsy, fill = biopsy)) + 
    geom_boxplot(outlier.alpha = 0) + 
    geom_text(data=summary, aes(x=biopsy, y=-0.8, label=paste0("N=",n)),color="black", size = 3)+
    theme_bw()+
    scale_fill_manual(values = c("#64CCC9", "#ff4c00", "darkgreen")) +
    theme(legend.position="none",
          strip.background = element_blank(),
          axis.title=element_text(size=14, face = "bold",colour = "black"),
          axis.text.x = element_text(size = 10,colour = "black"),
          axis.text.y  = element_text(size = 10,colour = "black"),
          plot.title = element_text(face = "bold", hjust=0.5)) +
    coord_cartesian(ylim=c(0,30))+
    xlab("Biopsy location") +
    ylab("TMB")+
    ggtitle(disease) +
    geom_signif(test = "wilcox.test", comparisons = combn(levels(data_4B$biopsy),2, simplify = F),textsize = 4,tip_length = 0.0001, y_position = c(20,24,28), margin_top = 0, map_signif_level = T, hjust = 0.99)
  
  print(p)
  dev.off()
}

#Figure 4C
data_4C = read.table("./data_4C.csv", header=TRUE, sep=",")

data_4C$value <- as.numeric(data_4C$value)

pdf('Figure_4C.pdf', width = 4, height = 3)
p <- ggplot(data_4C, aes(x = variable, y = value)) +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  geom_line(aes(group = Patient_number), alpha = 0.4, color = "gray50") +     
  geom_point(size = 2, aes(color = newDG_local)) + 
  theme_minimal() +
  labs(x = "Biopsy type", y = "TMB", color = "Tumor type") +
  scale_x_discrete(labels = c("TMB_local" = "Local", "TMB_brainMets"  = "Brain")) +
  stat_compare_means (method = "wilcox.test") + 
  stat_summary(fun = median, geom = "text", aes(label = sprintf("%.2f", ..y..)), vjust = -0.7, color = "black", size = 4) +
  scale_color_manual(values = c("#0072B2", "#E69F00", "#ffd166", "#ef476f","#073b4c", "#009E73"))

print(p)
dev.off()

#Figure 4D
bm_local = read_excel( "./fig2_data.xlsx", sheet = "Brain vs. Local")
nbm_local = read_excel( "./fig2_data.xlsx", sheet = "Brain vs. NBM")

bm_local_msi <- bm_local %>%
  filter(`Gene/Biomarker` == 'MSI-High')

nbm_local_msi <- nbm_local %>%
  filter(`Gene/Biomarker` == 'MSI-High')

bm_nbm_local_msi <- bm_local_msi %>%
  bind_cols(select(nbm_local_msi, -any_of(names(bm_local_msi)))) %>%
  mutate(`Brain Mets` = `% in Brain`) %>%
  mutate(`Local Biopsies` = `% in Local`) %>%
  mutate(`Non-Brain Mets` = `% in NBM`) 

bm_nbm_local_msi <- bm_nbm_local_msi %>%
  mutate(BM_total = `BM|Gene/Biomarker(+)` + `BM|Gene/Biomarker(-)`) %>%
  mutate(Local_total = `Local|Gene/Biomarker(+)` + `Local|Gene/Biomarker(-)`) %>%
  mutate(NBM_total =  `NBM|Gene/Biomarker(+)`+ `NBM|Gene/Biomarker(-)`) %>%
  mutate(BM_msiH = `BM|Gene/Biomarker(+)`) %>%
  mutate(Local_msiH = `Local|Gene/Biomarker(+)`) %>%
  mutate(NBM_msiH = `NBM|Gene/Biomarker(+)`)

bm_nbm_local_msi_percent <- bm_nbm_local_msi %>%
  select(Disease, `Brain Mets`, `Local Biopsies`, `Non-Brain Mets`) %>%
  melt()

bm_nbm_local_msi_percent$Disease <- factor(bm_nbm_local_msi_percent$Disease, levels = c("CRC", "Esophagus", "Breast", "NSCLC", "Melanoma","Neuroendocrine"))

bm_nbm_local_msi_numbers <- bm_nbm_local_msi %>%
  select(Disease, BM_msiH, BM_total, Local_msiH, Local_total, NBM_msiH, NBM_total)

bm_nbm_local_msi_numbers <- bm_nbm_local_msi_numbers %>%
  pivot_longer(cols = -Disease,names_to = c("variable", ".value"),names_sep = "_") %>%
  mutate(variable = recode(variable, BM = "Brain Mets", Local = "Local Biopsies", NBM = "Non-Brain Mets"), 
         variable = factor(variable, levels = c("Brain Mets", "Local Biopsies", "Non-Brain Mets")))

bm_nbm_local_msi_numbers$Disease <- factor(bm_nbm_local_msi_numbers$Disease, levels = c("CRC", "Esophagus", "Breast", "NSCLC", "Melanoma","Neuroendocrine"))

bm_nbm_local_msi_numbers <- bm_nbm_local_msi_numbers %>%
  left_join(bm_nbm_local_msi_percent %>% select(Disease, variable, value), by = c("Disease", "variable"))

pdf("Figure4D.pdf", height = 7, width=9)
p <- ggplot(bm_nbm_local_msi_percent, aes(x = Disease, y = value*100, fill = variable)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(data = bm_nbm_local_msi_numbers, aes(x= Disease, y=value*100, label = sprintf("frac(%d,%d)", msiH, total), group = variable),
            position = position_dodge(width = 0.8),
            vjust = -0.2,
            size = 3,
            lineheight = 0.9,
            inherit.aes = FALSE,
            parse = TRUE) +
  scale_fill_manual(values = c("Brain Mets" = "#FF4C00", "Non-Brain Mets" = "darkgreen", "Local Biopsies" = "#64CCC9")) +
  coord_cartesian(ylim = c(0, 9)) +
  labs(y = "Prevalence (%)", x = NULL, fill = NULL) +
  theme_classic() +
  theme(legend.position = "bottom")

print(p)
dev.off()

            