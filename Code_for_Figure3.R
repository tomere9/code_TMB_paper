#####Required input format:
# Sheet name: Brain vs. Local
# Disease	Gene/Biomarker	% in Brain	l95 in Brain	u95 in Brain	% in Local	l95 in Local	u95 in Local	BM|Gene/Biomarker(+)	BM|Gene/Biomarker(-)	Local|Gene/Biomarker(+)	Local|Gene/Biomarker(-)	OR	pval (raw)	pval (fdr)
# NSCLC	TMB-High	53.36%	51.49%	55.21%	32.46%	31.94%	32.98%	1502	1313	10113	21044	2.380419463	3.2303E-105	7.1066E-103
# NSCLC	TP53	77.73%	76.14%	79.25%	68.05%	67.53%	68.57%	2188	627	21204	9956	1.638501598	1.2071E-27	5.31123E-26
# NSCLC	KRAS	37.19%	35.40%	39.01%	30.53%	30.02%	31.04%	1047	1768	9512	21648	1.34775316	5.30216E-13	5.83237E-12
# NSCLC	CDKN2A	32.65%	30.92%	34.41%	28.35%	27.85%	28.85%	919	1896	8833	22327	1.225178368	1.87041E-06	1.05511E-05

#Sheetname: Brain vs. NBM
# Disease	Gene/Biomarker	% in Brain	l95 in Brain	u95 in Brain	% in NBM	l95 in NBM	u95 in NBM	BM|Gene/Biomarker(+)	BM|Gene/Biomarker(-)	NBM|Gene/Biomarker(+)	NBM|Gene/Biomarker(-)	OR	pval (raw)	pval (fdr)
# NSCLC	TMB-High	53.36%	51.49%	55.21%	36.24%	35.29%	37.20%	1502	1313	3548	6242	2.012543887	6.06895E-59	1.33517E-56
# NSCLC	MSI-High	0.98%	0.65%	1.42%	0.63%	0.48%	0.81%	27	2728	60	9506	1.568071848	0.069322075	0.206092656
# NSCLC	PD-L1(+)	57.26%	54.73%	59.76%	61.46%	60.18%	62.73%	872	651	3472	2177	0.839874139	0.003129228	0.016790978
# NSCLC	PD-L1 High(+)	30.07%	27.78%	32.44%	31.65%	30.44%	32.88%	458	1065	1788	3861	0.928641649	0.249435729	0.496848383
  
#Sheetname: ERK-MAPK BM vs. Local
# Disease	Gene/Biomarker	% in Brain	l95 in Brain	u95 in Brain	% in Local	l95 in Local	u95 in Local	BM|Gene/Biomarker(+)	BM|Gene/Biomarker(-)	Local|Gene/Biomarker(+)	Local|Gene/Biomarker(-)	OR	pval (raw)	pval (fdr)
# NSCLC	ERK-MAPK path(+)	67.32%	65.55%	69.05%	58.45%	57.90%	59.00%	1895	920	18213	12947	1.464229146	1.87907E-20	1.69116E-19
# NSCLC	BRAF	5.54%	4.73%	6.45%	4.54%	4.31%	4.77%	156	2659	1414	29746	1.234199668	0.016787085	0.021583395
# NSCLC	EGFR	13.64%	12.39%	14.96%	14.25%	13.87%	14.65%	384	2431	4441	26719	0.950354624	0.382273632	0.382273632
# NSCLC	ERBB2	4.58%	3.84%	5.42%	3.48%	3.27%	3.68%	129	2686	1083	30077	1.333797076	0.003459326	0.006226788

#Sheetname: ERK-MAPK BM vs. NBM
# Disease	Gene/Biomarker	% in Brain	l95 in Brain	u95 in Brain	% in NBM	l95 in NBM	u95 in NBM	BM|Gene/Biomarker(+)	BM|Gene/Biomarker(-)	NBM|Gene/Biomarker(+)	NBM|Gene/Biomarker(-)	OR	pval (raw)	pval (fdr)
# NSCLC	ERK-MAPK path(+)	67.32%	65.55%	69.05%	60.40%	59.42%	61.37%	1895	920	5915	3878	1.350437355	2.07699E-11	9.34646E-11
# NSCLC	KRAS	37.19%	35.40%	39.01%	30.38%	29.47%	31.30%	1047	1768	2975	6818	1.357170615	1.25554E-11	9.34646E-11
# NSCLC	NF1	10.12%	9.03%	11.30%	7.54%	7.02%	8.08%	285	2530	738	9055	1.38215399	1.59122E-05	4.77366E-05
# NSCLC	EGFR	13.64%	12.39%	14.96%	15.48%	14.77%	16.21%	384	2431	1516	8277	0.862422383	0.016767468	0.037726803

library(ggplot2)
library(gdata)
library(readxl)
library(ggsignif)
library(plyr)
library(dplyr)
library(RColorBrewer)
library(ggrepel)

fig_data_file_path <- "fig3_data.xlsx"

#assigning colors
Revised_color = "darkgreen"
FMI_Fire = "#FF4C00"
FMI_Sea = "#64CCC9"

color_map = c('Brain Mets'=FMI_Fire, 'Local Biopsies'=FMI_Sea, "Non-Brain Mets"=Revised_color)

diseases = c('Breast', 'CRC', 'Esophagus', 'NSCLC', 'Melanoma', 'Neuroendocrine')

for (comparator in c('Local', 'NBM')) {
  for (sheetype in c('ERK-MAPK BM vs. ')) {
    sheetname = paste0(sheetype, comparator)
    for (disease in diseases) {
      threshold = 0
      
      if (sheetype == 'ERK-MAPK BM vs. ') {
        title1 = paste('ERK-MAPK Pathway(+) Prevalence', 'in', disease)
        outname = paste0(disease, '.', 'Brain_vs_', comparator, '.', 'ERK-MAPK', '.paired_genes.bar_plot.pdf')
      } else {
        title1 = paste(disease , 'Comparative Prevalence')
        outname = paste0(disease, '.', 'Brain_vs_', comparator, '.', 'all_genes_biomarkers', '.paired_genes.bar_plot.pdf')
        threshold = 5
      }
      outname = gsub('\\*', '', outname)    
      #Gene level data
      if (comparator == 'Local') {
        df=read_excel( fig_data_file_path, sheet = sheetname)
        title2 = "Brain Metastases vs. Local Biopsies"
      } else if (comparator == 'NBM') {
        df=read_excel( fig_data_file_path, sheet = sheetname)
        title2 = "Brain Metastases vs. Non-Brain Metastases"
      } else {
        df = NA
      }

      df <- subset(df, Disease == disease)
      df$sig_status[T]='not significant'
      df$sig_status[df$`pval (fdr)`<0.05]='significant (FDR <= 0.05)'
      df$sig_status[df$`pval (fdr)`>=0.05]='not significant'
      
      if (comparator == 'Local') {
        df$`Gene/Biomarker Prevalence (%)`=(df$`BM|Gene/Biomarker(+)`+df$`Local|Gene/Biomarker(+)`)*100/(df$`BM|Gene/Biomarker(+)`+df$`Local|Gene/Biomarker(+)`+df$`BM|Gene/Biomarker(-)`+df$`Local|Gene/Biomarker(-)`)
      } else if (comparator == 'NBM') {
        df$`Gene/Biomarker Prevalence (%)`=(df$`BM|Gene/Biomarker(+)`+df$`NBM|Gene/Biomarker(+)`)*100/(df$`BM|Gene/Biomarker(+)`+df$`NBM|Gene/Biomarker(+)`+df$`BM|Gene/Biomarker(-)`+df$`NBM|Gene/Biomarker(-)`)
      } else {
        df$`Gene/Biomarker Prevalence (%)`= NA
      }
      
      df = df %>% mutate(orderer=ifelse(grepl('path', `Gene/Biomarker`) & sheetype != 'Brain vs. ', 100 + `Gene/Biomarker Prevalence (%)`, `Gene/Biomarker Prevalence (%)`))
      df <- subset(df, `Gene/Biomarker Prevalence (%)` >= threshold)
      df_in_brain = subset(df, T)
      df_in_brain = df_in_brain %>% mutate(plot_attr="Brain Mets")
      df_in_brain = df_in_brain %>% mutate(plot_prev=`% in Brain`)
      df_in_brain = df_in_brain %>% mutate(l95=`l95 in Brain`)
      df_in_brain = df_in_brain %>% mutate(u95=`u95 in Brain`)
      df_in_brain = subset(df_in_brain, select=c('Gene/Biomarker', 'plot_attr', 'plot_prev', 'Gene/Biomarker Prevalence (%)', 'sig_status', 'orderer', 'l95', 'u95'))
      if (comparator == 'Local') {
        df_in_local = subset(df, T)
        df_in_local = df_in_local %>% mutate(plot_attr="Local Biopsies")
        df_in_local = df_in_local %>% mutate(plot_prev=`% in Local`)
        df_in_local = df_in_local %>% mutate(l95=`l95 in Local`)
        df_in_local = df_in_local %>% mutate(u95=`u95 in Local`)
        df_in_local = subset(df_in_local, select=c('Gene/Biomarker', 'plot_attr', 'plot_prev', 'Gene/Biomarker Prevalence (%)', 'sig_status', 'orderer', 'l95', 'u95'))
        df_long <- rbind(df_in_brain, df_in_local)
      } else if (comparator == 'NBM') {
        df_in_nbs = subset(df, T)
        df_in_nbs = df_in_nbs %>% mutate(plot_attr="Non-Brain Mets")
        df_in_nbs = df_in_nbs %>% mutate(plot_prev=`% in NBM`)
        df_in_nbs = df_in_nbs %>% mutate(l95=`l95 in NBM`)
        df_in_nbs = df_in_nbs %>% mutate(u95=`u95 in NBM`)
        df_in_nbs = subset(df_in_nbs, select=c('Gene/Biomarker', 'plot_attr', 'plot_prev', 'Gene/Biomarker Prevalence (%)', 'sig_status', 'orderer', 'l95', 'u95'))
        df_long <- rbind(df_in_brain, df_in_nbs)
      }
      
      df_long$sig_status <- factor(df_long$sig_status, levels = c('significant (FDR <= 0.05)', 'not significant'))
      title = paste0(title1, ': ', title2)
      ggplot(df_long, aes(x = reorder(`Gene/Biomarker`, -orderer), y = plot_prev, fill = plot_attr, alpha=sig_status, color=sig_status)) +
        geom_bar(stat = 'identity', position = 'dodge') +
        geom_errorbar(aes(x=`Gene/Biomarker`, ymin=`l95`, ymax=`u95`), width=0.45,position = position_dodge(0.9)) +
        scale_y_continuous(labels = scales::percent_format())+
        scale_alpha_discrete(name="Statistical Comparison", range = c(1, 0.5), drop = FALSE)+
        scale_fill_manual(name='Sample Type', values=color_map) +
        scale_color_manual(name="Statistical Comparison", values=c('significant (FDR <= 0.05)'='black', 'not significant'=NA), drop=FALSE) + 
        theme_bw() + theme(
          plot.title = element_text(face='bold', hjust=0.5),
          axis.text.x = element_text(size = 12, colour = "black", face = "italic"),
          legend.position="bottom"
        ) +
        xlab('') +
        ylab('Prevalence (%)') +
        labs(fill='Sample Status')
      ggsave(outname, device='pdf', width=10.5, height=5.5, units = "in")
    }
  }
}
