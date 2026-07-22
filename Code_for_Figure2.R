#####Required input format:
# Sheet name: Brain vs. Local
# Disease	Gene/Biomarker	% in Brain	l95 in Brain	u95 in Brain	% in Local	l95 in Local	u95 in Local	BM|Gene/Biomarker(+)	BM|Gene/Biomarker(-)	Local|Gene/Biomarker(+)	Local|Gene/Biomarker(-)	OR	pval (raw)	pval (fdr)
# NSCLC	TMB-High	53.36%	51.49%	55.21%	32.46%	31.94%	32.98%	1502	1313	10113	21044	2.380419463	3.2303E-105	7.1066E-103
# NSCLC	TP53	77.73%	76.14%	79.25%	68.05%	67.53%	68.57%	2188	627	21204	9956	1.638501598	1.2071E-27	5.31123E-26
# NSCLC	KRAS	37.19%	35.40%	39.01%	30.53%	30.02%	31.04%	1047	1768	9512	21648	1.34775316	5.30216E-13	5.83237E-12
# NSCLC	CDKN2A	32.65%	30.92%	34.41%	28.35%	27.85%	28.85%	919	1896	8833	22327	1.225178368	1.87041E-06	1.05511E-05

# Sheet name: Brain vs. NBM
# Disease	Gene/Biomarker	% in Brain	l95 in Brain	u95 in Brain	% in NBM	l95 in NBM	u95 in NBM	BM|Gene/Biomarker(+)	BM|Gene/Biomarker(-)	NBM|Gene/Biomarker(+)	NBM|Gene/Biomarker(-)	OR	pval (raw)	pval (fdr)
# NSCLC	TMB-High	53.36%	51.49%	55.21%	36.24%	35.29%	37.20%	1502	1313	3548	6242	2.012543887	6.06895E-59	1.33517E-56
# NSCLC	MSI-High	0.98%	0.65%	1.42%	0.63%	0.48%	0.81%	27	2728	60	9506	1.568071848	0.069322075	0.206092656
# NSCLC	PD-L1(+)	57.26%	54.73%	59.76%	61.46%	60.18%	62.73%	872	651	3472	2177	0.839874139	0.003129228	0.016790978
# NSCLC	PD-L1 High(+)	30.07%	27.78%	32.44%	31.65%	30.44%	32.88%	458	1065	1788	3861	0.928641649	0.249435729	0.496848383

library(ggplot2)
library(gdata)
library(readxl)
library(ggsignif)
library(plyr)
library(dplyr)
library(RColorBrewer)
library(ggrepel)

fig_data_file_path <- "fig2_data.xlsx"

#assigning colors
Revised_color = "darkgreen"
FMI_Fire = "#FF4C00"
FMI_Sea = "#64CCC9"
fmi_slate_dark = "black"

diseases = c('Breast', 'CRC', 'Esophagus', 'NSCLC', 'Melanoma', 'Neuroendocrine')
for (comparator in c('Local', 'NBM')) {
  for (disease in diseases) {
    title1 = disease
    #Gene level data
    if (comparator == 'Local') {
      volcano_data_gene=read_excel( fig_data_file_path, sheet = "Brain vs. Local")
      title2 = "Brain Metastases vs. Local Biopsies"
    } else if (comparator == 'NBM') {
      volcano_data_gene=read_excel( fig_data_file_path, sheet = "Brain vs. NBM")
      title2 = "Brain Metastases vs. Non-Brain Metastases"
    } else {
      volcano_data_gene = NA
    }
    
    volcano_data_gene=subset(volcano_data_gene, !is.na(OR))
    volcano_data_gene=subset(volcano_data_gene, Disease == disease)
    volcano_data_gene$OR=as.numeric(volcano_data_gene$OR)
    volcano_data_gene$OR_log=log2(volcano_data_gene$OR)
    volcano_data_gene$pval_log=-log10(volcano_data_gene$`pval (fdr)`)
    
    volcano_data_gene$Category <- ""
    volcano_data_gene$Category[volcano_data_gene$OR_log>0 & volcano_data_gene$`pval (fdr)`<0.05]='Enriched in Brain'
    volcano_data_gene$Category[volcano_data_gene$OR_log<0 & volcano_data_gene$`pval (fdr)`<0.05]='Depleted in Brain'
    volcano_data_gene$Category[volcano_data_gene$`pval (fdr)`>=0.05]='Not significant'
    
    
    if (comparator == 'Local') {
      volcano_data_gene$`Gene/Biomarker Prevalence (%)`=(volcano_data_gene$`BM|Gene/Biomarker(+)`+volcano_data_gene$`Local|Gene/Biomarker(+)`)*100/(volcano_data_gene$`BM|Gene/Biomarker(+)`+volcano_data_gene$`Local|Gene/Biomarker(+)`+volcano_data_gene$`BM|Gene/Biomarker(-)`+volcano_data_gene$`Local|Gene/Biomarker(-)`)
      color_map = c('Enriched in Brain'=FMI_Fire, 'Depleted in Brain'=FMI_Sea, 'Not significant'='khaki')
    } else if (comparator == 'NBM') {
      volcano_data_gene$`Gene/Biomarker Prevalence (%)`=(volcano_data_gene$`BM|Gene/Biomarker(+)`+volcano_data_gene$`NBM|Gene/Biomarker(+)`)*100/(volcano_data_gene$`BM|Gene/Biomarker(+)`+volcano_data_gene$`NBM|Gene/Biomarker(+)`+volcano_data_gene$`BM|Gene/Biomarker(-)`+volcano_data_gene$`NBM|Gene/Biomarker(-)`)
      color_map = c('Enriched in Brain'=FMI_Fire, 'Depleted in Brain'=Revised_color,  'Not significant'='khaki')
    } else {
      volcano_data_gene$`Gene/Biomarker Prevalence (%)`= NA
    }
    
    outname = paste0(disease, '.Brain_vs_', comparator, '.volcano.pdf')
    outname = gsub('\\*', '', outname)
    pdf(file = outname, height=6, width=8.5)
    p=ggplot(volcano_data_gene,aes(x=OR_log,y=pval_log))+
      geom_point(aes(size=`Gene/Biomarker Prevalence (%)`,color=Category),alpha=0.7)+
      scale_size_continuous(limits = c(0, 100)) +
      geom_text_repel(data=subset(volcano_data_gene, pval_log>-log10(0.05)),aes(label=`Gene/Biomarker`, color=Category),vjust=-0.1,size=6, 
                      nudge_y = 0.3,
                      nudge_x = 0.3, 
                      show.legend = FALSE,
                      fontface='italic',
                      max.overlaps = 6)+
      theme_classic()+
      ggtitle(title1, subtitle=title2) +
      geom_hline(yintercept = -log10(0.05),color=fmi_slate_dark, linetype="dashed")+
      geom_vline(xintercept = 0)+
      xlab(bquote('log' [2] ~ 'OR'))+
      ylab(bquote('-log' [10] ~ 'P-value'))+
      scale_color_manual(name='Association Category', values=color_map)+
      theme(axis.title = element_text(face="bold", size=20),
            legend.title = element_text(size = 14),
            legend.text = element_text(size = 12),
            axis.text = element_text(size=15),
            plot.title = element_text(hjust = 0.5, size=20, face="bold"),
            plot.subtitle = element_text(hjust = 0.5, size=15, face="bold"))+
      labs(size='Prevalence (%)') +
      coord_cartesian(clip='off')
    print(p)
    dev.off()
  }
}
