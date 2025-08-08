
library(tidyverse)
library(DESeq2)
library(GSVA)
library(ggpubr)

# load in data ------------------------------------------------------------

##expr mat
tpm_mat <- read_tsv('samples_mixture.tsv')

##gene signature
EFGR_genes <- read.csv('EGFR_wiki_GO.csv')
gene_sig_list <- split(EFGR_genes$gene_symbol, EFGR_genes$gene_set_name)

hypoxia_genes <- read.csv('hypoxia_genes.csv')

gene_sig_list$Hypoixa <- hypoxia_genes$Hypoxia
names(gene_sig_list) <- c('EGFR_GO', 'EGFR_wiki', 'Hypoxia')

## select the same samples from patients with multiple samples, as used in the CIBERSORT analysis
metadata <- read_tsv('cibersort_samples.tsv')


# compute gene set scores -------------------------------------------------


tpm_mat2 <- as.matrix(tpm_mat[,-1])
rownames(tpm_mat2) <- tpm_mat$gene

tpm_mat_log <- log2(tpm_mat2+1)
ssgsea_results <- gsva(tpm_mat_log, gene_sig_list, method = 'ssgsea')

selected_genes <- c('EGFR', 'AREG', 'TGFA', 'EREG', 'EPGN', 'HBEGF', 'KDR')

expr_selected_genes <- tpm_mat_log[selected_genes,]

plotting_df <- t(rbind(ssgsea_results[,metadata$analysis_id], expr_selected_genes[,metadata$analysis_id])) %>%
  as.data.frame() %>%
  rownames_to_column(var = 'analysis_id') %>%
  left_join(metadata)

save(plotting_df, file = 'plotting_df_250430.Rdata')


# plotting ----------------------------------------------------------------

library(cowplot)
library(ggplot2)
library(pathwork)

plotting_df <- plotting_df %>%
  mutate(hpv_status = if_else(hpv_status == 'Positive', 'p16 Positive', 'p16 Negative'))

plotting_df

hpv_color <- c(`p16 Positive` = '#e58604', `p16 Negative` = '#99c945')

p1 <- ggplot(plotting_df, aes(x = hpv_status, y = Hypoxia, fill = hpv_status)) +
  geom_boxplot(outlier.shape = NA) +  
  scale_fill_manual(values = hpv_color) +
  geom_jitter() + 
  theme_classic2() +
  stat_compare_means() +
  ggtitle('Hypoxia score') +
  theme(axis.title.x = element_blank())+
  labs(y = 'ssgsea score') + 
  theme(legend.position = "none")


p2 <- ggplot(plotting_df, aes(x = hpv_status, y = EGFR_GO, fill = hpv_status)) +
  geom_boxplot(outlier.shape = NA) +  
  geom_jitter() +   
  scale_fill_manual(values = hpv_color) +
  theme_classic2() +
  stat_compare_means() +
  ggtitle('EGFR-GO score') +
  theme(axis.title.x = element_blank())+
  labs(y = 'ssgsea score')+ 
  theme(legend.position = "none")


p3 <- ggplot(plotting_df, aes(x = hpv_status, y = EGFR_wiki, fill = hpv_status)) +
  geom_boxplot(outlier.shape = NA) +  
  scale_fill_manual(values = hpv_color) +
  geom_jitter() + 
  theme_classic2() +
  stat_compare_means() +
  ggtitle('EGFR-wiki score') +
  theme(axis.title.x = element_blank()) +
  labs(y = 'ssgsea score')+ 
  theme(legend.position = "none")



p4 <- ggplot(plotting_df, aes(x = hpv_status, y = EGFR, fill = hpv_status)) +
  geom_boxplot(outlier.shape = NA) +  
  scale_fill_manual(values = hpv_color) +
  geom_jitter() + 
  theme_classic2() +
  stat_compare_means() +
  ggtitle('EGFR') +
  theme(axis.title.x = element_blank()) +
  labs(y = 'log2(TPM+1)')+ 
  theme(legend.position = "none")



p5 <- ggplot(plotting_df, aes(x = hpv_status, y = AREG, fill = hpv_status)) +
  geom_boxplot(outlier.shape = NA) +  
  scale_fill_manual(values = hpv_color) +
  geom_jitter() + 
  theme_classic2() +
  stat_compare_means() +
  ggtitle('AREG') +
  theme(axis.title.x = element_blank())+
  labs(y = 'log2(TPM+1)')+ 
  theme(legend.position = "none")


p6 <- ggplot(plotting_df, aes(x = hpv_status, y = TGFA, fill = hpv_status)) +
  geom_boxplot(outlier.shape = NA) +  
  scale_fill_manual(values = hpv_color) +
  geom_jitter() + 
  theme_classic2() +
  stat_compare_means() +
  ggtitle('TGFA') +
  theme(axis.title.x = element_blank())+
  labs(y = 'log2(TPM+1)')+ 
  theme(legend.position = "none")


p7 <- ggplot(plotting_df, aes(x = hpv_status, y = EREG, fill = hpv_status)) +
  geom_boxplot(outlier.shape = NA) +  
  scale_fill_manual(values = hpv_color) +
  geom_jitter() + 
  theme_classic2() +
  stat_compare_means() +
  ggtitle('EREG') +
  theme(axis.title.x = element_blank())+
  labs(y = 'log2(TPM+1)')+ 
  theme(legend.position = "none")


p8 <- ggplot(plotting_df, aes(x = hpv_status, y = EPGN, fill = hpv_status)) +
  geom_boxplot(outlier.shape = NA) +  
  scale_fill_manual(values = hpv_color) +
  geom_jitter() + 
  theme_classic2() +
  stat_compare_means() +
  ggtitle('EPGN') +
  theme(axis.title.x = element_blank())+
  labs(y = 'log2(TPM+1)')+ 
  theme(legend.position = "none")


p9 <- ggplot(plotting_df, aes(x = hpv_status, y = HBEGF, fill = hpv_status)) +
  geom_boxplot(outlier.shape = NA) +  
  scale_fill_manual(values = hpv_color) +
  geom_jitter() + 
  theme_classic2() +
  stat_compare_means() +
  ggtitle('HBEGF') +
  theme(axis.title.x = element_blank())+
  labs(y = 'log2(TPM+1)')+ 
  theme(legend.position = "none")


plot_grid(p1,NULL,p2,NULL,p3,
          p4,NULL,p5,NULL,p6,
          p7,NULL,p8,NULL,p9,
          rel_widths = c(1, 0.05, 1, 0.05, 1), ncol = 5)

top <- plot_grid(p1,NULL,p2,NULL,p3,
          rel_widths = c(1, 0.05, 1, 0.05, 1), ncol = 5)
mid <- plot_grid(p4,NULL,p5,NULL,p6,
          rel_widths = c(1, 0.05, 1, 0.05, 1), ncol = 5)
bottom <- plot_grid(p7,NULL,p8,NULL,p9,
          rel_widths = c(1, 0.05, 1, 0.05, 1), ncol = 5)

plot_grid(top,
          NULL,
          mid,
          NULL,
          bottom,
          rel_heights = c(1, 0.05, 1, 0.05, 1), ncol = 1)

ggsave('Figure 3.pdf', width = 8, height = 9)












