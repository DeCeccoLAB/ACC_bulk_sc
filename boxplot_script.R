library(ggplot2)
library(tidyverse)
library(dplyr)
library(ggpubr)
library(rstatix)

data <- read.table('../ACC292_genes.txt',sep = '\t', header = T)
genes <- c("VTCN1","MYB","NECTIN4","TACSTD2","STING1","KIT","FAP","RBPJ","CDK4","CDK6",
            "EGFR","MET")

### Tumor vs Normal
long_df <- data %>%
  pivot_longer(cols = all_of(genes),
               names_to = "Gene",
               values_to = "Expression")

long_df$Covariate.1 <- ifelse(long_df$Covariate.1 == 'tum','Tumor','Normal')

long_df$Expression <- as.numeric(long_df$Expression)
long_df$Gene <- factor(long_df$Gene, levels = as.vector(genes))

my_comp <- list(c("Tumor","Normal"))

ggplot(long_df,aes(x = Covariate.1, y = Expression)) +
  geom_boxplot(aes(color = Covariate.1), outlier.shape = NA, width = 0.8, lwd = 0.8, notch = T, notchwidth = 0.3) +
  ggpubr::stat_compare_means(comparisons = my_comp,size = 5.5, label.x.npc = 0.3, label.y.npc = 0.9,
                             paired = F, method = 'wilcox.test',aes(label = custom_p)) +
  scale_color_manual(values = c("Tumor"  = "#EF3B2C","Normal" = "#009E73")) +
  facet_wrap(vars(Gene), scales = "free", ncol = 6) +
  theme_bw() + ylab('log2(expression)') + ggtitle('') + ylim(0,13) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size=14),
        axis.text.x = element_text(angle = 55, vjust = 1, hjust=1,size = 14),
        axis.text.y = element_text(size=14),
        legend.title = element_blank(),
        legend.text = element_text(size=14),
        strip.text = element_text(size = 14),
  )

### Normal, ACCI, ACCII

long_df$ClustersB <- gsub('norm','Normal',long_df$ClustersB)

my_comp <- list(c('Normal','ACC1'),
                c('Normal','ACC2'))

ggplot(long_df,aes(x = ClustersB, y = Expression)) +
  geom_boxplot(aes(color = ClustersB), outlier.shape = NA, width = 0.8, lwd = 0.8, notch = T, notchwidth = 0.3) +
  ggpubr::stat_compare_means(comparisons = my_comp,size = 5.5, label.x.npc = 0.3, label.y.npc = 0.9,
                             paired = F, method = 'wilcox.test',aes(label = custom_p)) +
  scale_color_manual(values = c("ACC1"  = "royalblue",ACC2 = 'darkorange1',"Normal" = "#009E73")) +
  facet_wrap(vars(Gene), scales = "free", ncol = 6) +
  theme_bw() + ylab('log2(expression)') + ggtitle('') + ylim(0,13) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size=14),
        axis.text.x = element_text(angle = 55, vjust = 1, hjust=1,size = 14),
        axis.text.y = element_text(size=14),
        legend.title = element_blank(),
        legend.text = element_text(size=14),
        strip.text = element_text(size = 14),
  )
