library("ggplot2")
library("ggrepel")
library("pheatmap")

exprs <- read.csv("D:/McGill/projects/project BRAFi/lighthouse/TCGA/TCGA-SKCM.htseq_fpkm gene name.csv")
BRAF <- read.csv("D:/McGill/projects/project BRAFi/lighthouse/TCGA/BRAF_mutation.csv")
BRAFp <- BRAF$sample[which(BRAF$gene != "BRAF")]
BRAFp <- gsub("-", ".", BRAFp)

unique(exprs$X)
nrow(exprs)
df <- as.data.frame(t(exprs[exprs$X %in% c("DDIT3", "ASNS", "TRIB3", "ATF3", "VEGFA", "MTHFD2", "SLC7A11", "AARS", "CEBPB", "CHAC1",
                     "DDIT4", "GPT2", "MAP1LC3B", "PPP1R15A", "PSAT1", "WARS", "ALDH18A1", "ATG7", "EIF2S2",
                     "EPRS", "FGF19", "GARS", "GDF15", "HERPUD1", "HSPA5", "IARS", "JDP2", "KDM7A", "LARS", "MKNK2", "NARS",
                     "PTGS2", "SARS", "SQSTM1", "VARS", "VLDLR", "YARS", "ATF4"), ]))
colnames(df) <- df[1, ]
df <- df[-(1 : 2),]
df <- df[BRAFp,] # exclude BRAF
ATF4 <- as.numeric(df$ATF4)
df <- subset(df, select = -ATF4)
cor_data_df <- data.frame(colnames(df))
for(i in 1 : ncol(df)){
  test <- cor.test(as.numeric(df[, i]), ATF4, type = "spearman")
  cor_data_df[i, 2] <- test$estimate
  cor_data_df[i, 3] <- test$p.value
}
colnames(cor_data_df) <- c("symbol", "correlation", "pvalue")
cor_data_df$pvalue <- cor_data_df$pvalue * nrow(cor_data_df)

cor_data_df[, 4] <- ifelse(cor_data_df$pvalue < 0.05 & cor_data_df$correlation > 0.33, "Y", "N")
cor_data_df[, 5] <- rank(cor_data_df[, 2])

colnames(cor_data_df) <- c("symbol", "correlation", "pvalue", "sig", "rank")


cor_data_df$gene <- rep(NA, nrow(cor_data_df))
cor_data_df$gene[cor_data_df$sig == "Y"] <- cor_data_df$symbol[cor_data_df$sig == "Y"]

ggplot(data = cor_data_df, aes(x = rank, y = correlation, colour = sig, size = -log10(pvalue))) +
  xlab("Gene rank position") +
  ylab("Pearson correlation coefficient") +
  geom_point() +
  ggtitle("TCGA (SKCM)") +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1, colour = "grey") +
#  geom_vline(xintercept = 0.3, linetype = 4, size = 1) +
  theme_classic(base_size = 15) +
  scale_color_manual(values = c("grey", "#E41A1C")) +
  geom_text_repel(aes(label = gene), size = 4, color = 'black')
  theme(# legend.title = element_blank(),
    legend.position = "none",
    legend.background = element_rect(fill = 'transparent'),
    axis.text = element_text(colour = "black"),
    panel.border = element_rect(size = 1, fill = 'transparent'))






# heatmap

df <- as.data.frame(exprs[exprs$X %in% c("ATF4", "DDIT3", "ASNS", "TRIB3", "ATF3", "VEGFA", "MTHFD2", "SLC7A11", "AARS", "CEBPB", "CHAC1",
                                           "DDIT4", "GPT2", "MAP1LC3B", "PPP1R15A", "PSAT1", "WARS", "ALDH18A1", "ATG7", "EIF2S2",
                                           "EPRS", "FGF19", "GARS", "GDF15", "HERPUD1", "HSPA5", "IARS", "JDP2", "KDM7A", "LARS", "MKNK2", "NARS",
                                           "PTGS2", "SARS", "SQSTM1", "VARS", "VLDLR", "YARS"),])
rownames(df) <- df[, 1]
df <- df[, -(1 : 2)]

dat <- df[, order(as.numeric(df["ATF4",]))]

n <- t(scale(t(dat)))
n[n > 2] <- 2
n[n < -2] <- -2
pheatmap(n , show_colnames = F, show_rownames = T, cluster_cols = F, clustering_distance_rows = "euclidean",
         clustering_method = "complete", fontsize = 12)






# MTHFD2
df <- as.data.frame(t(exprs[which(exprs$X %in% c("ATF4", "MTHFD2")), -(1:2)]))
colnames(df) <- c("MTHFD2", "ATF4")
df <- df[BRAFp,]

cor = cor.test(df$MTHFD2, df$ATF4)$estimate
p_value = cor.test(df$MTHFD2, df$ATF4)$p.value

labels = data.frame(
  x = rep(5, 2),
  y = c(6, 5.5),
  labels = rep( c( paste0('r = ', round(cor, 2)),
                   paste0('p.value = ', formatC(p_value))), 4)
)

add_label = geom_text(data = labels, aes(x = x, y = y, label = labels), size = 4, hjust = 0)

ggplot(data = df, aes(x = ATF4, y = MTHFD2)) +
  geom_point(size = 3.5, color = '#9B70F9', shape = 1) +
  theme_bw() +
  xlab('ATF4') +
  ylab('MTHFD2') +
  geom_smooth(method = 'lm',color = '#FD6F71',fill = '#FD6F71',alpha = 0.25)+
  add_label +
  ggtitle("TCGA (SKCM)")

# GPT2
df <- as.data.frame(t(exprs[which(exprs$X == "ATF4" | exprs$X == "GPT2"), -(1:2)]))
colnames(df) <- c("GPT2", "ATF4")

cor = cor.test(df$GPT2, df$ATF4)$estimate
p_value = cor.test(df$GPT2, df$ATF4)$p.value

labels = data.frame(
  x = rep(5, 2),
  y = c(5, 4.5),
  labels = rep( c( paste0('r = ', round(cor, 2)),
                   paste0('p.value = ', formatC(p_value))), 4)
)

add_label = geom_text(data = labels, aes(x = x, y = y, label = labels), size = 4, hjust = 0)

ggplot(data = df, aes(x = ATF4, y = GPT2)) +
  geom_point(size = 3.5, color = '#9B70F9', shape = 1) +
  theme_bw() +
  xlab('ATF4') +
  ylab('GPT2') +
  geom_smooth(method = 'lm',color = '#FD6F71',fill = '#FD6F71',alpha = 0.25)+
  add_label +
  ggtitle("TCGA (SKCM)")

# PSAT1
df <- as.data.frame(t(exprs[which(exprs$X == "ATF4" | exprs$X == "PSAT1"), -(1:2)]))
colnames(df) <- c("ATF4", "PSAT1")

cor = cor.test(df$PSAT1, df$ATF4)$estimate
p_value = cor.test(df$PSAT1, df$ATF4)$p.value

labels = data.frame(
  x = rep(5, 2),
  y = c(8, 7.5),
  labels = rep( c( paste0('r = ', round(cor, 2)),
                   paste0('p.value = ', formatC(p_value))), 4)
)

add_label = geom_text(data = labels, aes(x = x, y = y, label = labels), size = 4, hjust = 0)

ggplot(data = df, aes(x = ATF4, y = PSAT1)) +
  geom_point(size = 3.5, color = '#9B70F9', shape = 1) +
  theme_bw() +
  xlab('ATF4') +
  ylab('PSAT1') +
  geom_smooth(method = 'lm',color = '#FD6F71',fill = '#FD6F71',alpha = 0.25)+
  add_label +
  ggtitle("TCGA (SKCM)")
