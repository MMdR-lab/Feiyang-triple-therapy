BiocManager::install("GSVA")

library("DESeq2")
library("ggplot2")
library("ggrepel")
library("ggpubr")
library("org.Hs.eg.db")
library("clusterProfiler")
library("dplyr")
library("pheatmap")
library("enrichplot")
library("msigdbr")
library("fgsea")
library("GSVA")
library("GSEABase")

#expression matrix
rawdat <- read.table("D:/McGill/projects/project BRAFi/RNAseq/WM3406-main/salmon.merged.gene_counts.tsv", header = T, sep = "\t")

dat <- rawdat[, 3 : 14]
dat <- apply(dat, 2, as.numeric)
dat <- apply(dat, 2, round)
rownames(dat) <- rawdat$gene_id
dat <- dat[, c(7 : 12, 4 : 6, 1 : 3)]

write.csv(dat, "D:/McGill/projects/project BRAFi/RNAseq/output/data.csv")

boxplot(dat, las = 2)  #check batch effect / las means layout direction


condition <- factor(c(rep("DMSO", 3), rep("INK", 3), rep("BC", 3), rep("BCI", 3)), levels = c("DMSO", "INK", "BC", "BCI"))
colData <- data.frame(row.names = colnames(dat), condition)
colData$condition <- relevel(colData$condition, ref = "DMSO")

dds <- DESeqDataSetFromMatrix(dat, colData, design = ~condition)
dds <- DESeq(dds)




# ---------------------------------------------------------------------------------------

# DMSO vs BC
res <- as.data.frame(results(dds, name = "condition_BC_vs_DMSO"))
res <- res[order(res$padj),]

head(res)
table(res$padj < 0.05)

diff_gene_deseq2 <- subset(res, padj < 0.05 & abs(log2FoldChange) > 1)
head(diff_gene_deseq2)

#threshold for differential expression
res$threshold = ifelse(res$padj > 0.05 | is.na(res$padj) | is.na(res$log2FoldChange) , "stable" ,
                       ifelse(res$log2FoldChange > 1, "up" ,
                              ifelse(res$log2FoldChange < -1, "down" , "stable")
                       )
)

res$gene <- rep(NA, nrow(res))
res[c("LDHA", "PDK1", "ALDOC", "ABAT", "PGK1", "EGLN3", "PGAM1", "BNIP3", "FAM162A"),]$gene <-
  c("LDHA", "PDK1", "ALDOC", "ABAT", "PGK1", "EGLN3", "PGAM1", "BNIP3", "FAM162A")

DEG <- res[which(res$threshold == "up" | res$threshold == "down"),]

write.csv(DEG, "D:/McGill/projects/project BRAFi/RNAseq/output/DMSO_BC_DEG1.csv")

# volcano
ggplot(res, aes(x = log2FoldChange, y = -log10(padj), colour = threshold)) +
  xlab("log2 fold-change") +
  ylab("-log10 padj") +
  geom_point() +
  geom_hline(yintercept = 1.3, linetype = "dashed", size = 1) +
  geom_vline(xintercept = c(-1, 1), linetype = 4, size = 1) +
  xlim(-10, 10) + ylim(0, 270) +
  theme_classic(base_size = 15) +
  scale_color_manual(values = c("#377EB8", "grey", "#E41A1C")) +
  geom_text_repel(aes(label = gene), size = 3, color = 'black') +
  theme(legend.title = element_blank(),
        legend.position = c(0.5, 0.9),
        legend.background = element_rect(fill = 'transparent'),
        axis.text = element_text(colour = "black"),
        panel.border = element_rect(size = 1, fill = 'transparent'))

# ---------------------------------------------------------------------------------------

# DMSO vs INK
res <- as.data.frame(results(dds, name = "condition_INK_vs_DMSO"))
res <- res[order(res$padj),]

head(res)
table(res$padj < 0.05)

diff_gene_deseq2 <- subset(res, padj < 0.05 & abs(log2FoldChange) > 1)
head(diff_gene_deseq2)

#threshold for differential expression
res$threshold = ifelse(res$padj > 0.05 | is.na(res$padj) | is.na(res$log2FoldChange) , "stable" ,
                       ifelse(res$log2FoldChange > 1, "up" ,
                              ifelse(res$log2FoldChange < -1, "down" , "stable")
                       )
)

res$gene <- rep(NA, nrow(res))
res[1 : 15,]$gene <- rownames(res[1 : 15,])

DEG <- res[which(res$threshold == "up" | res$threshold == "down"),]

write.csv(DEG, "D:/McGill/projects/project BRAFi/RNAseq/output/DMSO_INK_DEG.csv")


ggplot(res, aes(x = log2FoldChange, y = -log10(padj), colour = threshold)) +
  xlab("log2 fold-change") +
  ylab("-log10 padj") +
  geom_point() +
  geom_hline(yintercept = 1.3, linetype = "dashed", size = 1) +
  geom_vline(xintercept = c(-1, 1), linetype = 4, size = 1) +
  xlim(-3, 3) + ylim(0, 30) +
  theme_classic(base_size = 15) +
  scale_color_manual(values = c("#377EB8", "grey", "#E41A1C")) +
  geom_text_repel(aes(label = gene), size = 3, color = 'black') +
  theme(legend.title = element_blank(),
        legend.position = c(0.5, 0.9),
        legend.background = element_rect(fill = 'transparent'),
        axis.text = element_text(colour = "black"),
        panel.border = element_rect(size = 1, fill = 'transparent'))

# ---------------------------------------------------------------------------------------

# DMSO vs BCI
res <- as.data.frame(results(dds, name = "condition_BCI_vs_DMSO"))
res <- res[order(res$padj),]

head(res)
table(res$padj < 0.05)

diff_gene_deseq2 <- subset(res, padj < 0.05 & abs(log2FoldChange) > 1)
head(diff_gene_deseq2)

#threshold for differential expression
res$threshold = ifelse(res$padj > 0.05 | is.na(res$padj) | is.na(res$log2FoldChange) , "stable" ,
                       ifelse(res$log2FoldChange > 1.5, "up" ,
                              ifelse(res$log2FoldChange < -1.5, "down" , "stable")
                       )
)

res$gene <- rep(NA, nrow(res))
res[1 : 10,]$gene <- rownames(res[1 : 10,])

DEG <- res[which(res$threshold == "up" | res$threshold == "down"),]

write.csv(DEG, "D:/McGill/projects/project BRAFi/RNAseq/output/DMSO_BCI_DEG.csv")


ggplot(res, aes(x = log2FoldChange, y = -log10(padj), colour = threshold)) +
  xlab("log2 fold-change") +
  ylab("-log10 padj") +
  geom_point() +
  geom_hline(yintercept = 1.3, linetype = "dashed", size = 1) +
  geom_vline(xintercept = c(-1.5, 1.5), linetype = 4, size = 1) +
  xlim(-10, 10) + ylim(0, 230) +
  theme_classic(base_size = 15) +
  scale_color_manual(values = c("#377EB8", "grey", "#E41A1C")) +
  geom_text_repel(aes(label = gene), size = 3, color = 'black') +
  theme(legend.title = element_blank(),
        legend.position = c(0.5, 0.9),
        legend.background = element_rect(fill = 'transparent'),
        axis.text = element_text(colour = "black"),
        panel.border = element_rect(size = 1, fill = 'transparent'))




# ---------------------------------------------------------------------------------------

# BC vs BCI
res <- as.data.frame(results(dds, c("condition", "BCI", "BC")))
res <- res[order(res$padj),]

head(res)
table(res$padj < 0.05)

diff_gene_deseq2 <- subset(res, padj < 0.05 & abs(log2FoldChange) > 1)
head(diff_gene_deseq2)

#threshold for differential expression
res$threshold = ifelse(res$padj > 0.05 | is.na(res$padj) | is.na(res$log2FoldChange) , "stable" ,
                       ifelse(res$log2FoldChange > 1, "up" ,
                              ifelse(res$log2FoldChange < -1, "down" , "stable")
                       )
)

res$gene <- rep(NA, nrow(res))
res$gene1 <- rep(NA, nrow(res))
res[c("LDHA", "PDK1", "ALDOC", "ABAT", "PGK1", "EGLN3", "PGAM1", "BNIP3", "FAM162A"),]$gene <-
  c("LDHA", "PDK1", "ALDOC", "ABAT", "PGK1", "EGLN3", "PGAM1", "BNIP3", "FAM162A")
res[c("CDT1", "ANGPTL7", "GAS5", "PIK3R3", "RHCG", "IL1B", "SLC16A9", "SLC16A6", "IL24"),]$gene1 <-
  c("CDT1", "ANGPTL7", "GAS5", "PIK3R3", "RHCG", "IL1B", "SLC16A9", "SLC16A6", "IL24")

DEG <- res[which(res$threshold == "up" | res$threshold == "down"),]

write.csv(DEG, "D:/McGill/projects/project BRAFi/RNAseq/output/BC_BCI_DEG.csv")


ggplot(res, aes(x = log2FoldChange, y = -log10(padj), colour = threshold)) +
  xlab("log2 fold-change") +
  ylab("-log10 padj") +
  geom_point() +
  geom_hline(yintercept = 1.3, linetype = "dashed", size = 1) +
  geom_vline(xintercept = c(-1, 1), linetype = 4, size = 1) +
  xlim(-3, 3) + ylim(0, 30) +
  theme_classic(base_size = 15) +
  scale_color_manual(values = c("#377EB8", "grey", "#E41A1C")) +
  geom_text_repel(aes(label = gene), size = 4, color = 'blue') +
  geom_text_repel(aes(label = gene1), size = 3, color = 'black') +
  theme(legend.title = element_blank(),
        legend.position = c(0.5, 0.9),
        legend.background = element_rect(fill = 'transparent'),
        axis.text = element_text(colour = "black"),
        panel.border = element_rect(size = 1, fill = 'transparent'))

# -------------------------------------------------------------------------

# heatmap

BC <- read.csv("D:/McGill/projects/project BRAFi/RNAseq/output/DMSO_BC_DEG1.csv")
INK <- read.csv("D:/McGill/projects/project BRAFi/RNAseq/output/DMSO_INK_DEG.csv")
BCI <- read.csv("D:/McGill/projects/project BRAFi/RNAseq/output/BC_BCI_DEG.csv")
BCI1 <- read.csv("D:/McGill/projects/project BRAFi/RNAseq/output/DMSO_BCI_DEG.csv")

length(which(log(BC$padj) < -40))

BC_3 <- BC[which(log(BC$padj) < -40 & abs(BC$log2FoldChange) > 1.5),]

cg <- intersect(BCI$X[which(BCI$threshold == "up")], BCI1$X[which(BCI1$threshold == "up")]) # BCI vs BC
cg <- setdiff(BCI$X[which(BCI$threshold == "down")], BCI1$X[which(BCI1$threshold == "up")]) # BCI vs BC

cg <- union(BC_3$X, union(INK$X, BCI$X)) # DEG
cg <- c("MITF", "MLANA", "ALCAM", "GFRA2", "PMEL", "SOX9", "SOX10" ,"NGFR", "CD36", "AXL", "ZEB1", "ZEB2") # pheno swi
cg <- c("ELOVL6", "PLD4", "HMGCS1", "LDHA", "HKDC1", "TK1", "PGAM1", "SNCA",
        "LDLR", "SCD", "SCD5", "LIPH", "PDK1", "AK4", "ACSBG1", "ALDOC", "DHCR24",
        "ERG28", "DHCR7", "ACAT2", "SQLE", "MSMO1", "VLDLR") # metabolism
cg <- c("ELOVL6", "ELOVL7", "PLD1", "PLD4", "HMGCS1", "CD36",
        "LDLR", "SCD", "SCD5", "LIPH", "ACSBG1", "DHCR24", "TM7SF2",
        "DHCR7", "ACAT2", "SQLE", "MSMO1", "VLDLR", "VLDLR-AS1") # lipid
cg <- c("PLD1", "PLD4", "CD36", "SCD", "SCD5", "LIPH", "ACSBG1", "VLDLR", "VLDLR-AS1") # TG
cg <- c("LDHA", "HKDC1", "PDK1", "ALDOC", "ABAT", "PGK1", "SLC2A1", "SLC2A1-AS1", "SLC2A8", "SLC2A9",
        "SLC2A12", "EGLN3", "E2F1", "CD24", "PGAM1", "HIF1A", "HIF1A-AS1", "MYC") # hypoxia

cg <- c("LDHA", "PDK1", "ALDOA", "ALDOC", "ABAT", "PGK1", "SLC16A1", "SLC16A2", "SLC16A3", "SLC16A4",
        "EGLN3", "CD24", "PGAM1", "HIF1A", "HIF1A-AS1", "EPAS1", "BNIP3", "BNIP3L", "FAM162A") # hypoxia

cg <- c("RPA3", "LIG1", "POLD1", "POLA1", "PRIM1", "POLD3", "RFC2", "RFC4", "NME1", "RFC5", "TYMS",
        "POLA2", "RFC3", "RAD51", "FEN1", "PNP", "ZWINT", "ALYREF", "UMPS", "HPRT1", "POLR1C", "POLR2H",
        "SSRP1", "DUT", "POLR2D", "PCNA", "RPA2", "ERCC1", "BRCA1", "BRCA2", "TP53", "TOP2A", "TOP2B") # repair
cg <- c("RFC2", "RFC4", "RFC5", "TYMS", "RFC3", "RAD51", "FEN1", "PCNA", "RPA2", "BRCA1", "BRCA2", "TP53") # repair hypo

cg <- c("RPA3", "LIG1", "POLD1", "POLA1", "PRIM1", "POLD3", "RFC2", "RFC4", "RFC5",
        "POLA2", "RFC3", "RAD51", "FEN1", "ZWINT", "POLR1C", "POLR2H",
        "POLR2D", "PCNA", "RPA2", "ERCC1", "BRCA1", "BRCA2", "TOP2A") # HR related
cg <- c("TYMS", "PNP", "UMPS", "HPRT1", "DUT", "PNPT1", "AK4", "AK1", "AK2", "TK1", "SLC2A9", "SLC16A9") # nucleotide synthesis related
cg <- c("NME1", "ALYREF", "SSRP1") # others

cg <- c("PNP", "PNPT1", "AK4", "AK8", "DBI", "TK1", "SLC2A9", "SLC16A9")
cg <- c("CDC6", "CDT1", "PLK1", "NDC80", "TTK", "CENPF", "CENPE") # cell cycle..

cg <- c("PPARGC1A", "PPARGC1B", "IDH3A", "NDUFS3", "CYCS", "COX5A", "ESRRA", "ACADM", "ACADL") # PGC1a/b

cg <- c("CCL26", "CXCL1", "CXCL11", "CXCL12", "STAT2", "IRF5", "IRF7", "IRF9", "IFNLR1", "GZMM")

cg <- c("BCL2", "BCL2L1", "MCL1", "BAD", "BAX", "BIK", "BCL2L11", "BID", "BAK1", "BBC3")

cg <- c("FOS", "FOSL1", "FOSL2", "JUN", "JUNB", "JUND", "ATF2", "ATF3", "ATF4", "ATF5", "ATF6") # AP-1 transcription factor

cg <- c("FKBP2", "YWHAH", "PPP2R1B", "MARS", "GARS", "AUTS2", "BCAT1", "PTPRG", "SLC38A1", "SLC38A2", "XBP1", "HERPUD1",
        "ATF4", "DDIT3", "SLC1A4", "BNIP3L", "EIF1", "FOXD3", "MAP1LC3A", "NUPR1", "PYCR1", "TPM1", "FTH1", "RAB32",
        "MYL12A", "DDIT4", "ATF3", "CEBPB", "CEBPG", "ERN1", "HMOX1", "PPP1R15A") # ATF4
cg <- c("CEBPB", "CEBPG", "PPP1R15A", "ATF4", "EIF1", "DDIT3", "ATF3")
cg <- c("DDIT3", "ASNS", "TRIB3", "ATF3", "VEGFA", "MTHFD2", "SLC7A11", "AARS", "CEBPB", "CHAC1",
        "DDIT4", "GPT2", "MAP1LC3B", "PPP1R15A", "PSAT1", "WARS", "ALDH18A1", "ATG7", "EIF2S2",
        "EPRS", "FGF19", "GARS", "GDF15", "HERPUD1", "HSPA5", "IARS", "JDP2", "KDM7A", "LARS", "MKNK2", "NARS",
        "PTGS2", "SARS", "SQSTM1", "VARS", "VLDLR", "YARS") # ATF4 target
cg <- c("DDIT3", "ASNS", "TRIB3", "ATF3", "VEGFA", "SLC7A11", "AARS", "CEBPB", "CHAC1", "GPT2",
        "MAP1LC3B", "PPP1R15A", "PSAT1", "WARS", "MTHFD2", "DDIT4") # ATF4 target high confidence
cg <- c("DDIT3", "TRIB3", "CEBPB", "CEBPD", "CEBPG", "ATF3", "JDP2", "NFE2L1") # ATF4 target & interactor
cg <- c("BECN1", "DPF2", "G0S2", "GHITM", "MCL1", "NLRP1", "PMAIP1", "BBC3", "SNAI2", "TP53BP2", "DNAJA3") # ATF4 target & apoptosis
cg <- c("MAP1LC3B", "ATG3", "ATG7", "SQSTM1", "BECN1", "PPM1D") # ATF4 target & autophagy
cg <- c("MTHFD2", "GPT2", "ALDH18A1", "ALDH1L2", "ALDH2", "DNAJA3", "GHITM", "LONP1", "PCK2", "TMEM11") # ATF4 target & mito
cg <- c("XPOT", "AARS", "WARS", "EPRS", "GARS", "IARS", "LARS", "NARS", "SARS", "VARS", "YARS", "CARS", "FARSB",
        "HARS", "NARS", "TARS") # ATF4 target & tRNA
cg <- c("SLC7A11", "SLC3A2", "SLC7A1", "SLC7A5", "SLC38A2") #ATF4 target & amino acid transporter

cg <- c("IDH1", "IDH2", "ME1", "G6PD", "PGD", "GLS", "SLC1A5", "SLC1A3", "GSS", "SLC7A11",
        "MTHFD2", "PHGDH", "PSAT1", "HMOX1") # NADPH GSH one-carbon

cg <- c("PPAT", "GART", "PRPS1", "IMPDH1", "PRPS2", "ADSL", "PAICS", "ATIC", "IMPDH2") # purine synthesis
cg <- c("RPIA", "TALDO1", "G6PD", "TKT", "RPE") # pentose phosphate
cg <- c("MTHFD2", "GLDC", "PSPH", "SHMT2", "PSAT1", "DHFR", "SHMT1", "MTHFD1L", "MTHFD1", "PHGDH") # THF pathway

cg <- c("PRKDC", "XRCC5", "XRCC6", "MRE11A", "NHEJ1", "XRCC4", "LIG4", "RAD50", "NBN", "DCLRE1C") # NHEJ

cg <- c("KMT2A", "KMT2B", "KMT2C", "KMT2D", "KMT2E", "KDM1A", "KDM1B", "KDM5A", "KDM5B", "KDM5C", "SETD1A") # H3K4me3
cg <- c("SUV39H1", "SUV39H2", "SETDB1", "SETDB2", "KDM4A", "KDM4B", "KDM4C", "KDM4D") # H3K9me3

# cg <- intersect(BC$X[which(BC$threshold == "up")], BCI$X[which(BCI$threshold == "down")])
# cg <- intersect(BC$X[which(BC$threshold == "down")], BCI$X[which(BCI$threshold == "up")])

cg <- c("ATG3", "ATG7", "BECN1") # ATF4 target & autophagy
cg <- c("BECN1", "G0S2", "MCL1", "PMAIP1", "DNAJA3") # ATF4 target & apoptosis
cg <- c("MTHFD2", "GPT2", "DNAJA3", "LONP1", "TMEM11") # ATF4 target & mito
cg <- c("XPOT", "GARS", "IARS", "LARS", "NARS", "VARS", "YARS", "FARSB", "HARS", "NARS", "TARS") # ATF4 target & tRNA
cg <- c("SLC7A11", "SLC3A2", "SLC7A1", "SLC7A5") #ATF4 target & amino acid transporter




n <- t(scale(t(dat[cg,])))
n[n > 2] <- 2
n[n < -2] <- -2
condition <- factor(c(rep("DMSO", 3), rep("INK", 3), rep("BC", 3), rep("BCI", 3)), levels = c("DMSO", "INK", "BC", "BCI"))
ac <- data.frame(groupList = condition)  #choose subject to analyze
rownames(ac) <- colnames(n)
pheatmap(n , show_colnames = F, show_rownames = T, cluster_cols = F, clustering_distance_rows = "euclidean",
         cellwidth = 25, cellheight = 15, clustering_method = "complete", annotation_col = ac, fontsize = 12)

x <- pheatmap(n , show_colnames = F, show_rownames = T, cluster_cols = F, clustering_distance_rows = "euclidean",
              clustering_method = "complete", annotation_col = ac)

cluster <- x$tree_row
plot(cluster,hang = -1,cex=0.6,axes=FALSE,ann=FALSE)

cut <- cutree(cluster, 3)
cg <- names(cut)[cut == 3]

cut <- cutree(cluster, 5)
cg <- names(cut)[cut == 5]

# ------------------------------------------------------------------------

BC <- read.csv("D:/McGill/projects/project BRAFi/RNAseq/output/DMSO_BC_DEG1.csv")
INK <- read.csv("D:/McGill/projects/project BRAFi/RNAseq/output/DMSO_INK_DEG.csv")
BCI <- read.csv("D:/McGill/projects/project BRAFi/RNAseq/output/BC_BCI_DEG.csv")
BCI1 <- read.csv("D:/McGill/projects/project BRAFi/RNAseq/output/DMSO_BCI_DEG.csv")

rownames(BC) <- BC$X
BC_repair <- BC[cg,]
write.csv(BC_repair, "D:/McGill/projects/project BRAFi/RNAseq/output/DMSO_BC_repair.csv", row.names = F)

# -------------------------------------------------------------

BC <- read.csv("D:/McGill/projects/project BRAFi/RNAseq/output/DMSO_BC_DEG1.csv")
INK <- read.csv("D:/McGill/projects/project BRAFi/RNAseq/output/DMSO_INK_DEG.csv")
BCI <- read.csv("D:/McGill/projects/project BRAFi/RNAseq/output/BC_BCI_DEG.csv")
BCI1 <- read.csv("D:/McGill/projects/project BRAFi/RNAseq/output/DMSO_BCI_DEG.csv")

cg <- intersect(BC$X[which(BC$threshold == "up")], BCI$X[which(BCI$threshold == "up")])
cg <- intersect(BC$X[which(BC$threshold == "down")], BCI$X[which(BCI$threshold == "down")])
cg <- intersect(BC$X[which(BC$threshold == "up")], BCI$X[which(BCI$threshold == "down")])
cg <- intersect(BC$X[which(BC$threshold == "down")], BCI$X[which(BCI$threshold == "up")])

cg <- intersect(BC$X[which(BC$threshold == "up")], BCI1$X[which(BCI1$threshold == "up")])
cg <- intersect(BC$X[which(BC$threshold == "down")], BCI1$X[which(BCI1$threshold == "down")])

cg <- setdiff(BCI1$X[which(BCI1$threshold == "up")], BC$X[which(BC$threshold == "down")])
cg <- setdiff(BCI1$X[which(BCI1$threshold == "down")], BC$X[which(BC$threshold == "up")])

id2symbol = toTable(org.Hs.egSYMBOL)
id <- id2symbol$gene_id[na.omit(match(cg, id2symbol$symbol))]

# -----------------------------------------------------------------------------

#GO enrichment
symbol_up = res$X[which(res$threshold == "up")]
symbol_down = res$X[which(res$threshold == "down")]
id2symbol = toTable(org.Hs.egSYMBOL)
id_up = id2symbol$gene_id[na.omit(match(symbol_up, id2symbol$symbol))]
id_down = id2symbol$gene_id[na.omit(match(symbol_down, id2symbol$symbol))]
id <- union(id_up, id_down)

# ---------------------------------------------------------------------

#choose up or down
ego <- enrichGO(OrgDb = "org.Hs.eg.db", gene = id, ont = "ALL", pvalueCutoff = 0.05, readable = T) # BP MF CC

saveRDS(ego, "D:/McGill/projects/project BRAFi/RNAseq/output/DMSO_BC_GO.rds")
saveRDS(ego, "D:/McGill/projects/project BRAFi/RNAseq/output/DMSO_INK_GO.rds")
saveRDS(ego, "D:/McGill/projects/project BRAFi/RNAseq/output/BC_BCI_GO.rds")
saveRDS(ego, "D:/McGill/projects/project BRAFi/RNAseq/output/DMSO_BCI_GO.rds")

saveRDS(ego, "D:/McGill/projects/project BRAFi/RNAseq/output/BC_BCI_down_GO.rds")

# a <- as.data.frame(ego@result$Description)

ego <- readRDS("D:/McGill/projects/project BRAFi/RNAseq/output/DMSO_BCI_GO.rds")

ego_plot <- ego

ego_plot@result <- ego@result[c(2, 3, 4, 5, 6, 7, 8, 28, 29, 33),] # DMSO vs BC
ego_plot@result <- ego@result[c(3, 4, 6, 12, 28, 35, 44, 57, 59, 69),] # DMSO vs INK
ego_plot@result <- ego@result[c(1, 2, 18, 20, 23, 26, 30, 37, 39, 40),] # BC vs BCI


dotplot(ego_plot, showCategory = 10, title = "Enrichment GO")  # bubble plot
dotplot(ego, showCategory = 10, title = "Enrichment GO")  # bubble plot
barplot(ego, showCategory = 10, title = "EnrichmentGO")

write.csv(ego@result, "D:/McGill/projects/project BRAFi/RNAseq/output/BC_BCI_down_GO.csv")

ego <- readRDS("D:/McGill/projects/project BRAFi/RNAseq/output/down_GO.rds")

#GSEA
d <- res[which(res$threshold == "up" | res$threshold == "down"),]
d <- BC[which(BC$threshold %in% c("up", "down")),]
d <- BCI1[which(BCI1$threshold == "up" | BCI1$threshold == "down"),]

names(d)[1] <- "SYMBOL"
rownames(d) <- d$SYMBOL

id <- bitr(rownames(d), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")

geneList <- merge(d, id, by = "SYMBOL", all = FALSE)
geneList = geneList[order(geneList$log2FoldChange, decreasing = TRUE),]

gene.expr <- geneList$log2FoldChange
names(gene.expr) <- geneList$ENTREZID

go <- gseGO(gene.expr, ont = "ALL", OrgDb = org.Hs.eg.db)
sortgo <- go[order(go$enrichmentScore, decreasing = TRUE),]
gseaplot2(go, 
          go@result$ID[125],
          title = go@result$Description[125],
          base_size = 15,
          color = "green",
          pvalue_table = T,
          ES_geom = "line")
write.csv(sortgo, "D:/McGill/projects/project BRAFi/RNAseq/output/BC_BCI_GSEA.csv")
write.csv(go@result, "D:/McGill/projects/project BRAFi/RNAseq/output/DMSO_BC_GSEA.csv")


# hallmarks
msigdbr_df <- msigdbr(species = "Homo sapiens", category = "H")
pathwaysH <- msigdbr_df[, c("gs_name", "entrez_gene")]
gsea <- GSEA(gene.expr, TERM2GENE =  pathwaysH)

gseaplot2(gsea, 
          gsea@result$ID[11],
          title = gsea@result$Description[11],
          base_size = 15,
          color = "green",
          rel_heights = c(1.5, 0.5, 1),
          pvalue_table = F,
          ES_geom = "line")


go.gene <- go@result[["core_enrichment"]][272] # 48 PGC1a
go.gene <- unlist(strsplit(go.gene, split = "/"))
go.gene <- bitr(go.gene, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = "org.Hs.eg.db")
cg <- go.gene$SYMBOL

# ssGSEA
gs <- getGmt("D:/McGill/projects/project BRAFi/RNAseq/h.all.v2023.2.Hs.symbols.gmt")
set <- names(gs)
geneset <- list()
for (i in set){
  geneset[[i]] <- gs[[i]]@geneIds
}
ssgsea_score <- gsva(as.matrix(FPKM), geneset, method = "ssgsea", ssgsea.norm = T)


# ----------------------------------------------------------------------------------------------

#FPKM
GRCh38 <- read.csv("D:/McGill/projects/project BRAFi/RNAseq/coding_genes_GRCH38.csv" , head = T)

reads_sum <- apply(dat, 2, sum)
#exon length
loc <- match(rownames(dat), GRCh38$gene_name)
lens <- rep(1, nrow(dat))
for (n in 1 : nrow(dat)){
  if(!is.na(loc[n])){
    lens[n] <- GRCh38$end[loc[n]] - GRCh38$start[loc[n]]
  }
}
reads_sum <- reads_sum / 10 ^ 6
lens <- lens / 10 ^ 3
FPKM <- t(apply(dat, 1, function(x)x / reads_sum))
FPKM <- apply(FPKM, 2, function(x)x / lens)
FPKM <- FPKM[-which(lens == 1),]  #filter NA gene

write.csv(FPKM, "D:/McGill/projects/project BRAFi/RNAseq/output/FPKM.csv")

# -----------------------------------------------------------------------------------------------
FPKM <- read.csv("D:/McGill/projects/project BRAFi/RNAseq/output/FPKM.csv")
rownames(FPKM) <- FPKM$X
FPKM <- FPKM[, -1]
FPKM_repair <- FPKM[c("RPA3", "LIG1", "POLD1", "POLA1", "PRIM1", "POLD3", "RFC2", "RFC4", "NME1", "RFC5", "TYMS",
                      "POLA2", "RFC3", "RAD51", "FEN1", "PNP", "ZWINT", "ALYREF", "UMPS", "HPRT1", "POLR1C", "POLR2H",
                      "SSRP1", "DUT", "POLR2D", "PCNA", "RPA2", "ERCC1"),]
write.csv(FPKM_repair, "D:/McGill/projects/project BRAFi/RNAseq/output/FPKM_repair.csv")

# ATF4
ATF4 <- read.csv("D:/McGill/projects/project BRAFi/lighthouse/ATF4 target genes/mmc6.csv")
FPKM_ATF4 <- FPKM[unique(ATF4$Gene.Symbol),]
FPKM_ATF4 <- FPKM_ATF4[-110,]

# ATF4 target

ATF4_target <- res[cg, c(2, 6)]
ATF4_target <- ATF4_target[order(ATF4_target[,1]),]
write.csv(ATF4_target, "D:/McGill/projects/project BRAFi/RNAseq/output/ATF4_target.csv")
