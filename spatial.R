library(dplyr)
library(devtools)
library(Seurat)
library(ggplot2)
library(DoubletFinder)
library(GSEABase)
library(GSVA)

sample02 <- Load10X_Spatial(data.dir = "D:/McGill/projects/project BRAFi/lighthouse/GSE207592/sample02",
                            filename = "filtered_feature_bc_matrix.h5")
sample03 <- Load10X_Spatial(data.dir = "D:/McGill/projects/project BRAFi/lighthouse/GSE207592/sample03",
                            filename = "filtered_feature_bc_matrix.h5")
sample04 <- Load10X_Spatial(data.dir = "D:/McGill/projects/project BRAFi/lighthouse/GSE207592/sample04",
                            filename = "filtered_feature_bc_matrix.h5")

sample <- merge(x = sample02, y = c(sample04, sample03))
sample <- SCTransform(sample, assay = "Spatial", verbose = T)

SpatialFeaturePlot(sample, features = "Atf4", pt.size.factor = 2, alpha = c(0.1, 1))

ATF4 <- as.data.frame(sample@assays$SCT$data[c("Atf4", "Atf3", "Eif1", "Ddit3", "Vegfa",
                                               "Aars", "Map1lc3b", "Gdf15", "Herpud1", "Hspa5",
                                               "Sars", "Sqstm1"),])
ATF4_score <- apply(ATF4, 2, sum) / sqrt(12)
sample$ATF4 <- ATF4_score





mapk <- as.data.frame(sample@assays$SCT$data[c("Phlda1", "Spry2", "Spry4", "Dusp4", "Dusp6",
                                                 "Ccnd1", "Epha2", "Epha4", "Etv4", "Etv5"),])
MPAS <- apply(mapk, 2, sum) / sqrt(10)
sample$MPAS <- MPAS
SpatialFeaturePlot(sample, features = "MPAS", pt.size.factor = 2, alpha = c(0.1, 1),
                   min.cutoff = 1.1)

# ssGSEA
gs <- getGmt("D:/McGill/projects/project BRAFi/RNAseq/mh.all.v2023.2.Mm.symbols.gmt")
set <- names(gs)
geneset <- list()
for (i in set){
  geneset[[i]] <- gs[[i]]@geneIds
}

ssgsea_score <- gsva(as.matrix(sample@assays$SCT$data), geneset, method = "ssgsea", ssgsea.norm = T)

sample@meta.data <- cbind(sample@meta.data, t(ssgsea_score))

# sample$pathway <- ifelse(sample$MPAS > 2, 
#                         ifelse(sample$HALLMARK_MTORC1_SIGNALING > 0.33, "hh", "hl"),
#                         ifelse(sample$HALLMARK_MTORC1_SIGNALING > 0.33, "lh", "ll")
#                         )
# 
# VlnPlot(sample, group.by = "pathway",
#         features = c("MPAS", "HALLMARK_MTORC1_SIGNALING", "Atf4"), pt.size = 0)


# sample$pathway <- ifelse(sample$MPAS > 2.5, 
#                          ifelse(sample$HALLMARK_MTORC1_SIGNALING > 0.27, "hh", "hl"),
#                          ifelse(sample$HALLMARK_MTORC1_SIGNALING > 0.27, "lh", "ll")
# )
# 
# VlnPlot(sample, group.by = "pathway",
#         features = c("MPAS", "HALLMARK_MTORC1_SIGNALING", "Atf4"), pt.size = 0)

sample$pathway <- ifelse(sample$MPAS > quantile(sample$MPAS, 0.8), 
                         ifelse(sample$HALLMARK_MTORC1_SIGNALING > quantile(sample$HALLMARK_MTORC1_SIGNALING, 0.8), "hh", "hl"),
                         ifelse(sample$HALLMARK_MTORC1_SIGNALING > quantile(sample$HALLMARK_MTORC1_SIGNALING, 0.8), "lh", "ll")
)

sample$ATF4 <- exp(sample@assays$SCT$data["Atf4",])

VlnPlot(sample, group.by = "pathway",
        features = c("MPAS", "HALLMARK_MTORC1_SIGNALING", "Atf4"), pt.size = 0)



saveRDS(sample, "D:/McGill/projects/project BRAFi/lighthouse/GSE207592/output/sample.rds")


sample_mpas <- subset(sample, subset = pathway == "hl")
mean(sample_mpas@assays$SCT$data["Atf4",])













# MPAS
df <- cbind(sample@assays$SCT$data["Atf4",], sample$MPAS)

df <- as.data.frame(df)
colnames(df) <- c("Atf4", "MPAS")

cor = cor.test(df[, 2], df[, 1])$estimate
p_value = cor.test(df[, 2], df[, 1])$p.value

labels = data.frame(
  x = rep(0, 2),
  y = c(5, 4.5),
  labels = rep( c( paste0('r = ', round(cor, 2)),
                   paste0('p.value = ', formatC(p_value))), 4)
)

add_label = geom_text(data = labels, aes(x = x, y = y, label = labels), size = 4, hjust = 0)

ggplot(data = df, aes(x = Atf4, y = MPAS)) +
  geom_point(size = 3.5, color = '#9B70F9', shape = 1) +
  theme_bw() +
  xlab('Atf4 Expression Level') +
  ylab('MAPK Pathway Activity Score') +
  geom_smooth(method = 'lm',color = '#FD6F71',fill = '#FD6F71',alpha = 0.25)+
  add_label

# mTOR
df <- cbind(sample@assays$SCT$data["Atf4",], sample$HALLMARK_MTORC1_SIGNALING)

df <- as.data.frame(df)
colnames(df) <- c("Atf4", "mTOR")

cor = cor.test(df[, 2], df[, 1])$estimate
p_value = cor.test(df[, 2], df[, 1])$p.value

labels = data.frame(
  x = rep(0, 2),
  y = c(0.6, 0.5),
  labels = rep( c( paste0('r = ', round(cor, 2)),
                   paste0('p.value = ', formatC(p_value))), 4)
)

add_label = geom_text(data = labels, aes(x = x, y = y, label = labels), size = 4, hjust = 0)

ggplot(data = df, aes(x = Atf4, y = mTOR)) +
  geom_point(size = 3.5, color = '#9B70F9', shape = 1) +
  theme_bw() +
  xlab('Atf4 Expression Level') +
  ylab('mTOR Pathway Activity Score') +
#  geom_smooth(method = 'lm',color = '#FD6F71',fill = '#FD6F71',alpha = 0.25)+
  add_label

# Pmel_Mlana
df <- cbind(sample@assays$SCT$data["Pmel",], sample@assays$SCT$data["Mlana",])

df <- as.data.frame(df)
colnames(df) <- c("Pmel", "Mlana")

cor = cor.test(df[, 2], df[, 1])$estimate
p_value = cor.test(df[, 2], df[, 1])$p.value

labels = data.frame(
  x = rep(0, 2),
  y = c(0.6, 0.5),
  labels = rep( c( paste0('r = ', round(cor, 2)),
                   paste0('p.value = ', formatC(p_value))), 4)
)

add_label = geom_text(data = labels, aes(x = x, y = y, label = labels), size = 4, hjust = 0)

ggplot(data = df, aes(x = Pmel, y = Mlana)) +
  geom_point(size = 2, color = '#9B70F9', shape = 1) +
  theme_bw() +
  xlim(0, 4) +
  ylim(0, 5) +
  xlab('Pmel Expression Level') +
  ylab('Mlana Expression Level')
  #  geom_smooth(method = 'lm',color = '#FD6F71',fill = '#FD6F71',alpha = 0.25)





# MAPK_mTOR
df <- cbind(sample$MPAS, sample$HALLMARK_MTORC1_SIGNALING, sample$pathway)

df <- as.data.frame(df)
colnames(df) <- c("MAPK", "mTOR", "pathway")
df$MAPK <- as.numeric(df$MAPK)
df$mTOR <- as.numeric(df$mTOR)
df$pathway <- as.factor(df$pathway)

cor = cor.test(df$MAPK, df$mTOR)$estimate
p_value = cor.test(df[, 2], df[, 1])$p.value

labels = data.frame(
  x = rep(0, 2),
  y = c(0.6, 0.5),
  labels = rep( c( paste0('r = ', round(cor, 2)),
                   paste0('p.value = ', formatC(p_value))), 4)
)

pathway <- df$pathway
ggplot(data = df, aes(x = MAPK, y = mTOR, colour = pathway)) +
  geom_point(size = 3.5, shape = 1) +
  theme_bw() +
  xlab('MAPK Pathway Activity Score') +
  ylab('mTOR Pathway Activity Score') +
  #  geom_smooth(method = 'lm',color = '#FD6F71',fill = '#FD6F71',alpha = 0.25)+
  theme_classic(base_size = 15)+
  scale_color_manual(values = c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF"))





# MPAS_low
sample_mpas <- subset(sample, subset = pathway %in% c("lh", "ll"))
df <- cbind(sample_mpas$ATF4, sample_mpas$HALLMARK_MTORC1_SIGNALING)

df <- as.data.frame(df)
colnames(df) <- c("Atf4", "mTOR")

cor = cor.test(df[, 2], df[, 1])$estimate
p_value = cor.test(df[, 2], df[, 1])$p.value

labels = data.frame(
  x = rep(0.1, 2),
  y = c(5, 4.7),
  labels = rep( c( paste0('r = ', round(cor, 2)),
                   paste0('p.value = ', formatC(p_value))), 4)
)

add_label = geom_text(data = labels, aes(x = x, y = y, label = labels), size = 4, hjust = 0)

ggplot(data = df, aes(x = mTOR, y = Atf4)) +
  geom_point(size = 3.5, color = '#9B70F9', shape = 1) +
  theme_bw() +
  xlab('mTOR Pathway Activity Score') +
  ylab('ISR score') +
  geom_smooth(method = 'lm',color = '#FD6F71',fill = '#FD6F71',alpha = 0.25)+
  add_label

# MPAS_high
sample_mpas <- subset(sample, subset = pathway %in% c("hh", "hl"))
df <- cbind(sample_mpas$ATF4, sample_mpas$HALLMARK_MTORC1_SIGNALING)

df <- as.data.frame(df)
colnames(df) <- c("Atf4", "mTOR")

cor = cor.test(df[, 2], df[, 1])$estimate
p_value = cor.test(df[, 2], df[, 1])$p.value

labels = data.frame(
  x = rep(0.1, 2),
  y = c(6, 5.7),
  labels = rep( c( paste0('r = ', round(cor, 2)),
                   paste0('p.value = ', formatC(p_value))), 4)
)

add_label = geom_text(data = labels, aes(x = x, y = y, label = labels), size = 4, hjust = 0)

ggplot(data = df, aes(x = mTOR, y = Atf4)) +
  geom_point(size = 3.5, color = '#9B70F9', shape = 1) +
  theme_bw() +
  xlab('mTOR Pathway Activity Score') +
  ylab('ISR score') +
  geom_smooth(method = 'lm',color = '#FD6F71',fill = '#FD6F71',alpha = 0.25)+
  add_label




df <- cbind(sample@assays$SCT$data["Atf4",], sample$MPAS, sample$HALLMARK_MTORC1_SIGNALING, sample$pathway)
colnames(df) <- c("Atf4", "MAPK", "mTOR", "class")
write.csv(df, "D:/McGill/projects/project BRAFi/lighthouse/GSE207592/output/df.csv")






# test
df <- cbind(sample@assays$SCT$data["Atf4",], sample@assays$SCT$data["Yars",])

df <- as.data.frame(df)

cor = cor.test(df[, 2], df[, 1])$estimate
p_value = cor.test(df[, 2], df[, 1])$p.value

