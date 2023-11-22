setwd('D:/My Drive/PostDoc UCI/New Project PCOS Leandro UCI 2021/BxD PCOS Project 2021/Results/2022_2023_Data table and  R environment/temporal Leandro')

library(ggplot2)
library(dplyr)
library(ggrepel)
library(forcats)
library(reshape2)
library(RColorBrewer)
library(ggpubr)
library(WGCNA)
traits_data = all_traits_BxdPcos
table(traits_data$trait_name)
head(traits_data)
traits_data$value = as.numeric(traits_data$value)
traits_data = traits_data[!is.na(traits_data$value),]

##let
pcos_data = traits_data[traits_data$treatment=='let',]
avg_table = pcos_data


#Make and example scatterplot, example with 2 traits who shows differential correlation in control and let
trait1 = 'AUC'
trait2 = 'lean/fat_ratio'

tr1 = avg_table[avg_table$trait_name %in% trait1,]
tr1$trait = paste0(trait1)
tr2 = avg_table[avg_table$trait_name %in% trait2,]
tr2$trait = paste0(trait2)
tr1$tr2 = tr2$value[match(tr1$mouse_number, tr2$mouse_number)]
head(tr1)

sum_sf = tr1 %>%  dplyr::group_by(strain) %>% dplyr::summarise(t1_mean = mean(value, na.rm=T), tr1_sd = sd(value, na.rm=T), tr2_mean = mean(tr2, na.rm=T) ,tr2_sd = sd(tr2, na.rm=T) )
sum_sf
# pdf(file = paste0('trait scatter ', trait1, ' - ', trait2, '.pdf'))
ggplot(data = sum_sf ,aes(x = tr2_mean,y = t1_mean, color = strain, fill = strain)) + 
  geom_point() + 
  geom_errorbar(aes(ymin = t1_mean-tr1_sd,ymax = t1_mean+tr1_sd)) + 
  geom_errorbarh(aes(xmin = tr2_mean-tr2_sd,xmax = tr2_mean+tr2_sd)) +
  theme_pubr() + theme(legend.position = "none") + geom_text_repel(label=sum_sf$strain) +  geom_smooth(aes(group=1), 
  method = "lm", se=T, formula = 'y ~ x') +
 ggtitle('AUC x lean/fat_ratio, strain means (+/- sd) pcos-m BxD PCOS')
# dev.off()



##ctrl
pcos_data = traits_data[traits_data$treatment=='ctrl',]
avg_table = pcos_data

#Make and example scatterplot
trait1 = 'AUC'
trait2 = 'lean/fat_ratio'

tr1 = avg_table[avg_table$trait_name %in% trait1,]
tr1$trait = paste0(trait1)
tr2 = avg_table[avg_table$trait_name %in% trait2,]
tr2$trait = paste0(trait2)
tr1$tr2 = tr2$value[match(tr1$mouse_number, tr2$mouse_number)]
head(tr1)

sum_sf = tr1 %>%  dplyr::group_by(strain) %>% dplyr::summarise(t1_mean = mean(value, na.rm=T), tr1_sd = sd(value, na.rm=T), tr2_mean = mean(tr2, na.rm=T) ,tr2_sd = sd(tr2, na.rm=T) )
sum_sf
# pdf(file = paste0('trait scatter ', trait1, ' - ', trait2, '.pdf'))
ggplot(data = sum_sf ,aes(x = tr2_mean,y = t1_mean, color = strain, fill = strain)) + 
  geom_point() + 
  geom_errorbar(aes(ymin = t1_mean-tr1_sd,ymax = t1_mean+tr1_sd)) + 
  geom_errorbarh(aes(xmin = tr2_mean-tr2_sd,xmax = tr2_mean+tr2_sd)) +
  theme_pubr() + theme(legend.position = "none") + geom_text_repel(label=sum_sf$strain) +  geom_smooth(aes(group=1), 
  method = "lm", se=T, formula = 'y ~ x') +
  ggtitle('AUC x lean/fat_ratio, (+/- sd) controls BxD PCOS')
# dev.off()




tr1$mouse_number[1:10]

resid_table_extract = function(trait1, trait2){
pcos_data = traits_data[traits_data$treatment=='let',]
avg_table = pcos_data
tr1 = avg_table[avg_table$trait_name %in% trait1,]
tr2 = avg_table[avg_table$trait_name %in% trait2,]
tr1$mean_trt2 = tr2$value[match(tr1$mouse_number, tr2$mouse_number)]
fit <- lm(value ~ mean_trt2, data=tr1)  
resids1 = as.data.frame(fit$residuals)
colnames(resids1) = 'residuals'
#row.names(resids1) = row.names(tr1)
resids1$strain = tr1$strain[match(row.names(resids1), row.names(tr1))]
resids1$category = paste0('PCOS')
full_table = resids1

pcos_data = traits_data[traits_data$treatment=='ctrl',]
avg_table = pcos_data
tr1 = avg_table[avg_table$trait_name %in% trait1,]

tr2 = avg_table[avg_table$trait_name %in% trait2,]
tr1$mean_trt2 = tr2$value[match(tr1$mouse_number, tr2$mouse_number)]

fit <- lm(value ~ mean_trt2, data=tr1)  
resids1 = as.data.frame(fit$residuals)
colnames(resids1) = 'residuals'
#row.names(resids1) = row.names(tr1)
resids1$strain = tr1$strain[match(row.names(resids1), row.names(tr1))]
resids1$category = paste0('cntl')

full_table = as.data.frame(rbind(full_table, resids1))
full_table$trait_comboID = paste0(trait1, '_', trait2, '-', full_table$category)
return(full_table)
}

g1 = resid_table_extract('AUC', 'lean/fat_ratio')
g2 = resid_table_extract('Testosterone', 'Insulin')
g3 = resid_table_extract('E2/T', 'Insulin')
g4 = resid_table_extract('estradiol', 'stroke_vol')
g5 = resid_table_extract('LH_FSH_ratio', 'eject_frac')

full_combined_traits = as.data.frame(rbind(g1, g2, g3, g4, g5))
full_combined_traits$strain
nn1 = dcast(full_combined_traits, strain ~ trait_comboID,  value.var = 'residuals', fun.aggregate = mean, na.rm=T)
nn1 = nn1[!is.na(nn1$strain),]
row.names(nn1) = nn1$strain
nn1$strain=NULL
nn1[is.na(nn1)] <- 0

j <- sapply(nn1, is.numeric)
nn1[j] <- scale(nn1[j])
#pdf(file = 'strain residuals fro trait x trait correlations.pdf')
breaksList = seq(min(nn1), max(nn1), by = .1)
pheatmap::pheatmap(nn1,  cluster_cols = F, cluster_rows = T, color = colorRampPalette(rev(brewer.pal(n = 7, name = "PuOr")))(length(breaksList)), main = 'strain residuals for trait x trait correlations' )
#dev.off()


####Correlate genes with residuals

pcos_data = traits_data[traits_data$treatment=='let',]
avg_table = pcos_data
tr1 = avg_table[avg_table$trait_name %in% trait1,]
tr2 = avg_table[avg_table$trait_name %in% trait2,]
tr1$mean_trt2 = tr2$value[match(tr1$mouse_number, tr2$mouse_number)]
fit <- lm(value ~ mean_trt2, data=tr1)  
resids1 = as.data.frame(fit$residuals)
colnames(resids1) = 'residuals'
row.names(resids1) = tr1$mouse_number

counts_ovary = cnts_mat_ov

counts_ovary$sampleID = rownames(counts_ovary)


ov1 = counts_ovary[as.character(counts_ovary$sampleID) %in% row.names(resids1),]



rr1 = resids1[row.names(resids1) %in% row.names(ov1),]
cc1 = bicorAndPvalue(rr1, ov1, use = 'p')
full_gene_list = reshape2::melt(cc1$bicor)
full_gene_list$Var1=NULL
colnames(full_gene_list) = c('gene_symbol', 'bicor')
full_gene_list$pvalue = reshape2::melt(cc1$p)$value
full_gene_list = full_gene_list[order(full_gene_list$pvalue, decreasing = F),]
ff1 = full_gene_list[!grepl('Gm', full_gene_list$gene_symbol),]

ov2 = counts_ovary
trt1 = traits_data[traits_data$trait_name==trait1,]
trt2 = traits_data[traits_data$trait_name==trait2,]
ov2$cat1 = trt1$value[match(ov2$sampleID, trt1$mouse_number)]
ov2$treatment = trt1$treatment[match(ov2$sampleID, trt1$mouse_number)]
ov2$cat2 = trt2$value[match(ov2$sampleID, trt2$mouse_number)]


ov11 = ov2[ov2$gene_symbol ==ff1$gene_symbol[1],]
g1 = ggscatter(ov11, x = "tpm", y = "cat1",
          color = "treatment", palette = "jco", 
          add = "reg.line",
          shape = "treatment", 
          fullrange = TRUE,         
          rug = TRUE
)+ stat_cor(aes(color = treatment), label.x = 3) + ggtitle(paste0(ff1$gene_symbol[1], ' ~ ',  trait1)) + xlab(paste0('Normalized ', ff1$gene_symbol[1], ' (tpm)')) + ylab(paste0(trait1))
g2 = ggscatter(ov11, x = "tpm", y = "cat2",
               color = "treatment", palette = "jco", 
               add = "reg.line",
               shape = "treatment", 
               fullrange = TRUE,         
               rug = TRUE
)+ stat_cor(aes(color = treatment), label.x = 3) + ggtitle(paste0(ff1$gene_symbol[1], ' ~ ',  trait2)) + xlab(paste0('Normalized ', ff1$gene_symbol[1], ' (tpm)')) + ylab(paste0(trait2))
gridExtra::grid.arrange(g1, g2)

#maybe include a pathway enrichment for the residuals?