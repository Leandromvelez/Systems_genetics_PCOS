############################################################################
#                                                                          #
#                        WGCNA Ovary and Fat PCOS model                    #
#                                                                          #
############################################################################




here::i_am("code/WGCNA on cntl or let mice and correlation with traits.R")
load(here::here('data/raw/Sequencing_and_traits_dataBxDpcosAugust2023.RData'))

# load('Sequencing_and_traits_dataBxDpcosAugust2023 - LV.RData')

rm(immune_data, ovarian_morph_counts, secreted_and_steroid_factors, steroid_genes)

library(WGCNA)
library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(reshape2)
library(ggplot2)
library(cowplot)
library(here)

### 0.0 - inspect data ####
colnames(all_traits_BxdPcos)
unique(all_traits_BxdPcos$treatment)
unique(all_traits_BxdPcos$mouse_number)

rownames(cnts_mat_ov) # n 43 for Ovary rnaseq
rownames(cnts_mat_ov_fat)# n 93 for ovary fat rnaseq

# we want to investigate cntrl, let, ov and ov_fat
# let ov
# let ov_fat
# let ov + ov_fat

# cntrl ov
# cntrl ov_fat
# cntrl ov + ov_fat



### 1.0 - Subset for Treatment and Traits ####


let_mice = unique(subset(all_traits_BxdPcos, treatment == 'let')$mouse_number)
ctrl_mice = unique(subset(all_traits_BxdPcos, treatment == 'ctrl')$mouse_number)

# 43 ov
cnts_mat1 = cnts_mat_ov %>%
  rename_all(~paste0(., "_ov")) %>%
  mutate(Var.2_ov = NULL) %>%
  rownames_to_column(var = "mouse") %>%
  arrange(as.integer(mouse)) %>%
  column_to_rownames(var = "mouse")

# 93 ov_fat
cnts_mat2 = cnts_mat_ov_fat %>%
  rename_all(~paste0(., "_ovfat")) %>%
  mutate(Var.2_ovfat = NULL) %>%
  rownames_to_column(var = "mouse") %>%
  arrange(as.integer(mouse)) %>%
  column_to_rownames(var = "mouse")

# thus 43 with both
cnts_mat_all = merge(cnts_mat1, cnts_mat2, by = 'row.names', all = F) %>%
  arrange(as.integer(Row.names)) %>%
  column_to_rownames(var = "Row.names")



let_ov_cor = cnts_mat1[row.names(cnts_mat1) %in% let_mice,] # 21
let_ovfat_cor = cnts_mat2[row.names(cnts_mat2) %in% let_mice,] # 
let_all_cor = cnts_mat_all[row.names(cnts_mat_all) %in% let_mice,] # 20

ctrl_ov_cor = cnts_mat1[row.names(cnts_mat1) %in% ctrl_mice,] # 
ctrl_ovfat_cor = cnts_mat2[row.names(cnts_mat2) %in% ctrl_mice,] # 
ctrl_all_cor= cnts_mat_all[row.names(cnts_mat_all) %in% ctrl_mice,] #


# change this one
datExpr0 = let_ovfat_cor
data = "let ovfat expr"

# remove objects
rm( list = ls()[grep("_cor$", ls())] )
rm( cnts_mat1, cnts_mat2, cnts_mat_all )


select_traits = c(
  "weight-6_minus_weight-0",
  "fat_mass",
  "AUC",
  "Mean_vel_Asc_Aorta",
  "eject_frac",
  "Testosterone",
  "estradiol",
  "LH_FSH_ratio",
  "ovary_area.mm2.",
  "corpus_luteum",
  "antral"
)


# subset mice in select trt group for select traits
traits1 = all_traits_BxdPcos[all_traits_BxdPcos$trait_name %in% select_traits,]
traits1 =  subset(traits1, mouse_number %in% row.names(datExpr0))

write.csv(traits1, here("data/traits1 let.csv"), row.names = FALSE)




### 2.0 - Check for outlier/bad genes and samples ####

## ... in the future we could also examine outliers with PCA

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)

# check if all genes/samples are good
gsg = goodSamplesGenes(datExpr0, verbose = 3)
gsg$allOK

# remove bad genes/samples
if (!gsg$allOK){
  
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste0("Removing genes: ", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", "), "\n"));
  
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste0("Removing samples: ", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", "), "\n"));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

# Cluster samples and determine cuttoff for outliers
sampleTree = hclust(dist(datExpr0), method = "average")
plot(sampleTree, main = paste0(data,": Sample clustering to detect outliers"), sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
abline(h = 500, col = "red")
# ggsave('plots/outliers and TOM/hclust on let ov and ovfat samples.pdf', width = 12, height = 9)

clust = cutreeStatic(sampleTree, cutHeight = 500, minSize = 4)
table(clust)

keepSamples = (clust==1)
datExpr = datExpr0[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
collectGarbage()

# Re-cluster samples
sampleTree2 = hclust(dist(datExpr), method = "average")
plot(sampleTree2, main = paste0(data, ": Sample clustering to detect outliers"), sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2) 



### 3.0 - Soft thresholding and network topology analysis ####

powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

# Plot the results:
sft.data = sft$fitIndices

# Scale-free topology fit index as a function of the soft-thresholding power
g1 = ggplot(sft.data, aes(Power, SFT.R.sq, label = Power) ) +
  geom_text(nudge_y = 0.05) +
  geom_point() +
  geom_hline(yintercept = 0.8, color = "red") +
  labs(
    x="Soft Threshold (power)",y="Scale Free Topology Model Fit,signed R^2",
  ) +
  theme_classic()

# Mean connectivity as a function of the soft-thresholding power
g2 = ggplot(sft.data, aes(Power, mean.k., label = Power) ) +
  geom_text(nudge_y = 100) +
  geom_point() +
  geom_hline(yintercept = 0.0, color = "red") +
  labs(
    x="Soft Threshold (power)", y="Mean Connectivity",
  ) +
  theme_classic()

plot_grid(g1, g2)


file = paste0('plots/outliers and TOM/', data ,' - scale independence and mean connectivity.pdf')
ggsave( here(file), height = 5, width = 9)


### 4.0 - Module Generation and heatmap (Marcus) ####

# **** this is where we fiddle for module generation:
# net colors is n modules
# amongst the traits that leandro is interested in, how do modules cluster at the end?
# we want to end up with a decent number of modules that correlate with the traits
# ... at the end look at which of these Let modules have representation from DEGs between let and cntrl

allowWGCNAThreads()

# set power decided from soft thresholding and min genes per modules
pwr = 5
ngenes = 50

# takes several hours
net = blockwiseModules(datExpr,
                       power = pwr, #determined above
                       maxBlockSize = 25000, # we can set this high because we have a lot of ram. We have 20000 genes so this will do them all at onec
                       TOMType = "unsigned",
                       minModuleSize = ngenes,
                       reassignThreshold = 0.3,
                       mergeCutHeight = 0.3,
                       numericLabels = TRUE, # names modules as numbers isntead of colors
                       pamRespectsDendro = FALSE,
                       saveTOMs = F,
                       nThreads = 48,
                       verbose = 3)

#save(net, file =here(paste0("data/", data, " - pwr12 min50 net blockwisemodules .Rdata")))
save(net, file =here(paste0("data/", data, " - pwr", pwr, " min", ngenes," net blockwisemodules .Rdata")))


# *** can pickup here ***
#load( "let ov expr net blockwisemodules power 7 - min 200.Rdata")

table(net$colors)
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);

all_tog = MEs

# ME0 is error module
all_tog$ME0=NULL
#colnames(all_tog) = paste0('Gene_', colnames(all_tog))

gene_modules = all_tog

module_membership = as.data.frame(net$colors)
colnames(module_membership) = 'module'
module_membership$ID = row.names(module_membership)
head(module_membership)


#Ife we want to look at prportions per module of a given pathway this is an example
#gene_list1 = c(some_vector_of_genes)
#module_membership$data_ID = ifelse(module_membership$ID %in%gene_list1, 'Pathway', 'Other')
#perc_table <- module_membership %>%
#  dplyr::group_by(data_type, module) %>%
#  dplyr::summarise(cnt = n()) %>%
#  dplyr::mutate(freq = round(cnt / sum(cnt), 3)) %>% 
#  dplyr::arrange(desc(freq))

#file = paste0("data/",data, " - pwr12 min50 module membership full list .csv")
file = paste0("data/",data, " - pwr", pwr," min", ngenes, " module membership full list .csv")
write.csv(module_membership, file = here(file), row.names = F)

#pdf(file = 'module proportions by datatype.pdf')
#ggplot(perc_table, aes(x=module, y=freq, fill=data_type)) + geom_col() + theme_minimal() + ggtitle('module proportions by datatype')
#dev.off()

cor_traits = dcast(traits1, mouse_number ~ trait_name, value.var = 'value', fun.aggregate = mean)

row.names(cor_traits) = cor_traits$mouse_number
cor_traits$mouse_number=NULL
row.names(cor_traits)[1:20]
row.names(all_tog)[1:20]

length(intersect(row.names(cor_traits)[1:20], row.names(all_tog)[1:20]))
setdiff(row.names(cor_traits)[1:20], row.names(all_tog)[1:20])
setdiff(row.names(all_tog)[1:20], row.names(cor_traits)[1:20])
# 43 and 177 are not overlapping

#same




all_tog = all_tog[row.names(all_tog) %in% row.names(cor_traits),]
cor_traits = cor_traits[row.names(cor_traits) %in% row.names(all_tog),]
pcor = corAndPvalue(cor_traits, all_tog, use = 'p', method='pearson')

df = melt(pcor$cor)
df$pval = melt(pcor$p)$value
df$obs = melt(pcor$nObs)$value
df$sig = ifelse(df$pval<0.05, '*', '')





traitXmod = dcast(df, Var1 ~ Var2, value.var = 'value', fun.aggregate = mean) %>%
  #traitXmod = traitXmod %>%
  mutate(
    
    Var1 = case_when(
      Var1 == "estradiol" ~ "Serum Estradiol",
      Var1 == "Testosterone" ~ "Serum Testosterone",
      Var1 == "weight-6_minus_weight-0" ~ "Delta Bodyweight",
      Var1 == "Mean_vel_Asc_Aorta"~ "Mean Velocity Asc. Aorta",
      Var1 == "fat_mass"~ "Fat Mass by MRI",
      Var1 == "LH_FSH_ratio"~ "LH:FSH Ratio",
      Var1 == "eject_frac"~ "Ejection Fraction",
      Var1 == "AUC" ~ "Glucose Tolerance AUC"
      
    )
  )


rownames(traitXmod) = traitXmod$Var1
traitXmod$Var1 = NULL

sig_table = dcast(df, Var1 ~ Var2, value.var = 'sig') %>%
  mutate(
    
    Var1 = case_when(
      Var1 == "estradiol" ~ "Serum Estradiol",
      Var1 == "Testosterone" ~ "Serum Testosterone",
      Var1 == "weight-6_minus_weight-0" ~ "Delta Bodyweight",
      Var1 == "Mean_vel_Asc_Aorta"~ "Mean Velocity Asc. Aorta",
      Var1 == "fat_mass"~ "Fat Mass by MRI",
      Var1 == "LH_FSH_ratio"~ "LH:FSH Ratio",
      Var1 == "eject_frac"~ "Ejection Fraction",
      Var1 == "AUC" ~ "Glcucose Tolerance AUC"
      
    )
  )
rownames(sig_table) = sig_table$Var1
sig_table$Var1 = NULL

range(traitXmod)


breaksList = seq(-0.5, 0.5, by = 0.1)
pheatmap(traitXmod,
         fontsize_number = 30,
         display_numbers = sig_table, show_rownames = T, show_colnames = T,
         number_color = "black",
         color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(length(breaksList)),
         breaks = breaksList, main='WGCNA Modules vs Letrozole Traits',
         fontsize_row = 8, fontsize_col = 5,
         cluster_rows = T, cluster_cols = T)

ggsave("plots/WGCNA Modules vs Letrozole Traits.pdf", height = 6, width = 12)


### 6.0 - build vertical version of heatmap (mods on Y axis, traits on X) ####

traitXmod2 = as.data.frame(t(traitXmod)) %>%
  rownames_to_column(var = "module") %>%
  mutate( module = as.numeric(gsub("ME", "", module))) %>%
  column_to_rownames(var = "module" )
traitXmod2 = traitXmod2[order(as.integer(rownames(traitXmod2))),]

sig_table2 = as.data.frame(t(sig_table)) %>%
  rownames_to_column(var = "module") %>%
  mutate( module = as.numeric(gsub("ME", "", module))) %>%
  column_to_rownames(var = "module" )
sig_table2 =sig_table2[order(as.integer(rownames(sig_table2))),]

breaksList = seq(-0.5, 0.5, by = 0.1)

heatmap = pheatmap(traitXmod2,
                   fontsize_number = 20,
                   display_numbers = sig_table2, # note the change here
                   show_rownames = T, show_colnames = T,
                   number_color = "black",
                   color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(length(breaksList)),
                   breaks = breaksList, main='Module Correlations with Letrozole Traits',
                   fontsize_row = 8, fontsize_col = 12,
                   cluster_rows = F, cluster_cols = T,
                   angle_col = 315)




### 5.0 - barplot of Ov vs Gwat tissue representation ####


t = module_membership %>%
  filter( grepl("_ov$", ID) ) %>%
  group_by(module) %>%
  summarise(
    ov_count = n()
  )

t2 = module_membership %>%
  filter( grepl("_ovfat$", ID) ) %>%
  group_by(module) %>%
  summarise(
    ovfat_count = n()
  ) %>%
  left_join(
    t
  ) %>%
  pivot_longer(
    cols = c(ov_count, ovfat_count), values_to = "counts"
  ) %>%
  filter(module != 0)

t2 = t2 %>%
  group_by( module) %>%
  mutate(
    relab = counts/sum(counts)
  )

p1 = ggplot(t2, aes(fill=name, y=relab, x= reorder(module, -module))) + 
  geom_bar(position="stack", stat="identity") +
  geom_text(aes(label = counts), position = position_stack(), size = 3, hjust = 1.5, fontface = "bold") +
  #scale_x_discrete( labels = function(module) as.character(as.integer(module)) ) +
  coord_flip() + scale_y_reverse() +
  theme_bw() +
  labs(title = "Relative Tissue Abundance", x = "Module") +
  scale_fill_manual(values = c("#F1A340","#998EC3"), labels=c('Ovary', 'GWAT')) +
  theme(
    legend.position = "none", axis.title.x = element_blank(),
    plot.title = element_text(hjust = 0.5)
  )




### 6.0 - Barplot of DE genes ####

# ov
de_ov = read.csv(here("data/DE genes/results from limma on let over cntl, ovary BxD PCOS.csv"), fileEncoding = "UTF-8-BOM") %>%
  filter( P.Value <= 0.05 )

mod_ov = module_membership %>%
  filter( grepl("_ov$", ID) ) %>%
  mutate( ID = gsub("_ov$","",ID) )

de_ov = left_join(de_ov, mod_ov) %>%
  filter(is.na(module) == FALSE) %>%
  group_by(module) %>%
  summarise(
    counts_ov = n()
  )

# ovfat
de_ovfat = read.csv(here("data/DE genes/results from limma on let over cntl, periov fat BxD PCOS.csv"), fileEncoding = "UTF-8-BOM") %>%
  filter( P.Value <= 0.05 )

mod_ovfat = module_membership %>%
  filter( grepl("_ovfat$", ID) ) %>%
  mutate( ID = gsub("_ovfat$","",ID) )

de_ovfat = left_join(de_ovfat, mod_ovfat) %>%
  filter(is.na(module) == FALSE) %>%
  group_by(module) %>%
  summarise(
    counts_ovfat = n()
  )


# join and melt
de_mod = merge(de_ov, de_ovfat) %>%
  filter(module != 0) %>%
  pivot_longer( cols = c(counts_ov, counts_ovfat), values_to = "counts" ) %>%
  group_by( module) %>%
  mutate(
    relab = counts/sum(counts)
  )


p2 = ggplot(de_mod, aes(fill = name, y=relab, x= reorder(module, -module))) + 
  geom_bar(position="stack", stat="identity") +
  geom_text(aes(label = counts), position = position_stack(), size = 3, hjust = 1.5, fontface = "bold") +
  #scale_x_discrete( labels = function(module) as.character(as.integer(module)) ) +
  coord_flip() + scale_y_reverse() +
  scale_fill_manual(values = c("#F1A340","#998EC3")) +
  theme_bw() +
  labs(title = "Proportions of DEGs \n Letr. vs Cntrl (p < 0.05)", x = "Module") +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5),
    axis.title.x = element_blank(), axis.title.y = element_blank()
  ) 




### 7.0 - Plot together DEGs, tissue representation and the heatmap for the combined modules ####

# legend
legend = get_legend( p1 + theme(legend.position = "top", legend.title = element_blank()) )

# title
title = ggdraw() + draw_label("Module Compositions", fontface = "bold")

# barplots
gg = plot_grid(p1, p2, align = 'h', axis = 'l')
barplot = plot_grid(  title, legend, gg, ncol = 1, rel_heights = c(.05, .05, 1)  ) + theme(plot.margin = margin( b =90))

# heatmap
heatobj = plot_grid( ggplot()+theme_minimal(), heatmap$gtable, ncol =1, rel_heights = c(0.02, 1) )

# together
plot_grid(barplot, heatobj, rel_widths = c(.9,1) )

ggsave(here("plots/letrezole modules traits and tissue representation3.pdf"), height = 8, width = 10,)











here::i_am("code/WGCNA on cntl or let mice and correlation with traits.R")
load(here::here('data/raw/Sequencing_and_traits_dataBxDpcosAugust2023.RData'))

# load('Sequencing_and_traits_dataBxDpcosAugust2023 - LV.RData')

rm(immune_data, ovarian_morph_counts, secreted_and_steroid_factors, steroid_genes)

library(WGCNA)
library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(reshape2)
library(ggplot2)
library(cowplot)
library(here)

# The appropriate data objects for this analysis were already made in 'WGCNA on cntrl or let mice and correlation with traits.Rdata'
# here we can just load the WGCNA objects for ovary and ovfat


# *** can pickup here ***
load( "../data/let ovary expr - pwr12 min50 net blockwisemodules.Rdata")
net_ov = net
rm(net)
mods_ov = read.csv("../data/let ovary expr - pwr12 min50 module membership full list.csv")


load( "../data/let ovfat expr - pwr5 min50 net blockwisemodules .Rdata")
net_fat = net
rm(net)
mods_fat = read.csv("../data/let ovfat expr - pwr5 min50 module membership full list .csv")

traits1 = read.csv(here("data/traits1 let.csv"), fileEncoding = "UTF-8-BOM")


### 1.0 - Build ovary heatmap objects ####
table(net_ov$colors)
moduleLabels = net_ov$colors
moduleColors = labels2colors(net_ov$colors)
MEs = net_ov$MEs;


all_tog = MEs
# ME0 is error module
all_tog$ME0=NULL


cor_traits = dcast(traits1, mouse_number ~ trait_name, value.var = 'value', fun.aggregate = mean)

row.names(cor_traits) = cor_traits$mouse_number
cor_traits$mouse_number=NULL
row.names(cor_traits)[1:20]
row.names(all_tog)[1:20]

length(intersect(row.names(cor_traits)[1:20], row.names(all_tog)[1:20]))
setdiff(row.names(cor_traits)[1:20], row.names(all_tog)[1:20])
setdiff(row.names(all_tog)[1:20], row.names(cor_traits)[1:20])
# 43 and 177 are not overlapping

all_tog = all_tog[row.names(all_tog) %in% row.names(cor_traits),]
cor_traits = cor_traits[row.names(cor_traits) %in% row.names(all_tog),]
pcor = corAndPvalue(cor_traits, all_tog, use = 'p', method='pearson')

df = melt(pcor$cor)
df$pval = melt(pcor$p)$value
df$obs = melt(pcor$nObs)$value
df$sig = ifelse(df$pval<0.05, '*', '')

traitXmod_ov = dcast(df, Var1 ~ Var2, value.var = 'value', fun.aggregate = mean) %>%
  #traitXmod = traitXmod %>%
  mutate(
    
    Var1 = case_when(
      Var1 == "estradiol" ~ "Serum Estradiol",
      Var1 == "Testosterone" ~ "Serum Testosterone",
      Var1 == "weight-6_minus_weight-0" ~ "Delta Bodyweight",
      Var1 == "Mean_vel_Asc_Aorta"~ "Mean Velocity Asc. Aorta",
      Var1 == "fat_mass"~ "Fat Mass by MRI",
      Var1 == "LH_FSH_ratio"~ "LH:FSH Ratio",
      Var1 == "eject_frac"~ "Ejection Fraction",
      Var1 == "AUC" ~ "Glucose Tolerance AUC"
      
    )
  )


rownames(traitXmod_ov) = traitXmod_ov$Var1
traitXmod_ov$Var1 = NULL

sig_table_ov = dcast(df, Var1 ~ Var2, value.var = 'sig') %>%
  mutate(
    
    Var1 = case_when(
      Var1 == "estradiol" ~ "Serum Estradiol",
      Var1 == "Testosterone" ~ "Serum Testosterone",
      Var1 == "weight-6_minus_weight-0" ~ "Delta Bodyweight",
      Var1 == "Mean_vel_Asc_Aorta"~ "Mean Velocity Asc. Aorta",
      Var1 == "fat_mass"~ "Fat Mass by MRI",
      Var1 == "LH_FSH_ratio"~ "LH:FSH Ratio",
      Var1 == "eject_frac"~ "Ejection Fraction",
      Var1 == "AUC" ~ "Glucose Tolerance AUC"
      
    )
  )
rownames(sig_table_ov) = sig_table_ov$Var1
sig_table_ov$Var1 = NULL

range(traitXmod_ov)

breaksList = seq(-0.5, 0.5, by = 0.1)
pheatmap(traitXmod_ov,
         fontsize_number = 30,
         display_numbers = sig_table_ov,
         show_rownames = T, show_colnames = T,
         number_color = "black",
         color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(length(breaksList)),
         breaks = breaksList, main='WGCNA Modules vs Letrozole Traits',
         fontsize_row = 8, fontsize_col = 5,
         cluster_rows = T, cluster_cols = T)







### 2.0 - Build ovary GWAT heatmap objects ####
table(net_fat$colors)
moduleLabels = net_fat$colors
moduleColors = labels2colors(net_fat$colors)
MEs = net_fat$MEs


all_tog = MEs
# ME0 is error module
all_tog$ME0=NULL


cor_traits = dcast(traits1, mouse_number ~ trait_name, value.var = 'value', fun.aggregate = mean)

row.names(cor_traits) = cor_traits$mouse_number
cor_traits$mouse_number=NULL
row.names(cor_traits)[1:20]
row.names(all_tog)[1:20]

length(intersect(row.names(cor_traits)[1:20], row.names(all_tog)[1:20]))
setdiff(row.names(cor_traits)[1:20], row.names(all_tog)[1:20])
setdiff(row.names(all_tog)[1:20], row.names(cor_traits)[1:20])
# 43 and 177 are not overlapping

all_tog = all_tog[row.names(all_tog) %in% row.names(cor_traits),]
cor_traits = cor_traits[row.names(cor_traits) %in% row.names(all_tog),]
pcor = corAndPvalue(cor_traits, all_tog, use = 'p', method='pearson')

df = melt(pcor$cor)
df$pval = melt(pcor$p)$value
df$obs = melt(pcor$nObs)$value
df$sig = ifelse(df$pval<0.05, '*', '')

traitXmod_fat = dcast(df, Var1 ~ Var2, value.var = 'value', fun.aggregate = mean) %>%
  #traitXmod = traitXmod %>%
  mutate(
    
    Var1 = case_when(
      Var1 == "estradiol" ~ "Serum Estradiol",
      Var1 == "Testosterone" ~ "Serum Testosterone",
      Var1 == "weight-6_minus_weight-0" ~ "Delta Bodyweight",
      Var1 == "Mean_vel_Asc_Aorta"~ "Mean Velocity Asc. Aorta",
      Var1 == "fat_mass"~ "Fat Mass by MRI",
      Var1 == "LH_FSH_ratio"~ "LH:FSH Ratio",
      Var1 == "eject_frac"~ "Ejection Fraction",
      Var1 == "AUC" ~ "Glucose Tolerance AUC"
      
    )
  )


rownames(traitXmod_fat) = traitXmod_fat$Var1
traitXmod_fat$Var1 = NULL

sig_table_fat = dcast(df, Var1 ~ Var2, value.var = 'sig') %>%
  mutate(
    
    Var1 = case_when(
      Var1 == "estradiol" ~ "Serum Estradiol",
      Var1 == "Testosterone" ~ "Serum Testosterone",
      Var1 == "weight-6_minus_weight-0" ~ "Delta Bodyweight",
      Var1 == "Mean_vel_Asc_Aorta"~ "Mean Velocity Asc. Aorta",
      Var1 == "fat_mass"~ "Fat Mass by MRI",
      Var1 == "LH_FSH_ratio"~ "LH:FSH Ratio",
      Var1 == "eject_frac"~ "Ejection Fraction",
      Var1 == "AUC" ~ "Glucose Tolerance AUC"
      
    )
  )

rownames(sig_table_fat) = sig_table_fat$Var1
sig_table_fat$Var1 = NULL

range(traitXmod_fat)

breaksList = seq(-0.5, 0.5, by = 0.1)
pheatmap(traitXmod_fat,
         fontsize_number = 30,
         display_numbers = sig_table_fat,
         show_rownames = T, show_colnames = T,
         number_color = "black",
         color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(length(breaksList)),
         breaks = breaksList, main='WGCNA Modules vs Letrozole Traits',
         fontsize_row = 8, fontsize_col = 5,
         cluster_rows = T, cluster_cols = T)


### 6.0 - build vertical version of heatmaps (mods on Y axis, traits on X) ####


## 6.1 - Ovary
traitXmod_ov2 = as.data.frame(t(traitXmod_ov)) %>%
  rownames_to_column(var = "module") %>%
  mutate( module = as.numeric(gsub("ME", "", module))) %>%
  column_to_rownames(var = "module" )
traitXmod_ov2 = traitXmod_ov2[order(as.integer(rownames(traitXmod_ov2))),]

sig_table_ov2 = as.data.frame(t(sig_table_ov)) %>%
  rownames_to_column(var = "module") %>%
  mutate( module = as.numeric(gsub("ME", "", module))) %>%
  column_to_rownames(var = "module" )
sig_table_ov2 =sig_table_ov2[order(as.integer(rownames(sig_table_ov2))),]

breaksList = seq(-0.5, 0.5, by = 0.1)

heatmap_ov = pheatmap(
  traitXmod_ov2,
  fontsize_number = 20,
  display_numbers = sig_table_ov2, # note the change here
  show_rownames = T, show_colnames = T,
  number_color = "black",
  color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(length(breaksList)),
  breaks = breaksList, main='Module Correlations with Letrozole Traits',
  fontsize_row = 8, fontsize_col = 12,
  cluster_rows = F, cluster_cols = T,
  angle_col = 315
)


## 6.2 - Ovary fat
traitXmod_fat2 = as.data.frame(t(traitXmod_fat)) %>%
  rownames_to_column(var = "module") %>%
  mutate( module = as.numeric(gsub("ME", "", module))) %>%
  column_to_rownames(var = "module" )
traitXmod_fat2 = traitXmod_fat2[order(as.integer(rownames(traitXmod_fat2))),]

sig_table_fat2 = as.data.frame(t(sig_table_fat)) %>%
  rownames_to_column(var = "module") %>%
  mutate( module = as.numeric(gsub("ME", "", module))) %>%
  column_to_rownames(var = "module" )
sig_table_fat2 =sig_table_fat2[order(as.integer(rownames(sig_table_fat2))),]

breaksList = seq(-0.5, 0.5, by = 0.1)

heatmap_gwat = pheatmap(traitXmod_fat2,
                        fontsize_number = 20,
                        display_numbers = sig_table_fat2, # note the change here
                        show_rownames = T, show_colnames = T,
                        number_color = "black",
                        color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(length(breaksList)),
                        breaks = breaksList, main='Module Correlations with Letrozole Traits',
                        fontsize_row = 8, fontsize_col = 12,
                        cluster_rows = F, cluster_cols = T,
                        angle_col = 315)



plot_grid(heatmap_ov$gtable, heatmap_gwat$gtable)


### 7.0 Plot teh heatmaps with module on x and trait on y; bind the ov and fat matrices and cluster by traits ####

t_ov =  traitXmod_ov %>% rename_all(~ paste0(., "_OV"))
t_fat =  traitXmod_fat %>% rename_all(~ paste0(., "_GWAT"))
s_ov = sig_table_ov %>% rename_all(~ paste0(., "_OV"))
s_fat = sig_table_fat %>% rename_all(~ paste0(., "_GWAT"))



t = cbind(t_ov, t_fat)
numeric_part <- as.numeric(gsub("ME([0-9]+)_(OV|GWAT)", "\\1", names(t)))
OV_GWAT_part <- gsub("ME([0-9]+)_(OV|GWAT)", "\\2", names(t))
sorted_columns <- names(t)[order(OV_GWAT_part, numeric_part)]
t <- t[, sorted_columns]

s =  cbind(s_ov, s_fat)
numeric_part <- as.numeric(gsub("ME([0-9]+)_(OV|GWAT)", "\\1", names(s)))
OV_GWAT_part <- gsub("ME([0-9]+)_(OV|GWAT)", "\\2", names(s))
sorted_columns <- names(s)[order(OV_GWAT_part, numeric_part)]
s <- s[, sorted_columns]

heatmap_join = pheatmap(t,
                        fontsize_number = 20,
                        display_numbers = s, # note the change here
                        show_rownames = T, show_colnames = T,
                        number_color = "black",
                        color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(length(breaksList)),
                        breaks = breaksList,
                        main='Module Correlations with Letrozole Traits',
                        fontsize_row = 8, fontsize_col = 12,
                        cluster_rows = T, cluster_cols = F,
                        gaps_col = 15,
                        angle_col = 315)


ggsave("../data/plots/ov and gwat wgcna with traits.pdf", height = 9, width = 12)

### 8.0 - Evaluate module compositions ####

library(org.Hs.eg.db) # bioconductor # or mouse, Org.Mm.eg.db
library(clusterProfiler) #bioconductor package
library(enrichplot)

# enrichment for Ovary wgcna modules
pdf(here("plots/ovary wgcna clusterprofiler and enrichr.pdf"), width = 10, height = 8)
for(m in c("1", "3", "4", "6", "8", "12") ){
  
  #m = "1"
  modov = mods_ov %>%
    mutate(ID = toupper(gsub("_ov$", "", ID))) %>%
    filter(module == m)
  
  pp1 = modov$ID
  
  organism = "org.Hs.eg.db"
  go = enrichGO(pp1, OrgDb = organism, pvalueCutoff = 1,ont = "BP", qvalueCutoff = 1, pAdjustMethod = "none", keyType = "SYMBOL")
  go = pairwise_termsim(go)
  p1 = enrichplot::dotplot(go,showCategory = 20, label_format = 80) +
    ggtitle(paste0("Letrezole Ovary WGCNA Module ", m))
  
  p2 = emapplot(go,showCategory = 20, cluster.params = list(label_format = 80)) +
    ggtitle(paste0("Letrezole Ovary WGCNA Module ", m))
  
  
  print(p1)
  print(p2)
  
}
dev.off()


# enrichment for Ovary FAT wgcna modules
pdf(here("plots/ov fat wgcna clusterprofiler and enrichr.pdf"), width = 10, height = 8)
for(m in c("1", "3", "5", "8", "12", "15") ){
  
  #m = "1"
  modfat = mods_fat %>%
    mutate(ID = toupper(gsub("_ovfat$", "", ID))) %>%
    filter(module == m)
  
  pp1 = modfat$ID
  
  organism = "org.Hs.eg.db"
  go = enrichGO(pp1, OrgDb = organism, pvalueCutoff = 1,ont = "BP", qvalueCutoff = 1, pAdjustMethod = "none", keyType = "SYMBOL")
  go = pairwise_termsim(go)
  p1 = enrichplot::dotplot(go,showCategory = 20, label_format = 80) +
    ggtitle(paste0("Letrezole Ov Fat WGCNA Module ", m))
  
  p2 = emapplot(go,showCategory = 20, cluster.params = list(label_format = 80)) +
    ggtitle(paste0("Letrezole Ov Fat WGCNA Module ", m))
  
  
  print(p1)
  print(p2)
  
}
dev.off()





## z - old

pp1 = modov$ID

pp1 = modfat$ID


#pp1 = as.character(shared_genes$human_orth) # if needed, cut number of genes at 500
organism = "org.Hs.eg.db"
go = enrichGO(pp1, OrgDb = organism, pvalueCutoff = 1,ont = "BP", qvalueCutoff = 1, pAdjustMethod = "none", keyType = "SYMBOL")
go = pairwise_termsim(go)
enrichplot::dotplot(go,showCategory = 20, label_format = 80) + ggtitle('')
emapplot(go,showCategory = 20, cluster.params = list(label_format = 80)) + ggtitle('ORA Shared genes Mouse gwat BxD PCOS - Human Adipose obese vs ctrls GSE203371')































setwd("G:/My Drive/lab files/PCOS screen/12-4-23 working figure analyses/WGCNA on PCOS and trait inspections")
library(enrichplot)
library(enrichR)
library(WGCNA)
library(qvalue)
load('G:/My Drive/lab files/PCOS screen/12-4-23 working figure analyses/Sequencing_and_traits_dataBxDpcosAugust2023.RData')

mod_meblership = read.csv('G:/My Drive/lab files/PCOS screen/12-4-23 working figure analyses/WGCNA on PCOS and trait inspections/data/let both expr - pwr7 min200 module membership full list.csv')
head(mod_meblership)
mod_meblership$tissue =gsub(".*_","",mod_meblership$ID)
mod_meblership$gene_symbol =gsub("\\_.*","",mod_meblership$ID)

dbs <- listEnrichrDbs()
dbs = dbs[order(dbs$geneCoverage, decreasing = T),]
setEnrichrSite("Enrichr")
dbs1 <- c("GO_Biological_Process_2021", "HDSigDB_Mouse_2021", "Reactome_2022", "Enrichr_Users_Contributed_Lists_2020")

run_tissue_paths = function(tissue2){
  deg_set2 = mod_meblership[mod_meblership$tissue %in% tissue2 & mod_meblership$module==18,]
  gg1 = deg_set2$gene_symbol
  enriched <- enrichr(gg1, dbs1)
  g1 = plotEnrich(enriched[[1]], showTerms = 10, numChar = 40, y = "Ratio", orderBy = "P.value") + ggtitle(names(enriched[1])) + theme(axis.text=element_text(size=5),
                                                                                                                                       axis.title=element_text(size=6,face="bold"))
  g2 = plotEnrich(enriched[[2]], showTerms = 10, numChar = 40, y = "Ratio", orderBy = "P.value") + ggtitle(names(enriched[2]))+ theme(axis.text=element_text(size=5),
                                                                                                                                      axis.title=element_text(size=6,face="bold"))
  
  g3 = plotEnrich(enriched[[3]], showTerms = 10, numChar = 40, y = "Ratio", orderBy = "P.value") + ggtitle(names(enriched[3]))+ theme(axis.text=element_text(size=5),
                                                                                                                                      axis.title=element_text(size=6,face="bold"))
  
  g4 = plotEnrich(enriched[[4]], showTerms = 10, numChar = 40, y = "Ratio", orderBy = "P.value") + ggtitle(names(enriched[4]))+ theme(axis.text=element_text(size=5),
                                                                                                                                      axis.title=element_text(size=6,face="bold"))
  
  pdf(file = paste0( tissue2, ' Top pathways based on normalized Ssec.pdf'))
  gg2 = gridExtra::grid.arrange(g1, g2, g3, g4)
  print(gg2)
  dev.off()
}

run_tissue_paths('ov')
run_tissue_paths('ovfat')

#looks like tnf signaling in OV and TCA cycle in OVFAT

table(mod_meblership$tissue)
full_mod = deg_set2 = mod_meblership[ mod_meblership$module==18,]
#correlate genes and traits
colnames(all_traits_BxdPcos)
table(all_traits_BxdPcos$trait_name)
ms_use = all_traits_BxdPcos[all_traits_BxdPcos$treatment=='let',]
use_traits = c('Testosterone', 'eject_frac', 'estradiol')
ms_use = ms_use[ms_use$trait_name %in% use_traits,]

trt_table = dcast(ms_use, mouse_number ~ trait_name, value.var = 'value')
row.names(trt_table) = trt_table$mouse_number
trt_table$mouse_number=NULL
table(full_mod$tissue)
ov_seq = cnts_mat_ov
colnames(ov_seq) = paste0(colnames(ov_seq), '_ov')
ovfat_seq = cnts_mat_ov_fat
colnames(ovfat_seq) = paste0(colnames(ovfat_seq), '_ovfat')
ov_seq = ov_seq[,colnames(ov_seq) %in% full_mod$ID]
ovfat_seq = ovfat_seq[,colnames(ovfat_seq) %in% full_mod$ID]
l = list(ov_seq, ovfat_seq, trt_table)
combined_metrt = Reduce(merge, lapply(l, function(x) data.frame(x, rn = row.names(x))))
row.names(combined_metrt) = combined_metrt$rn
combined_metrt$rn=NULL


cc1 = bicorAndPvalue(combined_metrt, combined_metrt, use = 'p')
mm1 = reshape2::melt(cc1$bicor)
colnames(mm1) = c('ID1', 'ID2', 'bicor')
mm1$pvalue = reshape2::melt(cc1$p)$value
head(mm1)
mm1 = mm1[order(abs(mm1$bicor), decreasing = T),]
est1 = mm1[mm1$ID1=='estradiol',]
head(est1)

est1 = est1[1:100,]
full_set = c(as.vector(est1$ID2), 'Testosterone', 'eject_frac')

paths_table = read.delim('G:/My Drive/Datasets/Mouse/genome files/uniprot-mouse-with goterms.tab')
sec_prots = read.delim('G:/My Drive/Datasets/Mouse/genome files/secreted_proteins.txt')


net_genes = mm1[mm1$ID1 %in% full_set,]

net_genes = net_genes[net_genes$ID2 %in% full_set,]
net_genes1 = net_genes[!net_genes$ID1==net_genes$ID2,]
cc4 = net_genes %>% group_by(ID1) %>% summarise(meanp = mean(-log(pvalue)), meanabsbic = mean(abs(bicor)))
cc4 = cc4[order(cc4$meanabsbic, decreasing = T),]
head(cc4)


net_genes$gene_symbol1 =gsub("\\_.*","",net_genes$ID1)
head(net_genes)
library(forcats)
top_tissue = net_genes[net_genes$gene_symbol1 %in% sec_prots$Gene.names...primary..,]
nn1 = top_tissue
nn1$meanbics = cc4$meanabsbic[match(nn1$ID1, cc4$ID1)]
nn1$tissue1 = gsub(".*_","",nn1$ID1)
#nn1$gene_symbol =gsub("\\_.*","",nn1$ID1)
head(nn1)

nn1 = nn1[!nn1$ID1==nn1$ID2,]
pdf(file = paste0( timept, 'top-ranked secreted candidates mean connectivity.pdf'))
table(as.vector(nn1$gene_tissue1))
library(forcats)

nn1 = net_genes
nn1$bicor[nn1$bicor>0.9999] = 0
map2 = reshape2::dcast(nn1, ID1 ~ ID2, value.var = 'bicor')
row.names(map2) = map2$ID1
map2$ID1=NULL
cols_set = as.data.frame(row.names(map2))
row.names(cols_set) = row.names(map2)
cols_set$cols = ifelse(grepl('_ov', cols_set$`row.names(map2)`), 'darkorange3', 'gray5')
cols_set$cols = ifelse(grepl('_ovfat', cols_set$`row.names(map2)`), 'darkorchid2', paste0(cols_set$cols))

head(map2)
map2 = as.data.frame(map2)
cols_set$labs = ifelse('estradiol|Testosterone|eject_frac' %in% cols_set$`row.names(map2)`, paste0(cols_set$`row.names(map2)`), '')
library(qgraph)
pdf(file = 'undirected network - With labels.pdf')
qgraph(map2, minimum = 0.1, cut = 0.4, vsize = 4, color=cols_set$cols, legend = F, borders = F, layout='spring', posCol = "dodgerblue3", negCol = "firebrick3", label.cex=0,  directed=F, labels = gsub("\\_.*","",row.names(map2))) 
dev.off()

head(gene_set3)
























######################################################################################
#                                                                                    #
#    Integration with GTEX data and correlation analysis through relevant tissues    #
#                                                                                    #
######################################################################################









load('G:/My Drive/lab files/PCOS screen/12-4-23 working figure analyses/Sequencing_and_traits_dataBxDpcosAugust2023.RData')
setwd('G:/My Drive/lab files/PCOS screen/12-4-23 working figure analyses/WGCNA on PCOS and trait inspections/')
library(dplyr)
library(ggplot2)
library(ggupset)
library(WGCNA)
library(reshape2)
pathways_table = read.delim('G:/My Drive/Datasets/Human/genome files/uniprot-human-genes and goterms mapping.tab')

load('G:/My Drive/Datasets/Human/GTEx NA included env.RData')


origin_gene_set = c('ESR1', 'GPER', 'AR', 'TNF', 'IL1B', 'CPT1A', 'ACAT')
origin_tissues= c('Ovary', 'Adipose - Subcutaneous', 'Adipose - Visceral (Omentum)')

origin_set1 = apply(expand.grid(origin_gene_set, origin_tissues), 1, paste, collapse="_")


target_tissues = c('Thyroid', 'Uterus', 'Stomach', 'Spleen', 'Small Intestine - Terminal Ileum', 'Pituitary', 'Pancreas', 'Muscle - Skeletal', 'Heart - Left Ventricle', 'Brain - Hypothalamus', 'Brain - Hippocampus',  'Adrenal Gland')

path_term = 'inflam'
run_pathcors = function(path_term){
gg1 = pathways_table[grepl(path_term, pathways_table$Gene.ontology..biological.process.) |  grepl(path_term, pathways_table$Gene.ontology..molecular.function.),]

path_gene_set = unique(gg1$Gene.names...primary..)
target_set = apply(expand.grid(path_gene_set, target_tissues), 1, paste, collapse="_")
origin1 = GTEx_full[,colnames(GTEx_full) %in% origin_set1]
target1 = GTEx_full[,colnames(GTEx_full) %in% target_set]

head(origin1)
new_origin = origin1[complete.cases(origin1), ]
new_target = target1[row.names(target1) %in% row.names(new_origin),]

cc1 = bicorAndPvalue(new_origin, new_target, use = 'p')
cc2 = melt(cc1$bicor)
colnames(cc2) = c('ori_genetissue', 'targ_genetissue', 'bicor')
cc2$pvallue = melt(cc1$p)$value
cc2$ori_genesymbol = gsub("\\_.*","",cc2$ori_genetissue)
cc2$targ_genesymbol = gsub("\\_.*","",cc2$targ_genetissue)
cc2$ori_tissue = gsub(".*_","",cc2$ori_genetissue)
cc2$targ_tissue = gsub(".*_","",cc2$targ_genetissue)
cc2$pathterm = paste0(path_term)
sum_table = cc2 %>% dplyr::group_by(ori_genetissue, targ_tissue) %>% dplyr::summarise(logp = mean(-log10(pvallue), na.rm = T))
sum_table$pathterm = paste0(path_term)
return(cc2)
}


plots_tables = as.data.frame(rbind(run_pathcors('inflam'), run_pathcors('oxidation'), run_pathcors('stress')))
table(plots_tables$ori_genetissue)
plots_tables$logp = -log10(plots_tables$pvallue)

c434 = 'TNF_Ovary'
g1 = plots_tables[plots_tables$ori_genetissue==c434,]
g1 = g1[g1$logp>1.3,]
g1 = na.omit(g1)

gg1 = ggplot(g1, aes(x=targ_tissue,  fill=pathterm)) + theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + geom_bar(stat="count", position = 'dodge') + ggtitle(paste0(c434)) + scale_fill_manual(values = c('forestgreen', 'dodgerblue4', 'darkorchid4'))


c434 = 'IL1B_Ovary'
g1 = plots_tables[plots_tables$ori_genetissue==c434,]
g1 = g1[g1$logp>1.3,]
g1 = na.omit(g1)

gg2 = ggplot(g1, aes(x=targ_tissue,  fill=pathterm)) + theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + geom_bar(stat="count", position = 'dodge') + ggtitle(paste0(c434)) + scale_fill_manual(values = c('forestgreen', 'dodgerblue4', 'darkorchid4'))

c434 = 'ESR1_Ovary'
g1 = plots_tables[plots_tables$ori_genetissue==c434,]
g1 = g1[g1$logp>1.3,]
g1 = na.omit(g1)

gg3 = ggplot(g1, aes(x=targ_tissue,  fill=pathterm)) + theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + geom_bar(stat="count", position = 'dodge') + ggtitle(paste0(c434)) + scale_fill_manual(values = c('forestgreen', 'dodgerblue4', 'darkorchid4'))

c434 = 'AR_Ovary'
g1 = plots_tables[plots_tables$ori_genetissue==c434,]
g1 = g1[g1$logp>1.3,]
g1 = na.omit(g1)

gg4 = ggplot(g1, aes(x=targ_tissue,  fill=pathterm)) + theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + geom_bar(stat="count", position = 'dodge') + ggtitle(paste0(c434)) + scale_fill_manual(values = c('forestgreen', 'dodgerblue4', 'darkorchid4'))

#+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
pdf(file ='pathway focus pantissue.pdf')
gridExtra::grid.arrange(gg1, gg2, gg3,  gg4)
dev.off()
