setwd('E:/My Drive/PostDoc UCI/New Project PCOS Leandro UCI 2021/BxD PCOS Project 2021/Results/2022_2023_Data table, R environment and paper in preparation 2023/temporal Leandro')
load('E:/My Drive/PostDoc UCI/New Project PCOS Leandro UCI 2021/BxD PCOS Project 2021/Results/2022_2023_Data table, R environment and paper in preparation 2023/Sequencing_and_traits_dataBxDpcosAugust2023.RData')

all_traits_BxdPcos$strain = ifelse(all_traits_BxdPcos$strain == '129', paste0("	129X1/SvJ"), all_traits_BxdPcos$strain)


# BiocManager::install("GOSemSim")
library(stringr)
library(GOSemSim)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(forcats)
library(reshape2)
library(RColorBrewer)
library(ggpubr)
library(impute)
library(preprocessCore)
library(WGCNA)
# BiocManager::install("org.Mm.eg.db")
library(org.Mm.eg.db)
# BiocManager::install("clusterProfiler")
library(clusterProfiler)
library(pheatmap)
library(colormap)
library(cowplot)






### DEG genes ovary BxD PCOS, in shared folder
DEG_ovary = read.csv('E:/My Drive/PostDoc UCI/New Project PCOS Leandro UCI 2021/BxD PCOS Project 2021/Results/2022_2023_Data table, R environment and paper in preparation 2023/results from limma on let over cntl, ovary BxD PCOS.csv')
DEG_ov_fat = read.csv('E:/My Drive/PostDoc UCI/New Project PCOS Leandro UCI 2021/BxD PCOS Project 2021/Results/2022_2023_Data table, R environment and paper in preparation 2023/results from limma on let over cntl, periov fat BxD PCOS.csv')


### traits: load provided environment

## All traits x all traits Bi Corr Heatmaps

## traits matrices, let and ctrl

traits_cat=unique(all_traits_BxdPcos$trait_category)

trait_matrix = all_traits_BxdPcos
trait_matrix$trait_category <- factor(trait_matrix$trait_category, levels= traits_cat, ordered=TRUE)

trait_matrix <- trait_matrix[order(trait_matrix$trait_category, trait_matrix$trait_id),]
trait_matrix = trait_matrix[!trait_matrix$trait_category=='Age',]
trait_matrix$trait_name = str_replace_all(trait_matrix$trait_name, "/", "_")    # Applying str_replace_all
trait_matrix$trait_name = str_replace_all(trait_matrix$trait_name, "%", "Perc")
trait_matrix$trait_name = str_replace_all(trait_matrix$trait_name, " ", "")    # Applying str_replace_all


traits_name = unique(trait_matrix$trait_name)


trait_matrix= dcast(trait_matrix, mouse_number ~ trait_name, value.var = 'value', fun.aggregate = mean)
row.names(trait_matrix) = trait_matrix$mouse_number
trait_matrix$mouse_number=NULL
#reorder as all_traits

trait_matrix = trait_matrix[,traits_name]




####genes matrices

## ovary genes matrixes
# cnts_mat_ov)


## periovarian fat genes matrixes

# cnts_mat_ov_fat


### Corr Analysis, endocrine communication ovary - periovarian fat (origin-target)(worthwhile to try the inverse)
## ovary x periov fat
## for gene x gene corrs, I use the "secreted and steroid related" group of genes for the origin tissue. It is the list we use of secreted factors in
# mouse, and additionally I incorporated all steroid related factors from diff databases (GO paths lists, etc)

# ctrls
select_ctrl = sample_info$mouse_number[sample_info$treatment=='ctrl']
# let
select_let = sample_info$mouse_number[sample_info$treatment=='let']

# ctrl samples in ovary RNA seq
select_ctrl_ovarygenes = rownames(cnts_mat_ov)[rownames(cnts_mat_ov) %in% select_ctrl]
# let samples in ovary RNA seq
select_let_ovarygenes = rownames(cnts_mat_ov)[rownames(cnts_mat_ov) %in% select_let]
##
# ctrl samples in periov RNA seq
select_ctrl_ovfatgenes = rownames(cnts_mat_ov_fat)[rownames(cnts_mat_ov_fat) %in% select_ctrl]
# let samples in ovary RNA seq
select_let_ovfatgenes = rownames(cnts_mat_ov_fat)[rownames(cnts_mat_ov_fat) %in% select_let]
## ovary genes x periovarian fat genes
# as origin (ovary), I use the list of secreted + steroid related
# ctrls
cc1 = bicorAndPvalue(cnts_mat_ov[rownames(cnts_mat_ov) %in% select_ctrl_ovarygenes,
      colnames(cnts_mat_ov) %in% secreted_and_steroid_factors$secreted], 
      cnts_mat_ov_fat[rownames(cnts_mat_ov_fat) %in% select_ctrl_ovarygenes, ], 
      use = 'p')

cc2 = cc1$bicor
tt4 = cc1$p

set1 = melt(cc1$bicor)
head(set1)
colnames(set1) = c('gene_ov', 'gene_ovfat', 'bicor')
set2 = melt(tt4)
set1$pvalue = set2$value
set2=NULL
set1$gene_ov_gene_ovfat = paste0(set1$gene_ov,"_", set1$gene_ovfat)
### Ssec scores for ctrl, ovary to periovarian fat
## to know about the Ssec scores, see Marcus paper https://www.sciencedirect.com/science/article/pii/S1550413118301931?via%3Dihub


Ssec_ov_to_fat= set1[set1$pvalue<0.99999999 & !is.na(set1$gene_ov),]
Ssec_ov_to_fat = Ssec_ov_to_fat[!is.na(Ssec_ov_to_fat$pvalue),]
scores = Ssec_ov_to_fat %>%   dplyr::group_by(gene_ov) %>%  dplyr::summarise(Score = sum(-log(pvalue), na.rm=TRUE))
scores1 = dplyr::count(Ssec_ov_to_fat, gene_ov)
Ssec_ov_to_fat = NULL
scores$lenght = scores1$n[match(scores$gene_ov, scores1$gene_ov)]
scores1= NULL
scores$Ssec = scores$Score / scores$lenght
scores = scores[!scores$Score==Inf,]
scores = scores[order(-scores$Ssec),]

Ssec_ovTofat_ctrl = scores
corr_ovTofat_ctrl = set1[set1$pvalue<0.05 & set1$bicor<0.9999999,]
test1 = na.omit(corr_ovTofat_ctrl[corr_ovTofat_ctrl$gene_ov=='Angpt4',])

corr_ovTofat_ctrl = corr_ovTofat_ctrl[order(corr_ovTofat_ctrl$pvalue),]
corr_ovTofat_ctrl$qvalue = signif(p.adjust(corr_ovTofat_ctrl$pvalue, "BH"), 3)
corr_ovTofat_ctrl = corr_ovTofat_ctrl[!is.na(corr_ovTofat_ctrl$gene_ov),]
set1 = NULL
scores = NULL
cc1 = NULL
cc2 = NULL
tt4 = NULL
gc()




## let

cc1 = bicorAndPvalue(cnts_mat_ov[rownames(cnts_mat_ov) %in% select_let_ovarygenes & !rownames(cnts_mat_ov)=='68',
                     colnames(cnts_mat_ov) %in% secreted_and_steroid_factors$secreted], 
                     cnts_mat_ov_fat[rownames(cnts_mat_ov_fat) %in% select_let_ovarygenes, ], 
                     use = 'p')

cc2 = cc1$bicor
tt4 = cc1$p

set1 = melt(cc1$bicor)
head(set1)
colnames(set1) = c('gene_ov', 'gene_ovfat', 'bicor')
set2 = melt(tt4)
set1$pvalue = set2$value
set2=NULL
set1$gene_ov_gene_ovfat = paste0(set1$gene_ov,"_", set1$gene_ovfat)
### Ssec scores for ctrl, ovary to periovarian fat
## to know about the Ssec scores, see Marcus paper https://www.sciencedirect.com/science/article/pii/S1550413118301931?via%3Dihub

Ssec_ov_to_fat= set1[set1$pvalue<0.99999999 & !is.na(set1$gene_ov),]
Ssec_ov_to_fat = Ssec_ov_to_fat[!is.na(Ssec_ov_to_fat$pvalue),]
scores = Ssec_ov_to_fat %>%   dplyr::group_by(gene_ov) %>%  dplyr::summarise(Score = sum(-log(pvalue), na.rm=TRUE))
scores1 = dplyr::count(Ssec_ov_to_fat, gene_ov)
Ssec_ov_to_fat = NULL
scores$lenght = scores1$n[match(scores$gene_ov, scores1$gene_ov)]
scores1= NULL
scores$Ssec = scores$Score / scores$lenght
scores = scores[!scores$Score==Inf,]
scores = scores[order(-scores$Ssec),]

Ssec_ovTofat_let = scores
corr_ovTofat_let = set1[set1$pvalue<0.05 & set1$bicor<0.9999999,]
corr_ovTofat_let = corr_ovTofat_let[order(corr_ovTofat_let$pvalue),]
corr_ovTofat_let$qvalue = signif(p.adjust(corr_ovTofat_let$pvalue, "BH"), 3)
corr_ovTofat_let = corr_ovTofat_let[!is.na(corr_ovTofat_let$gene_ov),]
## It is interesting to see the highest Ssec scores in ctrls and lets, and how different are top Ssec genes in ctrl and let
# the pattern of ovarian signaling to ovary fat is different from ctrls and let (pcos-like)
## locate, i.e., Hsd17b3 (synthesis of potent androgens) and Hsd17b2 (described as antiandrogenic/antiestrogenic)
#
set1 = NULL
scores = NULL
cc1 = NULL
cc2 = NULL
tt4 = NULL
gc()



## histogram distributions of Ssec scores
hist(Ssec_ovTofat_ctrl$Ssec, breaks=120, ylim = c(0, 150),main="Ssec scores distribution secreted ovarian factors, ctrl, ovary x periov fat")
hist(Ssec_ovTofat_let$Ssec, breaks=120, ylim = c(0, 150),main="Ssec scores distribution secreted ovarian factors, let, ovary x periov fat")
mean(Ssec_ovTofat_ctrl$Ssec)
mean(Ssec_ovTofat_let$Ssec)


## paste DE info, noting if the ovarian gene is differentially expressed in the origin tissue

Ssec_ovTofat_ctrl$DE_ovary = DEG_ovary$P.Value[match(Ssec_ovTofat_ctrl$gene_ov, DEG_ovary$ID)]
Ssec_ovTofat_ctrl$DE_ovary = ifelse(Ssec_ovTofat_ctrl$DE_ovary<0.005, paste0('**'),
                              ifelse(Ssec_ovTofat_ctrl$DE_ovary<0.05, paste0('*'),
                                     paste0('NS'))      
)

Ssec_ovTofat_let$DE_ovary = DEG_ovary$P.Value[match(Ssec_ovTofat_let$gene_ov, DEG_ovary$ID)]
Ssec_ovTofat_let$DE_ovary = ifelse(Ssec_ovTofat_let$DE_ovary<0.005, paste0('**'),
                                      ifelse(Ssec_ovTofat_let$DE_ovary<0.05, paste0('*'),
                                             paste0('NS'))      
)

library(forcats)
Ssec_ovTofat_ctrl[1:10,] %>% mutate(gene = fct_reorder(gene_ov, Ssec))  %>%
  ggplot(  aes(x= gene, y= Ssec, fill= DE_ovary)) + geom_col() + 
  scale_fill_manual(values = c( 'lightblue','dodgerblue3')) +
  theme(axis.text.x = element_text(angle = 90)) + coord_flip() + 
  ggtitle('Top 20 Ssec ovary to fat, controls BxD PCOS')


Ssec_ovTofat_let[1:10,] %>% mutate(gene = fct_reorder(gene_ov, Ssec))  %>%
  ggplot(  aes(x= gene, y= Ssec, fill= DE_ovary)) + geom_col() + 
  scale_fill_manual(values = c( 'lightgoldenrod','dodgerblue3','darkgoldenrod2')) +
  theme(axis.text.x = element_text(angle = 90)) + coord_flip() + 
  ggtitle('Top 20 Ssec ovary to fat, pcos-m BxD PCOS')




### pheat map with top Ssec which are also top DEG in ovary, filtering by secreted or steroid related

## Ssec scores in let and ctrls, the DEGs that also are secreted_steroid related
Ssec_matrix = as.data.frame(Ssec_ovTofat_ctrl[Ssec_ovTofat_ctrl$gene_ov %in% DEG_ovary$ID[DEG_ovary$P.Value<0.05],])
# rownames(Ssec_matrix) = Ssec_matrix$gene_ov
# Ssec_matrix$gene_ov = NULL
Ssec_matrix$let = Ssec_ovTofat_let$Ssec[match(Ssec_matrix$gene_ov, Ssec_ovTofat_let$gene_ov)]
# Ssec_matrix = Ssec_matrix[,c('Ssec','let')]
# colnames(Ssec_matrix) = c('ctrl','let')
Ssec_matrix$pvalue = DEG_ovary$P.Value[match(Ssec_matrix$gene_ov, DEG_ovary$ID )]
Ssec_matrix$delta_Ssec = Ssec_matrix$let - Ssec_matrix$Ssec
Ssec_matrix = Ssec_matrix[order(abs(Ssec_matrix$delta_Ssec),decreasing=T),]




ggplot(Ssec_matrix[1:20,], aes(x=fct_reorder(gene_ov, abs(delta_Ssec), .desc=T), y=delta_Ssec, fill=gene_ov)) + 
  geom_col(position = position_dodge2(reverse = TRUE)) +  theme_classic() + theme(legend.position = "none") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))  + xlab('') +  
  ggtitle(paste0('Top DEG ovarian and delta Ssec (ovary to periovarian fat)')) + ylab('delta Ssec (PCOS Ssec - cntl Ssec')


 ## some interesting papers, this one related to Col6a5 (the biggest diff in Ssec in ctrl vs let) in a model of pcos, relating ovary to adipose
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8181728/

## if interested, look for Sult1e1, Fshr, Hsd17b3, Angpt4, Ptgds, Lhcgr,Fbln7 which are included in this heatmap
# I'm following with more interest to Sult1e1, Hsd17b3, because those are markers of leydig cells, which are the male counter-part
# of granulosa cells (females), and leydig cells sinthetize potent androgens (through hsd17b3) and inactivates estrogen (through Sult1e1)
# Here, Sult1e1 scores is almost the same in both conditions, but this gene is the number one DE gene in ovary (increased in the PCOS condition) and also in periov fat.


## what about signaling from periovarian fat to ovary???? (all we did so far, bur reversed in order)

# testing specific ovarian genes sign to ov fat in let
test1 = na.omit(corr_ovTofat_let[corr_ovTofat_let$gene_ov=='Col6a5',])
test1 = test1[order(test1$pvalue,decreasing=F),]

### trait x genes corrs


## ovary genes x all traits
## again, I filter the secreted and steroid related genes

##ctrl

cc1 = bicorAndPvalue(cnts_mat_ov[rownames(cnts_mat_ov) %in% select_ctrl_ovarygenes, ], 
      trait_matrix[rownames(trait_matrix) %in% select_ctrl_ovarygenes,], 
      use = 'p')

cc2 = cc1$bicor 
tt4 = cc1$p

set1 = melt(cc1$bicor)
head(set1)
colnames(set1) = c('gene_ov', 'trait', 'bicor')
set2 = melt(tt4)
set1$pvalue = set2$value
set2=NULL
set1$gene_ov_trait = paste0(set1$gene_ov,"_", set1$trait)

## I keep only signf corrs
corr_ovTotraits_ctrl = set1[set1$pvalue<0.05 & set1$bicor<0.9999999,]
set1= NULL
## I add the number of genes correlating with each trait

## It might be interesting to know how many ovarian genes correlates significantly with each trait
corr_ovTotraits_ctrl = corr_ovTotraits_ctrl %>%
  group_by(trait) %>%
  mutate(genes_count = n())

genes_trait_count_ctrl = corr_ovTotraits_ctrl[!duplicated(corr_ovTotraits_ctrl$trait), c('trait', 'genes_count')]
genes_trait_count_ctrl = genes_trait_count_ctrl[order(genes_trait_count_ctrl$genes_count,decreasing=T),]
cc1 = NULL
cc2 = NULL
tt4 = NULL
gc()
##let

cc1 = bicorAndPvalue(cnts_mat_ov[rownames(cnts_mat_ov) %in% select_let_ovarygenes,], 
     trait_matrix[rownames(trait_matrix) %in% select_let_ovarygenes,], 
     use = 'p')
cc2 = cc1$bicor
tt4 = cc1$p

set1 = melt(cc1$bicor)
head(set1)
colnames(set1) = c('gene_ov', 'trait', 'bicor')
set2 = melt(tt4)
set1$pvalue = set2$value
set2=NULL
set1$gene_ov_trait = paste0(set1$gene_ov,"_", set1$trait)

## I keep only signf corrs
corr_ovTotraits_let = set1[set1$pvalue<0.05 & set1$bicor<0.9999999,]
set1= NULL
cc1 = NULL
cc2 = NULL
tt4 = NULL
gc()
## I add the number of genes correlating with each trait

## It might be interesting to know how many ovarian genes correlates significantly with each trait
corr_ovTotraits_let = corr_ovTotraits_let %>%
  group_by(trait) %>%
  mutate(genes_count = n())

genes_trait_count_let = corr_ovTotraits_let[!duplicated(corr_ovTotraits_let$trait), c('trait', 'genes_count')]
genes_trait_count_let = genes_trait_count_let[order(genes_trait_count_let$genes_count,decreasing=T),]


### Compare number of genes and traits in ctrls and let!
genes_trait_count_ov = as.data.frame(genes_trait_count_ctrl)
genes_trait_count_ov$let = genes_trait_count_let$genes_count[match(genes_trait_count_ov$trait, genes_trait_count_let$trait)]
genes_trait_count_let = NULL
genes_trait_count_ctrl = NULL
genes_trait_count_ov = as.data.frame(na.omit(genes_trait_count_ov))
rownames(genes_trait_count_ov) = genes_trait_count_ov$trait
genes_trait_count_ov$trait = NULL
colnames(genes_trait_count_ov) = c('ctrl','let')


##get delta gene counts

genes_trait_count_ov$delta = genes_trait_count_ov$let - genes_trait_count_ov$ctrl

genes_trait_count_ov = genes_trait_count_ov[order(abs(genes_trait_count_ov$delta),decreasing=T),]

filtered_list_traits = unique(all_traits_BxdPcos$trait_name)
filtered_list_traits = c('LH_FSH_ratio', 'FSH','','estradiol','cyclicity', 'LH','Testosterone','E2_T','Insulin', 'AUC','glucose_0mins','Perc_final_increase_BW', 
                         'weight-6_minus_weight-0','weight-6', 'eject_frac','stroke_vol','lv_mass_corr','heart_rate', 'cardiac_output',
                         'fat_mass','lean_fat_ratio','lean_mass','Mean_vel_Asc_Aorta','Mean_vel_Desc_Aorta','peakVel_Asc_Desc'  )

genes_trait_count_ov = genes_trait_count_ov[rownames(genes_trait_count_ov) %in% filtered_list_traits,]
genes_trait_count_ov$trait = rownames(genes_trait_count_ov)
ggplot(genes_trait_count_ov, aes(x=fct_reorder(trait, abs(delta), .desc=T), y=delta, fill=trait)) + 
  geom_col(position = position_dodge2(reverse = TRUE)) +  theme_classic() + theme(legend.position = "none") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))  + xlab('') +  
  ggtitle(paste0('Selected traits and delta ovarian gene count correlations (pcos-m - ctrl) p<0.05')) + ylab('delta gene count (PCOS count - cntl count)')

test1 = na.omit(corr_ovTotraits_let[corr_ovTotraits_let$trait=='Perc_final_increase_BW',])
test1 = test1[order(test1$pvalue,decreasing=F),]

## periovarian fat genes x all traits
## again, I filter the secreted and steroid related genes

##ctrl

cc1 = bicorAndPvalue(cnts_mat_ov_fat[rownames(cnts_mat_ov_fat) %in% select_ctrl_ovfatgenes,], 
     trait_matrix[rownames(trait_matrix) %in% select_ctrl_ovfatgenes,], 
       use = 'p')

cc2 = cc1$bicor
tt4 = cc1$p

set1 = melt(cc1$bicor)
head(set1)
colnames(set1) = c('gene_ovFat', 'trait', 'bicor')
set2 = melt(tt4)
set1$pvalue = set2$value
set2=NULL
set1$gene_ovFat_trait = paste0(set1$gene_ovFat,"_", set1$trait)

## I keep only signf corrs
corr_ovfatTotraits_ctrl = set1[set1$pvalue<0.05 & set1$bicor<0.9999999,]
set1= NULL
cc1 = NULL
cc2 = NULL
tt4 = NULL
gc()
## I add the number of genes correlating with each trait

## It might be interesting to know how many periov fat genes correlates significantly with each trait
corr_ovfatTotraits_ctrl = corr_ovfatTotraits_ctrl %>%
  group_by(trait) %>%
  mutate(genes_count = n())

genes_trait_count_ctrl = corr_ovfatTotraits_ctrl[!duplicated(corr_ovfatTotraits_ctrl$trait), c('trait', 'genes_count')]
genes_trait_count_ctrl = genes_trait_count_ctrl[order(genes_trait_count_ctrl$genes_count,decreasing=T),]
genes_trait_count_ctrl = na.omit(genes_trait_count_ctrl)

##let

cc1 = bicorAndPvalue(cnts_mat_ov_fat[rownames(cnts_mat_ov_fat) %in% select_let_ovfatgenes,], 
                     trait_matrix[rownames(trait_matrix) %in% select_let_ovfatgenes & !rownames(trait_matrix)=='68',], 
                     use = 'p')

cc2 = cc1$bicor
tt4 = cc1$p

set1 = melt(cc1$bicor)
head(set1)
colnames(set1) = c('gene_ovFat', 'trait', 'bicor')
set2 = melt(tt4)
set1$pvalue = set2$value
set2=NULL
set1$gene_ovFat_trait = paste0(set1$gene_ovFat,"_", set1$trait)


cc1 = NULL
cc2 = NULL
tt4 = NULL

## I keep only signf corrs
corr_ovfatTotraits_let = set1[set1$pvalue<0.05 & set1$bicor<0.9999999,]
set1= NULL
gc()
## I add the number of genes correlating with each trait

## It might be interesting to know how many ovarian genes correlates significantly with each trait
corr_ovfatTotraits_let = corr_ovfatTotraits_let %>%
  group_by(trait) %>%
  mutate(genes_count = n())

genes_trait_count_let = corr_ovfatTotraits_let[!duplicated(corr_ovfatTotraits_let$trait), c('trait', 'genes_count')]
genes_trait_count_let = genes_trait_count_let[order(genes_trait_count_let$genes_count,decreasing=T),]

### Compare number of genes and traits in ctrls and let!
genes_trait_count_ovfat = genes_trait_count_ctrl
genes_trait_count_ovfat$let = genes_trait_count_let$genes_count[match(genes_trait_count_ovfat$trait, genes_trait_count_let$trait)]
genes_trait_count_let = NULL
genes_trait_count_ctrl = NULL
genes_trait_count_ovfat = as.data.frame(na.omit(genes_trait_count_ovfat))
rownames(genes_trait_count_ovfat) = genes_trait_count_ovfat$trait
genes_trait_count_ovfat$trait = NULL
colnames(genes_trait_count_ovfat) = c('ctrl','let')
genes_trait_count_ovfat$delta = genes_trait_count_ovfat$let - genes_trait_count_ovfat$ctrl
genes_trait_count_ovfat = genes_trait_count_ovfat[order(genes_trait_count_ovfat$delta,decreasing=T),]

genes_trait_count_ovfat = genes_trait_count_ovfat[rownames(genes_trait_count_ovfat) %in% filtered_list_traits,]


genes_trait_count_ovfat$trait = rownames(genes_trait_count_ovfat)
ggplot(genes_trait_count_ovfat, aes(x=fct_reorder(trait, abs(delta), .desc=T), y=delta, fill=trait)) + 
  geom_col(position = position_dodge2(reverse = TRUE)) +  theme_classic() + theme(legend.position = "none") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))  + xlab('') +  
  ggtitle(paste0('Selected traits and delta periov fat gene count correlations (pcos-m - ctrl) p<0.05')) + ylab('delta gene count (PCOS count - cntl count)')

test1 = na.omit(corr_ovfatTotraits_let[corr_ovfatTotraits_let$trait=='fat_mass',])
test1 = test1[order(test1$pvalue,decreasing=F),]

### Once I identified some interesting corrs between ovary x traits, ov fat x traits, I can plot that corr differentiating in ctrls and let
# and plotting the let vs ctrl for the gene
## 



### bind ovarian genes matrix and periovarian genes matrix
#### ggscatter plots

cnts_mat_all_ov_ovfat =  cnts_mat_ov #editing
colnames(cnts_mat_all_ov_ovfat) = paste0(colnames(cnts_mat_all_ov_ovfat),'_ovary')

cnts_mat_ov_fat1 = cnts_mat_ov_fat
colnames(cnts_mat_ov_fat1) = paste0(colnames(cnts_mat_ov_fat1),'_ovfat')
cnts_mat_ov_fat1 = cnts_mat_ov_fat1[rownames(cnts_mat_ov_fat1) %in% rownames(cnts_mat_all_ov_ovfat),]
cnts_mat_all_ov_ovfat = cbind(cnts_mat_all_ov_ovfat[!rownames(cnts_mat_all_ov_ovfat)=='68',], cnts_mat_ov_fat1)
cnts_mat_ov_fat1 = NULL

cnts_mat_all_ov_ovfat$dm = sample_info$treatment[match(rownames(cnts_mat_all_ov_ovfat), sample_info$mouse_number)]
table(cnts_mat_all_ov_ovfat$Cyp3a13_ovary, cnts_mat_all_ov_ovfat$dm)
cnts_mat_all_ov_ovfat$dm <- factor(cnts_mat_all_ov_ovfat$dm, levels=c('ctrl', 'let'))

# run this function which will make corr plots and let vs ctrls, when ovarian genes and periov genes are provided

# control and pcos like conditions
test1 = na.omit(corr_ovTofat_let[corr_ovTofat_let$gene_ov=='Cyp17a1',])

## run function, then pick ovarian and periov fat genes

plot_corrs_ovarytofat = function(ov_gene,ovfat_gene ){
g1 = ggscatter(cnts_mat_all_ov_ovfat[ cnts_mat_all_ov_ovfat$dm=='ctrl',], 
               x = paste0(ov_gene,'_ovary'), y = paste0(ovfat_gene,'_ovfat'),  color = 'dodgerblue3',
               xlab=paste0(ov_gene,'_ovary'),
               ylab = paste0(ovfat_gene,'_ovfat'),
               ellipse = F, mean.point = TRUE,
               # add = "reg.line",  # Add regressin line
               add.params = list(color = "black", fill = "lightgray"), # Customize reg. line
               conf.int = TRUE) + # stat_cor(method = "pearson") + 
               ggtitle(paste0(ov_gene,'_ovary'," x ", ovfat_gene,'_ovfat',"_ctrl BxD PCOS"))

g2 = ggscatter(cnts_mat_all_ov_ovfat[ cnts_mat_all_ov_ovfat$dm=='let',], 
               x = paste0(ov_gene,'_ovary'), y = paste0(ovfat_gene,'_ovfat'),  color = "darkgoldenrod2", # shape = "dm",
               # palette = c('dodgerblue3', 'darkgoldenrod2'),
               xlab= paste0(ov_gene,'_ovary'),
               ylab = paste0(ovfat_gene,'_ovfat'),
               ellipse = F, mean.point = TRUE,
               # add = "reg.line",  # Add regressin line
               add.params = list(color = "black", fill = "lightgray"), # Customize reg. line
               conf.int = TRUE) + # stat_cor(method = "pearson") + 
               ggtitle(paste0(ov_gene,'_ovary'," x ", ovfat_gene,'_ovfat',"_let BxD PCOS"))                                                                    

g3= ggboxplot(cnts_mat_all_ov_ovfat, x = "dm", y = paste0(ov_gene,'_ovary'),
               color = "dm", palette = c('dodgerblue3', 'darkgoldenrod2'), outlier.shape = NA,
               add = "jitter") + stat_compare_means(aes(label = after_stat(p.format)), 
               method = "wilcox.test", label.x.npc = "center", label.y.npc = "center" ) + 
               labs(title= paste0(ov_gene,'_ovary', '_let vs ctrl, BxD PCOS'), x="treatment") 

g4= ggboxplot(cnts_mat_all_ov_ovfat, x = "dm", y = paste0(ovfat_gene,'_ovfat'),
               color = "dm", palette = c('dodgerblue3', 'darkgoldenrod2'), outlier.shape = NA,
               add = "jitter") + stat_compare_means(aes(label = after_stat(p.format)), 
               method = "wilcox.test", label.x.npc = "center", label.y.npc = "center" ) + 
               labs(title= paste0(ovfat_gene,'_ovfat', '_let vs ctrl, BxD PCOS'), x="treatment") 

g5 = cowplot::plot_grid(g1, g2, ncol=2) 
g6 = cowplot::plot_grid(g3, g4, ncol=2) 
title <- cowplot::ggdraw() + cowplot::draw_label(paste0(ov_gene,'_ovary', ' x ', paste0(ovfat_gene,'_ovfat'), '_BxD PCOS 2023'), fontface='bold')
pdf(file = paste0("Corr ",ov_gene,'_ovary_x_', ovfat_gene,'_ovfat'," PCOS-like vs Ctrl", ".pdf"))
g7 = cowplot::plot_grid(title, g5, g6, ncol=1, rel_heights=c(0.1, 1, 1)) # rel_heights values control title margins for each element (title, g5, g6)
print(g7)
dev.off()
}
# need to check which genes to plot, in "corr ovtofat ctrl" or "corr ovtofat let", and next line will plot both genes corr as pdf

plot_corrs_ovarytofat('Hsd17b3','Ar') # check in the working directory

## I can add the bicor and p to the scatter plots after, since ggscatter does not offer the biweight method we apply with WGCNA package


## trait x gene scatter plots and pcos vs ctrl plot



#all ovarian genes x traits matrix

trait_matrix$dm = sample_info$treatment[match(rownames(trait_matrix), sample_info$mouse_number)]
trait_matrix$dm = factor(trait_matrix$dm, levels=c('ctrl', 'let'))

matrix_all_ovary_trait = cbind(cnts_mat_ov,trait_matrix[rownames(trait_matrix) %in% rownames(cnts_mat_ov),])
matrix_all_ovary_trait$dm = factor(matrix_all_ovary_trait$dm, levels=c('ctrl', 'let'))

cnts_mat_ov$dm = sample_info$treatment[match(rownames(cnts_mat_ov), sample_info$mouse_number)]
cnts_mat_ov$dm = factor(cnts_mat_ov$dm, levels=c('ctrl', 'let'))
# cnts_mat_ov$dm = factor(cnts_mat_ov$dm, levels=c('ctrl', 'let'))

 
plot_corrs_ovarytotrait = function(ov_gene,trait1 ){
g1=  ggscatter(matrix_all_ovary_trait[ rownames(matrix_all_ovary_trait) %in% select_ctrl_ovarygenes,], 
               x = paste0(ov_gene), y = paste0(trait1),  color = 'dodgerblue3',
               xlab=paste0(ov_gene),
               ylab = paste0(trait1),
               ellipse = F, mean.point = TRUE,
               # add = "reg.line",  # Add regressin line
               # add.params = list(linetype = "dm") ,
               # palette = c('dodgerblue3', 'darkgoldenrod2'),
               # add.params = list(color = "black", fill = "lightgray"), # Customize reg. line
               conf.int = TRUE) + # + stat_cor(method = "pearson") 
  ggtitle(paste0(ov_gene," x ", trait1,"_BxD PCOS"))

 g2 = ggscatter(matrix_all_ovary_trait[ rownames(matrix_all_ovary_trait) %in% select_let_ovarygenes,], 
               x = paste0(ov_gene), y = paste0(trait1),  color = "darkgoldenrod2", # shape = "dm",
              # palette = c('dodgerblue3', 'darkgoldenrod2'),
               xlab= paste0(ov_gene),
               ylab = paste0(trait1),
              ellipse = F, mean.point = TRUE,
              # add = "reg.line",  # Add regressin line
              # add.params = list(color = "black", fill = "lightgray"), # Customize reg. line
              conf.int = TRUE) + # stat_cor(method = "pearson") + ggtitle(paste0(ov_gene," x ", trait1,"_let BxD PCOS"))                                                                    
   ggtitle(paste0(ov_gene," x ", trait1,"_BxD PCOS"))
g3 = ggboxplot(cnts_mat_ov, x = "dm", y = paste0(ov_gene),
               color = "dm", palette = c('dodgerblue3', 'darkgoldenrod2'), outlier.shape = NA,
               add = "jitter") + stat_compare_means(aes(label = after_stat(p.format)), method = "t.test", label.x.npc = "left", label.y.npc = "top" ) + 
  labs(title= paste0('let vs ctrl, BxD PCOS'), x="treatment") 


g4 = ggboxplot(trait_matrix, x = "dm", y = paste0(trait1),
               color = "dm", palette = c('dodgerblue3', 'darkgoldenrod2'), outlier.shape = NA,
               add = "jitter") + stat_compare_means(aes(label = after_stat(p.format)), method = "t.test", label.x.npc = "left", label.y.npc = "top" ) + 
  labs(title= paste0('let vs ctrl, BxD PCOS'), x="treatment") 

g4 = plot_grid(g1, g2, g3, g4, ncol=2) 
# g5 = plot_grid(g4, g3, ncol=1) 
title <- ggdraw() + draw_label(paste0('Ovarian_',ov_gene, ' x ', paste0(trait1), '_BxD PCOS 2023'), fontface='bold')
pdf(file = paste0("Correlation ",ov_gene, ' x ', trait1," and PCOS-like vs Ctrl", ".pdf"))
g5 = plot_grid(title, g4, ncol=1, rel_heights=c(0.1, 1)) # rel_heights values control title margins
print(g5)
dev.off()
}

plot_corrs_ovarytotrait("Ptgds","glucose_0mins")

### all periovarian genes x traits

matrix_all_ovfat_trait = cbind(cnts_mat_ov_fat,trait_matrix[rownames(trait_matrix) %in% rownames(cnts_mat_ov_fat),])
matrix_all_ovfat_trait$dm = factor(matrix_all_ovfat_trait$dm, levels=c('ctrl', 'let'))

cnts_mat_ov_fat$dm = sample_info$treatment[match(rownames(cnts_mat_ov_fat), sample_info$mouse_number)]
cnts_mat_ov_fat$dm = factor(cnts_mat_ov_fat$dm, levels=c('ctrl', 'let'))
# cnts_mat_ov_fat$dm = factor(cnts_mat_ov_fat$dm, levels=c('ctrl', 'let'))


plot_corrs_ovfattotrait = function(ov_gene,trait1 ){
  g1=  ggscatter(matrix_all_ovfat_trait[rownames(matrix_all_ovfat_trait) %in% select_ctrl_ovfatgenes,], 
                 x = paste0(ov_gene), y = paste0(trait1),  color = 'dodgerblue3',
                 xlab=paste0(ov_gene),
                 ylab = paste0(trait1),
                 ellipse = F, mean.point = TRUE,
                 # add = "reg.line",  # Add regressin line
                 # add.params = list(linetype = "dm") ,
                 # palette = c('dodgerblue3', 'darkgoldenrod2'),
                 # add.params = list(color = "black", fill = "lightgray"), # Customize reg. line
                 conf.int = TRUE) + # + stat_cor(method = "pearson") 
    ggtitle(paste0(ov_gene," x ", trait1,"_BxD PCOS"))
  
  g2 = ggscatter(matrix_all_ovfat_trait[rownames(matrix_all_ovfat_trait) %in% select_let_ovfatgenes,], 
                 x = paste0(ov_gene), y = paste0(trait1),  color = "darkgoldenrod2", # shape = "dm",
                 # palette = c('dodgerblue3', 'darkgoldenrod2'),
                 xlab= paste0(ov_gene),
                 ylab = paste0(trait1),
                 ellipse = F, mean.point = TRUE,
                 # add = "reg.line",  # Add regressin line
                 # add.params = list(color = "black", fill = "lightgray"), # Customize reg. line
                 conf.int = TRUE) + # stat_cor(method = "pearson") + ggtitle(paste0(ov_gene," x ", trait1,"_let BxD PCOS"))                                                                    
    ggtitle(paste0(ov_gene," x ", trait1,"_BxD PCOS"))
  g3 = ggboxplot(cnts_mat_ov_fat, x = "dm", y = paste0(ov_gene),
                 color = "dm", palette = c('dodgerblue3', 'darkgoldenrod2'), outlier.shape = NA,
                 add = "jitter") + stat_compare_means(aes(label = after_stat(p.format)), method = "t.test", label.x.npc = "left", label.y.npc = "top" ) + 
    labs(title= paste0('let vs ctrl, BxD PCOS'), x="treatment") 
  
  
  g4 = ggboxplot(trait_matrix, x = "dm", y = paste0(trait1),
                 color = "dm", palette = c('dodgerblue3', 'darkgoldenrod2'), outlier.shape = NA,
                 add = "jitter") + stat_compare_means(aes(label = after_stat(p.format)), method = "t.test", label.x.npc = "left", label.y.npc = "top" ) + 
    labs(title= paste0('let vs ctrl, BxD PCOS'), x="treatment") 
  
  g4 = plot_grid(g1, g2, g3, g4, ncol=2) 
  # g5 = plot_grid(g4, g3, ncol=1) 
  title <- ggdraw() + draw_label(paste0('Ovfat_',ov_gene, ' x ', paste0(trait1), '_BxD PCOS 2023'), fontface='bold')
  pdf(file = paste0("Correlation ",ov_gene, ' x ', trait1," and PCOS-like vs Ctrl", ".pdf"))
  g5 = plot_grid(title, g4, ncol=1, rel_heights=c(0.1, 1)) # rel_heights values control title margins
  print(g5)
  dev.off()
}

plot_corrs_ovfattotrait("Jak2","lean_fat_ratio")

