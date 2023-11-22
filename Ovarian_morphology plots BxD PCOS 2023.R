setwd('D:/My Drive/PostDoc UCI/New Project PCOS Leandro UCI 2021/BxD PCOS Project 2021/Results/2022_2023_Data table, R environment and paper in preparation 2023/temporal Leandro')
library(reshape2)
library(ggplot2)
library(Rmisc)
library(rstatix)
library(dplyr)

### Follicular and ovarian morphology

## Dont run lines 13 to 29, these are to show the formula I applied to estimate total number of follicles, cysts, and corp luteums

# follicle_counts = read.csv('G:/My Drive/PostDoc UCI/New Project PCOS Leandro UCI 2021/BxD PCOS Project 2021/Results/Histology/BXD PCOS Ovary analysis 2023_csv.csv')  # full melted counts unfiltered 
# Here I applied a formula to estimate a total number of each follicular, corpus luteum, and cysts population in the whole ovary. 
# It is applied to all the samples, so its not going to change the comparisons, but it will give a more expected final number.

# follicles_melt = melt(follicle_counts, id = c("mouse_number", "slide_number")) 
# follicles_melt$value=as.numeric(follicles_melt$value)
# follicles_melt = follicles_melt[!is.na(follicles_melt$value),]
# follicles_melt$value = ifelse(follicles_melt$variable=='primordial', (follicles_melt$value + 10)*50, follicles_melt$value)
# follicles_melt$value = ifelse(follicles_melt$variable=='primary.secon', (follicles_melt$value + 5)*10, follicles_melt$value)
# follicles_melt$value = ifelse(follicles_melt$variable=='antral', (follicles_melt$value + 3)*10, follicles_melt$value)
# follicles_melt$value = ifelse(follicles_melt$variable=='atretic_antral', (follicles_melt$value + 3)*10, follicles_melt$value)
# follicles_melt$value = ifelse(follicles_melt$variable=='corpus_luteum', (follicles_melt$value + 1)*10, follicles_melt$value)
# follicles_melt$value = ifelse(follicles_melt$variable=='cysts', (follicles_melt$value)*3, follicles_melt$value)
# follicles_melt1 = follicles_melt %>% dplyr::group_by(mouse_number, variable) %>% dplyr::summarise(mean_value=mean(value, na.rm=T), sd=sd(value, na.rm=T), sem=sd(value, na.rm=T)/n())                      

# follicles_melt1$treatment = all_traits_BxdPcos$treatment[match(follicles_melt1$mouse_number, all_traits_BxdPcos$mouse_number)]
# ovarian_morph_counts = follicles_melt1  



tests <- summarySE(ovarian_morph_counts, measurevar="mean_value", groupvars=c("treatment", "variable"))
tests$treatment <- factor(tests$treatment, levels = c("ctrl", "let"))
ggplot(tests, aes(x=treatment, y=mean_value, fill= treatment, group = variable)) + 
  scale_fill_manual(values = c( 'dodgerblue3','darkgoldenrod2')) +
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=mean_value-se, ymax= mean_value +se), width=.2,
                position=position_dodge(.9)) + labs(title="", x="variable", y = "Estimated number") + facet_wrap(~ variable, scales="free")


stat.test <- ovarian_morph_counts %>%
  group_by(variable) %>%
  t_test(mean_value ~ treatment) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()

# write.csv(stat.test, file = 'G:/My Drive/PostDoc UCI/New Project PCOS Leandro UCI 2021/BxD PCOS Project 2021/Results/Histology/stats_follicleCounts_bxdPCOS.csv', row.names = F)
# write.csv(tests, file = 'G:/My Drive/PostDoc UCI/New Project PCOS Leandro UCI 2021/BxD PCOS Project 2021/Results/Histology/mean_sd_se_follicleCounts_bxdPCOS.csv', row.names = F)



## combine in a matrix with all the traits

ovarian_morph_counts_matrix = dcast(ovarian_morph_counts, mouse_number ~ variable, value.var = 'mean_value', fun.aggregate = mean)

## I order the traits (except for now the ovarian morph) by category
traits_cat=unique(all_traits_BxdPcos$trait_category)

trait_matrix = all_traits_BxdPcos
trait_matrix$trait_category <- factor(trait_matrix$trait_category, levels= traits_cat, ordered=TRUE)

trait_matrix <- trait_matrix[order(trait_matrix$trait_category, trait_matrix$trait_id),]
 # trait_matrix = trait_matrix[!trait_matrix$trait_category=='Age',]
traits_name = unique(trait_matrix$trait_name)
## matrix
trait_matrix= dcast(trait_matrix, mouse_number ~ trait_name, value.var = 'value', fun.aggregate = mean)
#reorder as all_traits
trait_matrix = trait_matrix[,traits_name]

trait_matrix$mouse_number = as.numeric( rownames(trait_matrix))

# row.names(trait_matrix) = trait_matrix$mouse_number
#trait_matrix$mouse_number=NULL




## Merge the two dataframes
trait_matrix_all = merge(trait_matrix, ovarian_morph_counts_matrix, by="mouse_number", all = T)
trait_matrix_all = trait_matrix_all[order(trait_matrix_all$mouse_number),] 
trait_matrix_all[trait_matrix_all=='NaN']  = NA
rownames(trait_matrix_all) = trait_matrix_all$mouse_number
trait_matrix_all$mouse_number = NULL

library(WGCNA)
cc1 = bicorAndPvalue(trait_matrix_all[,colnames(trait_matrix_all) %in% traits_name],
                     trait_matrix_all[,!colnames(trait_matrix_all) %in% traits_name & !colnames(trait_matrix_all)=='mouse_number'], 
                     use = 'p')
cc2 = cc1$bicor
tt4 = cc1$p
tt3= ifelse(tt4 < 0.05,"*","")

trait_matrix_all$treatment = sample_info$treatment[match(rownames(trait_matrix_all),sample_info$mouse_number)]
trait_matrix_all$treatment
library(pheatmap)

pheatmap(cc2, fontsize_number = 20, display_numbers = tt3, 
         number_color = "black", 
         color = colorRampPalette(c("black", "darkred","red","orange","gold","white","lightcyan","powderblue","skyblue","lightskyblue","dodgerblue"))(30), 
         cluster_rows = F, cluster_cols = T, na_col="white", main='Trait x ovarian morph traits corrs in BxD PCOS MICE', fontsize_row = 10, fontsize_col = 10)

#### Corr scatter plots and boxplots, selected traits
plot_corrs_ovarmorphTotraits = function(ov_trait,trait1 ){
  g1 = ggscatter(trait_matrix_all, 
                 x = paste0(ov_trait), y = paste0(trait1),  
                 color = "treatment", palette = c('dodgerblue3', 'darkgoldenrod2'),
                 xlab=paste0(ov_trait),
                 ylab = paste0(trait1),
                 point = T,
                 ellipse = F, mean.point = F, star.plot = F, 
                 # add = "reg.line",  # Add regressin line
                 add.params = list(color = "black", fill = "lightgray"), # Customize reg. line
                 conf.int = TRUE) + # stat_cor(method = "pearson") + 
    ggtitle(paste0(ov_trait," x ", trait1,"_ctrl BxD PCOS"))
  
g2 = ggboxplot(trait_matrix_all, x = "treatment", y = paste0(ov_trait),
               color = "treatment", palette = c('dodgerblue3', 'darkgoldenrod2'), outlier.shape = NA,
               add = "jitter") + stat_compare_means(aes(label = after_stat(p.format)), 
                                                    method = "wilcox.test", label.x.npc = "center", label.y.npc = "center" ) + 
  labs(title= paste0(ov_trait, '_let vs ctrl, BxD PCOS'), x="treatment") 

  
  g3= ggboxplot(trait_matrix_all, x = "treatment", y = paste0(trait1),
                color = "treatment", palette = c('dodgerblue3', 'darkgoldenrod2'), outlier.shape = NA,
                add = "jitter") + stat_compare_means(aes(label = after_stat(p.format)), 
                                                     method = "wilcox.test", label.x.npc = "center", label.y.npc = "center" ) + 
    labs(title= paste0(trait1, '_let vs ctrl, BxD PCOS'), x="treatment") 
  
#  g4= ggboxplot(cnts_mat_all_ov_ovfat, x = "treatment", y = paste0(trait1),
#                color = "treatment", palette = c('dodgerblue3', 'darkgoldenrod2'), outlier.shape = NA,
#                add = "jitter") + stat_compare_means(aes(label = after_stat(p.format)), 
#                                                     method = "wilcox.test", label.x.npc = "center", label.y.npc = "center" ) + 
    labs(title= paste0(trait1, '_let vs ctrl, BxD PCOS'), x="treatment") 
  
  g5 = cowplot::plot_grid(g2, g3, ncol=2) 
 # g6 = cowplot::plot_grid(g3, g4, ncol=2) 
  title <- cowplot::ggdraw() + cowplot::draw_label(paste0(ov_trait, ' x ', paste0(trait1), '_BxD PCOS 2023'), fontface='bold')
# pdf(file = paste0("Corr ",ov_trait,'_x_', trait1," PCOS-like vs Ctrl", ".pdf"))
cowplot::plot_grid(title, g1, g5, ncol=1, rel_heights=c(0.1, 1, 1)) # rel_heights values control title margins for each element (title, g5, g6)
#  print(g7)
#  dev.off()
}
# need to check which genes to plot, in "corr ovtofat ctrl" or "corr ovtofat let", and next line will plot both genes corr as pdf
library(ggpubr)
plot_corrs_ovarmorphTotraits('primordial','glucose_90mins')  # check in the working directory

