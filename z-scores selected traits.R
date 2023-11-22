setwd('D:/My Drive/PostDoc UCI/New Project PCOS Leandro UCI 2021/BxD PCOS Project 2021/Results/2022_2023_Data table and  R environment/temporal Leandro')
library(dplyr)
library(reshape2)
library(pheatmap)
library(scales)
#nstall.packages('RColorBrewer')
library(RColorBrewer)
full_traits = all_traits_BxdPcos
colnames(full_traits)
table(full_traits$treatment)
#ctrl  let 
#4961 5181 
table(full_traits$trait_name)
# traits_to_plot = c('Perc_final_increase_BW', 'lean_mass',  'fat_mass', 'lean/fat_ratio','stroke_vol',  'eject_frac', 'lv_mass','heart_rate','Mean_vel_Desc_Aorta','glucose_0mins','glucose_90mins','AUC','Insulin','estradiol','Testosterone','LH_FSH_ratio','last_cycle_stage')
#2nd# traits_to_plot = c('Perc_final_increase_BW', 'lean/fat_ratio', 'eject_frac', 'lv_mass','Mean_vel_Desc_Aorta','glucose_0mins', 'AUC','Insulin','estradiol','Testosterone','LH_FSH_ratio')
traits_to_plot = c('Perc_final_increase_BW', 'lean/fat_ratio', 'eject_frac', 'lv_mass','Mean_vel_Desc_Aorta','glucose_0mins', 'AUC','Insulin','estradiol','Testosterone','LH_FSH_ratio')

full_traits$value = as.numeric(as.character(full_traits$value))
full_traits = full_traits[!full_traits$strain=='BXD124*',]
filtered_traits = full_traits[full_traits$trait_name %in% traits_to_plot,]
filtered_traits = filtered_traits[!is.na(filtered_traits$value),]
filtered_traits$scaled_value  = ave(filtered_traits$value, filtered_traits$trait_name, FUN=scale)

#summary_z_score = filtered_traits %>% dplyr::select(trait_name, #strain, treatment, value) %>% 
# dplyr::group_by(trait_name, strain, treatment) %>% 
#  dplyr::mutate(z_score = scale(value))

summary_z_score = as.data.frame(filtered_traits)
## exclude the strains without lets or ctrls
summary_z_score = summary_z_score %>% 
  group_by(trait_name, strain) %>% 
  filter(length(unique(treatment)) == 2) 

head(summary_z_score)
summary_z_score$strain_treatment = paste0(summary_z_score$strain, '_', summary_z_score$treatment)
mheat = reshape2::dcast(summary_z_score,  trait_name~ strain_treatment, value.var = 'scaled_value', fun.aggregate = mean, na.rm=T)
## specifi order of traits to plot
mheat = mheat %>% arrange(factor(trait_name, levels = c('Perc_final_increase_BW', 'lean/fat_ratio', 'eject_frac', 'lv_mass','Mean_vel_Desc_Aorta','glucose_0mins', 'AUC','Insulin','estradiol','Testosterone','LH_FSH_ratio')))

row.names(mheat) = mheat$trait_name
mheat$trait_name=NULL
m2 = as.matrix(mheat)
m2[is.na(m2)] = 0
range(m2)

rg <- max(abs(m2))
# breaksList = seq(min(m2), max(m2), by = (max(m2)-min(m2))/100)
breaksList = seq(-rg, rg, length.out = 100)

m2 = t(na.omit(m2))


pdf(file = 'Row zscore for select traits.pdf')
pheatmap(m2, color = colorRampPalette(c("black","red","orange","gold","white","lightcyan","powderblue","skyblue","dodgerblue"))(length(breaksList)), breaks = breaksList,  
         cluster_rows = T, cluster_cols = F, main='Row zscore for select traits', fontsize_row = 10, fontsize_col = 12)
dev.off()

## Really interesting to compare how strains fall in this clustered heatmaps, some strains seems to be really affected by let
# and their relative scores, while some others, particularly Cast and Akr/j, seems unaffected by the treatment
# locate the clustered strains that showed high LH/FSh ratio and testosterone, two of the key markers in PCOS, then this finding
# does not correlate for every trait, meaning that responses are very variable.
# will discuss in the presentation











######## different way (needs editing and cleaning)



setwd('G:/My Drive/PostDoc UCI/New Project PCOS Leandro UCI 2021/BxD PCOS Project 2021/Results/2022_2023_Data table and  R environment/temporal Leandro')
library(dplyr)
library(reshape2)
library(pheatmap)
library(scales)
#nstall.packages('RColorBrewer')
library(RColorBrewer)
traits_data = all_traits_BxdPcos
colnames(traits_data)
table(traits_data$treatment)
#ctrl  let 
#4961 5181 
table(traits_data$trait_name)
traits_to_plot = c('Testosterone','weight-6_minus_weight-0', 'heart_rate', 'stroke_vol', 'AUC', 'eject_frac', 'Perc_final_increase_BW', 'glucose_0mins', 'lv_mass',  'fat_mass')

traits_data$value = as.numeric(as.character(traits_data$value))
traits_data = traits_data[!traits_data$strain=='BXD124*',]
filtered_traits = traits_data[traits_data$trait_name %in% traits_to_plot,]
filtered_traits = filtered_traits[!is.na(filtered_traits$value),]
filtered_traits$scaled_value  = ave(filtered_traits$value, filtered_traits$trait_name, FUN=scale)

#summary_z_score = filtered_traits %>% dplyr::select(trait_name, #strain, treatment, value) %>% 
# dplyr::group_by(trait_name, strain, treatment) %>% 
#  dplyr::mutate(z_score = scale(value))

summary_z_score = as.data.frame(filtered_traits)
head(summary_z_score)
summary_z_score$strain_treatment = paste0(summary_z_score$strain, '_', summary_z_score$treatment)
mheat = reshape2::dcast(summary_z_score,  trait_name~ strain_treatment, value.var = 'scaled_value', fun.aggregate = mean, na.rm=T)
row.names(mheat) = mheat$trait_name
mheat$trait_name=NULL
m2 = as.matrix(mheat)
m2[is.na(m2)] = 0
range(m2)
breaksList = seq(min(m2), max(m2), by = (max(m2)-min(m2))/1000)
m2 = na.omit(m2)
pdf(file = 'Row zscore for select traits.pdf')
pheatmap(m2, color = colorRampPalette(brewer.pal(n = 7, name = "PuOr"))(length(breaksList)), breaks = breaksList,  cluster_rows = T, cluster_cols = F, main='Row zscore for select traits', fontsize_row = 5, fontsize_col = 5)
dev.off()


m3 = as.data.frame(m2)

write.csv(m3, "z_scores_table.csv")

## I want to make a table with relative changes ("increase", "decrease", "no change") of key traits through the strains

library(broom)
library(Rmisc)
ttests = summary_z_score[!summary_z_score$trait_category=='Age',] %>% 
  group_by(trait_name, strain) %>% 
  filter(length(unique(treatment)) == 2) %>% 
  summarise(t_test = tidy(t.test(scaled_value ~ treatment, data = tibble(treatment, strain, scaled_value)))) %>% 
  unnest(t_test)
ttests = ttests[order(ttests$p.value),]

ttests$let_minus_ctrl = ttests$estimate2 - ttests$estimate1

### this probably needs to be rethinked, I just assigned 0.4 as threshold for assigning increase (>0.4), decrease(< -0.4), or no change (0.4 > z_score > -0.4)
# I'm not using the p values, because a lot of these would be not significant, because the n number is really low when separating by strain and trait
# So I consider the global p value for most of the traits (grouping by trait and comparing let vs ctrl)
ttests$change = ifelse(ttests$let_minus_ctrl> 0.4, paste0('Increase'),
                       ifelse(ttests$let_minus_ctrl< -0.4, paste0('Decrease'),
                              paste0('No_change') ))
ttests$change = as.factor(ttests$change)
summary_stats_z = t(dcast(ttests, trait_name ~ strain, value.var = 'change', fun.aggregate = NULL))
colnames(summary_stats_z) = summary_stats_z[1,]
## table with relative changes in key traits
summary_stats_z = as.data.frame(summary_stats_z[2:25,])
write.csv(summary_stats_z, "z_scores_change_ranking_table.csv")


