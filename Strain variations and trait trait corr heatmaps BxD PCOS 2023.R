setwd('E:/My Drive/PostDoc UCI/New Project PCOS Leandro UCI 2021/BxD PCOS Project 2021/Results/2022_2023_Data table, R environment and paper in preparation 2023/temporal Leandro')
load('E:/My Drive/PostDoc UCI/New Project PCOS Leandro UCI 2021/BxD PCOS Project 2021/Results/2022_2023_Data table, R environment and paper in preparation 2023/Sequencing_and_traits_dataBxDpcosAugust2023.RData')
library(ggplot2)
library(dplyr)
library(forcats)
library(ggpubr)
library(stringr)
library(tidyr)
library(purrr)
library(broom)
library(reshape2)
library(WGCNA)
library(pheatmap)
library(rstatix)

## Load R Environment 

### Strain variation function and plots (plots will be saved in working directory)
traits_data = all_traits_BxdPcos
traits_data$treatment= relevel(as.factor(traits_data$treatment), ref = 'ctrl')

##replace characters that will cause issues
traits_data$trait_name = str_replace_all(traits_data$trait_name, "/", "_")    # Applying str_replace_all
traits_data$trait_name = str_replace_all(traits_data$trait_name, "%", "Perc")
traits_data$trait_name = str_replace_all(traits_data$trait_name, " ", "")    # Applying str_replace_all


plot_strain_variance = function(select_trait){
  pp1 = traits_data[traits_data$trait_name %in% select_trait,]
pp1$strain = gsub(' ', '', pp1$strain)
sum_traits = pp1 %>% dplyr::group_by(strain, treatment) %>% dplyr::summarise(avg=mean(value, na.rm=T), sd=sd(value, na.rm=T), sem=sd(value, na.rm=T)/sqrt(n()))
str_list = sum_traits %>%
  filter(duplicated(strain))
sum_traits = sum_traits[sum_traits$strain %in% str_list$strain,]
sum_traits = na.omit(sum_traits)
sel_order <- 
  sum_traits %>% 
  filter(treatment == "ctrl") %>% 
  arrange(plyr::desc(-avg)) %>% 
  mutate(labels = factor(strain))

plot_data = sum_traits %>% 
  mutate(labels = factor(strain, levels = sel_order$labels, ordered = TRUE))
pdf(file = paste0('PCOS strain variance of ', select_trait, '.pdf'))
g1 = ggplot(plot_data, aes(x=labels, y=avg, fill=treatment)) +  geom_bar(stat="identity", color="black", 
  position=position_dodge()) + theme_classic() + theme(axis.text.x = element_text(angle = 60, vjust = 0.5, hjust=1)) +   ylab(paste0(select_trait)) + xlab('') + ggtitle(paste0(select_trait, ' strain variance')) +
  geom_errorbar(aes(ymin=avg-sem, ymax=avg+sem), width=.2,
  position=position_dodge(.9)) + scale_fill_manual(values=c('dodgerblue3', 'darkgoldenrod2'))
print(g1)
dev.off()
}


### Check for different variables

traits_name = as.data.frame(list(unique(traits_data$trait_name)))
## plots outcomes (check working directory)
# examples
plot_strain_variance('Testosterone')
plot_strain_variance('LH_FSH_ratio')
plot_strain_variance('lean_mass')
plot_strain_variance('fat_mass')
plot_strain_variance('cardiac_output')
plot_strain_variance('lean_fat_ratio')
plot_strain_variance('eject_frac')
plot_strain_variance('AUC')

summary(levels(traits_data$treatment))

### This will generate global let vs ctrls comparison for each trait, in a pdf in the working directory, 
## check the traits dataframe and run the comparison you want

str(traits_data$treatment)

all_ctrlVslet = function(select_trait){
my_comparisons <- rev(list(c("ctrl","let")))
levels(traits_data$treatment) <- c('ctrl','let')
pdf(file = paste0('PCOS vs Ctrl, ', select_trait, '.pdf'))
g1 = ggboxplot(traits_data[traits_data$trait_name %in% select_trait,], x = "treatment", y = "value",
     color = "treatment", palette = "jco",
     add = "jitter", title = paste0(select_trait, "_pcos vs ctrls, Bxd PCOS 2023")) + 
     stat_compare_means(comparisons = my_comparisons, label = "p.signif") 
print(g1)
dev.off()

}

##call the function for diff traits
all_ctrlVslet('Testosterone')
all_ctrlVslet('AUC')

## t tests all samples together, ctrl vs let
tests2 = traits_data[!traits_data$trait_category=='Age',] %>% 
  group_by(trait_name) %>% 
  filter(length(unique(treatment)) == 2) %>% 
  dplyr::summarise(t_test = tidy(t.test(value ~ treatment, data = tibble(treatment, value)))) %>% 
  unnest(t_test)


##### t tests all traits, by strains,  ctrl vs let
## delete constant data
traits_data1 = traits_data %>% 
  group_by(trait_name, strain, treatment) %>% 
  filter(!all(duplicated(value)[-1L]))

stats <- compare_means(value ~ treatment, group.by = c("strain", "trait_name"), data = traits_data1, method = "t.test")
stats_wilc <- compare_means(value ~ treatment, group.by = c("strain", "trait_name"), data = traits_data1, method = "wilcox.test")





#### Multiple ANovas
# split and apply function, result is a named list:

res.aov2= lapply(split(traits_data1, traits_data1$trait_name), function(i){
  anova(lm(value ~ strain + treatment , data = i))
})


res.aov3 = lapply(split(traits_data1, traits_data1$trait_name), function(i){
  aov(value ~ strain * treatment , data = i
)
})

# print results
summary(res.aov3$AUC)

tukeys = TukeyHSD(res.aov3$fat_mass, ordered = FALSE, conf.level = 0.95)
tukeys = as.data.frame(tukeys[3])
tukeys$strain_strain = rownames(tukeys)
plot(res.aov3$AUC, 1)
plot(res.aov3$AUC, 2)

## for unbalanced desings (unequal number of mice per group)
library(car)
res.aov4 <- lapply(split(traits_data1, traits_data1$trait_name), function(i){
  aov(value ~ strain * treatment , data = i)
}) 
Anova(res.aov4$AUC, type = "III",singular.ok = TRUE)
tukeys = TukeyHSD(res.aov4$fat_mass, ordered = FALSE, conf.level = 0.95)
tukeys = as.data.frame(tukeys[3])
tukeys$strain_strain = rownames(tukeys)
plot(res.aov3$AUC, 1)
plot(res.aov3$AUC, 2)


## All traits x all traits Bi Corr Heatmaps


traits_cat=unique(traits_data$trait_category)

traits_data$trait_category <- factor(traits_data$trait_category, levels= traits_cat, ordered=TRUE)
traits_data$trait_id = as.factor(traits_data$trait_id)
traits_data <- traits_data[order(traits_data$trait_category, traits_data$trait_id, decreasing = F),]


## ctrls and let samples
select_ctrl= sample_info$mouse_number[sample_info$treatment=='ctrl']
select_let= sample_info$mouse_number[sample_info$treatment=='let']

trait_matrix= dcast(traits_data, mouse_number ~ trait_name, value.var = 'value', fun.aggregate = mean)
row.names(trait_matrix) = trait_matrix$mouse_number
trait_matrix$mouse_number=NULL
#reorder as all_traits
#traits_name = unique(traits_data$trait_name)

traits_name = c ("weight-0"      ,          "weight-1"       ,         "weight-2"     ,           "weight-3"      ,         
 "weight-4"        ,        "weight-5"         ,       "weight-6"          ,      "weight-6_minus_weight-0",
 "Perc_final_increase_BW",  "glucose_0mins"     ,      "glucose_15mins"     ,     "glucose_30mins"       ,  
 "glucose_45mins"         , "glucose_60mins"       ,   "glucose_90mins"       ,   "AUC"      ,              
 "VTI_Asc_Aorta"         ,  "Mean_vel_Asc_Aorta"   ,   "Mean_grad_Asc_Aorta"   ,  "peak_vel_Asc_Aorta"   ,  
 "peak_grad_Asc_Aorta"   ,  "VTI_Desc_Aorta"       ,  "VTI_Asc/Desc"   , "Mean_vel_Desc_Aorta"   ,  "Mean_grad_Desc_Aorta"   ,
 "peak_vel_Desc_Aorta"   , "peakVel_Asc/Desc" , "peak_grad_Desc_Aorta"  ,  "heart_rate"         ,     "Diameter;s"       ,      
"Diameter;d"           ,   "Volume;s"              ,  "Volume;d"           ,     "stroke_vol"     ,        
 "eject_frac"          ,    "frac_short"           ,   "cardiac_output"    ,      "lv_mass"     ,           
 "lv_mass_corr"         ,   "LVAW;s"               ,   "LVAW;d"            ,      "LVPW;s"      ,           
"LVPW;d"    ,  "res_index_Asc_Aorta"    , "res_index_Desc_Aorta" ,   "pul_index_Asc_Aorta"   ,  "pul_index_Desc_Aorta" ,   "cyclicity"   ,           
"last_cycle_stage"     ,   "lean_mass"             ,  "fat_mass"            ,    "free_water"  ,           
 "total_water"         ,    "lean/fat_ratio"       ,   "fat/weight-6"        ,    "total_water/weight-6" ,  
 "weeks_at_sac"        ,    "weeks_at_start"       ,   "Testosterone"       ,     "FSH"     ,               
 "Insulin"             ,    "estradiol"            ,   "E2/T"              ,      "LH"   ,                  
 "LH_FSH_ratio")  
trait_matrix = trait_matrix[,traits_name]

## exclude non informative traits, redundant or difficult to interpret traits
trait_matrix = trait_matrix[,!colnames(trait_matrix) %in% c('free_water','lv_mass','cyclicity','res_index_Asc_Aorta','weeks_at_start',
 'weeks_at_sac', 'total_water', 'pul_index_Desc_Aorta','pul_index_Asc_Aorta','res_index_Desc_Aorta', 'res_index_Asc_Aorta', 
 'peak_grad_Desc_Aorta','Mean_grad_Desc_Aorta',"peak_grad_Asc_Aorta" ,"Mean_grad_Asc_Aorta" 
                                                           )]



### ctrl


cc1 = bicorAndPvalue(trait_matrix[rownames(trait_matrix) %in% select_ctrl,], 
      trait_matrix[rownames(trait_matrix) %in% select_ctrl,], use = 'p')
cc2 = cc1$bicor
tt4 = cc1$p

## let (PCOS like)

cc1 = bicorAndPvalue(trait_matrix[rownames(trait_matrix) %in% select_let,], 
      trait_matrix[rownames(trait_matrix) %in% select_let,], use = 'p')
cc2_l = cc1$bicor
tt4_l = cc1$p
cc1=NULL

## ctrl
tt3_c = ifelse(tt4 < 0.05 & tt4_l > 0.05,"*","")
tt3_c[lower.tri(tt3_c)] = ""

cc2[lower.tri(cc2)] <- NA

pheatmap(cc2, fontsize_number = 20, display_numbers = tt3_c, 
         number_color = "black", 
         color = colorRampPalette(c("black", "darkred","red","orange","gold","white","lightcyan","powderblue","skyblue","lightskyblue","dodgerblue"))(30), 
         cluster_rows = F, cluster_cols = F, na_col="white", main='Trait x Trait Correlations in Ctrl BxD PCOS MICE', fontsize_row = 5, fontsize_col = 5)


## let (PCOS like)
tt3_l = ifelse(tt4_l < 0.05 & tt4 > 0.05,"*","")
tt3_l[upper.tri(tt3_l)] = ""

cc2_l[upper.tri(cc2_l)] <- NA

pheatmap(cc2_l, fontsize_number = 20, display_numbers = tt3_l, 
         number_color = "black", 
         color = colorRampPalette(c("black", "darkred","red","orange","gold","white","lightcyan","powderblue","skyblue","lightskyblue","dodgerblue"))(30), 
         cluster_rows = F, cluster_cols = F, na_col="white", main='Trait x Trait Correlations in Let BxD PCOS MICE', fontsize_row = 5, fontsize_col = 5)




# sum_stats_table = traits_data %>% dplyr::group_by(strain, trait_category, trait_name, treatment) %>% dplyr::summarise(avg=mean(value), sem = sd(value)/2, n=n())

## Undirected networkds (needs editing and improving! Ideas welcome!)
## controls
map1_ctrl = melt(as.matrix(cc2))
map1_ctrl$value[map1_ctrl$value > 0.999999] <- 0
map1_ctrl = na.omit(map1_ctrl)
map2 = dcast(map1_ctrl, Var1 ~ Var2, value.var = 'value')
row.names(map2) = map2$Var1
map2$Var1 = NULL
#map2 = na.omit(map2)
#dev.off()
col_key = colnames(map2)
new_cols = ifelse(grepl('weight | BW', col_key), 'grey70', 
           ifelse(grepl('glucose | AUC', col_key), 'orange',       
           ifelse(grepl('Aorta', col_key), 'red',
                  
                  met.brewer("Redon"))))
#new_labs = as.data.frame(colnames(map2))
#new_labs$new_lab = goID$ME[match(new_labs$`colnames(map2)`, goID$ME)]
#new_labs1 = ifelse(is.na(new_labs$new_lab), paste0(new_labs$`colnames(map2)`), paste0(new_labs$new_lab))
#pdf(file = 'E:/My Drive/PostDoc UCI/Collaboration Dequina Lab/Networks from lps pituitary with only females gtex august 2022/Undirected network 18 modules.pdf')
qgraph(map2, minimum = 0.3, cut = 0, # vsize = 5, 
legend = F, title ='Undirected Network traits BxD PCOS controls mice' , 
borders = FALSE, layout='spring', posCol = "dodgerblue3", negCol = "firebrick3", label.cex=3, 
directed=F, color=new_cols, labels = colnames(map2))
#dev.off()

## pcos like
map1_let = melt(as.matrix(cc2_l))
map1_let$value[map1_let$value > 0.999999] <- 0
map1_let = na.omit(map1_let)
map2 = dcast(map1_let, Var1 ~ Var2, value.var = 'value')
row.names(map2) = map2$Var1
map2$Var1 = NULL
#map2 = na.omit(map2)
#dev.off()
col_key = colnames(map2)
new_cols = ifelse(grepl('weight | BW', col_key), 'grey70', 
                  ifelse(grepl('glucose | AUC', col_key), 'orange',       
                         ifelse(grepl('Aorta', col_key), 'red',
                                
                                met.brewer("Redon"))))
# new_labs = as.data.frame(colnames(map2))
# new_labs$new_lab = goID$ME[match(new_labs$`colnames(map2)`, goID$ME)]
# new_labs1 = ifelse(is.na(new_labs$new_lab), paste0(new_labs$`colnames(map2)`), paste0(new_labs$new_lab))
#pdf(file = 'E:/My Drive/PostDoc UCI/Collaboration Dequina Lab/Networks from lps pituitary with only females gtex august 2022/Undirected network 18 modules.pdf')
qgraph(map2, minimum = 0.3, cut = 0, # vsize = 5, 
legend = F, title ='Undirected Network traits BxD PCOS. pcos-m mice' , 
borders = FALSE, layout='spring', posCol = "dodgerblue3", negCol = "firebrick3", label.cex=3, 
directed=F, color=new_cols, labels = colnames(map2))
#dev.off()