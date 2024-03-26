

## All traits x all traits Bi Corr Heatmaps

trait_matrix = all_traits_and_morph_matrix
trait_matrix$mouse_number = rownames(trait_matrix)

trait_matrix$treatment = as.factor(sample_info$treatment[match(trait_matrix$mouse_number,sample_info$mouse_number)])
trait_matrix$strain = as.factor(sample_info$strain[match(trait_matrix$mouse_number,sample_info$mouse_number)])


##replace characters that will cause issues
colnames(trait_matrix) = str_replace_all(colnames(trait_matrix), "/", "_")    # Applying str_replace_all
colnames(trait_matrix) = str_replace_all(colnames(trait_matrix), "%", "Perc")
colnames(trait_matrix) = str_replace_all(colnames(trait_matrix), " ", "")    # Applying str_replace_all
colnames(trait_matrix) = str_replace_all(colnames(trait_matrix), ";", "_") 
colnames(trait_matrix) = str_replace_all(colnames(trait_matrix), "-", "_")

## Exclude strain BxD124* (reasons explained in metadata)
trait_matrix = trait_matrix[!trait_matrix$strain == 'BXD124*',]




## ctrls and let samples
select_ctrl= sample_info$mouse_number[sample_info$treatment=='ctrl']
select_let= sample_info$mouse_number[sample_info$treatment=='let']


#reorder as all_traits

## order of traits
traits_name = c ("weight_0"      ,          "weight_1"       ,         "weight_2"     ,           "weight_3"      ,         
                 "weight_4"        ,        "weight_5"         ,       "weight_6"          ,      "weight_6_minus_weight_0",
                 "Perc_final_increase_BW",  "glucose_0mins"     ,      "glucose_15mins"     ,     "glucose_30mins"       ,  
                 "glucose_45mins"         , "glucose_60mins"       ,   "glucose_90mins"       ,   "AUC"      ,              
                 "VTI_Asc_Aorta"         ,  "Mean_vel_Asc_Aorta"   ,   "Mean_grad_Asc_Aorta"   ,  "peak_vel_Asc_Aorta"   ,  
                 "peak_grad_Asc_Aorta"   ,  "VTI_Desc_Aorta"       ,  "VTI_Asc_Desc"   , "Mean_vel_Desc_Aorta"   ,  "Mean_grad_Desc_Aorta"   ,
                 "peak_vel_Desc_Aorta"   , "peakVel_Asc_Desc" , "peak_grad_Desc_Aorta"  ,  "heart_rate"         ,     "Diameter_s"       ,      
                 "Diameter_d"           ,   "Volume_s"              ,  "Volume_d"           ,     "stroke_vol"     ,        
                 "eject_frac"          ,    "frac_short"           ,   "cardiac_output"    ,      "lv_mass"     ,           
                 "lv_mass_corr"         ,   "LVAW_s"               ,   "LVAW_d"            ,      "LVPW_s"      ,           
                 "LVPW_d"    ,  "res_index_Asc_Aorta"    , "res_index_Desc_Aorta" ,   "pul_index_Asc_Aorta"   ,  "pul_index_Desc_Aorta" ,   "cyclicity"   ,           
                 "last_cycle_stage"     ,   "lean_mass"             ,  "fat_mass"            ,    "free_water"  ,           
                 "total_water"         ,    "lean_fat_ratio"       ,   "fat_weight_6"        ,    "total_water_weight_6" ,  
                 "weeks_at_sac"        ,    "weeks_at_start"       ,   "Testosterone"       ,     "FSH"     ,               
                 "Insulin"             ,    "estradiol"            ,   "E2_T"              ,      "LH"   ,                  
                 "LH_FSH_ratio", "primary_to_smallantral" , "ovary_area.mm2."   ,      "perimeter..mm."   ,      
                  "primordial"       ,       "antral"           ,       "atretic_antral"    ,      "cysts"   ,               
                  "corpus_luteum"     ,      "blood_infilt"   ,         "fatty_ovary"      ,       "debris_atretic_follic" , 
                  "high_vasc"     )  
trait_matrix2 = trait_matrix[,traits_name]

## exclude non informative traits, redundant or difficult to interpret traits
trait_matrix2 = trait_matrix2[,!colnames(trait_matrix2) %in% c('free_water','lv_mass','cyclicity','res_index_Asc_Aorta','weeks_at_start',
                                                            'weeks_at_sac', 'total_water', 'pul_index_Desc_Aorta','pul_index_Asc_Aorta','res_index_Desc_Aorta', 'res_index_Asc_Aorta', 
                                                            'peak_grad_Desc_Aorta','Mean_grad_Desc_Aorta',"peak_grad_Asc_Aorta" ,"Mean_grad_Asc_Aorta" ,
                                                            "blood_infilt"     ,      "fatty_ovary"      ,      "debris_atretic_follic" ,
                                                            "high_vasc"                                    )]



### ctrl


cc1 = bicorAndPvalue(trait_matrix2[rownames(trait_matrix2) %in% select_ctrl,], 
                     trait_matrix2[rownames(trait_matrix2) %in% select_ctrl,], use = 'p')
cc2 = cc1$bicor
tt4 = cc1$p


set1 = melt(cc2)
head(set1)
colnames(set1) = c('pcos_trait1', 'pcos_trait2', 'bicor')
set2 = melt(tt4)
set1$pvalue = set2$value
set2=NULL
set1$trait_trait = paste0(set1$pcos_trait1,"_", set1$pcos_trait2)
set1 = set1[!(set1$pcos_trait1 == set1$pcos_trait2),]
## save trait x trait bicor and p values
set1_ctrl = set1
set1=  NULL
cc1= NULL




## let (PCOS like)

cc1 = bicorAndPvalue(trait_matrix2[rownames(trait_matrix2) %in% select_let,], 
                     trait_matrix2[rownames(trait_matrix2) %in% select_let,], use = 'p')
cc2_l = cc1$bicor
tt4_l = cc1$p

set1 = melt(cc2_l)
head(set1)
colnames(set1) = c('pcos_trait1', 'pcos_trait2', 'bicor')
set2 = melt(tt4_l)
set1$pvalue = set2$value
set2=NULL
set1$trait_trait = paste0(set1$pcos_trait1,"_", set1$pcos_trait2)
set1 = set1[!(set1$pcos_trait1 == set1$pcos_trait2),]
## save trait x trait bicor and p values
set1_let = set1
set1=  NULL
cc1= NULL


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
cc2=NULL
cc2_l = NULL
tt3_c = NULL
tt3_l = NULL
cc1 = NULL


######## PCOS impact network

 ## in PCOS

tt1 = trait_matrix2[rownames(trait_matrix2) %in% select_let,]

cors1 = bicorAndPvalue(tt1, tt1, use = 'p')
cors2 = melt(cors1$bicor)
colnames(cors2) = c('trait_1', 'trait_2', 'bicor')
cors2$pvalue = melt(cors1$p)$value
cors2 = cors2[!cors2$trait_1==cors2$trait_2,]
cors2 = cors2[order(cors2$pvalue, decreasing = F),]
cors2$trait1_trait2 = paste0(cors2$trait_1, ' ~ ', cors2$trait_2)

let_table = cors2  

 ## In ctrls

tt1 = trait_matrix2[rownames(trait_matrix2) %in% select_ctrl,]

cors1 = bicorAndPvalue(tt1, tt1, use = 'p')
cors2 = melt(cors1$bicor)
colnames(cors2) = c('trait_1', 'trait_2', 'bicor')
cors2$pvalue = melt(cors1$p)$value
cors2 = cors2[!cors2$trait_1==cors2$trait_2,]
cors2 = cors2[order(cors2$pvalue, decreasing = F),]
cors2$trait1_trait2 = paste0(cors2$trait_1, ' ~ ', cors2$trait_2)

cntl_table = cors2 

cntl_table$let_bicor = let_table$bicor[match(cntl_table$trait1_trait2, let_table$trait1_trait2)]
cntl_table$let_pvalue = let_table$pvalue[match(cntl_table$trait1_trait2, let_table$trait1_trait2)]

cntl_table$delta_bicor = cntl_table$let_bicor - cntl_table$bicor
cntl_table = cntl_table[order(abs(cntl_table$delta_bicor), decreasing = T),]
head(cntl_table)

# cntl_table = cntl_table[!duplicated(abs(cntl_table$delta_bicor)),]
write.csv(cntl_table, file = 'trait delta table.csv', row.names = F)

table(cntl_table$trait_1)
###build network for change - preselect traits
traits_set = c('LH_FSH_ratio', 'FSH','','estradiol','last_cycle_stage', 'LH','Testosterone','E2_T','Insulin', 'AUC','glucose_0mins', 
               'weight_6_minus_weight_0','weight_6', 'weight_1', 'eject_frac','stroke_vol','lv_mass_corr','heart_rate', 'cardiac_output',
               'fat_mass','lean_fat_ratio','lean_mass','Mean_vel_Asc_Aorta','Mean_vel_Desc_Aorta', "cysts", "corpus_luteum", "ovary_area.mm2.")


map1 = cntl_table[cntl_table$trait_1 %in% traits_set & cntl_table$trait_2 %in% traits_set,]

map2 = reshape2::dcast(map1, trait_1 ~ trait_2, value.var = 'delta_bicor', fun.aggregate = mean)
row.names(map2) = map2$trait_1
map2$trait_1 = NULL
map2 = as.matrix(map2)

#add color scheme
#install.packages("igraph")
library(igraph)
dput(map2)
names<-c("1","2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26")

network_graph = qgraph(map2, minimum = 0.1, cut = 0.25, # vsize = 5, 
                       legend = T, title ='Undirected network for PCOS impact', 
                       borders = T, layout='spring', posCol = "dodgerblue3", negCol = "firebrick3", label.cex=2, 
                       directed=F, nodeNames = colnames(map2), labels = names)

plot(network_graph)




## centrality estimates

g1 = centralityPlot(network_graph, include = 
                      c("Strength","Betweenness","Closeness"),
                    orderBy ="Strength")
title <- cowplot::ggdraw() + cowplot::draw_label('Centrality estimates for delta bicor undirected Network, BxD PCOS', fontface='bold')


cowplot::plot_grid(title, g1, ncol=1, rel_heights=c(0.05, 1)) # rel_heights values control title margins for each element (title, g5, g6)






############waterfallplot
plot_table = cntl_table[cntl_table$trait_1 %in% traits_set & cntl_table$trait_2 %in% traits_set,]
plot_table = plot_table[!duplicated(plot_table$bicor),]

ggplot(plot_table[1:20,], aes(x=fct_reorder(trait1_trait2, abs(delta_bicor), .desc=T), y=delta_bicor, fill=trait1_trait2)) + 
  geom_col(position = position_dodge2(reverse = TRUE)) +  theme_classic() + theme(legend.position = "none") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))  + xlab('') + ylab('bicor coefficent') + 
  ggtitle(paste0('Top trait correlations shifting in PCOS over cntl ')) + ylab('PCOS bicor - cntl bicor')

## count freq of traits in the top 20 trait-trait corr shifting

trait1_n = plot_table[1:20,] %>%   dplyr::group_by(trait_1) %>%  dplyr::summarise(n())
trait2_n = plot_table[1:20,] %>%   dplyr::group_by(trait_2) %>%  dplyr::summarise(n())
colnames(trait2_n) = colnames(trait1_n)  
trait_n = rbind(trait1_n, trait2_n)
colnames(trait_n)= c("trait", "freq")
trait_n$trait = as.character(trait_n$trait)
trait_n = group_by(trait_n, trait) %>%
  summarise(Freq = sum(freq))
## plot
trait_n %>% mutate(trait_1 = fct_reorder(trait, Freq))  %>%
  ggplot(  aes(x= trait_1, y= Freq, fill= trait_1)   ) + geom_col(show.legend = FALSE) + 
  
  theme(axis.text.x = element_text(angle = 90)) + coord_flip() + 
  ggtitle('frequency of traits from top shifting corrs PCOS - ctrls')

trait_matrix2 = NULL




##### Scatter plots for unique trait corrs 


## define unique in PCOS and look for interesting correlations to make scattter plots
set1_unique_pcos = set1_let[set1_let$pvalue<0.05,]
set1_unique_pcos = set1_unique_pcos[!set1_unique_pcos$trait_trait %in% set1_ctrl$trait_trait[set1_ctrl$pvalue <0.05],]


### make function for plot trait trait scatter plots

# summarize mean, sd and sem per trait
sum_stats_table = traits_data %>% dplyr::group_by(trait_name, treatment) %>% dplyr::summarise(avg=mean(value, na.rm=TRUE), sem = sd(value, na.rm=TRUE)/sqrt(length(value)), n=n())

## function
plot_corrs_traits = function(trait1,trait2 ){
  g1 = ggscatter(trait_matrix, 
                 x = paste0(trait1), y = paste0(trait2),  color = 'treatment', shape = 'treatment',
                 xlab=paste0(trait1),
                 ylab = paste0(trait2),
                 ellipse = F, mean.point = TRUE,
                 # add = "reg.line",  # Add regressin line
                 # add.params = list(linetype = "treatment") ,
                 palette = c('dodgerblue3', 'darkgoldenrod2'),
                 # add.params = list(color = "black", fill = "lightgray"), # Customize reg. line
                 conf.int = TRUE) + # stat_cor(method = "pearson") +
    ggtitle(paste0(trait1,"_x_", trait2,"_BxD PCOS Project All samples"))
  
  g3= ggplot(sum_stats_table[sum_stats_table$trait_name == trait1,], aes(x=treatment, y=avg, fill=treatment)) +  geom_bar(stat="identity", color="black", 
                                                                                                                          osition=position_dodge()) + theme_classic() + theme(axis.text.x = element_text(angle = 60, vjust = 0.5, hjust=1)) +   ylab(paste0(trait1)) + xlab('') + ggtitle(paste0(trait1, ' BxD PCOS')) +
    geom_errorbar(aes(ymin=avg-sem, ymax=avg+sem), width=.2,
                  position=position_dodge(.9)) + scale_fill_manual(values=c('dodgerblue3', 'darkgoldenrod2')) 
  
  g4= ggplot(sum_stats_table[sum_stats_table$trait_name == trait2,], aes(x=treatment, y=avg, fill=treatment)) +  geom_bar(stat="identity", color="black", 
                                                                                                                          position=position_dodge()) + theme_classic() + theme(axis.text.x = element_text(angle = 60, vjust = 0.5, hjust=1)) +   ylab(paste0(trait2)) + xlab('') + ggtitle(paste0(trait2, ' BxD PCOS')) +
    geom_errorbar(aes(ymin=avg-sem, ymax=avg+sem), width=.2,
                  position=position_dodge(.9)) + scale_fill_manual(values=c('dodgerblue3', 'darkgoldenrod2'))
  g5 = cowplot::plot_grid(g3, g4, ncol=1) 
  g2 = cowplot::plot_grid(g1, g5, ncol=2) 
  title <- cowplot::ggdraw() + cowplot::draw_label(paste0(trait1, ' x ', paste0(trait2), '_BxD PCOS'), fontface='bold')
  pdf(file = paste0("Corr ",trait1,'_x_', trait2, ".pdf"), width = 8, height = 6)
  g7 = cowplot::plot_grid(title, g2, ncol=1, rel_heights=c(0.1, 1)) # rel_heights values control title margins for each element (title, g5, g6)
  print(g7)
  dev.off()
}
 ## call function 
plot_corrs_traits("Testosterone","LH_FSH_ratio") ## check working directory for saved plots
plot_corrs_traits("Diameter_d","weight_6")
## check t tests from figure 1 scripts,p value for the traits and paste in the plot 
