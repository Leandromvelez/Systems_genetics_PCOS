## set working directory
setwd('F:/My Drive/PostDoc UCI/New Project PCOS Leandro UCI 2021/BxD PCOS Project 2021/Results/2022_2023_Data table, R environment and paper in preparation 2023/temporal Leandro')
## Load required data
load('F:/My Drive/PostDoc UCI/New Project PCOS Leandro UCI 2021/BxD PCOS Project 2021/Results/2022_2023_Data table, R environment and paper in preparation 2023/Sequencing_and_traits_dataBxDpcosAugust2023.RData')

## Load required packages (install from CRAN or Bioconductor if necessary)
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
library(MetBrewer)
library(qgraph)


### melt traits data
traits_data = all_traits_and_morph_matrix
traits_data$mouse_number = rownames(traits_data)
traits_data = melt(traits_data, id = 'mouse_number')
traits_data = traits_data %>% rename(trait_name = variable)
traits_data$treatment = sample_info$treatment[match(traits_data$mouse_number, sample_info$mouse_number)]
traits_data$treatment= relevel(as.factor(traits_data$treatment), ref = 'ctrl')
traits_data$strain = sample_info$strain[match(sample_info$mouse_number, traits_data$mouse_number)]

##replace characters that will cause issues
traits_data$trait_name = str_replace_all(traits_data$trait_name, "/", "_")    # Applying str_replace_all
traits_data$trait_name = str_replace_all(traits_data$trait_name, "%", "Perc")
traits_data$trait_name = str_replace_all(traits_data$trait_name, " ", "")    # Applying str_replace_all
traits_data$trait_name = str_replace_all(traits_data$trait_name, ";", "_") 
traits_data$trait_name = str_replace_all(traits_data$trait_name, "-", "_")


## For Model Validation

### This will generate global pcos vs ctrls comparison for each trait, in a pdf in the working directory, 

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
all_ctrlVslet('cysts')

## All t tests for traits, pcos vs ctrls


t_tests_traits = traits_data %>% 
  group_by(trait_name) %>% 
  filter(length(unique(treatment)) == 2) %>% 
  dplyr::summarise(t_test = tidy(t.test(value ~ treatment, data = tibble(treatment, value)))) %>% 
  unnest(t_test)

write.csv(t_tests_traits, 'all_t_tests_traitsBXDPCOS.csv')


#### Phenotypic Variation plots

### This will generate variation plots accounting for strain and treatment, for each trait

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

## call the function (plots will be save in working directory)

plot_strain_variance('fat_mass')
plot_strain_variance('eject_frac')
plot_strain_variance('AUC')



##### Heritability estimates

## Load packages
library(tidyverse)
library(ggrepel)
library(plyr)
library(ggbeeswarm)
library(lme4)
library(ggridges)
# install.packages("nlme")
# install.packages("unglue")
library(unglue)
library(nlme)

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

# Extract the target variable names, excluding "strain", "mouse_number", and "treatment" (grouping variables)

target_vars <- colnames(trait_matrix) [!colnames(trait_matrix) %in% c("mouse_number", "strain", "treatment")]
no.traits <- length(target_vars)


##Estimate the total non-strain variance e.g. LMM1<-lme(expression ~ 1, random=~1|strain, data=data)
# create a named list to hold the fitted models
fitlist1 <- as.list(1:no.traits)
names(fitlist1) <- target_vars
Var1list <- as.list(1:no.traits)  

# loop over trait names
for(i in target_vars){ 
  
  # create temporary data matrix and model formula
  tmp1 <- trait_matrix[, c(i,"strain","treatment")]
  fml1 <- as.formula( paste( i, "~", paste(c("1") ) ))
  
  # assign fit to list by name
  fitlist1[[i]] <- lme(fml1, random=~1|strain, method="REML", data=tmp1, na.action=na.omit)
  
  summ1 <- summary(fitlist1[[i]])
  output1 <- capture.output(summ1)
  ranef1 <- output1[8:9]
  ranef1c <- paste(ranef1, collapse = "\n")
  ranef1m <- readr::read_fwf(ranef1c, readr::fwf_empty(ranef1c), )
  ranef1d <- ranef1m[-1, -1]
  ranef1df <- as.data.frame(ranef1d)
  ranef1df <- ranef1df %>% 
    dplyr::rename(
      strain_mod1 = X2,
      Residual_mod1 = X3
    )
  Var1list[[i]] <- ranef1df
} 

Var1df <- do.call(rbind, Var1list)
Var1df = Var1df[78:154,]

write.csv(Var1df, "mod1_variance_total_nonstrain.csv")


##Estimate the among-strain variance after accounting for treatment e.g. LMM2<-lme(expression ~ treatment, random=~1|strain, data=data)
# create a named list to hold the fitted models
fitlist2 <- as.list(1:no.traits)
names(fitlist2) <- target_vars
Var2list <- as.list(1:no.traits)

# loop over traits names
for(i in target_vars){ 
  
  # create temporary data matrix and model formula
  tmp2 <- trait_matrix[, c(i,"strain","treatment")]
  fml2 <- as.formula( paste( i, "~", paste(c("treatment") ) ))
  
  # assign fit to list by name
  fitlist2[[i]] <- lme(fml2, random=~1|strain, method="REML", data=tmp2, na.action=na.omit)
  
  summ2 <- summary(fitlist2[[i]])
  output2 <- capture.output(summ2)
  ranef2 <- output2[8:9]
  ranef2c <- paste(ranef2, collapse = "\n")
  ranef2m <- readr::read_fwf(ranef2c, readr::fwf_empty(ranef2c), )
  ranef2d <- ranef2m[-1, -1]
  ranef2df <- as.data.frame(ranef2d)
  ranef2df <- ranef2df %>% 
    dplyr::rename(
      strain_mod2 = X2,
      Residual_mod2 = X3
    )
  Var2list[[i]] <- ranef2df
}
Var2df <- do.call(rbind, Var2list)
Var2df = Var2df[78:154,]

write.csv(Var2df, "mod2_variance_amongstrain_treatmentcontrolled.csv")


##Fixed effects
fitlist3 <- as.list(1:no.traits)
names(fitlist3) <- target_vars
FixVar3list <- as.list(1:no.traits)

# loop over traits names
for(i in target_vars){ 
  
  # create temporary data matrix and model formula
  tmp3 <- trait_matrix[, c(i,"strain","treatment")]
  fml3 <- as.formula( paste( i, "~", paste(c("treatment") ) ))
  
  # assign fit to list by name
  fitlist3[[i]] <- lme(fml3, random=~1|strain, method="REML", data=tmp3, na.action=na.omit)
  
  fixsumm3 <- summary(fitlist3[[i]])
  fixoutput3 <- capture.output(fixsumm3)
  fixef3 <- fixoutput3[11:20]
  fixef3c <- paste(fixef3, collapse = "\n")
  fixef3m <- readr::read_fwf(fixef3c, readr::fwf_empty(fixef3c), )
  fixef3d <- fixef3m%>%dplyr::slice(4)
  fixef3d <- fixef3d[1]
  fixef3df <- as.data.frame(fixef3d)
  fixef3df2 <- unglue_unnest(fixef3df, X1, "{strain=\\D+}{Fixed}")
  fixef3df3 <- fixef3df2 %>%
    tidyr::separate(Fixed, c("Fixed_mod3", "col2", "col3", "col4", "col5", "col6"), " ")
  fixef3df4 <- fixef3df3[2]
  FixVar3list[[i]] <- fixef3df4
}

FixVar3df <- do.call(rbind, FixVar3list)
FixVar3df$traits = rownames(FixVar3df)
FixVar3df = as.data.frame(FixVar3df)
FixVar3df = FixVar3df[!rownames(FixVar3df) == FixVar3df$Fixed_mod3,]
FixVar3df$traits = NULL

colnames(FixVar3df) = c("pcos")

write.csv(FixVar3df, "mod3_variance_pcoseffect_straincontrolled.csv")


###Generate delta mice
table(unique(trait_matrix$strain))

trait_matrix$strain <- factor(trait_matrix$strain,
                              levels = c("C57BL/6J","129", "A/J", "AKR/J", "BALBc/j", "Cast/EiJ","CBA/CaJ", "CH3/Hej", "DBA/2J","FVB/NJ", "NOD", "NZBWF1/J", "BXD124", "BXD125", "BXD128", "BXD128a","BXD44", "BXD48a", "BXD51", "BXD60", "BXD70", "BXD71", "BXD73","BXD73b", "BXD75", "BXD77", "BXD78", "BXD79", "BXD89"))

summary(trait_matrix$strain)
str(trait_matrix$strain)

Ctrl_data <- subset(trait_matrix, treatment == "ctrl")

PCOS_data <- subset(trait_matrix, treatment == "let")

trait.names.long <- rep(target_vars, each = 2)
trait.no.long <- length(trait.names.long)
deltalist <- as.list(1:trait.no.long)

for(i in target_vars){
  Ctrl_df <- Ctrl_data[, c("strain","treatment", i)]
  
  Ctrl_seq <- ddply(Ctrl_df,"strain",transform,Number=seq(from=1,by=1,length.out=length(strain)))
  
  PCOS_df <- PCOS_data[, c("strain","treatment", i)]
  
  PCOS_seq <- ddply(PCOS_df,"strain",transform,Number=seq(from=1,by=1,length.out=length(strain)))
  
  Ctrl_seq$strainnumber = paste0(Ctrl_seq$strain,"_", Ctrl_seq$Number)
  
  PCOS_seq$strainnumber = paste0(PCOS_seq$strain,"_", PCOS_seq$Number)
  
  joined_df = Ctrl_seq
  joined_df$pcos = NA
  joined_df$pcos = ifelse(match(joined_df$strainnumber, PCOS_seq$strainnumber), PCOS_seq[,3], joined_df$pcos)

  
  joined_df2 <- joined_df %>% mutate(delta = joined_df[,6]-joined_df[,3])
    sample_df <- joined_df2 %>% group_by(strain, Number) %>% sample_n(4, replace = TRUE)
  sample_df2 <- sample_df[, c("strain","delta")]
    deltalist[[i]] <- sample_df2$delta
} 


deltalistdf <- do.call(rbind, deltalist)
deltalistdf =  deltalistdf[155:231,]
deltalistdf <- data.frame(deltalistdf)


write.csv(deltalistdf, "delta_data.csv")
write.csv(sample_df$strain, "strain_names.csv")

# Manually transpose delta_data and add strain as a column 

delta_data <- read.csv("delta_data.csv", header = TRUE)
delta_data_t = t(delta_data)
colnames(delta_data_t) = delta_data_t[1,]
delta_data_t = as.data.frame(delta_data_t)
delta_data_t = delta_data_t[2:357,]
delta_data_t <- mutate_all(delta_data_t, function(x) as.numeric(as.character(x)))
delta_data_t$strain = sample_df$strain  # I did it this way cause the rownames are the strain names modified (to comply with the unique rowname condition)
# I checked visually that "strain" column coincides with rownames
# Of course there is a better way to do it, this just worked


##Estimate the GxE (gene-by-environment) using delta data e.g. LMM3<-lme(expression ~ 1, random=~1|strain, data=data)
# create a named list to hold the fitted models
fitlist4 <- as.list(1:no.traits)
names(fitlist4) <- target_vars
Var4list <- as.list(1:no.traits)


# loop over traits names
for(i in target_vars){ 
  
  # create temporary data matrix and model formula
  tmp4 <- delta_data_t[, c(i,"strain")]
  fml4 <- as.formula( paste( i, "~", paste(c("1") ) ))
  
  # assign fit to list by name
  fitlist4[[i]] <- lme(fml4, random=~1|strain, method="REML", data=tmp4, na.action=na.omit)
  
  summ4 <- summary(fitlist4[[i]])
  output4 <- capture.output(summ4)
  ranef4 <- output4[8:9]
  ranef4c <- paste(ranef4, collapse = "\n")
  ranef4m <- readr::read_fwf(ranef4c, readr::fwf_empty(ranef4c), )
  ranef4d <- ranef4m[-1, -1]
  ranef4df <- as.data.frame(ranef4d)
  ranef4df <- ranef4df %>% 
    dplyr::rename(
      GxE_mod4 = X2,
      Residual_mod4 = X3
    )
  Var4list[[i]] <- ranef4df
}

Var4df <- do.call(rbind, Var4list)
Var4df =  Var4df[78:154,]

write.csv(Var4df, "mod4_variance_GxE.csv")

###now put it all together to calculate partitioned variances

Mod1 = Var1df
Mod1$trait = rownames(Mod1)
Mod2 = Var2df
Mod2$trait = rownames(Mod2)
Mod3 = FixVar3df
Mod3$trait = rownames(Mod3)
colnames(Mod3) = c("PCOS_mod3", "trait")
Mod4 = Var4df
Mod4$trait = rownames(Mod4)

Join1 <- merge(x = Mod1, y = Mod2, by = "trait", all.x = TRUE)
Join2 <- merge (x = Join1, y = Mod3, by = "trait", all.x = TRUE)
Join3 <- merge (x = Join2, y = Mod4, by = "trait", all.x = TRUE)
All_variances <- Join3 


i <- c(2, 3, 4, 5, 6,7,8)                                  # Specify columns you want to change
# We can now use the apply function to change columns to numeric:

All_variances[ , i] <- apply(All_variances[ , i], 2,            # Specify own function within apply
                             function(x) as.numeric(x))
str(All_variances)

Total <- 1

All_variances <- All_variances %>% 
  mutate(Total_var = strain_mod2+Residual_mod2+PCOS_mod3) %>% 
  mutate(Part_var_strain = (strain_mod2/Total_var)*Total) %>% 
  mutate(Part_var_diet = (PCOS_mod3/Total_var)*Total) %>% 
  mutate(Part_var_res = (Residual_mod2/Total_var)*Total) %>%
  mutate(Sum_part_var = Part_var_strain + Part_var_diet + Part_var_res) %>%
  mutate(Total_var_mod4 = (GxE_mod4 + Residual_mod4)) %>%
  mutate(Part_var_GxE_unadj = (GxE_mod4/Total_var_mod4)*Total) %>%
  mutate(GxE_corr_factor = Total_var/Total_var_mod4) %>%
  mutate(Part_var_GxE = Part_var_GxE_unadj*Part_var_res) %>%
  mutate(Part_var_res_adj = Total-(Part_var_strain+Part_var_diet+Part_var_GxE)) %>%
  mutate(Sum_part_var_2= Part_var_strain + Part_var_diet + Part_var_res_adj + Part_var_GxE)

Partitioned_variances <- All_variances[, c("trait", "Part_var_strain", "Part_var_diet", "Part_var_GxE", "Part_var_res_adj", "Total_var")]
colnames(Partitioned_variances) <- c("trait", "Strain (H2)", "PCOS", "Strain (H2)_x_PCOS", "Residual", "Total_variance")

write.csv(Partitioned_variances, "Partitioned_variances_traits.csv")

part_var_melt = melt(Partitioned_variances, id='trait' )  
part_var_melt$trait = as.factor(part_var_melt$trait)

traits = c('LH_FSH_ratio', 'FSH','estradiol','last_cycle_stage', 'LH','Testosterone','E2_T','Insulin', 'AUC','glucose_0mins',
           'Perc_final_increase_BW', 'weight_6_minus_weight_0','weight_6','weight_1','eject_frac','stroke_vol','lv_mass_corr',
           'heart_rate', 'cardiac_output','fat_mass','lean_fat_ratio','lean_mass','Mean_vel_Asc_Aorta',
           'Mean_vel_Desc_Aorta','peakVel_Asc_Desc')

data = part_var_melt[part_var_melt$trait %in% traits & !part_var_melt$variable == 'Total_variance',]
data$variable <- factor(data$variable, levels = c("Strain (H2)", "PCOS", "Strain (H2)_x_PCOS", 'Residual'))

data %>%
  mutate(name = fct_relevel(trait, 
                            "weight_1", 'weight_6',"lean_mass",'lean_fat_ratio',"fat_mass",'stroke_vol',"AUC","Mean_vel_Desc_Aorta",
                            'cardiac_output',"glucose_0mins",'weight_6_minus_weight_0','Insulin','estradiol','heart_rate',
                            "lv_mass_corr",'peakVel_Asc_Desc', "Mean_vel_Asc_Aorta","eject_frac" ,"Perc_final_increase_BW", 
                            'E2_T', "last_cycle_stage","Testosterone", "LH", "LH_FSH_ratio","FSH")) %>%
  ggplot(aes(x = name, y = value, fill = variable)) +  theme(text = element_text(size = 10)) +
  geom_col(position = "fill") + xlab("traits") + ylab("Variance proportion") + 
  ggtitle('Heritability estimates key traits, BxD PCOS 2023')  + coord_flip() +  scale_y_reverse()





##### Strain personalization figures


library(pheatmap)
library(scales)
#nstall.packages('RColorBrewer')
library(RColorBrewer)
full_traits = traits_data
#ctrl  let 
#4961 5181 
table(full_traits$trait_name)
traits_to_plot = c('Perc_final_increase_BW', 'lean_fat_ratio', 'eject_frac', 'lv_mass','Mean_vel_Desc_Aorta','glucose_0mins', 'AUC','Insulin','estradiol','Testosterone','LH_FSH_ratio')

full_traits$value = as.numeric(as.character(full_traits$value))
full_traits = full_traits[!full_traits$strain=='BXD124*',]
filtered_traits = full_traits[full_traits$trait_name %in% traits_to_plot,]
filtered_traits = filtered_traits[!is.na(filtered_traits$value),]
filtered_traits$scaled_value  = ave(filtered_traits$value, filtered_traits$trait_name, FUN=scale)



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


### Sum of absolute delta Z scores (pcos - ctrl), to get estimation per strain

m2_t = t(m2)

# Initialize an empty dataframe to store the differences
num_rows <- nrow(m2_t)
num_cols <- ncol(m2_t)

# Initialize the resulting dataframe
diff_m2_t <- data.frame(matrix(NA, nrow = num_rows, ncol = num_cols / 2))

# Loop through the rows
for (i in 1:num_rows) {
  # Loop through the columns, taking differences
  for (j in seq(2, num_cols, by = 2)) {
    diff_m2_t[i, j/2] <- m2_t[i, j] - m2_t[i, j - 1]
  }
}

# Rename the columns of the resulting dataframe
colnames(diff_m2_t) <- paste0("Diff_", seq(1, ncol(diff_m2_t)))

# View the resulting dataframe
print(diff_m2_t)

## sum absolute values
abs_sum <- colSums(abs(diff_m2_t))

# View the absolute sum for each column
print(abs_sum)

cleaned_colnames <- sub("_ctrl$|_let$", "", colnames(m2_t))
unique_colnames <- unique(cleaned_colnames)
names(abs_sum) = unique_colnames

abs_sum = as.data.frame(abs_sum)
abs_sum$abs_sum = as.numeric(abs_sum$abs_sum)
abs_sum$strain = rownames(abs_sum)
abs_sum = abs_sum[order(abs_sum$abs_sum),]

abs_sum$strain = NULL

head(abs_sum)
summary(abs_sum)
rg <- max(abs(abs_sum))
# breaksList = seq(min(m2), max(m2), by = (max(m2)-min(m2))/100)
breaksList = seq(2.797, rg, length.out = 100)


pheatmap( abs_sum, color = colorRampPalette(c("black","red","orange","gold","white","lightcyan","powderblue","skyblue","dodgerblue"))(length(breaksList)), breaks = breaksList,  
          cluster_rows =F, cluster_cols = F, main='Row zscore for select traits', fontsize_row = 10, fontsize_col = 12)

## Select the most impacted strains, the less, and black six as the common strain to compare

selected_strains = c('AKR/J', 'CBA/CaJ', 'BXD125', '129', 'BXD60', 'BXD124', 'C57BL/6J')
selected_traits = c('Perc_final_increase_BW', 'lean/fat_ratio', 'glucose_0mins', 'AUC','Insulin','E2/T','Testosterone','LH_FSH_ratio', 'eject_frac', 'lv_mass','cysts', 'corpus_luteum')

filtered_traits = full_traits[full_traits$trait_name %in% selected_traits & full_traits$strain %in% selected_strains,]
filtered_traits = filtered_traits[!is.na(filtered_traits$value),]
filtered_traits$scaled_value  = ave(filtered_traits$value, filtered_traits$trait_name, FUN=scale)


summary_z_score = as.data.frame(filtered_traits)
## exclude the strains without lets or ctrls
summary_z_score = summary_z_score %>% 
  group_by(trait_name, strain) %>% 
  filter(length(unique(treatment)) == 2) 

head(summary_z_score)
summary_z_score$strain_treatment = paste0(summary_z_score$strain, '_', summary_z_score$treatment)
mheat = reshape2::dcast(summary_z_score,  trait_name~ strain_treatment, value.var = 'scaled_value', fun.aggregate = mean, na.rm=T)
## specifi order of traits to plot
mheat = mheat %>% arrange(factor(trait_name, levels = c('Testosterone','E2/T','LH_FSH_ratio', 'cysts', 'corpus_luteum','Perc_final_increase_BW', 'lean/fat_ratio', 'eject_frac', 'lv_mass','glucose_0mins', 'AUC','Insulin')))

row.names(mheat) = mheat$trait_name
mheat$trait_name=NULL
m2 = as.matrix(mheat)
m2[is.na(m2)] = 0
range(m2)

rg <- max(abs(m2))
# breaksList = seq(min(m2), max(m2), by = (max(m2)-min(m2))/100)
breaksList = seq(-rg, rg, length.out = 100)

m2 = t(na.omit(m2))


# pdf(file = 'Row zscore for select traits.pdf')
pheatmap(t(m2), color = colorRampPalette(c("black","red","orange","gold","white","lightcyan","powderblue","skyblue","dodgerblue"))(length(breaksList)), breaks = breaksList,  
         cluster_rows = F, cluster_cols = F, main='Row zscore for select traits', fontsize_row = 10, fontsize_col = 12)
# dev.off()

