setwd ('D:/My Drive/PostDoc UCI/New Project PCOS Leandro UCI 2021/BxD PCOS Project 2021/Results/2022_2023_Data table and  R environment/temporal Leandro')
library(ggplot2)
library(dplyr)
library(forcats)
library(tidyverse)
library(ggpubr)
library(ggrepel)
library(ggplot2)
library(plyr)
library(ggbeeswarm)
library(lme4)
library(ggridges)
library(reshape2)
# install.packages("nlme")
# install.packages("unglue")
library(unglue)
library(nlme)
library(stringr)
# traits_data = read.csv('BxD PCOS Screen data table all traits FINAL Seldin Lab May 2023_csv.csv')
#replace "/" character and introduce "_"

#edit data to avoid caracter issues
traits_data = all_traits_BxdPcos
traits_data$trait_name = str_replace_all(traits_data$trait_name, "/", "_")    # Applying str_replace_all
traits_data$trait_name = str_replace_all(traits_data$trait_name, "%", "Perc")
traits_data$trait_name = str_replace_all(traits_data$trait_name, " ", "")    # Applying str_replace_all

traits_data$value = as.numeric(traits_data$value)

traits_data = traits_data[!traits_data$strain=='BXD124*',]
table(traits_data$trait_name)

traits_name=unique(traits_data$trait_category)

traits_data$trait_category <- factor(traits_data$trait_category, levels= traits_name, ordered=TRUE)



traits_data <- traits_data[order(traits_data$trait_category, traits_data$trait_id),]
traits_data <- traits_data[!traits_data$trait_category=='Age',]


# colnames(traits_data) = gsub("_","", colnames(traits_data))
colnames(traits_data) = gsub("-","_", colnames(traits_data))
colnames(traits_data) = gsub(";","_", colnames(traits_data))

# traits_data$trait_name = gsub("_","", traits_data$trait_name)
traits_data$trait_name = gsub("-","_", traits_data$trait_name)
traits_data$trait_name = gsub(";","_", traits_data$trait_name)

table(traits_data$trait_name,is.na(traits_data$value))

## data does not need to be normally distributed (residuals of the regressions do, ideally)

##anyway, we can check how variables are distributed
# ggplot(traits_data, aes(x= value )) + geom_histogram() + facet_wrap(~ trait_name, scales="free", strip.position = c( "bottom")) + ggtitle("Distribution of traits values")
# ggplot(traits_data, aes(x=log(abs(value) + 0.1))) + geom_histogram() + facet_wrap(~ trait_name, scales="free", strip.position = c( "bottom")) + ggtitle("log(Distribution of traits values)")


## I will use the original data, 

## make data matrix
trait_matrix= dcast(traits_data, mouse_number ~ trait_name, value.var = 'value', fun.aggregate = mean)
trait_matrix$treatment = as.factor(traits_data$treatment[match(trait_matrix$mouse_number,traits_data$mouse_number)])
trait_matrix$strain = as.factor(traits_data$strain[match(trait_matrix$mouse_number,traits_data$mouse_number)])
row.names(trait_matrix) = trait_matrix$mouse_number
trait_matrix$mouse_number=NULL

# trait_matrix <- trait_matrix[ , apply(data, 2, function(x) !any(is.na(x)))]


  # Extract the target variable names, excluding "strain", "mouse_number", and "treatment" (grouping variables)




  target_vars <- colnames(trait_matrix) [!colnames(trait_matrix) %in% c("mouse_number", "strain", "treatment")]
  no.traits <- length(target_vars)
  
 
  ##Estimate the total non-strain variance e.g. LMM1<-lme(expression ~ 1, random=~1|strain, data=data)
  # create a named list to hold the fitted models
  fitlist1 <- as.list(1:no.traits)
  names(fitlist1) <- target_vars
  Var1list <- as.list(1:no.traits)  
  
  # loop over protein names
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
  Var1df = Var1df[64:126,]
 ## Variance1 <- Var1df[-seq(from=1, length.out=57, by=1),]
  
  write.csv(Var1df, "mod1_variance_total_nonstrain.csv")
  
  
  ##Estimate the among-strain variance after accounting for treatment e.g. LMM2<-lme(expression ~ treatment, random=~1|strain, data=data)
  # create a named list to hold the fitted models
  fitlist2 <- as.list(1:no.traits)
  names(fitlist2) <- target_vars
  Var2list <- as.list(1:no.traits)
  
  # loop over protein names
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
  Var2df = Var2df[64:126,]
 ## Variance2 <- Var2df[-seq(from=1, length.out=2950, by=1),]
  
  write.csv(Var2df, "mod2_variance_amongstrain_treatmentcontrolled.csv")
  
  
  ##Fixed effects
  fitlist3 <- as.list(1:no.traits)
  names(fitlist3) <- target_vars
  FixVar3list <- as.list(1:no.traits)
  
  # loop over protein names
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
#    FixVariance3 <- FixVar3df[-seq(from=1, length.out=2950, by=1),]
#  FixVariance3 <- as.data.frame(FixVariance3)
#  FixVariance3_named <- FixVariance3 %>% dplyr::rename(treatment = FixVariance3)
#  rownames(FixVariance3_named) <- target_vars
  colnames(FixVar3df) = c("pcos")
  
  write.csv(FixVar3df, "mod3_variance_pcoseffect_straincontrolled.csv")
  
  
  
  ###Generate delta mice
  table(unique(trait_matrix$strain))
  
  trait_matrix$strain <- factor(trait_matrix$strain,
                             levels = c("C57BL/6J","129", "A/J", "AKR/J", "BALBc/j", "Cast/EiJ","CBA/CaJ", "CH3/Hej", "DBA/2J","FVB/NJ", "NOD", "NZBWF1/J", "BXD124", "BXD125", "BXD128", "BXD128a","BXD44", "BXD48a", "BXD51", "BXD60", "BXD70", "BXD71", "BXD73","BXD73b", "BXD75", "BXD77", "BXD78", "BXD79", "BXD89"))
  
  summary(trait_matrix$strain)
  str(trait_matrix$strain)
  # ,labels = c("C57Bl/6","BXH9/TyJ", "Balb/C", "129X1/SvJ", "BXD34/TyJ", "ILSXISS97", "A/J")
  
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
  
  # joined_df <- left_join(Ctrl_seq, PCOS_seq, by = "strain")
  
  
  joined_df2 <- joined_df %>% mutate(delta = joined_df[,6]-joined_df[,3])
  
  sample_df <- joined_df2 %>% group_by(strain, Number) %>% sample_n(4, replace = TRUE)
  sample_df2 <- sample_df[, c("strain","delta")]
  
  deltalist[[i]] <- sample_df2$delta
  } 
  
  
  
  
  deltalistdf <- do.call(rbind, deltalist)
  deltalistdf =  deltalistdf[127:189,]
  
  deltalistdf <- data.frame(deltalistdf)
  
  
  # deltalistdf2 <- deltalistdf%>%slice(5901:8850)
  write.csv(deltalistdf, "delta_data.csv")
  # strain.names.df <- as.data.frame(strain.names)
  write.csv(sample_df$strain, "strain_names.csv")
  
  #Manually transpose delta_data and add strain as a column
  
  delta_data <- read.csv("delta_data.csv", header = TRUE)
  delta_data_t = t(delta_data)
  colnames(delta_data_t) = delta_data_t[1,]
  delta_data_t = as.data.frame(delta_data_t)
  delta_data_t = delta_data_t[2:357,]
  delta_data_t <- mutate_all(delta_data_t, function(x) as.numeric(as.character(x)))
  delta_data_t$strain = sample_df$strain  # I did it this way cause the rownames are the strain names modified (to comply with the unique rowname condition)
  # I checked visually that "strain" column coincides with rownames
  # Of course there is a better way to do it, this just worked
  
  #delta_data_t = as.matrix(delta_data_t)
 
  
  ##Estimate the GxE using delta data e.g. LMM3<-lme(expression ~ 1, random=~1|strain, data=data)
  # create a named list to hold the fitted models
  fitlist4 <- as.list(1:no.traits)
  names(fitlist4) <- target_vars
  Var4list <- as.list(1:no.traits)
  
 
  # loop over protein names
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
  Var4df =  Var4df[58:114,]
 # Variance4 <- Var4df[-seq(from=1, length.out=2950, by=1),]
  
  write.csv(Var4df, "mod4_variance_GxE.csv")
  
  
  ###now put it all together to calculate partitioned variances
  
  # Mod1 <- read.csv("C:/Users/marin/Dropbox (Sydney Uni)/Big Mouse Study/Soleus proteomics/Heritability/mod1_variance_total_nonstrain.csv", header = TRUE)
  Mod1 = Var1df
  Mod1$trait = rownames(Mod1)
  # Mod2 <- read.csv("C:/Users/marin/Dropbox (Sydney Uni)/Big Mouse Study/Soleus proteomics/Heritability/mod2_variance_amongstrain_dietcontrolled.csv", header = TRUE)
 Mod2 = Var2df
 Mod2$trait = rownames(Mod2)
  # Mod3 <- read.csv("C:/Users/marin/Dropbox (Sydney Uni)/Big Mouse Study/Soleus proteomics/Heritability/mod3_variance_dieteffect_straincontrolled.csv", header = TRUE)
  Mod3 = FixVar3df
   Mod3$trait = rownames(Mod3)
   colnames(Mod3) = c("PCOS_mod3", "trait")
  # Mod4 <- read.csv("C:/Users/marin/Dropbox (Sydney Uni)/Big Mouse Study/Soleus proteomics/Heritability/mod4_variance_GxE.csv", header = TRUE)
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
# levels(part_var_melt$trait)[levels(part_var_melt$trait)=='tT'] <- 'Testosterone'

traits = c('LH_FSH_ratio', 'FSH','estradiol','last_cycle_stage', 'LH','Testosterone','E2_T','Insulin', 'AUC','glucose_0mins',
           'Perc_final_increase_BW', 'weight_6_minus_weight_0','weight_6','weight_1','eject_frac','stroke_vol','lv_mass_corr',
           'heart_rate', 'cardiac_output','fat_mass','lean_fat_ratio','lean_mass','Mean_vel_Asc_Aorta',
           'Mean_vel_Desc_Aorta','peakVel_Asc_Desc')

data = part_var_melt[part_var_melt$trait %in% traits & !part_var_melt$variable == 'Total_variance',]
data$variable <- factor(data$variable, levels = c("Strain (H2)", "PCOS", "Strain (H2)_x_PCOS", 'Residual'))

data %>%
  mutate(name = fct_relevel(trait, 
   "weight_1", 'weight_6',"lean_mass",'lean_fat_ratio',"fat_mass",'stroke_vol',
   "AUC","Mean_vel_Desc_Aorta",'cardiac_output',"glucose_0mins",'weight_6_minus_weight_0','Insulin','estradiol','heart_rate',
   "lv_mass_corr",
   'peakVel_Asc_Desc',
   "Mean_vel_Asc_Aorta","eject_frac" ,"Perc_final_increase_BW", 
    'E2_T', "last_cycle_stage","Testosterone", "LH", "LH_FSH_ratio","FSH")) %>%
ggplot(aes(x = name, y = value, fill = variable)) +  theme(text = element_text(size = 10)) +
  geom_col(position = "fill") + xlab("traits") + ylab("Variance proportion") + 
  ggtitle('Heritability estimates key traits, BxD PCOS 2023')  + coord_flip() +  scale_y_reverse()
