## bxd all female data from genenetworks
bxd_traits = read.csv('E:/My Drive/PostDoc UCI/New Project PCOS Leandro UCI 2021/BxD PCOS Project 2021/BxD data from genenetwork2/BXD_traits_forR.csv')
bxd_traits = bxd_traits[,grepl("Symbol|PubMed|Description|C57|DBA|BXD", colnames(bxd_traits))]
bxd_traits =  melt(bxd_traits, id = c("Symbol","Description","PubMed_ID")) 
bxd_traits = bxd_traits[bxd_traits$value != 'x',]
# Applying str_replace_all
bxd_traits$variable = str_replace_all(bxd_traits$variable, "C57BL.6J", "C57BL/6J")
bxd_traits$variable = str_replace_all(bxd_traits$variable, "DBA.2J", "DBA/2J")

bxd_traits = bxd_traits[bxd_traits$variable %in% all_traits_BxdPcos$strain,]

## MDP all female data from genenetworks
mdp_traits = read.csv('E:/My Drive/PostDoc UCI/New Project PCOS Leandro UCI 2021/BxD PCOS Project 2021/BxD data from genenetwork2/MDP_traits_forR.csv')

mdp_traits =  melt(mdp_traits, id = c("Symbol","Description")) 
mdp_traits = mdp_traits[mdp_traits$value != 'x',]
# Applying str_replace_all
mdp_traits$variable <- as.factor(gsub("\\.", "/", as.character(mdp_traits$variable)))
mdp_traits$variable <- as.factor(gsub("X129X1/SvJ", "129X1/SvJ", as.character(mdp_traits$variable)))
## bind all
bxd_mdp_traits = bind_rows(bxd_traits, mdp_traits)
## matrix bxd-mdp traits
bxd_mdp_traits$value = as.numeric(bxd_mdp_traits$value)
bxdmdp_matrix= dcast(bxd_mdp_traits, Description ~ variable, value.var = 'value', fun.aggregate = mean)
#add numbers to traits
bxdmdp_matrix$trait_n = c(1:3880)
trait_bxdmdp = bxdmdp_matrix$Description
names(trait_bxdmdp) = bxdmdp_matrix$trait_n
head(trait_bxdmdp)
trait_bxdmdp = as.data.frame(trait_bxdmdp)
## delete traits description (too long), replace by numbers
bxdmdp_matrix$Description = NULL
rownames(bxdmdp_matrix) = bxdmdp_matrix$trait_n
bxdmdp_matrix$trait_n=NULL
bxdmdp_matrix = t(bxdmdp_matrix)
## keep only columns with data
bxdmdp_matrix = bxdmdp_matrix[,colSums(is.na(bxdmdp_matrix)) < 19]


## cut very large characters tratis
## bxdmdp_matrix$Description =  strtrim(bxdmdp_matrix$Description, 80)

## function to cut strings from the end
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

## make matrix with bxd pcos traits, averaging per strain so we can correlate with mean from genenetworks

trait_matrix_strain= dcast(trait_matrix, strain ~ trait_name, value.var = 'value', fun.aggregate = mean)
row.names(trait_matrix_strain) = trait_matrix_strain$strain
trait_matrix_strain$strain=NULL
## keep common strains
trait_matrix_strain = trait_matrix_strain[rownames(trait_matrix_strain) %in% rownames(bxdmdp_matrix),]
bxdmdp_matrix = bxdmdp_matrix[rownames(bxdmdp_matrix) %in% rownames(trait_matrix_strain),]

## corrs

cc1 = bicorAndPvalue(trait_matrix_strain, 
                     bxdmdp_matrix, 
                     use = 'p')

cc2 = cc1$bicor
tt4 = cc1$p

set1 = melt(cc1$bicor)
head(set1)
colnames(set1) = c('pcos_trait', 'genenetworks_trait', 'bicor')
set2 = melt(tt4)
set1$pvalue = set2$value
set2=NULL
set1$pcos_genenetworks_traits = paste0(set1$pcos_trait,"_", set1$genenetworks_trait)
set1_sign = set1[set1$pvalue<0.05,]
set1_sign = set1_sign[order(set1_sign$pvalue),]
set1_sign = set1_sign[set1_sign$bicor > -1,]

set1_sign$genenetworks_traitname = trait_bxdmdp$trait_bxdmdp[match(set1_sign$genenetworks_trait, rownames(trait_bxdmdp))]
LH_FSH_ratio_sign = set1_sign[set1_sign$pcos_trait=='LH_FSH_ratio',]
LH_FSH_ratio_sign =na.omit(LH_FSH_ratio_sign)

testost_sign = set1_sign[set1_sign$pcos_trait=='Testosterone',]
testost_sign =na.omit(testost_sign)
