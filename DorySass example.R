library(fastDiversity)
library(dplyr)
library(ggplot2)
library(ggthemes)
data(example_data)

gt <- example_gt
genetic_group_variable <- example_meta$population
site_variable <- example_meta$site
allele_list_population <- make_allele_list(example_gt, example_meta$population)
ggvenn::ggvenn(allele_list_population)
private_total_alleles_population  <- calculate_private_alleles(allele_list_population)

basicstats <- faststats(example_gt, example_meta$population,
                      example_meta$site, minimum_n=6,
                       minimum_loci=50, maf=0.05, max_missingness=0.3)


pca <- prcomp(dist(example_gt), center=TRUE, scale. = TRUE)
s_pca <- summary(pca)
importance <- s_pca[["importance"]][2,]
axis_names <- paste0(names(importance),' (', round(importance*100,1),'%)')
pcs <- pca$x %>% data.frame()
pcs <- merge(pcs, example_meta, by.x=0, by.y='sample')
ggplot(pcs, aes(x=PC1, y=PC2, color=factor(population)))+
  geom_point()+
  theme_bw()+
  labs(x=axis_names[1], y=axis_names[2])
