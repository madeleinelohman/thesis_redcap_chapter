
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Libraries and set up
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(clusterCrit)
library(colorspace)
library(corrplot)
library(ggfortify)
library(ggbiplot)
library(factoextra)
library(devtools)
library(rgeoda)
library(spdep)
library(sf)
library(tidyverse)
library(terra)
library(units)

setwd("/Users/madelienelohman/Desktop/thesis_redcap")

source("ppr_n_clust.R")
source("order_clusts.R")

load("pca_res.RData")

importance <- as.data.frame(t(summ[["importance"]]))
importance$`Cumulative Proportion` <- round(importance$`Cumulative Proportion`, 3)
importance$PC <- as.character(1:9)


p <- ggplot(importance, aes(x=PC, y=`Cumulative Proportion`, fill=PC)) +
  geom_bar(stat = "identity") +
  scale_fill_hue(c = 40) +
  geom_text(aes(label=`Cumulative Proportion`, vjust = 1.5)) +
  theme_classic() +
  labs(title="Cumulative proportion of variance explained") +
  theme(legend.position='none')

png("plots/pca_cumvar.png", 7, 6, unit="in", res=400)
print(p)
dev.off()
