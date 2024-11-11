
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

load("map_dat.RData")
load("combined_ppr.RData")

### Prep data
all.new <- st_drop_geometry(all)
for(i in 1:ncol(all.new)){
  all.new[,i] <- as.numeric(unlist(all.new[,i]))
}
colnames(all.new)

need <- cor(all.new)
corrplot(need)

#all.new <- all.new[,-c(1,3,8,13:15)]
all.new <- all.new[,-c(1,3,8:9,14:17)]
colnames(all.new)

states <- unique(counties.ppr$STATE_NAME)
counties.ppr$state_id <- NA
for(i in 1:length(states)){
  counties.ppr$state_id[which(counties.ppr$STATE_NAME == states[i])] <- i
}

counties.ppr$country <- "US"
can.states <- c("Alberta", "Manitoba", "Saskatchewan")
counties.ppr$country[which(counties.ppr$STATE_NAME%in%can.states)] <- "CAN"

counties.ppr$country_id <- 1
counties.ppr$country_id[which(counties.ppr$STATE_NAME%in%can.states)] <- 2

states2 <- matrix(0, nrow(counties.ppr), length(states))
for(i in 1:length(states)){
  states2[which(counties.ppr$STATE_NAME == states[i]), i] <- 1
}
# states2[which(states2[,8:9]==1, arr.ind=T)[,1], 7] <- 1
# states2 <- states2[,1:7]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Run PCA
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### PCA
pc.res <- prcomp(~., data=all.new, scale.=T, center=T)

### Look at summary
summ <- summary(pc.res)
dim(summ)



### Initial plots
fviz_eig(pc.res, addlabels = TRUE) # Eigenvalues (What explains the most variance)
# fviz_pca_var(pc.res, axes=c(1,2)) # See how variables are distributed along PC axes
# fviz_pca_biplot(pc.res) # Biplot (How do the eigenvector looked compared to data points)

### Correlation of variables and PCs
# var <- get_pca_var(pc.res)
# corrplot(var$cor, is.corr=FALSE)
# 
# ### Contribution of variables to each PC
# fviz_contrib(pc.res, choice = "var", axes = 1, top = 10)
# fviz_contrib(pc.res, choice = "var", axes = 2, top = 10)

# fviz_pca_var(pc.res, col.var = "contrib",
#              gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07")
# )

### Can the variables be grouped?
# Create a grouping variable using kmeans
# Create 3 groups of variables (centers = 3)
# res.km <- kmeans(var$coord, centers = 3, nstart = 25)
# grp <- as.factor(res.km$cluster)
# # Color variables by groups
# fviz_pca_var(pc.res, col.var = grp,
#              palette = c("#0073C2FF", "#EFC000FF", "#868686FF", "#FC4E07",
#                          "green", "purple"),
#              legend.title = "Cluster")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Run Geary's C
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pc.want = 1:5
dat <- pc.res$x[,pc.want]
weights <- nb2listw(poly2nb(as(counties.ppr, "Spatial")))

geary_c <- data.frame(PC=1:ncol(dat), C=NA)
for(i in 1:ncol(dat)){
  geary_c[i,2] <- round(geary(dat[,i], weights, length(dat[,i]), length(dat[,i])-1, Szero(weights))$C, 3)
}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Run Redcap
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### How to bound the clusters
counties.ppr$area <- st_area(counties.ppr)
counties.ppr$area <- as.numeric(set_units(counties.ppr$area, km^2))
bound_vals <- counties.ppr["area"]
min_bound <- 8500

grd2$area <- st_area(grd2)
grd2$area <- as.numeric(set_units(grd2$area, km^2))
min_bound/grd2$area


### How many clusters to use?
min.clust = 20
max.clust = 60
CH=c("Calinski_Harabasz", "Calinski-Harabasz")
SD=c("SD_Scat", "SD")
DB=c("Davies_Bouldin", "Davies-Bouldin")
S=c("Silhouette", "Silhouette")
D=c("Dunn", "Dunn")

indices <- as.data.frame(rbind(CH, SD, DB, S, D))
colnames(indices) <- c("index", "name")
indices$clust <- NA


### Put first # of important PCs into data frame and combine with geographic data
dat <- pc.res$x[,pc.want]
dat <- cbind(counties.ppr, dat)

### Neighborhood weight matrices
rook_w <- rook_weights(dat) 

### Specific data values we want
data <- dat[,grep("PC",colnames(dat))] 
data <- cbind(states2, data)
data <- st_drop_geometry(data)

# How correlated are the states with the PCs?
states.cor <- cor(as.matrix(data))
corrplot(states.cor)
mean(states.cor)

for(i in 1:4){
  index = indices[i,]
  
  test <- num.clust(counties.ppr, pc.res, min.clust, max.clust, pc.want,
                    index[1], bound_vals, min_bound, states2)
  # print(test["rs.plot"])
  # print(test["index.plot"])
  n.clust = test$want
  indices$clust[i] <- n.clust
  
  ### Run Redcap!!
  cr <- redcap(n.clust, rook_w, data, "fullorder-completelinkage", scale_method="raw",
               bound_vals, min_bound) 
  cr$`The ratio of between to total sum of squares`

  #~~~~~~~~~~~~~~~~~~
  # Plotting
  #~~~~~~~~~~~~~~~~~~
  ### Get ready to plot
  # Put clusters into data frame 
  dat$pred <- as.factor(cr$Clusters)
  
  # Put together counties in the same cluster
  new <- dat %>%
    group_by(pred) %>%
    summarize(cat = first(pred))
  
  # Plot colors for each cluster
  col.samples <- sample(1:n.clust, n.clust) # Randomize color assignments
  # divergingx_palettes(plot = TRUE)
  cols <- divergingx_hcl(n.clust, "Spectral")[col.samples]
  new$cols <- cols[match(new$pred, col.samples)]
  
  ### Plot!!
  p <- ggplot(new) +
    geom_sf(aes(fill=pred))+
    xlim(st_bbox(grd2)[c(1,3)]) + ylim(st_bbox(grd2)[c(2, 4)]) +
    geom_sf(data=ppr_states, fill=NA, color="grey50", size=0.25) +
    scale_fill_manual(values=cols, aesthetics="fill") +
    theme_classic() +
    theme(legend.position="none") +
    labs(title="REDCAP regionalization", 
         subtitle=paste(index[2], "cluster optimization index"))
  
  png(paste0("plots/",index[1],"_regions.png"), 7, 6, unit="in", res=400)
  print(p)
  dev.off()
  
  assign(paste0(index[1], "_cr"), cr)
}





save.image("pca_res.RData")
