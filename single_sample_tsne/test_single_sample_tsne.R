########## Single Sample T-SNE ###########
### This is a basic semi-temp script just for one T-SNE
## Need to work on it to get it to work on lots of samples at once

library(FlowSOM)
library(flowCore)
library(Biobase)
library(ggplot2)
library(MEM)
library(tidyverse)
library(Rtsne)
library(uwot)
library(RColorBrewer)

# Read FCS 
setwd("C:/Users/rbuch/DataShare/T-SNE/FCS/single_sample_tsne_copies/CD3+CD4+/stim/")
files <-  dir(pattern = "*.fcs")
print(files)

data <- lapply(lapply(files, read.FCS), exprs)
ID <- c(1:length(data))
combined.data <- as.data.frame(do.call(rbind, mapply(
  cbind, data, "File ID" = ID, SIMPLIFY = F
)))
combined.data$`File ID` <- as.numeric(combined.data$`File ID`)

# Subset
data_subset <- combined.data[,c(5, 8, 10:11)] # Change to lines in the .fcs that want to subset
chosen.markers <- data_subset

# Transform as part of the function

choose.markers <- function(exp_data) {
  print("Numbered column names, in order they appear in file: ")
  print(paste(c(1:(ncol(exp_data))), ": ", 
              colnames(exp_data[, c(1:(ncol(exp_data)))]), sep = ""))
  markers = readline("Enter column numbers to include (e.g. 1:5,6,8:10).\n")
  sep_vals = unlist(strsplit(markers, ","))
  list_vals = vector()
  for (i in 1:length(sep_vals)) {
    val = sep_vals[i]
    if (length(unlist(strsplit(val, ":"))) > 1) {
      new_val = as.numeric(unlist(strsplit(val, ":"))[1]):
        as.numeric(unlist(strsplit(val, ":"))[2])
    } else{
      new_val = as.numeric(sep_vals[i])
    }
    list_vals = c(list_vals, new_val)
  }
  markerList = c(list_vals)
  return(markerList)
}

transformed.chosen.markers <- chosen.markers %>%
  mutate_all(function(x)
    asinh(x / 15))

# Set seed
overall_seed = 7
set.seed(overall_seed)

# Perplexity 10
my_tsne10 <- Rtsne(transformed.chosen.markers,
                 dims = 2,
                 initial_dims = 12,
                 perplexity = 10,
                 check_duplicates = F,
                 max_iter = 5000
)

tSNE.data <- as.data.frame(my_tsne10$Y)

range <- apply(apply(tSNE.data, 2, range), 2, diff)
graphical.ratio.t <- (range[1] / range[2])

# t-SNE flat dot plot (1 dot = 1 patient)
tSNE.plot <- data.frame(x = tSNE.data[, c(1)], y = tSNE.data[, c(2)])

tnse10p <- ggplot(tSNE.plot) + coord_fixed(ratio = graphical.ratio.t) + 
  geom_point(aes(x = x, y = y), cex = 1) + 
  labs(x = "t-SNE 1", y = "t-SNE 2", title = "t-SNE on PBMC Data Th17 Staining Perplexity 10") + 
  theme_bw()


# Perplexity 30 
my_tsne30 <- Rtsne(transformed.chosen.markers,
                   dims = 2,
                   initial_dims = 12,
                   perplexity = 30,
                   check_duplicates = F,
                   max_iter = 5000
)

tSNE.data <- as.data.frame(my_tsne30$Y)

range <- apply(apply(tSNE.data, 2, range), 2, diff)
graphical.ratio.t <- (range[1] / range[2])

# t-SNE flat dot plot (1 dot = 1 patient)
tSNE.plot <- data.frame(x = tSNE.data[, c(1)], y = tSNE.data[, c(2)])

tsne30p <- ggplot(tSNE.plot) + coord_fixed(ratio = graphical.ratio.t) + 
  geom_point(aes(x = x, y = y), cex = 1) + 
  labs(x = "t-SNE 1", y = "t-SNE 2", title = "t-SNE on PBMC Data Th17 Staining Perplexity 30") + 
  theme_bw()

# Perplexity 45 
my_tsne45 <- Rtsne(transformed.chosen.markers,
                   dims = 2,
                   initial_dims = 12,
                   perplexity = 45,
                   check_duplicates = F,
                   max_iter = 5000
)

tSNE.data <- as.data.frame(my_tsne45$Y)

range <- apply(apply(tSNE.data, 2, range), 2, diff)
graphical.ratio.t <- (range[1] / range[2])

# t-SNE flat dot plot (1 dot = 1 patient)
tSNE.plot <- data.frame(x = tSNE.data[, c(1)], y = tSNE.data[, c(2)])

tsne45p <- ggplot(tSNE.plot) + coord_fixed(ratio = graphical.ratio.t) + 
  geom_point(aes(x = x, y = y), cex = 1) + 
  labs(x = "t-SNE 1", y = "t-SNE 2", title = "t-SNE on PBMC Data Th17 Staining Perplexity 45") + 
  theme_bw()

# UMAP 
myumap <- umap(transformed.chosen.markers,  # input scaled data
       ret_model = TRUE,
       n_threads = 1, 
       verbose = TRUE)
umap.data = as.data.frame(myumap$embedding)

range <- apply(apply(umap.data, 2, range), 2, diff)
graphical.ratio <- (range[1] / range[2])

UMAP.plot <- data.frame(x = umap.data[, 1], y = umap.data[, 2])

UMAP_plot <- ggplot(UMAP.plot, aes(x = x, y = y)) + coord_fixed(ratio = graphical.ratio)  + 
  geom_bin2d(bins = 128) +
  scale_fill_viridis_c(option = "A", trans = "sqrt") + 
  scale_x_continuous(expand = c(0.1, 0)) +
  scale_y_continuous(expand = c(0.1, 0)) + labs(x = "UMAP 1", y = "UMAP 2", 
                                                title = "UMAP on PBMC Data Th17 Staining") + 
  theme_bw()

# Run FlowSOM on the UMAP axes
umap.matrix <- as.matrix(umap.data)

metadata <-  data.frame(name = dimnames(umap.matrix)[[2]],
             desc = paste('t-SNE', dimnames(umap.matrix)[[2]]))
metadata$range <- apply(apply(umap.matrix , 2, range), 2, diff)
metadata$minRange <- apply(umap.matrix , 2, min)
metadata$maxRange <- apply(umap.matrix , 2, max)
umap.flowframe <- new("flowFrame",
                      exprs = umap.matrix ,
                      parameters = AnnotatedDataFrame(metadata))

# create flowFrame for FlowSOM input
UMAP.metadata <-  data.frame(name = dimnames(umap.matrix)[[2]],
             desc = paste('UMAP', dimnames(umap.matrix)[[2]]))
UMAP.metadata$range <- apply(apply(umap.matrix, 2, range), 2, diff)
UMAP.metadata$minRange <- apply(umap.matrix, 2, min)
UMAP.metadata$maxRange <- apply(umap.matrix, 2, max)
umap.flowframe <- new("flowFrame",
                      exprs = umap.matrix,
                      parameters = AnnotatedDataFrame(UMAP.metadata))

# implement the FlowSOM on the data
fsom <- FlowSOM(
    umap.flowframe,      # input flowframe 
    colsToUse = c(1:2),  # columns to use 
    nClus = 6,          # target number of clusters 
    seed = overall_seed  # set seed
  )
FlowSOM.clusters <- as.matrix(fsom[[2]][fsom[[1]]$map$mapping[, 1]])

# plot FlowSOM clusters on UMAP axes
umap_flowsom_plot <- ggplot(UMAP.plot) + coord_fixed(ratio=graphical.ratio) + 
  geom_point(aes(x=x, y=y, color=FlowSOM.clusters),cex = 1.5) + 
  labs(x = "UMAP 1", y = "UMAP 2",title = "FlowSOM Clustering on UMAP Axes", 
       color = "FlowSOM Cluster") + theme_bw() + 
  guides(colour = guide_legend(override.aes = list(size=5)))

# Run MEM on the FlowSOM clusters found using the t-SNE axes
cluster <- as.numeric(as.vector((FlowSOM.clusters)))
MEMdata <- cbind(transformed.chosen.markers, cluster) 

MEM.values.tf = MEM(
  MEMdata,                # input data (last column should contain cluster 
  # values)
  transform = FALSE,      
  cofactor = 0,
  choose.markers = FALSE, 
  markers = "all",
  choose.ref = FALSE,     # each cluster will be compared to all other patient
  # clusters 
  zero.ref = FALSE,
  rename.markers = FALSE,
  new.marker.names = "TNFA, IL-17A, RORgT, IFNG",
  file.is.clust = FALSE,
  add.fileID = FALSE,
  IQR.thresh = NULL
)

# build MEM heatmap and output enrichment scores
MEMheatmap <- build.heatmaps(
  MEM.values.tf,          # input MEM values
  cluster.MEM = "both",
  display.thresh = 1,
  newWindow.heatmaps = FALSE,
  output.files = FALSE,
  labels = TRUE,
  only.MEMheatmap = TRUE
)

