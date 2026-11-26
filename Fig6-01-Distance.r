library(Seurat);library(Matrix);suppressMessages(library(dplyr));library(ggplot2);#suppressMessages(library(monocle3));suppressMessages(library(monocle));
suppressMessages(library(cowplot));suppressMessages(library('pagoda2'));suppressMessages(library("velocyto.R"));suppressMessages(library("patchwork"))
library("RColorBrewer");library(scales);suppressMessages(library(SeuratDisk));suppressMessages(library(SeuratWrappers));suppressMessages(library(SeuratData));
options(repr.plot.width=7, repr.plot.height=6); 
setwd("/home/ahmu/YLL/Xenium/dist/")
getwd();.libPaths();#setwd('/sdb/xxjiang/data');getwd();

length(dir('/home/ahmu/R/x86_64-conda-linux-gnu-library/4.4'))
length(dir('/home/ahmu/anaconda3/envs/r443/lib/R/library'))
packageVersion("Seurat")

xenium.obj = readRDS( file = "rctd/ashleyint/B10-2s_24675.rds");xenium.obj



m1=read.csv("/home/ahmu/YLL/Xenium/dist/rctd/Selection_1_cells_stats.csv",check.names = FALSE, header=TRUE,row.names=1);dim(m1);head(m1,2)##,skip = 2,
m2=read.csv("/home/ahmu/YLL/Xenium/dist/rctd/Selection_2_cells_stats.csv",check.names = FALSE, header=TRUE,row.names=1);dim(m2);head(m2,2)
m3=read.csv("/home/ahmu/YLL/Xenium/dist/rctd/Selection_3_cells_stats.csv",check.names = FALSE, header=TRUE,row.names=1);dim(m3);head(m3,2)


rownames(m1) <- paste0("B10-2_", rownames(m1))
rownames(m2) <- paste0("B10-2_", rownames(m2))
rownames(m3) <- paste0("B10-2_", rownames(m3))

X1 <- subset(xenium.obj, cells = rownames(m1));X1
X2 <- subset(xenium.obj, cells = rownames(m2));X2
X3 <- subset(xenium.obj, cells = rownames(m3));X3


Idents(xenium.obj)='typex2'
Idents(xenium.obj, cells = rownames(m3)) <- "Decidua-deep zone"
Idents(xenium.obj, cells = rownames(m2)) <- "Decidua-superficial zone"
Idents(xenium.obj, cells = rownames(m1)) <- "Chorion-side zone"
xenium.obj$region=Idents(xenium.obj)
table(Idents(xenium.obj))

colnames(xenium.obj)

options(repr.plot.width=17, repr.plot.height=10)
ImageDimPlot(xenium.obj, fov = "fov.8",coord.fixed =T,flip_xy = F)

Idents(xenium.obj, cells =  WhichCells(xenium.obj,idents = c('1   dMac',"2   dNK","4   EVT" ,"5   other"))) <- "Decidua-superficial zone"
xenium.obj$region=Idents(xenium.obj)
table(Idents(xenium.obj))

Idents(xenium.obj, cells =  WhichCells(xenium.obj,idents = "Chorion-side zone")) <- "Chorion-side zone"
xenium.obj$region=Idents(xenium.obj)
table(Idents(xenium.obj))

saveRDS(xenium.obj,"rctd/ashleyint/B10-2s_24675.rds",compress = F);xenium.obj



library(plyr); library(dplyr);library(scales)
meta=as.data.frame(table(xenium.obj$region,xenium.obj$typex2)) colnames(meta)=c("Region","Type","Cells")
meta <- ddply(meta, "Region", transform,  Fraction=Cells/sum(Cells)*100,Fraction2= percent(Cells/sum(Cells), accuracy = 0.1))
#head(meta);dim(meta)
meta
write.csv(meta, "No.of cells in regions.csv")

Idents(xenium.obj)='typex2'
xenium.obj1=subset(xenium.obj,idents="5   other",invert=T)
library(plyr); library(dplyr);library(scales)
meta=as.data.frame(table(xenium.obj1$region,xenium.obj1$typex2))colnames(meta)=c("Region","Type","Cells")
meta <- ddply(meta, "Region", transform,  Fraction=Cells/sum(Cells)*100,Fraction2= percent(Cells/sum(Cells), accuracy = 0.1))#
#head(meta);dim(meta)
meta






Idents(xenium.obj)='region'
table(Idents(xenium.obj))

options(repr.plot.width=15, repr.plot.height=6)
ImageDimPlot(xenium.obj, 
             fov = "fov.8",
             coord.fixed = T,flip_xy=F,
             molecules = c("CD68"),
             crop = T,
             cols = c("Chorion-side zone" = "black",
                      "Decidua-superficial zone" = "gray50",
                      "Decidua-deep zone" = "gray80"
                     ),
             mols.alpha = 0.5, 
             size = 0.75, nmols = 5000,#,border.size = 0.1,
             dark.background =F) & 
  theme(legend.title = element_blank(),
        text = element_text(face = "bold", size = 20))



table(xenium.obj$typex2)

search()
detach("package:plyr", unload = TRUE)

library(ggplot2)
#library(dplyr)
library(reshape2)

calculate_detailed_distances <- function(data, target_region, evt_type) {
  evt_cells <- data[data$region == "Chorion-side zone" & data$cell_type == evt_type, ]
  
  if (nrow(evt_cells) == 0) {
    cat("jinggao：Chorion-side zone", evt_type, "\n")
    return(NULL)
  }
  
  target_cells <- data[data$region == target_region, ]
  
  if (nrow(target_cells) == 0) {
    cat("jinggao：", target_region, "\n")
    return(NULL)
  }
  
  distances <- c()
  cell_types <- c()
  
  for (i in 1:nrow(target_cells)) {
    dists <- sqrt((target_cells$x[i] - evt_cells$x)^2 + (target_cells$y[i] - evt_cells$y)^2)
    k <- min(5, length(dists))
    avg_distance <- mean(sort(dists)[1:k])
    distances <- c(distances, avg_distance)
    cell_types <- c(cell_types, target_cells$cell_type[i])
  }

  result_df <- data.frame(
    cell_type = cell_types,
    distance = distances,
    region = target_region,
    stringsAsFactors = FALSE
  )
  
  return(result_df)
}

target_regions <- c("Chorion-side zone","Decidua-superficial zone", "Decidua-deep zone")
all_detailed_distances <- data.frame()

for (reg in target_regions) {
  dist_data <- calculate_detailed_distances(cell_data, reg, "4   EVT")
  if (!is.null(dist_data)) {
    all_detailed_distances <- rbind(all_detailed_distances, dist_data)
  }
}


breaks <- c(0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100,1200,1300,1400,Inf)
labels <- c("0-100", "100-200", "200-300", "300-400", "400-500", "500-600", 
            "600-700", "700-800","800-900","900-1000","1000-1100","1100-1200",
            "1200-1300","1300-1400",">1400")

all_detailed_distances$distance_group <- cut(all_detailed_distances$distance, 
                                              breaks = breaks, 
                                              labels = labels, 
                                              right = FALSE)

distance_stats <- all_detailed_distances %>%
  group_by(region, cell_type, distance_group) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(region, cell_type) %>%
  mutate(percentage = count / sum(count) * 100) %>%
  ungroup()
）
overall_stats <- all_detailed_distances %>%
  group_by(region, distance_group) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(region) %>%
  mutate(percentage = count / sum(count) * 100) %>%
  ungroup()

write.csv(distance_stats, "detailed_distance_distribution_stats.csv", row.names = FALSE)
write.csv(overall_stats, "overall_distance_distribution_stats.csv", row.names = FALSE)


distance_stats





library(ggplot2)
#library(dplyr)
library(reshape2)

calculate_detailed_distances <- function(data, target_region, evt_type) {
  evt_cells <- data[data$region == "Chorion-side zone" & data$cell_type == evt_type, ]
  
  if (nrow(evt_cells) == 0) {
    cat("jinggao：Chorion-side zone no", evt_type, "cell\n")
    return(NULL)
  }
  
  target_cells <- data[data$region == target_region, ]
  
  if (nrow(target_cells) == 0) {
    cat("jinggao：", target_region, "no cell\n")
    return(NULL)
  }
  
  
    distances <- c()
  cell_types <- c()
  
  for (i in 1:nrow(target_cells)) {
       dists <- sqrt((target_cells$x[i] - evt_cells$x)^2 + (target_cells$y[i] - evt_cells$y)^2)
        k <- min(5, length(dists))
    avg_distance <- mean(sort(dists)[1:k])
    distances <- c(distances, avg_distance)
    cell_types <- c(cell_types, target_cells$cell_type[i])
  }
  
  result_df <- data.frame(
    cell_type = cell_types,
    distance = distances,
    region = target_region,
    stringsAsFactors = FALSE
  )
  
  return(result_df)
}

target_regions <- c("Chorion-side zone","Decidua-superficial zone", "Decidua-deep zone")
all_detailed_distances <- data.frame()

for (reg in target_regions) {
  dist_data <- calculate_detailed_distances(cell_data, reg, "4   EVT")
  if (!is.null(dist_data)) {
    all_detailed_distances <- rbind(all_detailed_distances, dist_data)
  }
}

print(summary(all_detailed_distances$distance))

breaks <- c(0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, Inf)
labels <- c("0-100", "100-200", "200-300", "300-400", "400-500", "500-600", 
            "600-700", "700-800", "800-900", "900-1000", "1000-1100", "1100-1200",
            "1200-1300", "1300-1400", ">1400")

all_detailed_distances$distance_group <- cut(all_detailed_distances$distance, 
                                             breaks = breaks, 
                                             labels = labels, 
                                             right = FALSE)

distance_stats <- all_detailed_distances %>%
  group_by(cell_type, distance_group) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(cell_type) %>%
  mutate(percentage = count / sum(count) * 100) %>%
  ungroup()

print(distance_stats)

write.csv(distance_stats, "celltype_distance_distribution_stats.csv", row.names = FALSE)

distance_stats

meta=distance_stats

meta <- meta[meta$cell_type != 5, ]
head(meta)

desired_order <- c("0-100", "100-200", "200-300", "300-400", "400-500", "500-600", 
            "600-700", "700-800","800-900","900-1000","1000-1100","1100-1200",
            "1200-1300","1300-1400",">1400")

meta$distance_group <- factor(meta$distance_group, 
                               levels = desired_order, 
                               ordered = TRUE)


library(plyr)
meta <- ddply(meta, "distance_group", transform, 
              Fraction = count / sum(count) * 100)

head(meta)

library(ggplot2)
options(repr.plot.width=12, repr.plot.height=6)

ggplot(meta, aes(x = distance_group, y = Fraction, fill = cell_type))&
  geom_col(position = "fill") &
  theme(  
    legend.position = "right",
      plot.title=element_text(size=20,hjust=-0.1),axis.text = element_text(size=20,face='bold'),legend.key.height=unit(0.45,'inch'),
      legend.text=element_text(size=20), legend.title=element_blank())& 
  scale_y_continuous(labels = scales::percent_format()) & 
  labs(
    title = "Cell Type Composition by Distance Group",
    x = "",
    y = "cell proportion",
    fill = "Cell Type"
  ) &
  theme_minimal() 

meta$cell_type <- as.factor(meta$cell_type)

polychrome = c('#E31A1C','#00CC52','#E67E22',"black")


options(repr.plot.width=9, repr.plot.height=6)

ggplot(meta, aes(x = distance_group, y = Fraction, fill = cell_type)) +
  geom_col(position = "fill") +         
  theme_minimal() +  
  theme(
    legend.position = "right",
    plot.title = element_text(size = 20, hjust = -0.1),
    axis.text = element_text(size = 20, face = 'bold'),
    axis.text.x = element_text(angle = 45, hjust = 1),  
       axis.title.y = element_text(size = 20, face = 'bold'),
    legend.key.height = unit(0.45, 'inch'),
    legend.text = element_text(size = 20),
    legend.title = element_blank()
  ) + 
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = polychrome) +  
  labs(
    title = "Cell Type Composition by Distance Group",
    x = "",
    y = "cell proportion"
  )





meta=distance_stats

desired_order <- c("0-100", "100-200", "200-300", "300-400", "400-500", "500-600", 
            "600-700", "700-800","800-900","900-1000","1000-1100","1100-1200",
            "1200-1300","1300-1400",">1400")

meta$distance_group <- factor(meta$distance_group, 
                               levels = desired_order, 
                               ordered = TRUE)

library(plyr)
meta <- ddply(meta, "distance_group", transform, 
              Fraction = count / sum(count) * 100)

head(meta)

library(ggplot2)
options(repr.plot.width=12, repr.plot.height=6)

ggplot(meta, aes(x = distance_group, y = Fraction, fill = cell_type))&
  geom_col(position = "fill") &        
  theme(  
    legend.position = "right",
      plot.title=element_text(size=20,hjust=-0.1),axis.text = element_text(size=20,face='bold'),legend.key.height=unit(0.45,'inch'),
      legend.text=element_text(size=20), legend.title=element_blank())& 
  scale_y_continuous(labels = scales::percent_format()) & 
  labs(
    title = "Cell Type Composition by Distance Group",
    x = "",
    y = "cell proportion",
    fill = "Cell Type"
  ) &
  theme_minimal() 

meta$cell_type <- as.factor(meta$cell_type)

polychrome = c('#E31A1C','#00CC52','#E67E22',"black","grey")


options(repr.plot.width=9, repr.plot.height=6)

ggplot(meta, aes(x = distance_group, y = Fraction, fill = cell_type)) +
  geom_col(position = "fill") +         
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(size = 20, hjust = -0.1),
    axis.text = element_text(size = 20, face = 'bold'),
    axis.text.x = element_text(angle = 45, hjust = 1), 
       axis.title.y = element_text(size = 20, face = 'bold'),
    legend.key.height = unit(0.45, 'inch'),
    legend.text = element_text(size = 20),
    legend.title = element_blank()
  ) + 
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = polychrome) + 
  labs(
    title = "Cell Type Composition by Distance Group",
    x = "",
    y = "cell proportion"
  )
