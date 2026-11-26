library(Seurat);library(Matrix);suppressMessages(library(dplyr));library(ggplot2);#suppressMessages(library(monocle3));suppressMessages(library(monocle));
suppressMessages(library(cowplot));suppressMessages(library('pagoda2'));suppressMessages(library("velocyto.R"));suppressMessages(library("patchwork"))
library("RColorBrewer");library(scales);suppressMessages(library(SeuratDisk));suppressMessages(library(SeuratWrappers));suppressMessages(library(SeuratData));
options(repr.plot.width=7, repr.plot.height=6); 
setwd("/home/ahmu/YLL/Xenium/dist/")
getwd();.libPaths();#setwd('/sdb/xxjiang/data');getwd();

length(dir('/home/ahmu/R/x86_64-conda-linux-gnu-library/4.4'))
length(dir('/home/ahmu/anaconda3/envs/r443/lib/R/library'))
packageVersion("Seurat")

xenium.obj = readRDS( file = "rctd/ashleyint/B10-2s_24675.rds")

options(repr.plot.width=16, repr.plot.height=5)
celltype.plot <- ImageDimPlot(xenium.obj, group.by = "typex", size = 1, cols = "polychrome",
    dark.background = F) + ggtitle("Cell type")+theme(legend.title=element_blank(),text =element_text(face = "bold",size=16))
celltype.plot2 <- ImageDimPlot(xenium.obj, group.by = "typex1", size = 1, cols = "polychrome",
    dark.background = F) + ggtitle("Cell type")+theme(legend.title=element_blank(),text =element_text(face = "bold",size=16))
celltype.plot | celltype.plot2

DefaultAssay(xenium.obj)='Xenium'

Idents(xenium.obj)='typex1'##
options(repr.plot.width=8, repr.plot.height=5)
ImageDimPlot(xenium.obj,dark.background = T,flip_xy = F,cols = c(rep("grey0",0),Seurat::DiscretePalette(n = 36, palette = "polychrome")[c(1:4)],rep("grey0",200)), size = 1, axes =T)& 
theme(legend.title=element_blank(),text =element_text(face = "bold",size=20))& ggtitle("B10-2 Cell type")


coords <- GetTissueCoordinates(xenium.obj[["fov.8"]])
range(coords$x)   # 
range(coords$y)   #


xenium.obj

Idents(xenium.obj)='typex1'
names(table(Idents(xenium.obj)))

Idents(xenium.obj)='typex1'
xenium.obj2=xenium.obj
Idents(xenium.obj2) = '5   other'
Idents(xenium.obj2, cells = WhichCells(xenium.obj, idents = c('4   EVT'))) <- '4   EVT'
Idents(xenium.obj2, cells = WhichCells(xenium.obj, idents = c('3   dT'))) <- '3   dT'
Idents(xenium.obj2, cells = WhichCells(xenium.obj, idents = c('2   dNK'))) <- '2   dNK'
Idents(xenium.obj2, cells = WhichCells(xenium.obj, idents = c('1   dMac'))) <- '1   dMac'
xenium.obj2$typex2=Idents(xenium.obj2)
table(Idents(xenium.obj2))



xenium.obj=xenium.obj2

options(future.globals.maxSize = 1.0 * 1e12)
cropped.coords <- Crop(xenium.obj[["fov.8"]], x = c(2500, 3000), y = c(1500, 2000), coords = "plot")
xenium.obj[["zoom1"]] <- cropped.coords
# visualize cropped area with cell segmentations & selected molecules
DefaultBoundary(xenium.obj[["zoom1"]]) <- "segmentation" 

DefaultBoundary(xenium.obj[["fov.8"]]) <- "segmentation"  
Idents(xenium.obj)='typex2'
options(repr.plot.width=8, repr.plot.height=6)
ImageDimPlot(xenium.obj, fov = "zoom1",flip_xy = F, size=2,mols.size = 0.3,axes =F, border.color = "black", border.size = 0.03,
             cols = c('#E31A1C','#00CC52','#E67E22',
                      "black","grey99"), dark.background =F,
    coord.fixed =T,na.value = "grey90", nmols = 20000) & 
theme(legend.title=element_blank(),text =element_text(face = "bold.italic",size=20)) #&scale_colour_manual(values=y)




xenium.obj

options(future.globals.maxSize = 1.0 * 1e12)
cropped.coords <- Crop(xenium.obj[["fov.8"]], x = c(2600, 3100), y = c(2200, 2700), coords = "plot")
xenium.obj[["zoom2"]] <- cropped.coords
# visualize cropped area with cell segmentations & selected molecules
DefaultBoundary(xenium.obj[["zoom2"]]) <- "segmentation"  

DefaultBoundary(xenium.obj[["fov.8"]]) <- "segmentation"  
Idents(xenium.obj)='typex2'
options(repr.plot.width=8, repr.plot.height=6)
ImageDimPlot(xenium.obj, fov = "zoom2",flip_xy = F, size=2,mols.size = 0.3,axes =F, border.color = "black", border.size = 0.03,
             cols = c('#E31A1C','#00CC52','#E67E22',
                      "black","grey99"), dark.background =F,
    coord.fixed =T,na.value = "grey90", nmols = 20000) & 
theme(legend.title=element_blank(),text =element_text(face = "bold.italic",size=20)) 






options(future.globals.maxSize = 1.0 * 1e12)
cropped.coords <- Crop(xenium.obj[["fov.8"]], x = c(2700, 3200), y = c(3700, 4200), coords = "plot")
xenium.obj[["zoom3"]] <- cropped.coords
# visualize cropped area with cell segmentations & selected molecules
DefaultBoundary(xenium.obj[["zoom3"]]) <- "segmentation"

DefaultBoundary(xenium.obj[["fov.8"]]) <- "segmentation"  
Idents(xenium.obj)='typex2'
options(repr.plot.width=8, repr.plot.height=6)
ImageDimPlot(xenium.obj, fov = "zoom3",flip_xy = F, size=2,mols.size = 0.3,axes =F, border.color = "black", border.size = 0.03,
             cols = c('#E31A1C','#00CC52','#E67E22',
                      "black","grey99"), dark.background =F,
    coord.fixed =T,na.value = "grey90", nmols = 20000) & 
theme(legend.title=element_blank(),text =element_text(face = "bold.italic",size=20)) #&scale_colour_manual(values=y)




options(future.globals.maxSize = 1.0 * 1e12)
cropped.coords <- Crop(xenium.obj[["fov.8"]], x = c(2000, 2500), y = c(1200, 1700), coords = "plot")
xenium.obj[["zoom4"]] <- cropped.coords
# visualize cropped area with cell segmentations & selected molecules
DefaultBoundary(xenium.obj[["zoom4"]]) <- "segmentation"

DefaultBoundary(xenium.obj[["fov.8"]]) <- "segmentation" 
Idents(xenium.obj)='typex2'
options(repr.plot.width=8, repr.plot.height=6)
ImageDimPlot(xenium.obj, fov = "zoom4",flip_xy = F, size=2,mols.size = 0.3,axes =F, border.color = "black", border.size = 0.03,
             cols = c('#E31A1C','#00CC52','#E67E22',
                      "black","grey99"), dark.background =F,
    coord.fixed =T,na.value = "grey90", nmols = 20000) & 
theme(legend.title=element_blank(),text =element_text(face = "bold.italic",size=20)) #&scale_colour_manual(values=y)




options(future.globals.maxSize = 1.0 * 1e12)
cropped.coords <- Crop(xenium.obj[["fov.8"]], x = c(1900, 2400), y = c(3000, 3500), coords = "plot")
xenium.obj[["zoom5"]] <- cropped.coords
# visualize cropped area with cell segmentations & selected molecules
DefaultBoundary(xenium.obj[["zoom5"]]) <- "segmentation"

DefaultBoundary(xenium.obj[["fov.8"]]) <- "segmentation" 
Idents(xenium.obj)='typex2'
options(repr.plot.width=8, repr.plot.height=6)
ImageDimPlot(xenium.obj, fov = "zoom5",flip_xy = F, size=2,mols.size = 0.3,axes =F, border.color = "black", border.size = 0.03,
             cols = c('#E31A1C','#00CC52','#E67E22',
                      "black","grey99"), dark.background =F,
    coord.fixed =T,na.value = "grey90", nmols = 20000) & 
theme(legend.title=element_blank(),text =element_text(face = "bold.italic",size=20)) 




options(future.globals.maxSize = 1.0 * 1e12)
cropped.coords <- Crop(xenium.obj[["fov.8"]], x = c(2000, 2500), y = c(4700, 5200), coords = "plot")
xenium.obj[["zoom6"]] <- cropped.coords
# visualize cropped area with cell segmentations & selected molecules
DefaultBoundary(xenium.obj[["zoom6"]]) <- "segmentation"

DefaultBoundary(xenium.obj[["fov.8"]]) <- "segmentation"  
Idents(xenium.obj)='typex2'
options(repr.plot.width=8, repr.plot.height=6)
ImageDimPlot(xenium.obj, fov = "zoom6",flip_xy = F, size=2,mols.size = 0.3,axes =F, border.color = "black", border.size = 0.03,
             cols = c('#E31A1C','#00CC52','#E67E22',
                      "black","grey99"), dark.background =F,
    coord.fixed =T,na.value = "grey90", nmols = 20000) & 
theme(legend.title=element_blank(),text =element_text(face = "bold.italic",size=20)) #&scale_colour_manual(values=y)




options(future.globals.maxSize = 1.0 * 1e12)
cropped.coords <- Crop(xenium.obj[["fov.8"]],  x = c(1500, 2000), y = c(600, 1100), coords = "plot")
xenium.obj[["zoom7"]] <- cropped.coords
# visualize cropped area with cell segmentations & selected molecules
DefaultBoundary(xenium.obj[["zoom7"]]) <- "segmentation"

DefaultBoundary(xenium.obj[["fov.8"]]) <- "segmentation" 
Idents(xenium.obj)='typex2'
options(repr.plot.width=8, repr.plot.height=6)
ImageDimPlot(xenium.obj, fov = "zoom7",flip_xy = F, size=2,mols.size = 0.3,axes =F, border.color = "black", border.size = 0.03,
             cols = c('#E31A1C','#00CC52','#E67E22',
                      "black","grey99"), dark.background =F,
    coord.fixed =T,na.value = "grey90", nmols = 20000) & 
theme(legend.title=element_blank(),text =element_text(face = "bold.italic",size=20)) 




options(future.globals.maxSize = 1.0 * 1e12)
cropped.coords <- Crop(xenium.obj[["fov.8"]], x = c(1050, 1550), y = c(3200, 3700), coords = "plot")
xenium.obj[["zoom8"]] <- cropped.coords
# visualize cropped area with cell segmentations & selected molecules
DefaultBoundary(xenium.obj[["zoom8"]]) <- "segmentation"

DefaultBoundary(xenium.obj[["fov.8"]]) <- "segmentation" 
Idents(xenium.obj)='typex2'
options(repr.plot.width=8, repr.plot.height=6)
ImageDimPlot(xenium.obj, fov = "zoom8",flip_xy = F, size=2,mols.size = 0.3,axes =F, border.color = "black", border.size = 0.03,
             cols = c('#E31A1C','#00CC52','#E67E22',
                      "black","grey99"), dark.background =F,
    coord.fixed =T,na.value = "grey90", nmols = 20000) & 
theme(legend.title=element_blank(),text =element_text(face = "bold.italic",size=20)) 




options(future.globals.maxSize = 1.0 * 1e12)
cropped.coords <- Crop(xenium.obj[["fov.8"]], x = c(1500, 2000), y = c(4900, 5400), coords = "plot")
xenium.obj[["zoom9"]] <- cropped.coords
# visualize cropped area with cell segmentations & selected molecules
DefaultBoundary(xenium.obj[["zoom9"]]) <- "segmentation"

DefaultBoundary(xenium.obj[["fov.8"]]) <- "segmentation"  
Idents(xenium.obj)='typex2'
options(repr.plot.width=8, repr.plot.height=6)
ImageDimPlot(xenium.obj, fov = "zoom9",flip_xy = F, size=2,mols.size = 0.3,axes =F, border.color = "black", border.size = 0.03,
             cols = c('#E31A1C','#00CC52','#E67E22',
                      "black","grey99"), dark.background =F,
    coord.fixed =T,na.value = "grey90", nmols = 20000) & 
theme(legend.title=element_blank(),text =element_text(face = "bold.italic",size=20)) 











