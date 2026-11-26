library(Seurat);library(Matrix);suppressMessages(library(dplyr));library(ggplot2);#suppressMessages(library(monocle3));suppressMessages(library(monocle));
suppressMessages(library(cowplot));suppressMessages(library('pagoda2'));suppressMessages(library("velocyto.R"));suppressMessages(library("patchwork"))
library("RColorBrewer");library(scales);suppressMessages(library(SeuratDisk));suppressMessages(library(SeuratWrappers));suppressMessages(library(SeuratData));
options(repr.plot.width=7, repr.plot.height=6); 
setwd("/home/ahmu/YLL/Xenium")
getwd();.libPaths();

length(dir('/home/ahmu/R/x86_64-conda-linux-gnu-library/4.4'))
length(dir('/home/ahmu/anaconda3/envs/r443/lib/R/library'))
packageVersion("Seurat")

packageVersion("Seurat")
LoadXenium5=function (data.dir, fov = "fov", assay = "Xenium") 
{
    data <- ReadXenium(data.dir = data.dir, type = c("centroids", 
        "segmentations"), )
    segmentations.data <- list(centroids = CreateCentroids(data$centroids), 
        segmentation = CreateSegmentation(data$segmentations))
    coords <- CreateFOV(coords = segmentations.data, type = c("segmentation", 
        "centroids"), molecules = data$microns, assay = assay)
    xenium.obj <- CreateSeuratObject(counts = data$matrix[["Gene Expression"]], 
        assay = assay)
    if ("Blank Codeword" %in% names(data$matrix)) 
        xenium.obj[["BlankCodeword"]] <- CreateAssayObject(counts = data$matrix[["Blank Codeword"]])
    else xenium.obj[["BlankCodeword"]] <- CreateAssayObject(counts = data$matrix[["Unassigned Codeword"]])
    xenium.obj[["ControlCodeword"]] <- CreateAssayObject(counts = data$matrix[["Negative Control Codeword"]])
    xenium.obj[["ControlProbe"]] <- CreateAssayObject(counts = data$matrix[["Negative Control Probe"]])
    xenium.obj[[fov]] <- coords
    return(xenium.obj)
}



# Load the Xenium data
M1 <- LoadXenium5("/home/ahmu/YLL/Xenium-DEC/AM-20251118-0500/output-XETG00166", fov = "fov")


m1=read.csv("/home/ahmu/YLL/Xenium-DEC/AM-20251118-0500/Section/B5/B5_cells_stats.csv",check.names = FALSE, header=TRUE,row.names=1);dim(m2);head(m2,2)
m2=read.csv("/home/ahmu/YLL/Xenium-DEC/AM-20251118-0500/Section/B10-1/B10-1_cells_stats.csv",check.names = FALSE, header=TRUE,row.names=1);dim(m6);head(m6,2)
m3=read.csv("/home/ahmu/YLL/Xenium-DEC/AM-20251118-0500/Section/B10-2/B10-2_cells_stats.csv",check.names = FALSE, header=TRUE,row.names=1);dim(m7);head(m7,2)


X1 <- subset(M1, cells = rownames(m1));X1
X2 <- subset(M2, cells = rownames(m2));X2
X3 <- subset(M3, cells = rownames(m3));X3

library(future)
options(future.globals.maxSize = 1.0 * 1e12)

options(repr.plot.width=20, repr.plot.height=6)
p1 <- ImageFeaturePlot(X1, fov ="fov", features = c("nCount_Xenium"), max.cutoff="q95")#+scale_y_reverse()+scale_x_reverse()
p2 <- ImageFeaturePlot(X2, fov ="fov", features = c("nCount_Xenium"), max.cutoff="q95")
p3 <- ImageFeaturePlot(X3, fov ="fov", features = c("nCount_Xenium"), max.cutoff="q95")

p1|p2|p3

X1;X2;X3

saveRDS(X1, file = "rds/B5_61591.rds",compress = F)
saveRDS(X2, file = "rds/B10-1_19042.rds",compress = F)
saveRDS(X3, file = "rds/B10-2_25802.rds",compress = F)

ls()

xenium.obj <- merge(X1,                     y = c(X2, X3),
  add.cell.ids = c("B5","B10-1", "B10-2"));xenium.obj 

saveRDS(xenium.obj , file = "rds/DEC.rds",compress = F)


# remove cells with 0 counts
xenium.obj <- subset(xenium.obj, subset = nCount_Xenium > 0);xenium.obj
median(xenium.obj$nFeature_Xenium)

names(xenium.obj@images)

DefaultFOV(xenium.obj) <- "fov"

median(xenium.obj$nFeature_Xenium)

DefaultBoundary(xenium.obj[["fov"]])

options(repr.plot.width=6.5, repr.plot.height=7)
VlnPlot(xenium.obj, features = c("nFeature_Xenium", "nCount_Xenium"), ncol = 2, pt.size = 0)& theme(axis.title=element_blank(),axis.text =element_text(face = "bold",size=20),text =element_text(face = "bold",size=20))

names(xenium.obj@meta.data)

xenium.obj$Section <- gsub("_.*", "", colnames(xenium.obj))

names(xenium.obj@meta.data)

Idents(xenium.obj)="Section"
levels(xenium.obj)=names(table(Idents(xenium.obj)))[order( as.numeric(substr(names(table(Idents(xenium.obj))), 1, 2)))]
names(table(Idents(xenium.obj)))
xenium.obj$Section=Idents(xenium.obj)

table(xenium.obj$Section)

options(repr.plot.width=12, repr.plot.height=7)
VlnPlot(xenium.obj, features = c("nFeature_Xenium", "nCount_Xenium"), group.by = "Section", ncol = 2, pt.size = 0)& theme(axis.title=element_blank(),axis.text =element_text(face = "bold",size=20),text =element_text(face = "bold",size=20))


saveRDS(xenium.obj , file = "rds/DEC.rds",compress = F)


median(xenium.obj$nFeature_Xenium)
mean(xenium.obj$nFeature_Xenium)

library(dplyr)
xenium.obj@meta.data %>%
  group_by(Section) %>%
  summarise(median_nFeature = median(nFeature_Xenium))

xenium.obj@meta.data %>%
  group_by(Section) %>%
  summarise(median_nFeature = mean(nFeature_Xenium))


table(xenium.obj$orig.ident)
#table(xenium.obj$nCount_Xenium)
#table(xenium.obj$nFeature_Xenium)
table(xenium.obj$nCount_BlankCodeword)
table(xenium.obj$nFeature_BlankCodeword)
table(xenium.obj$nCount_ControlCodeword)
table(xenium.obj$nFeature_ControlCodeword)
table(xenium.obj$nCount_ControlProbe)
table(xenium.obj$nFeature_ControlProbe)

ImageFeaturePlot(xenium.obj, fov =names(xenium.obj@images), features = c("nFeature_Xenium"), max.cutoff="q95")

options(repr.plot.width=10, repr.plot.height=15)
ImageFeaturePlot(xenium.obj,fov =names(xenium.obj@images),features = c("nFeature_Xenium"),combine = TRUE, max.cutoff="q95",size = 0.5#, cols = c("white", "red")
                )&  theme(legend.title=element_blank(),plot.title = element_text(face = "bold.italic",size=16),text= element_text(face = "plain",size=16))

names(table(xenium.obj$Section))

library(patchwork)
options(repr.plot.width = 12, repr.plot.height = 12)


fovs <- names(xenium.obj@images)


plot_list <- lapply(fovs, function(f) {
  cells_in_fov <- Cells(xenium.obj[[f]])
  sample_name <- unique(xenium.obj$Section[cells_in_fov])[1]  #
  
  ImageFeaturePlot(xenium.obj, fov = f, features = "nFeature_Xenium", 
                   max.cutoff = "q95") +
    ggtitle(sample_name) +
    theme(legend.title = element_blank(),
          plot.title = element_text(face = "bold.italic", size = 16),
          text = element_text(face = "plain", size = 16))
})

wrap_plots(plot_list, ncol = 3)

library(patchwork)
options(repr.plot.width = 12, repr.plot.height = 12)


fovs <- names(xenium.obj@images)

plot_list <- lapply(fovs, function(f) {
  cells_in_fov <- Cells(xenium.obj[[f]])
  sample_name <- unique(xenium.obj$Section[cells_in_fov])[1]  
  ImageFeaturePlot(xenium.obj, fov = f, features = "nFeature_Xenium", 
                   max.cutoff = "q95",dark.background = FALSE) +
    ggtitle(sample_name) +
    theme(legend.title = element_blank(),
          plot.title = element_text(face = "bold.italic", size = 16),
          text = element_text(face = "plain", size = 16))
})

wrap_plots(plot_list, ncol = 3)

xenium.obj@images

Idents(xenium.obj)="Section"
options(repr.plot.width=60, repr.plot.height=40)
P1=ImageDimPlot(xenium.obj , fov = names(xenium.obj@images),split.by = 'Section', molecules = c("HLA-G", "KRT8", "CD68", "NCAM1"), nmols = 20000,axes = TRUE,crop=TRUE)& theme(legend.title=element_blank(),text =element_text(face = "bold.italic",size=20))
P2=ImageDimPlot(xenium.obj , fov = names(xenium.obj@images), molecules = c("HLA-G", "KRT8", "CD68", "NCAM1"),split.by = 'Section', nmols = 20000,mols.size = 0.0001,crop=TRUE)& theme(legend.title=element_blank(),text =element_text(face = "bold.italic",size=20))
P1
P2

unique(xenium.obj$Section)[2]

#Idents(xenium.obj)="SeuratObject"
Idents(xenium.obj)="Section"
options(repr.plot.width=15, repr.plot.height=15)
ImageDimPlot(xenium.obj, fov = "fov", molecules = c("HLA-G", "KRT8", "PLAC8", "NCAM1"),mols.size = 0.001, nmols = 20000)& theme(legend.title=element_blank(),text =element_text(face = "bold.italic",size=20))

#Idents(xenium.obj)="SeuratObject"
Idents(xenium.obj)='Section'
options(repr.plot.width=15, repr.plot.height=12)
ImageDimPlot(xenium.obj, fov = "fov", molecules = c("CD68","CD3E","HLA-G","NCAM1"),split.by = 'Section', crop=TRUE, nmols = 20000,cols = "polychrome")& theme(legend.title=element_blank(),text =element_text(face = "bold.italic",size=20))
ImageDimPlot(xenium.obj, fov = "fov.2", molecules = c("CD68", "CD3E","HLA-G", "NCAM1"),split.by = 'Section',crop=TRUE, nmols = 20000,cols = "polychrome")& theme(legend.title=element_blank(),text =element_text(face = "bold.italic",size=20))
ImageDimPlot(xenium.obj, fov = "fov.3", molecules = c("CD68", "CD3E","HLA-G", "NCAM1"),split.by = 'Section',crop=TRUE, nmols = 20000,cols = "polychrome")& theme(legend.title=element_blank(),text =element_text(face = "bold.italic",size=20))


xenium.obj <- SCTransform(xenium.obj, assay = "Xenium")
xenium.obj <- RunPCA(xenium.obj, npcs = 30, verbose = F)


options(repr.plot.width=12, repr.plot.height=6)
ElbowPlot(xenium.obj,ndims = 50) 

sct_models <- xenium.obj[["SCT"]]@SCTModel.list
umi_assays <- sapply(sct_models, function(m) m@umi.assay)
table(umi_assays)

xenium.obj <- RunUMAP(xenium.obj, dims = 1:22, verbose = F)


options(repr.plot.width=8, repr.plot.height=6)
DimPlot(xenium.obj,label = T,label.size = 8,raster=FALSE)+ theme(aspect.ratio=1,axis.title=element_blank(),axis.line=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),legend.key.height=unit(0.35,"inch"),legend.text=element_text(size=20,face="plain"))

options(repr.plot.width=10, repr.plot.height=12)
DimPlot(xenium.obj,label = T,label.size = 8,split.by="Section",ncol=3,raster=FALSE)+ theme(aspect.ratio=1,axis.title=element_blank(),axis.line=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),legend.key.height=unit(0.3,"inch"),legend.text=element_text(size=20,face="plain"),text=element_text(size=20,face="plain"))



table(xenium.obj$Section)

xenium.obj

saveRDS(xenium.obj , file = "rds/DEC.rds",compress = F)


suppressMessages(library(harmony))
xenium.obj <- RunHarmony(xenium.obj, group.by.vars = "Section",max.iter.harmony = 10)#,max.iter.harmony = 3

options(repr.plot.width=12, repr.plot.height=6)
ElbowPlot(xenium.obj,reduction = "harmony", ndims = 50) 

xenium.obj <- RunUMAP(xenium.obj, reduction = "harmony", dims = 1:26,verbose=F)


options(repr.plot.width=8, repr.plot.height=6)
DimPlot(xenium.obj,label = T,label.size = 8,raster=FALSE)+ theme(aspect.ratio=1,axis.title=element_blank(),axis.line=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),legend.key.height=unit(0.35,"inch"),legend.text=element_text(size=20,face="plain"))

Idents(xenium.obj)='Section'
options(repr.plot.width=10, repr.plot.height=10)
DimPlot(xenium.obj,label = T,label.size = 8,split.by="Section",ncol=3,raster=FALSE)+ theme(aspect.ratio=1,axis.title=element_blank(),axis.line=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),legend.key.height=unit(0.3,"inch"),legend.text=element_text(size=20,face="plain"),text=element_text(size=20,face="plain"))

xenium.obj <- FindNeighbors(xenium.obj, reduction = "harmony",  dims = 1:26) 
xenium.obj <- FindClusters(object = xenium.obj, resolution =0.3) 
options(repr.plot.width=8, repr.plot.height=6)
DimPlot(xenium.obj,label = T,label.size = 8,raster=FALSE)+ theme(aspect.ratio=1,axis.title=element_blank(),axis.line=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),legend.key.height=unit(0.35,"inch"),legend.text=element_text(size=20,face="plain"))

options(repr.plot.width=8, repr.plot.height=6)
DimPlot(xenium.obj,label = T,label.size = 8,raster=FALSE)+ theme(aspect.ratio=1,axis.title=element_blank(),axis.line=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),legend.key.height=unit(0.3,"inch"),legend.text=element_text(size=20,face="plain"))

#Idents(xenium.obj)="seurat_clusters"
xenium.obj$C19=Idents(xenium.obj)

Idents(xenium.obj)="C19"
options(repr.plot.width=16, repr.plot.height=8)
DimPlot(xenium.obj,cells.highlight = CellsByIdentities(object = xenium.obj, idents =c(11)),raster=FALSE, cols.highlight =rev(c('red','green','blue')),pt.size=0.00001,sizes.highlight=0.5,reduction="umap", label=F,label.size = 8) & theme(legend.key.height=unit(0.3,"inch"),axis.line=element_blank(),axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),text=element_text(size=20),legend.text=element_text(size=20,face = "plain"))

#Idents(xenium.obj)="SeuratObject"
options(repr.plot.width=20, repr.plot.height=20)
ImageDimPlot(xenium.obj, fov = "fov", cols = "polychrome", size = 1)& theme(legend.title=element_blank(),text =element_text(face = "bold",size=20))
ImageDimPlot(xenium.obj, fov = "fov.2", cols = "polychrome", size = 1)& theme(legend.title=element_blank(),text =element_text(face = "bold",size=20))
ImageDimPlot(xenium.obj, fov = "fov.3", cols = "polychrome", size = 1)& theme(legend.title=element_blank(),text =element_text(face = "bold",size=20))




Idents(xenium.obj)="seurat_clusters"
options(repr.plot.width=8, repr.plot.height=6)
DimPlot(xenium.obj, cols = "polychrome", size = 0.75,label=T,label.size = 8)& theme(aspect.ratio=1,axis.text=element_blank(),axis.ticks=element_blank(),axis.title=element_blank(),axis.line=element_blank(),plot.title =element_text(face = "bold.italic",size=20),legend.text =element_text(size=20))

Idents(xenium.obj)="seurat_clusters"
options(repr.plot.width=8, repr.plot.height=6)
DimPlot(xenium.obj, cols = "polychrome", size = 0.75,label=T,label.size =8 ,raster=FALSE)& theme(aspect.ratio=1,axis.text=element_blank(),axis.ticks=element_blank(),axis.title=element_blank(),axis.line=element_blank(),plot.title =element_text(face = "bold.italic",size=20),legend.text =element_text(size=20))


options(repr.plot.width=10, repr.plot.height=10)
DimPlot(xenium.obj, cols = "polychrome", size = 0.75,label = T,label.size = 8,split.by="Section",ncol=3,raster=FALSE)+ theme(aspect.ratio=1,axis.title=element_blank(),axis.line=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),legend.key.height=unit(0.3,"inch"),legend.text=element_text(size=20,face="plain"),text=element_text(size=20,face="plain"))

#判断低质量细胞群是否存在
options(repr.plot.width=20, repr.plot.height=4.5)
FeaturePlot(xenium.obj, c("nFeature_Xenium", "HLA-G", "VIM", "PTPRC"), raster=FALSE,min.cutoff = "q01", max.cutoff = "q90",ncol=4) & theme(aspect.ratio=1,axis.text=element_blank(),axis.ticks=element_blank(),axis.title=element_blank(),axis.line=element_blank(),plot.title =element_text(face = "bold.italic",size=20),legend.text =element_text(size=16))

options(repr.plot.width=20, repr.plot.height=4.5)
FeaturePlot(xenium.obj, c("VWF", "CGB3", "EGFR", "PECAM1"),raster=FALSE, min.cutoff = "q01", max.cutoff = "q90",ncol=4) & theme(aspect.ratio=1,axis.text=element_blank(),axis.ticks=element_blank(),axis.title=element_blank(),axis.line=element_blank(),plot.title =element_text(face = "bold.italic",size=20),legend.text =element_text(size=16))

options(repr.plot.width=20, repr.plot.height=4.5)
FeaturePlot(xenium.obj, c("ALAS2", "KRT7", "ITGA5", "PSG1"),raster=FALSE, min.cutoff = "q01", max.cutoff = "q90",ncol=4) & theme(aspect.ratio=1,axis.text=element_blank(),axis.ticks=element_blank(),axis.title=element_blank(),axis.line=element_blank(),plot.title =element_text(face = "bold.italic",size=20),legend.text =element_text(size=16))

options(repr.plot.width=20, repr.plot.height=4.5)
FeaturePlot(xenium.obj, c("CLEC9A", "ACTA2", "PECAM1", "VWF"),raster=FALSE, min.cutoff = "q01", max.cutoff = "q90",ncol=4) & theme(aspect.ratio=1,axis.text=element_blank(),axis.ticks=element_blank(),axis.title=element_blank(),axis.line=element_blank(),plot.title =element_text(face = "bold.italic",size=20),legend.text =element_text(size=16))

options(repr.plot.width=20, repr.plot.height=4.5)
FeaturePlot(xenium.obj, c("FOLR2","CDH1","CD19","CD34"),raster=FALSE, min.cutoff = "q01", max.cutoff = "q90",ncol=4) & theme(aspect.ratio=1,axis.text=element_blank(),axis.ticks=element_blank(),axis.title=element_blank(),axis.line=element_blank(),plot.title =element_text(face = "bold.italic",size=20),legend.text =element_text(size=16))

options(repr.plot.width=20, repr.plot.height=4.5)
FeaturePlot(xenium.obj, c("NKG7", "CD68", "DCN", "CD3E"),raster=FALSE, min.cutoff = "q01", max.cutoff = "q90",ncol=4) & theme(aspect.ratio=1,axis.text=element_blank(),axis.ticks=element_blank(),axis.title=element_blank(),axis.line=element_blank(),plot.title =element_text(face = "bold.italic",size=20),legend.text =element_text(size=16))

options(repr.plot.width=20, repr.plot.height=4.5)
FeaturePlot(xenium.obj, c("VIM", "LUM", "IGFBP1", "PSG1"),raster=FALSE, min.cutoff = "q01", max.cutoff = "q90",ncol=4) & theme(aspect.ratio=1,axis.text=element_blank(),axis.ticks=element_blank(),axis.title=element_blank(),axis.line=element_blank(),plot.title =element_text(face = "bold.italic",size=20),legend.text =element_text(size=16))

options(repr.plot.width=20, repr.plot.height=4.5)
FeaturePlot(xenium.obj, c("EGFR", "TRPV6", "ERVFRD-1", "CGB3"),raster=FALSE, min.cutoff = "q01", max.cutoff = "q90",ncol=4) & theme(aspect.ratio=1,axis.text=element_blank(),axis.ticks=element_blank(),axis.title=element_blank(),axis.line=element_blank(),plot.title =element_text(face = "bold.italic",size=20),legend.text =element_text(size=16))

options(repr.plot.width=20, repr.plot.height=8)
FeaturePlot(xenium.obj, c("HLA-G", "VIM", "PTPRC", "KRT7"),split.by='Section',raster=FALSE, min.cutoff = "q01", max.cutoff = "q90",ncol=4) & theme(aspect.ratio=1,axis.text=element_blank(),axis.ticks=element_blank(),axis.title=element_blank(),axis.line=element_blank(),plot.title =element_text(face = "bold.italic",size=20),legend.text =element_text(size=16))

saveRDS(xenium.obj , file = "rds/DEC_C19.rds",compress = F)


xenium.obj =readRDS( file = "rds/DEC_240873_C19.rds")


