library(Seurat);library(Matrix);suppressMessages(library(dplyr));library(ggplot2);#suppressMessages(library(monocle3));suppressMessages(library(monocle));
suppressMessages(library(cowplot));suppressMessages(library('pagoda2'));suppressMessages(library("velocyto.R"));suppressMessages(library("patchwork"))
library("RColorBrewer");library(scales);suppressMessages(library(SeuratDisk));suppressMessages(library(SeuratWrappers));suppressMessages(library(SeuratData));
options(repr.plot.width=7, repr.plot.height=6); 
setwd("/home/ahmu/YLL/Xenium/dist")
getwd();.libPaths();#setwd('/sdb/xxjiang/data');getwd();

length(dir('/home/ahmu/R/x86_64-conda-linux-gnu-library/4.4'))
length(dir('/home/ahmu/anaconda3/envs/r443/lib/R/library'))
packageVersion("Seurat")

library(future)
options(future.globals.maxSize = 1.0 * 1e12)



xenium.obj = readRDS( file = "/home/ahmu/YLL/Xenium/rds/DEC_240873_C19.rds")

Idents(xenium.obj)="C19"
options(repr.plot.width=8, repr.plot.height=6)
DimPlot(xenium.obj, cols = "polychrome", size = 0.75,label=T,label.size = 8,raster=FALSE)& theme(aspect.ratio=1,axis.text=element_blank(),axis.ticks=element_blank(),axis.title=element_blank(),axis.line=element_blank(),plot.title =element_text(face = "bold.italic",size=20),legend.text =element_text(size=20))

DefaultAssay(xenium.obj)='Xenium'
xenium.obj

xenium.obj=JoinLayers(xenium.obj);xenium.obj

names(xenium.obj@images)

library(spacexr)

query.counts <- GetAssayData(xenium.obj, assay = "Xenium", slot = "counts")[, Cells(xenium.obj[["fov.8"]])]
coords <- GetTissueCoordinates(xenium.obj[["fov.8"]], which = "centroids")
saveRDS(coords,"rctd/ashleyint/B10-2coo.rds")
rownames(coords) <- coords$cell
coords$cell <- NULL
query <- SpatialRNA(coords, query.counts, colSums(query.counts))

# allen.corted.ref can be downloaded here:
# https://www.dropbox.com/s/cuowvm4vrf65pvq/allen_cortex.rds?dl=1
allen.cortex.ref <- readRDS("/home/ahmu/YLL/Xenium/dist/other/human-early-all-120761.rds")
#DefaultAssay(allen.cortex.ref)='RNA'
allen.cortex.ref <- UpdateSeuratObject(allen.cortex.ref)
Idents(allen.cortex.ref) <- "typex1"

allen.cortex.ref

#remove CR cells because there aren't enough of them for annotation
allen.cortex.ref <- subset(allen.cortex.ref)#, subset = subclass != "CR")
counts <- GetAssayData(allen.cortex.ref, assay = "RNA", slot = "counts")
head(counts)

class(counts)

cluster <- as.factor(allen.cortex.ref$typex1)
names(cluster) <- colnames(allen.cortex.ref)
nUMI <- allen.cortex.ref$nCount_RNA
names(nUMI) <- colnames(allen.cortex.ref)
nUMI <- colSums(counts)
#levels(cluster) <- gsub("/", "-", levels(cluster))
reference <- Reference(counts, cluster, nUMI)



common_genes <- intersect(rownames(query@counts), rownames(reference@counts))
length(common_genes)


counts_query <- query@counts[common_genes, ]   
coords_query <- query@coords                    
nUMI_query <- query@nUMI                      

query_sub <- SpatialRNA(coords_query, counts_query, nUMI_query)

counts_ref <- reference@counts[common_genes, ]   
cell_types_ref <- reference@cell_types           
nUMI_ref <- reference@nUMI                        

reference_sub <- Reference(counts_ref, cell_types_ref, nUMI_ref)

table(allen.cortex.ref$typex1)

Sys.setenv(OPENBLAS_NUM_THREADS = 1) 

RCTD <- create.RCTD(query_sub, reference_sub, max_cores = 10, CELL_MIN_INSTANCE=18, DOUBLET_THRESHOLD = 18)
RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")



 head(RCTD@results$results_df)

annotations.df <- RCTD@results$results_df
annotations <- annotations.df$first_type
names(annotations) <- rownames(annotations.df)
head(annotations)

xenium.obj$predicted.celltype <- annotations
keep.cells <- Cells(xenium.obj)[!is.na(xenium.obj$predicted.celltype)]
B2s <- subset(xenium.obj, cells = keep.cells)
xenium.obj
B2s

xenium.obj<- AddMetaData(xenium.obj, metadata = RCTD@results$results_df)

saveRDS(B2s,'rctd/ashleyint/B10-2s_24675.rds',compress = F)

saveRDS(xenium.obj,'rctd/ashleyint/ALL-B10-2.rds',compress = F)

names(xenium.obj@meta.data)

#Idents(pbmc)='Spatial_snn_res.1'
options(repr.plot.width=12, repr.plot.height=8)
DimPlot(B2s,label = T,label.size = 8,group.by = 'predicted.celltype',raster=FALSE)+ theme(aspect.ratio=1,axis.title=element_blank(),axis.line=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),legend.key.height=unit(0.3,"inch"),legend.text=element_text(size=20,face="plain"),text=element_text(size=20,face="plain"))


table(B2s$predicted.celltype)

Idents(B2s)="predicted.celltype"
options(repr.plot.width=12, repr.plot.height=7)
DimPlot(B2s, cols = "polychrome", size = 0.75,label=T,label.size = 8,repel = T,raster=FALSE)& theme(aspect.ratio=1,axis.text=element_blank(),axis.ticks=element_blank(),axis.title=element_blank(),axis.line=element_blank(),plot.title =element_text(face = "bold.italic",size=20),legend.text =element_text(size=20))


table(xenium.obj$Section)

Idents(B2s)="predicted.celltype"
options(repr.plot.width=18, repr.plot.height=12)
DimPlot(B2s, cols = "polychrome",split.by = 'Section',ncol = 4, size = 0.75,label=T,repel = T,label.size = 8,raster=FALSE)& theme(aspect.ratio=1,axis.text=element_blank(),axis.ticks=element_blank(),axis.title=element_blank(),axis.line=element_blank(),plot.title =element_text(face = "bold.italic",size=20),legend.text =element_text(size=20))

Idents(xenium.obj)='predicted.celltype'
options(repr.plot.width=14, repr.plot.height=12)
ImageDimPlot(xenium.obj,fov = 'fov.8',cols = c(rep("grey30",34),Seurat::DiscretePalette(n = 36, palette = "polychrome")[c(35)],rep("grey10",200)), size = 1)& theme(legend.title=element_blank(),text =element_text(face = "bold",size=20))


Idents(xenium.obj)='predicted.celltype'
options(repr.plot.width=14, repr.plot.height=12)
ImageDimPlot(xenium.obj,fov = 'fov.8',cols = c(rep("grey30",17),Seurat::DiscretePalette(n = 36, palette = "polychrome")[c(18:20)],rep("grey10",200)), size = 1)& theme(legend.title=element_blank(),text =element_text(face = "bold",size=20))


Idents(xenium.obj)='predicted.celltype'
options(repr.plot.width=14, repr.plot.height=12)
ImageDimPlot(xenium.obj,fov = 'fov.8',cols = c(rep("grey30",21),Seurat::DiscretePalette(n = 36, palette = "polychrome")[c(22)],rep("grey10",200)), size = 1)& theme(legend.title=element_blank(),text =element_text(face = "bold",size=20))


Idents(xenium.obj)='predicted.celltype'
options(repr.plot.width=14, repr.plot.height=12)
ImageDimPlot(xenium.obj,fov = 'fov.8',cols = c(rep("grey30",11),Seurat::DiscretePalette(n = 36, palette = "polychrome")[c(12:14)],rep("grey10",200)), size = 1)& theme(legend.title=element_blank(),text =element_text(face = "bold",size=20))


Idents(xenium.obj)='predicted.celltype'
options(repr.plot.width=14, repr.plot.height=12)
ImageDimPlot(xenium.obj,fov = 'fov.8',cols = c(rep("grey30",4),Seurat::DiscretePalette(n = 36, palette = "polychrome")[c(5,6)],rep("grey10",200)), size = 1)& theme(legend.title=element_blank(),text =element_text(face = "bold",size=20))


Idents(xenium.obj)='predicted.celltype'
options(repr.plot.width=14, repr.plot.height=12)
ImageDimPlot(xenium.obj,fov = 'fov.8',cols = c(rep("grey30",28),Seurat::DiscretePalette(n = 36, palette = "polychrome")[c(29,30)],rep("grey10",200)), size = 1)& theme(legend.title=element_blank(),text =element_text(face = "bold",size=20))


Idents(xenium.obj)='predicted.celltype'
options(repr.plot.width=14, repr.plot.height=12)
ImageDimPlot(xenium.obj,fov = 'fov.8',cols = c(rep("grey30",14),Seurat::DiscretePalette(n = 36, palette = "polychrome")[c(15:17)],rep("grey10",200)), size = 1)& theme(legend.title=element_blank(),text =element_text(face = "bold",size=20))


Idents(xenium.obj)='predicted.celltype'
options(repr.plot.width=14, repr.plot.height=12)
ImageDimPlot(xenium.obj,fov = 'fov.8',cols = c(rep("grey30",31),Seurat::DiscretePalette(n = 36, palette = "polychrome")[c(32:34)],rep("grey10",200)), size = 1)& theme(legend.title=element_blank(),text =element_text(face = "bold",size=20))


Idents(xenium.obj)='predicted.celltype'
options(repr.plot.width=14, repr.plot.height=12)
ImageDimPlot(xenium.obj,fov = 'fov.8',cols = c(rep("grey30",35),Seurat::DiscretePalette(n = 36, palette = "polychrome")[c(36)],rep("grey10",200)), size = 1)& theme(legend.title=element_blank(),text =element_text(face = "bold",size=20))


ls()

rm('allen.cortex.ref','annotations','annotations.df',
'B2s','cell_types_ref','counts',
'counts_query',
'counts_ref',
'K_val','keep.cells','N_X',
'nUMI',
'nUMI_query',
'nUMI_ref',
'Q_mat',
'query',
'query_sub',
'query.counts',
'RCTD',
'reference',
'S_mat',
'SQ_mat',
'X_vals')
gc()





library(Seurat);library(Matrix);suppressMessages(library(dplyr));library(ggplot2);#suppressMessages(library(monocle3));suppressMessages(library(monocle));
suppressMessages(library(cowplot));suppressMessages(library('pagoda2'));suppressMessages(library("velocyto.R"));suppressMessages(library("patchwork"))
library("RColorBrewer");library(scales);suppressMessages(library(SeuratDisk));suppressMessages(library(SeuratWrappers));suppressMessages(library(SeuratData));
options(repr.plot.width=7, repr.plot.height=6); 
setwd("/home/ahmu/YLL/Xenium/dist")
getwd();.libPaths();#setwd('/sdb/xxjiang/data');getwd();

length(dir('/home/ahmu/R/x86_64-conda-linux-gnu-library/4.4'))
length(dir('/home/ahmu/anaconda3/envs/r443/lib/R/library'))
packageVersion("Seurat")

xenium.obj = readRDS( file = "rctd/ashleyint/B10-2s_24675.rds")



Idents(xenium.obj)='predicted.celltype'
xenium.obj2=xenium.obj
Idents(xenium.obj2, cells = WhichCells(xenium.obj2, idents = c('12 dMac1','13 dMac2','14 dMac3'))) <- '1   dMac'
Idents(xenium.obj2, cells = WhichCells(xenium.obj2, idents = c('18 dNK1','19 dNK2','20 dNK3','21 dNK_p'))) <- '2   dNK'
Idents(xenium.obj2, cells = WhichCells(xenium.obj2, idents = c('22 dT'))) <- '3   dT'
Idents(xenium.obj2, cells = WhichCells(xenium.obj2, idents = c('34 CC','35 iEVT','36 eEVT','37 GC'))) <- '4   EVT'
Idents(xenium.obj2, cells = WhichCells(xenium.obj2, idents = c('31 CTB'))) <- '5   CTB'
Idents(xenium.obj2, cells = WhichCells(xenium.obj2, idents = c('32 FC','33 STB'))) <- '6   STB'
Idents(xenium.obj2, cells = WhichCells(xenium.obj2, idents = c('29 pMes1','30 pMes2'))) <- '7   pMes'
Idents(xenium.obj2, cells = WhichCells(xenium.obj2, idents = c('1   dS1','2   dS2','3   dS3','4   SMC'))) <- '8   dS'
Idents(xenium.obj2, cells = WhichCells(xenium.obj2, idents = c('5   dPV1','6   dPV2'))) <- '9   dPV'
Idents(xenium.obj2, cells = WhichCells(xenium.obj2, idents = c('7   Epi1','8   Epi2'))) <- '10 Epi'
Idents(xenium.obj2, cells = WhichCells(xenium.obj2, idents = c('9   Endo1','10 Endo2','11 pEndo'))) <- '11 Endo'
Idents(xenium.obj2, cells = WhichCells(xenium.obj2, idents = c('15 HBC'))) <- '12 HBC'
Idents(xenium.obj2, cells = WhichCells(xenium.obj2, idents = c('16 dDC1','17 dDC2','23 NK_CD16+','24 B',
                                                               '25 Plasma','26 MO','27 Granulocytes','28 ILC3'))) <- '13 other immune cell'

xenium.obj2$typex1=Idents(xenium.obj2)
table(Idents(xenium.obj2))

Idents(xenium.obj2)="typex1"
levels(xenium.obj2)=names(table(Idents(xenium.obj2)))[order( as.numeric(substr(names(table(Idents(xenium.obj2))), 1, 2)))]
names(table(Idents(xenium.obj2)))
xenium.obj2$typex1=Idents(xenium.obj2)

table(xenium.obj2$typex1)

xenium.obj=xenium.obj2;xenium.obj

new.cluster.ids <- c(1:13)
names(x = new.cluster.ids) <- levels(x = xenium.obj)
xenium.obj <- RenameIdents(object = xenium.obj, new.cluster.ids, reorder.numeric = T)
xenium.obj$cx1=Idents(xenium.obj)  



saveRDS(xenium.obj,'rctd/ashleyint/B10-2s_24675.rds',compress=F)

options(repr.plot.width=10, repr.plot.height=5)
celltype.plot <- ImageDimPlot(xenium.obj2, group.by = "typex1", size = 1.5, cols = "polychrome",
    dark.background = F) + ggtitle("Cell type")

celltype.plot

DefaultBoundary(xenium.obj[["fov.8"]])="centroid"

y <- c('#E31A1C','#00CC52','#E67E22',"black",
       "#F693BF","#C89999","#97BAD3","#7077D1",
       "#89C1F5","#95E1FF","#7F6874", "#B31D89","#BBCCBB")
options(repr.plot.width=16, repr.plot.height=3)
barplot(rep(2,13),col=y)

Idents(xenium.obj)="typex1"
options(repr.plot.width=12, repr.plot.height=6)
ImageDimPlot(xenium.obj, cols = y,border.size = 0, na.value = "grey90", dark.background =F, border.color = "black", 
             mols.size =10,size =0.5, nmols = 10000,flip_xy = F)& theme(legend.title=element_blank(),text =element_text(face = "bold",size=20))

Idents(xenium.obj)="typex1"
options(repr.plot.width=12, repr.plot.height=6)
ImageDimPlot(xenium.obj, cols = y,border.size = 0, na.value = "grey90", dark.background =F, border.color = "black", 
             mols.size =10,size =0.5, nmols = 10000,flip_xy = F)& theme(legend.title=element_blank(),text =element_text(face = "bold",size=20))

Idents(xenium.obj)="typex1"
options(repr.plot.width=12, repr.plot.height=8)
ImageDimPlot(xenium.obj, cols =y,  border.size = 0.01,dark.background = F,coord.fixed = T,flip_xy = F)+   theme(legend.title = element_blank(),legend.text =element_text(size=16), legend.key.height=unit(0.2,"inch"),   legend.position = "bottom") +
  guides(fill = guide_legend(ncol = 7, override.aes = list(size = 7)),color= guide_legend(ncol =7, override.aes = list(size = 7)))  



DefaultBoundary(xenium.obj[["fov.8"]]) <- "segmentation"

options(repr.plot.width=12, repr.plot.height=8)
ImageDimPlot(xenium.obj, fov = "fov.8", cols = c(y[1:4], rep("grey90", 100)),
             border.color = "black", border.size = 0.02,
             dark.background = FALSE, coord.fixed = TRUE,flip_xy = F)

options(repr.plot.width=10, repr.plot.height=6)
ImageDimPlot(xenium.obj, fov = "fov.8", cols = c(y[1:4], rep("grey90", 100)),
             border.color = "black", border.size = 0.02,
             dark.background = FALSE, coord.fixed = TRUE,flip_xy = F)

DefaultBoundary(xenium.obj[["fov.8"]]) <- "centroid"

options(repr.plot.width=12, repr.plot.height=6)
ImageDimPlot(xenium.obj, fov = "fov.8", cols = c(y[1:4], rep("grey90", 100)),
             border.color = "black", border.size = 0.02,
             dark.background = FALSE, coord.fixed = TRUE,flip_xy = FALSE)






