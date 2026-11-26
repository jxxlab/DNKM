library(Seurat);library(Matrix);suppressMessages(library(dplyr));library(ggplot2);#suppressMessages(library(monocle3));suppressMessages(library(monocle));
suppressMessages(library(cowplot));suppressMessages(library('pagoda2'));suppressMessages(library("velocyto.R"));suppressMessages(library("patchwork"))
library("RColorBrewer");library(scales);suppressMessages(library(SeuratDisk));suppressMessages(library(SeuratWrappers));suppressMessages(library(SeuratData));
options(repr.plot.width=7, repr.plot.height=6); 
setwd("/home/ahmu/YLL/Xenium/dist")
getwd();.libPaths();#setwd('/sdb/xxjiang/data');getwd();

length(dir('/home/ahmu/R/x86_64-conda-linux-gnu-library/4.4'))
length(dir('/home/ahmu/anaconda3/envs/r443/lib/R/library'))
packageVersion("Seurat")

xenium.obj1 = readRDS( file = "rctd/ashleyint/B10-2s_24675.rds")

xenium.obj2 = readRDS( file = "rctd/ashleyint/B5s_56872.rds")

xenium.obj3 = readRDS( file = "rctd/ashleyint/B10-1s_16830.rds")

xenium.obj=merge(xenium.obj1,c(xenium.obj2,xenium.obj3));xenium.obj

Idents(xenium.obj)="typex1"
levels(xenium.obj)=names(table(Idents(xenium.obj)))[order( as.numeric(substr(names(table(Idents(xenium.obj))), 1, 2)))]
names(table(Idents(xenium.obj)))
xenium.obj$typex1=Idents(xenium.obj)

Idents(xenium.obj)="cx1"
levels(xenium.obj)=names(table(Idents(xenium.obj)))[order( as.numeric(substr(names(table(Idents(xenium.obj))), 1, 2)))]
names(table(Idents(xenium.obj)))
xenium.obj$cx1=Idents(xenium.obj)

table(xenium.obj$typex1,xenium.obj$Section)

table(xenium.obj$cx1,xenium.obj$Section)

Idents(xenium.obj)='Section'
Idents(xenium.obj, cells = WhichCells(xenium.obj,idents = "B10-1")) <- "Individual 3"
Idents(xenium.obj, cells = WhichCells(xenium.obj,idents = "B5")) <- "Individual 2"
Idents(xenium.obj, cells = WhichCells(xenium.obj,idents = "B10-2")) <- "Individual 1"
xenium.obj$Section2=Idents(xenium.obj)
table(Idents(xenium.obj))

xenium.obj

saveRDS(xenium.obj,"rctd//ashleyint/B10-2_5_10-1_merge98377.rds",compress = F)



xenium.obj <- SCTransform(xenium.obj, assay = "Xenium")

xenium.obj <- RunPCA(xenium.obj, npcs = 30, verbose = F)

options(repr.plot.width=12, repr.plot.height=6)
ElbowPlot(xenium.obj,ndims = 50) 

xenium.obj <- RunUMAP(xenium.obj, dims = 1:11, verbose = F)

options(repr.plot.width=8, repr.plot.height=6)
DimPlot(xenium.obj,label = T,label.size = 8,raster=FALSE)+ theme(aspect.ratio=1,axis.title=element_blank(),axis.line=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),legend.key.height=unit(0.35,"inch"),legend.text=element_text(size=20,face="plain"))


#Idents(pbmc)='Spatial_snn_res.1'
options(repr.plot.width=10, repr.plot.height=6)
DimPlot(xenium.obj,label = T,label.size = 8,split.by="Section",ncol=3,raster=FALSE)+ theme(aspect.ratio=1,axis.title=element_blank(),axis.line=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),legend.key.height=unit(0.3,"inch"),legend.text=element_text(size=20,face="plain"),text=element_text(size=20,face="plain"))


suppressMessages(library(harmony))
xenium.obj <- RunHarmony(xenium.obj, group.by.vars = "Section",max.iter.harmony = 1)#,max.iter.harmony = 3

options(repr.plot.width=12, repr.plot.height=6)
ElbowPlot(xenium.obj,reduction = "harmony", ndims = 50) 

xenium.obj <- RunUMAP(xenium.obj, reduction = "harmony", dims = 1:11,verbose=F)


polychrome = c("#5A5156", 
               "#E4E1E3", "#F6222E", "#FE00FA", "#16FF32", "#3283FE", 
               "#FEAF16", "#B00068", "#1CFFCE", "#90AD1C", "#2ED9FF", 
               "#DEA0FD", "#AA0DFE", "#F8A19F", "#325A9B", "#C4451C", 
               "#1C8356", "#85660D", "#B10DA1", "#FBE426", "#1CBE4F", 
               "#FA0087", "#FC1CBF", "#F7E1A0", "#C075A6", "#782AB6", 
               "#AAF400", "#BDCDFF", "#822E1C", "#B5EFB5", "#7ED7D1", 
               "#1C7F93", "#D85FF7", "#683B79", "#66B0FF", "#3B00FB")


options(repr.plot.width=8, repr.plot.height=6)
DimPlot(xenium.obj,label = T,label.size = 8,group.by = 'typex1',raster=FALSE)+ theme(aspect.ratio=1,axis.title=element_blank(),axis.line=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),legend.key.height=unit(0.35,"inch"),legend.text=element_text(size=20,face="plain"))


polychrome = c('#E31A1C','#00CC52','#E67E22',"black",
       "#F693BF","#C89999","#97BAD3","#7077D1",
       "#89C1F5","#95E1FF","#7F6874", "#B31D89","#BBCCBB")

Idents(xenium.obj)='typex1'
new.cluster.ids <- c(1:13)
names(x = new.cluster.ids) <- levels(x = xenium.obj)
xenium.obj <- RenameIdents(object = xenium.obj, new.cluster.ids, reorder.numeric = T)
xenium.obj$cx1=Idents(xenium.obj)  

xenium.obj1=xenium.obj
Idents(xenium.obj1)="cx1"
options(repr.plot.width=9, repr.plot.height=6)
p1=DimPlot(xenium.obj1,label = T,label.size = 8,raster=FALSE,repel=F, label.color = "white")+ 
theme(aspect.ratio = 1,axis.line=element_blank(),axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),legend.text=element_text(size=22,face="plain"),legend.key.height=unit(0.35,"inch"),text=element_text(size=22,face="plain"),strip.text=element_text(size=22)) + 
scale_colour_manual(values=polychrome[as.numeric(names(table(Idents(xenium.obj1))))],labels =rownames(table(xenium.obj1$typex1))) +guides(col = guide_legend(override.aes=list(size=3),ncol =1)) 
p1


### xenium.obj1=xenium.obj
Idents(xenium.obj1)="cx1"
options(repr.plot.width=16, repr.plot.height=6)
p1=DimPlot(xenium.obj1,label = T,label.size = 8,split.by = 'Section2',ncol=3,raster=FALSE,repel=F,label.color = 'white')+ 
theme(aspect.ratio = 1,axis.line=element_blank(),axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),legend.text=element_text(size=22,face="plain"),legend.key.height=unit(0.35,"inch"),text=element_text(size=22,face="plain"),strip.text=element_text(size=22)) + 
scale_colour_manual(values=polychrome[as.numeric(names(table(Idents(xenium.obj1))))],labels =rownames(table(xenium.obj1$typex1))) +guides(col = guide_legend(override.aes=list(size=3),ncol =1)) 
p1


names(xenium.obj@meta.data)

DimPlot(xenium.obj1, label = T, label.size = 8, raster = FALSE, repel = F, label.color = "white") + 
  theme(legend.position = "bottom",
        aspect.ratio = 1,
        axis.line = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.text = element_text(size = 22, face = "bold"),
        legend.key.height = unit(0.35, "inch"),
        text = element_text(size = 22, face = "plain"),
        strip.text = element_text(size = 22)) + 
  scale_colour_manual(values = polychrome[as.numeric(names(table(Idents(xenium.obj1))))],
                      labels = rownames(table(xenium.obj1$typex1))) +
  guides(col = guide_legend(override.aes = list(size = 3),
                            nrow = 1))  



library(Seurat);library(Matrix);suppressMessages(library(dplyr));library(ggplot2);#suppressMessages(library(monocle3));suppressMessages(library(monocle));
suppressMessages(library(cowplot));suppressMessages(library('pagoda2'));suppressMessages(library("velocyto.R"));suppressMessages(library("patchwork"))
library("RColorBrewer");library(scales);suppressMessages(library(SeuratDisk));suppressMessages(library(SeuratWrappers));suppressMessages(library(SeuratData));
options(repr.plot.width=7, repr.plot.height=6); 
setwd("/home/ahmu/YLL/Xenium/dist")
getwd();.libPaths();#setwd('/sdb/xxjiang/data');getwd();
length(dir('/home/ahmu/R/x86_64-conda-linux-gnu-library/4.4'))
length(dir('/home/ahmu/anaconda3/envs/r443/lib/R/library'))
packageVersion("Seurat")

xenium.obj=readRDS('rctd//ashleyint/B10-2_5_10-1_merge98377.rds')

DefaultAssay(xenium.obj)='Xenium'

xenium.obj <- NormalizeData(xenium.obj)

xenium.obj

xenium.obj <- JoinLayers(xenium.obj)

Idents(xenium.obj)='cx1'
pbmc.markers <- FindAllMarkers(object = xenium.obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25);dim(pbmc.markers)
write.csv(pbmc.markers, "DEG-98377-typex1.csv")

pbmc.markers %>% group_by(cluster) %>% top_n(n = 10,  wt = avg_log2FC) -> top10

pbmcx <- ScaleData(xenium.obj,  features = rownames(xenium.obj))

options(repr.plot.width=4, repr.plot.height=6)  
DoHeatmap(pbmcx, features = as.character(top10$gene),cells=(sample(colnames(pbmcx),1000)),
          group.colors = c("#5A5156","#E4E1E3", "#F6222E", "#FE00FA")) + NoLegend() +theme(text=element_text( size=16, face="italic"))

library(RColorBrewer)

options(repr.plot.width=8, repr.plot.height=12)  
DoHeatmap(pbmcx, 
          features = as.character(top10$gene),
          cells = (sample(colnames(pbmcx),10000)),
          group.colors = c("#FF6666","#442288","#87CEEB","#90EE90")) + 
  NoLegend() +
  theme(text=element_text(size=16, face="italic")) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0)


cells_metadata <- data.frame(  cell_id = Cells(pbmcx),  cluster = Idents(pbmcx))


set.seed(42) #
selected_cells <- lapply(split(cells_metadata$cell_id, cells_metadata$cluster), function(x) {
  if (length(x) <= 500) {    return(x) # 
  } else {    return(sample(x, 500)) # 
  }})

selected_cells <- unlist(selected_cells)
pbmcx_subset <- subset(pbmcx, cells = selected_cells)
table(Idents(pbmcx_subset))

library(RColorBrewer)

options(repr.plot.width=6.5, repr.plot.height=20)  
DoHeatmap(pbmcx_subset, 
          features = as.character(top10$gene),
          cells = (sample(colnames(pbmcx_subset),1000)),
          group.colors = y) + 
  NoLegend() +
  theme(text=element_text(size=16, face="italic")) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0)

pbmc.markers %>% group_by(cluster) %>% top_n(n = 5,  wt = avg_log2FC) -> top10

library(RColorBrewer)

options(repr.plot.width=6.5, repr.plot.height=10)  
DoHeatmap(pbmcx_subset, 
          features = as.character(top10$gene),
          cells = (sample(colnames(pbmcx_subset),1000)),
          group.colors = y) + 
  NoLegend() +
  theme(text=element_text(size=16, face="italic")) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0)

y <- c('#E31A1C','#00CC52','#E67E22',"black",
       "#F693BF","#C89999","#97BAD3","#7077D1",
       "#89C1F5","#95E1FF","#7F6874", "#B31D89","#BBCCBB")
options(repr.plot.width=16, repr.plot.height=3)
barplot(rep(2,13),col=y)






library(plyr); library(dplyr);library(scales)
meta=as.data.frame(table(xenium.obj$Section2,xenium.obj$typex1)) #
colnames(meta)=c("Section","Type","Cells")
meta <- ddply(meta, "Section", transform,  Fraction=Cells/sum(Cells)*100,Fraction2= percent(Cells/sum(Cells), accuracy = 0.1))
head(meta);dim(meta)



Idents(xenium.obj)='cx1'
options(repr.plot.width=8.5, repr.plot.height=5.5)
p1=ggplot(data=meta,aes(x =Section, y = Fraction,fill=Type))+  geom_bar(stat="identity", position="stack",color="NA",width=0.8,size=0.25)+
geom_text(aes( label = ifelse(Fraction >= 1, round(Fraction,digits =1),"")),position=position_stack(vjust = 0.5),size=6)   +
theme(  panel.background = element_blank(),
       panel.border = element_rect(color = "black", fill = NA, size = 1),  
      plot.title=element_text(size=20,hjust=-0.1),axis.text = element_text(size=20,face='bold'),legend.key.height=unit(0.45,'inch'),
       axis.title.y = element_text(size = 18, face = 'bold',), 
      legend.text=element_text(size=20), legend.title=element_blank())+ labs( title ="",x="",y="cell proportion (%)")+
scale_fill_manual(values= polychrome)+
scale_y_continuous(expand = c(0, 1))#+scale_y_continuous(  labels =scales::percent,breaks = seq(0, 1, 0.25)).Proportion Niches[as.numeric(names(table(Idents(xenium.obj))))]
p1


options(repr.plot.width=7, repr.plot.height=6)
ggplot(data=meta, aes(x = Section, y = Fraction, fill = Type)) + 
  geom_bar(stat = "identity", position = "stack", color = "NA", width = 0.8, size = 0.25) +
  geom_text(aes(label = ifelse(Fraction >= 1, round(Fraction, digits = 1), "")), 
            position = position_stack(vjust = 0.5), size = 6,color='white') +
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    plot.title = element_text(size = 20, hjust = -0.1),
    axis.text.x = element_text(size = 20, face = 'bold', angle = 45, hjust = 1, vjust = 1),  
    axis.text.y = element_text(size = 20, face = 'bold'),                                  
    legend.key.height = unit(0.45, 'inch'),
    axis.title.y = element_text(size = 18, face = 'bold'),
    legend.text = element_text(size = 20),
    legend.title = element_blank()
  ) +
  labs(title = "", x = "", y = "cell proportion (%)") +
  scale_fill_manual(values = polychrome) +
  scale_y_continuous(expand = c(0, 1))









