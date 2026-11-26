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



xenium.obj <- BuildNicheAssay(object = xenium.obj2, fov = "fov.8", group.by = "typex1",
    niches.k = 6, neighbors.k = 5)

options(repr.plot.width=16, repr.plot.height=5)
celltype.plot <- ImageDimPlot(xenium.obj, group.by = "typex1", size = 1.5, cols = "polychrome",
    dark.background = F) + ggtitle("Cell type")
niche.plot <- ImageDimPlot(xenium.obj, group.by = "niches", size = 1.5, dark.background = F,flip_xy = F) + ggtitle("Niches") #+
    #scale_fill_manual(values = c("#442288", "#6CA2EA", "#B5D33D", "#FED23F", "#EB7D5B"))
celltype.plot | niche.plot

options(repr.plot.width=16, repr.plot.height=5)
celltype.plot <- ImageDimPlot(xenium.obj, group.by = "typex1", size = 1, cols = "polychrome",
    dark.background = F) + ggtitle("Cell type")+theme(legend.title=element_blank(),text =element_text(face = "bold",size=16))
niche.plot <- ImageDimPlot(xenium.obj, group.by = "niches", size = 1, dark.background = F,flip_xy = F) + ggtitle("Regions") +
    scale_fill_manual(values = c(  "black", "grey50" ,"grey80","#6CA2EA", "#FED23F","red"))+
theme(legend.title=element_blank(),text =element_text(face = "bold",size=16))
celltype.plot | niche.plot

table(xenium.obj$typex1, xenium.obj$niches)

Idents(xenium.obj)='niches'
xenium.obj$niches_6=Idents(xenium.obj)

table(xenium.obj$niches)

Idents(xenium.obj)='niches'
#Idents(xenium.obj, cells = WhichCells(xenium.obj,idents = "4")) <- "4x"
Idents(xenium.obj, cells = WhichCells(xenium.obj,idents = c(1,3,4,5))) <- "3x"
Idents(xenium.obj, cells = WhichCells(xenium.obj,idents = c(2))) <- "2x"
Idents(xenium.obj, cells = WhichCells(xenium.obj,idents = c(6))) <- "1x"
xenium.obj$niches_6=Idents(xenium.obj)
table(Idents(xenium.obj))

new.cluster.ids <- 1:3
names(x = new.cluster.ids) <- levels(x = xenium.obj)
xenium.obj <- RenameIdents(object = xenium.obj, new.cluster.ids, reorder.numeric = T)
xenium.obj$niches_6=Idents(xenium.obj)  

table(xenium.obj$typex1, xenium.obj$niches_6)

options(repr.plot.width=16, repr.plot.height=5)
celltype.plot <- ImageDimPlot(xenium.obj, group.by = "typex1", size = 1, cols = "polychrome",
    dark.background = F) + ggtitle("Cell type")+theme(legend.title=element_blank(),text =element_text(face = "bold",size=16))
niche.plot <- ImageDimPlot(xenium.obj, group.by = "niches_6", size = 1, dark.background = F,flip_xy = F) + ggtitle("Regions") +
    scale_fill_manual(values = c(  "black", "grey50" ,"grey80","#6CA2EA", "#FED23F","red"))+
theme(legend.title=element_blank(),text =element_text(face = "bold",size=16))
celltype.plot | niche.plot



saveRDS(xenium.obj,'rctd/ashleyint/B10-2s_24675.rds',compress=F)

xenium.obj=readRDS('rctd/ashleyint/B10-2s_24675.rds')





Idents(xenium.obj)='typex1'
xenium.obj1=subset(xenium.obj,idents=c("4   EVT","1   dMac","2   dNK","3   dT"));xenium.obj1

library(plyr); library(dplyr);library(scales)
meta=as.data.frame(table(xenium.obj1$niches_6,xenium.obj1$typex1)) #
colnames(meta)=c("Niches","Type","Cells")
meta <- ddply(meta, "Niches", transform,  Fraction=Cells/sum(Cells)*100,Fraction2= percent(Cells/sum(Cells), accuracy = 0.1))
head(meta);dim(meta)



y <- c('#E31A1C','#00CC52','#E67E22',"black",
       "#F693BF","#C89999","#97BAD3","#7077D1",
       "#89C1F5","#95E1FF","#7F6874", "#B31D89","#BBCCBB")
options(repr.plot.width=16, repr.plot.height=3)
barplot(rep(2,13),col=y)

options(repr.plot.width=10, repr.plot.height=5)
p1=ggplot(data=meta,aes(x =Niches, y = Fraction,fill=Type))+  geom_bar(stat="identity", position="stack",color="black",width=0.7,size=0.25)+
geom_text(aes( label = round(Fraction,digits =2)),position= position_stack(vjust = 0.55),vjust=0.6,size=6,colour='white')   +
theme(  panel.background = element_blank(),  #
      axis.title=element_text(size=35),axis.text = element_text(size=25,face='bold'),legend.key.height=unit(0,'inch'),
      legend.text=element_blank(), legend.title=element_blank())+ xlab("")+ ylab("Cells %")+
scale_fill_manual(values=y)

p1


Idents(xenium.obj)='cx1'
options(repr.plot.width=6, repr.plot.height=2.8)
p1=ggplot(data=meta,aes(x =Niches, y = Fraction,fill=Type))+  geom_bar(stat="identity", position="stack",color=NA,width=0.7,size=0.25)+
geom_text(aes( label = ifelse(Fraction >= 2, round(Fraction,digits =1),"")),position=position_stack(vjust = 0.5),size=6, color = "white")   +
coord_flip() +
theme(  panel.background = element_blank(),
       panel.border = element_rect(color = "black", fill = NA, size = 1),  
     # panel.grid.major = element_line(color = "gray90", size = 0.5),  
      plot.title=element_text(size=20,hjust=-0.1,face = 'bold'),axis.text = element_text(size=20,face='bold'),legend.key.height=unit(0.3,'inch'),
      axis.title.x = element_text(size = 18, face = 'bold'), 
      legend.text=element_text(size=20), legend.title=element_blank())+ labs( title ="Individual 1",x="",y="cell proportion (%)")+
scale_fill_manual(values= y)+
scale_x_discrete(limits = rev(levels(meta$Niches))) +  
scale_y_continuous(expand = c(0, 0))  












