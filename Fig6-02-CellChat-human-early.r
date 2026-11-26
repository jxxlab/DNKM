library(Seurat);library(Matrix);suppressMessages(library(dplyr));library(ggplot2);#suppressMessages(library(monocle3));suppressMessages(library(monocle));
suppressMessages(library(cowplot));suppressMessages(library('pagoda2'));suppressMessages(library("velocyto.R"));suppressMessages(library("patchwork"))
library("RColorBrewer");library(scales);suppressMessages(library(SeuratDisk));suppressMessages(library(SeuratWrappers));suppressMessages(library(SeuratData));
options(repr.plot.width=7, repr.plot.height=6); 
setwd("/home/ahmu/YLL/Xenium/dist/")
getwd();.libPaths();#setwd('/sdb/xxjiang/data');getwd();

packageVersion('CellChat')

library(CellChat)
library(ggplot2)
library(ggalluvial)
library(NMF)
library(circlize)

pbmc=readRDS("other/human-early-35277.rds");pbmc;names(pbmc@meta.data)

DefaultAssay(pbmc)='RNA'
pbmc

Idents(pbmc)="typex"
names(table(Idents(pbmc)))
Idents(pbmc, cells = WhichCells(pbmc,idents = '20 EVT')) <- '4 EVT'
Idents(pbmc, cells = WhichCells(pbmc,idents = '11 dT')) <- '3 dT'
Idents(pbmc, cells = WhichCells(pbmc,idents = '10 dNK')) <- '2 dNK'
Idents(pbmc, cells = WhichCells(pbmc,idents = '7   dMac')) <- '1 dMac'
pbmc$typex=Idents(pbmc)
table(Idents(pbmc))


Normal.cellchat=pbmc

Normal.cellchat$active.ident <- Normal.cellchat$typex

Normal.cellchat@active.ident <- factor(Normal.cellchat$typex)

normal.data.input <- Normal.cellchat@assays$RNA$data 

Normal.cellchat

normal.meta <- Normal.cellchat@meta.data

dim(normal.data.input);head(normal.meta)[,1:4]

Normal.cellchat <- createCellChat(object = Normal.cellchat, meta = normal.meta, group.by = "active.ident")

Normal.cellchat <- addMeta(Normal.cellchat, meta = normal.meta)

Normal.cellchat

Normal.cellchat <- setIdent(Normal.cellchat, ident.use = "active.ident")

levels(Normal.cellchat@idents)





showDatabaseCategory(CellChatDB.human)

#set the ligand-receptor interaction database
Normal.cellchatDB <- CellChatDB.human

dplyr::glimpse(Normal.cellchatDB$interaction)

Normal.cellchatDB$interaction$annotation


Normal.cellchat@DB <- Normal.cellchatDB

Normal.cellchat <- subsetData(Normal.cellchat)



Normal.cellchat <- identifyOverExpressedGenes(Normal.cellchat)

Normal.cellchat <- identifyOverExpressedInteractions(Normal.cellchat)

#1) compute the communication probability and infer the cellular communication network
Normal.cellchat <- computeCommunProb(Normal.cellchat)

#2) infer the cell-cell communication at a signaling pathway level
Normal.cellchat <- computeCommunProbPathway(Normal.cellchat)



#3) calculate the aggregated cell-cell communication network
Normal.cellchat <- aggregateNet(Normal.cellchat)

saveRDS(Normal.cellchat, file = "rctd/ashleyint/cc-scsn.rds")

Normal.cellchat <- readRDS("rctd/ashleyint/cc-scsn.rds")



sum(Normal.cellchat@net$count)  

table(Normal.cellchat@idents)

Normal.cellchat@netP$pathways

groupSize <- as.numeric(table(Normal.cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(Normal.cellchat@net$count,
                 vertex.weight=groupSize,
                 weight.scale=T,
                 label.edge=F,                  
                 title.name="Number_of_interactions"
                )
netVisual_circle(Normal.cellchat@net$weight,                  
                 vertex.weight = groupSize,                  
                 weight.scale = TRUE, 
                  label.edge = FALSE,                  
                 title.name = "Interaction weights/strength"
                )



x=c("#5A5156","#E4E1E3", "#F6222E", "#FE00FA"   )
options(repr.plot.width=16, repr.plot.height=3)
barplot(rep(2,4),col=x)


options(repr.plot.width=12, repr.plot.height=6)
netVisual_aggregate(
  Normal.cellchat,
  signaling = "TGFb", #pathway_name
  vertex.receiver = c(1,2,3), #
  layout = "chord" #
)


options(repr.plot.width=12, repr.plot.height=6)
extractEnrichedLR(Normal.cellchat, signaling = "WNT", geneLR.return = FALSE)
netVisual_individual(
  Normal.cellchat,
  signaling = "WNT",
  vertex.receiver = c(1, 3, 5),
  pairLR.use = 'WNT2_FZD6_LRP6', #
  layout = "chord"
)


options(repr.plot.width=12, repr.plot.height=6)
netVisual_aggregate(
  Normal.cellchat,
  signaling = "WNT", #pathway_name
  vertex.receiver = c(1,2,3), 
  layout = "chord",#
    color.use=x
)





x=c("#5A5156","#E4E1E3", "#F6222E", "#FE00FA"   )
options(repr.plot.width=16, repr.plot.height=3)
barplot(rep(2,4),col=x)


Normal.cellchat@netP$pathways

options(repr.plot.width=8, repr.plot.height=15)
netVisual_bubble(
  Normal.cellchat,
  #sources.use = 8,
  #targets.use = c(1,2,3),
  remove.isolate = FALSE
)


bubble_output <- netVisual_bubble(
  Normal.cellchat,
  # sources.use = 8,
  # targets.use = c(1,2,3),
  remove.isolate = FALSE
)

 names(bubble_output)

df_bubble <- bubble_output$data

write.csv(df_bubble, file = "cellchat_bubble_data.csv", row.names = FALSE)

options(repr.plot.width=6, repr.plot.height=5)
netVisual_bubble(
  Normal.cellchat,
  sources.use = 8,
  targets.use = c(1,2,3),
  signaling = "NOTCH",
  remove.isolate = FALSE
)


options(repr.plot.width=6, repr.plot.height=5)
netVisual_bubble(
  Normal.cellchat,
  # sources.use = 2,
  # targets.use = c(1, 3, 5),
  pairLR.use = extractEnrichedLR(Normal.cellchat, signaling = "NOTCH"),
  remove.isolate = FALSE
)





