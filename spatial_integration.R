
## Run with Seurat v3.2.0

#setwd("/Users/u0114327/Dropbox/")
setwd("/Volumes/Transcend/Pascal/")

samplesc = c("sc_S5", "sc_S6")
sampleslide = c("20208_2020_Preview",  "20208_slide011_Preview", "20208_slide016_preview", "20208_slide10_preview"  )

## Keep only the ones for which we have data:
list_files = list.files("ResolveDataToAnalyze/", pattern = "sc")

library("Seurat")

folders = paste0(samplesc)
data_read <- Read10X(data.dir = paste0("ResolveDataToAnalyze/",folders, "/"  ))

### read spatial data
filenames = list.files(pattern = "all.tiff_measurements.csv", recursive = TRUE)

filenames = setdiff( filenames[grepl("[O|o]riginal", filenames)] ,filenames[grepl("[O|o]ld", filenames)])

datalist = lapply(filenames,
                    function(x){read.csv(file = x,
                                         row.names = 1,
                                         stringsAsFactors = FALSE)})

names(datalist) <- c("slide000_Control_A1", "slide000_Control_A2", "slide000_TBI_A1", "slide000_TBI_A2",  "slide011_Control_D2", "slide011_TBI_D1", "slide016_TBI_C1", "slide010_Control_D2", "slide010_TBI_D1")

for(var in names(datalist) ){
 #datalist[[var]][, c("Gene.ID.1", "Gene.ID.2", "X.1", "X.2", "X") ] <- NULL
 colnames(datalist[[var]]) = paste0(var, "zz", colnames(datalist[[var]]))
 datalist[[var]] = datalist[[var]][ , colSums(is.na(datalist[[var]])) == 0]
 datalist[[var]] = datalist[[var]][Reduce(intersect, lapply(datalist, rownames)),]
 datalist[[var]] = datalist[[var]][,colSums(datalist[[var]]) < max(colSums(datalist[[var]])/2 )]
  
}


### read spatial coordinates
filenames = list.files(pattern = "all_coordinates.csv", recursive = TRUE)
data_coord = lapply(filenames,
                    function(x){read.csv(file = x,
                                         stringsAsFactors = FALSE)})
names(data_coord) <- names(datalist)

for(var in names(data_coord) ){

 rownames(data_coord[[var]]) = paste0(var, "zzPathCellObject", ifelse( rownames(data_coord[[var]]) == "0", "", paste0( ".", rownames(data_coord[[var]])) ))

 data_coord[[var]] = data_coord[[var]][colnames(datalist[[var]] ),]
 
}



seurat_obj_list <-  list()

for(var in names(datalist))
{
  
  obj <- CreateSeuratObject(datalist[[var]], min.cells = 0, min.features = 0)
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 5000)
  obj <- NormalizeData(obj)
  obj <- ScaleData(obj, features = unique(c( rownames(obj))) )
  obj = RunPCA(obj,do.print=FALSE, pc.genes = unique(VariableFeatures(obj) ) )
  PCnumber = 15
  res = 0.5
  obj <- FindNeighbors(obj, dims = 1:PCnumber)
  obj <- FindClusters(obj, resolution = res, print.output = 0, save.SNN = T )
  obj <- RunUMAP(obj, dims = 1:PCnumber)
 # obj <- RunTSNE(obj, dims = 1:PCnumber)
  DimPlot(obj)
  
  
  obj@reductions$coord = obj@reductions$umap
  obj@reductions$coord@key = "Coord"
  obj@reductions$coord@cell.embeddings[,1] = data_coord[[var]][rownames(obj@reductions$umap@cell.embeddings),]$X
  obj@reductions$coord@cell.embeddings[,2] = data_coord[[var]][rownames(obj@reductions$umap@cell.embeddings),]$Y
  colnames(obj@reductions$coord@cell.embeddings) <- c("Coord1", "Coord2")
  obj@meta.data$region = data_coord[[var]]$Granular_Layer
  seurat_obj_list[[var]] = obj
  
}


for(var in names(seurat_obj_list))
{
  
  obj <- seurat_obj_list[[var]]
  
  DimPlot(obj, reduction = "coord", group.by = "region")
  
}


commongenes =  Reduce(intersect, lapply(datalist, rownames))
datalistcommon = datalist
for(var in names(datalistcommon) ){
    datalistcommon[[var]] = datalistcommon[[var]][commongenes,]
}


###### integration

data_read = datalistcommon [[names(datalistcommon)[1]]]
for(var in names(datalistcommon)[2:length(names(datalistcommon))] ){
    data_read = cbind( data_read, datalistcommon[[var]])
}


objSpatial <- CreateSeuratObject(data_read, min.cells = 10, min.features = 2)
l=length(colnames(objSpatial))

objSpatial@meta.data$orig.ident = unlist(strsplit(colnames(objSpatial), "zz"))[seq(1,l*2,2)]

objSpatial@meta.data$slide = unlist(strsplit(colnames(data_read), "_"))[seq(1,l*3,3)]
objSpatial@meta.data$condition = unlist(strsplit(colnames(data_read), "_"))[seq(2,l*3,3)]

objSpatial <- FindVariableFeatures(objSpatial, selection.method = "vst", nfeatures = 5000)
objSpatial <- NormalizeData(objSpatial)
objSpatial <- ScaleData(objSpatial, features = unique(c( rownames(objSpatial))) )
objSpatial = RunPCA(objSpatial ,do.print=FALSE, pc.genes = unique(VariableFeatures(objSpatial) ) )
PCnumber = 8
res = 1.5
objSpatial <- FindNeighbors(objSpatial, dims = 1:PCnumber)
objSpatial <- FindClusters(objSpatial, resolution = res, print.output = 0, save.SNN = T )
objSpatial <- RunUMAP(objSpatial, dims = 1:PCnumber)
# obj <- RunTSNE(obj, dims = 1:PCnumber)
DimPlot(objSpatial, label = TRUE)

### removing spatial coordinates if needed
seurat_obj_list$slide000_Control_A2 = seurat_obj_list$slide000_Control_A2[, setdiff(colnames(seurat_obj_list$slide000_Control_A2), "slide000_Control_A2zzPathCellObject.1604")]

AstroNeurontestdata <- as.data.frame(AstroNeurontest@assays$RNA@counts)[intersect(rownames(AstroNeurontest), commongenes),]
obj10x <- CreateSeuratObject(AstroNeurontestdata, min.cells = 100, min.features = 1)

source("SeuratRun.R")

mergeclusters <-SeuratV3MergeSpatial( obj10x@assays$RNA@counts, objSpatial@assays$RNA@counts,  "10x", "Genexyz", number.cc = 10, dimensions.align = 10)
DimPlot(mergeclusters, group.by = "sample")
DimPlot(mergeclusters, label= TRUE)

FeaturePlot(mergeclusters, c("Ascl1", "Aldoc", "Neurod1", "Frzb"))
mergeclusters@meta.data$previousID = "spatial"
mergeclusters@meta.data[intersect(rownames(mergeclusters@meta.data), rownames(AstroNeurontest@meta.data)),]$previousID = AstroNeurontest@meta.data[intersect(rownames(mergeclusters@meta.data), rownames(AstroNeurontest@meta.data)),]$previousIdent
DimPlot(mergeclusters, group.by = "previousID", label = TRUE)
DimPlot(mergeclusters, group.by = "previousID", split.by = "sample", label=TRUE)

###subsetdata
objSpatial1 = SubsetData(objSpatial, cells = intersect( colnames(objSpatial), WhichCells(mergeclusters, idents = setdiff(c(0:11), c(0,1,7,11,5) ) )))

## re-integrate
mergeclusters <-SeuratV3MergeSpatial( obj10x@assays$RNA@counts, objSpatial1@assays$RNA@counts,  "10x", "Genexyz", number.cc = 10, dimensions.align = 10)
DimPlot(mergeclusters, group.by = "sample")
DimPlot(mergeclusters, label= TRUE)


###subsetdata
objSpatial2 = SubsetData(objSpatial, cells = intersect( colnames(objSpatial), WhichCells(mergeclusters, idents = setdiff(c(1:10), c(6,9) ) )))

## re-integrate
mergeclusters <-SeuratV3MergeSpatial( obj10x@assays$RNA@counts, objSpatial2@assays$RNA@counts,  "10x", "Genexyz", number.cc = 15, dimensions.align = 10, res = 2.6)

fdata =FetchData(mergeclusters, c("UMAP_1", "UMAP_2") )
cellstoremove = rownames(fdata[which(fdata$UMAP_1 < 2.7 & fdata$UMAP_1 > -1.5 & fdata$UMAP_2 <3 & fdata$UMAP_2>-0.3),])

DimPlot( mergeclusters, cells.highlight =  cellstoremove )
mergeclusters = SubsetData(mergeclusters, cells = setdiff(colnames(mergeclusters),cellstoremove) )

DimPlot(mergeclusters, group.by = "sample")
DimPlot(mergeclusters, label= TRUE)

fdata =FetchData(mergeclusters, c("UMAP_1", "UMAP_2") )
DimPlot(mergeclusters, label= TRUE)+geom_vline(xintercept = -2)+geom_vline(xintercept = 1.7)+geom_hline(yintercept = -1.7)+geom_hline(yintercept = 1.5)


DefaultAssay(mergeclusters) <- "integrated"

Astro <- SubsetData(mergeclusters, ident.use =  c(4,5,6,1,8))
Neuro<- SubsetData(mergeclusters, ident.use =  c( 4,3,2,7))


Astro.cds <- as.cell_data_set(Astro)
Astro.cds <- cluster_cells(cds = Astro.cds, reduction_method = "UMAP")
Astro.cds <- learn_graph(Astro.cds, use_partition = TRUE)

Neuro.cds <- as.cell_data_set(Neuro)
Neuro.cds <- cluster_cells(cds = Neuro.cds, reduction_method = "UMAP")
Neuro.cds <- learn_graph(Neuro.cds, use_partition = TRUE)

# order cells
Astro.cds <- order_cells(Astro.cds, reduction_method = "UMAP", root_cells = WhichCells(mergeclusters, idents = c(4) ) )
Neuro.cds <- order_cells(Neuro.cds, reduction_method = "UMAP", root_cells = WhichCells(mergeclusters, idents = c(4) ) )


plot_cells(
  cds = Astro.cds,
  color_cells_by = "pseudotime",
  show_trajectory_graph = TRUE
)

plot_cells(
  cds = Neuro.cds,
  color_cells_by = "pseudotime",
  show_trajectory_graph = TRUE
)


mergeclusters <- AddMetaData(
  object = mergeclusters,
  metadata = Astro.cds@principal_graph_aux@listData$UMAP$pseudotime,
  col.name = "AstroLineage"
)


mergeclusters <- AddMetaData(
  object = mergeclusters,
  metadata = Neuro.cds@principal_graph_aux@listData$UMAP$pseudotime,
  col.name = "NeuroLineage"
)



FeaturePlot(mergeclusters, c("AstroLineage", "NeuroLineage"), pt.size = 0.1, cols =  c("yellow2", "red2"))


mergeclusters@meta.data$Astro_lineage_discrete = discretize(mergeclusters@meta.data$AstroLineage, method = "interval")
mergeclusters@meta.data$Neuro_lineage_discrete = discretize(mergeclusters@meta.data$NeuroLineage, method = "interval")


for(var in names(seurat_obj_list))
{
  
  obj <- seurat_obj_list[[var]]
  obj@meta.data$Astro_lineage_discrete = mergeclusters@meta.data[colnames(obj),]$Astro_lineage_discrete
  obj@meta.data$Neuro_lineage_discrete = mergeclusters@meta.data[colnames(obj),]$Neuro_lineage_discrete
  obj@meta.data$AstroLineage = mergeclusters@meta.data[colnames(obj),]$AstroLineage
  obj@meta.data$NeuroLineage = mergeclusters@meta.data[colnames(obj),]$NeuroLineage
  seurat_obj_list[[var]] <- obj
}

i= 1
names(seurat_obj_list)[i]
FeaturePlot(seurat_obj_list[[i]], c("AstroLineage", "NeuroLineage"), reduction = "coord", cols = c("yellow2", "red2"))

#### classify by % in region

library(arules)
i=0
i=i+1
var = names(datalist)[i]
print(var)
obj = seurat_obj_list[[var]]
subset = obj@meta.data[c("Astro_lineage_discrete", "Neuro_lineage_discrete", "region")]

tt = table(subset[, c(1,3)])/ colSums(table(subset[, c(3,1)])) * 100
barplot(t(tt), col = c("red", "green", "blue"), legend.text = rownames(t(tt)), beside = TRUE, names.arg = c("young", "juvenile","adult"), ylim = c(0, 150), main = paste0( var, " AstroLineage") )

tt = table(subset[, c(1,2)])/ colSums(table(subset[, c(2,1)])) * 100
barplot(t(tt), col = c("red", "green", "blue"), legend.text = rownames(t(tt)), beside = TRUE, names.arg = c("young", "juvenile","adult"), ylim = c(0, 150), main = paste0( var, " NeuroLineage") )

### creating the large astrocytic table
datatableAstroLineage = data.frame(region=character(),
Astro_lineage_discrete=character(),
Freq=double(),
slide=character(),
stringsAsFactors=FALSE)

for(var in names(datalist)){
    obj = seurat_obj_list[[var]]
    subset = obj@meta.data[c("Astro_lineage_discrete", "Neuro_lineage_discrete", "region")]
    #subset$Astro_lineage_discrete = discretize(subset$AstroLineage, method = "interval")
   # subset$Neuro_lineage_discrete = discretize(subset$NeuroLineage, method = "interval")
    tt = table(subset[, c(1,3)])/ colSums(table(subset[, c(3,1)])) * 100
    dataframeobj = as.data.frame(tt)
    dataframeobj$slide = var
    datatableAstroLineage = rbind(datatableAstroLineage, dataframeobj)
}
datatableAstroLineage$condition = unlist(strsplit(datatableAstroLineage$slide, "_"))[seq(2,length(datatableAstroLineage$slide)*3,3)]

### plotting astro proportions

datatableAstroLineageGranular = datatableAstroLineage[which(datatableAstroLineage$region == "Granular_Layer" ),]

### creating the big neuronal table
datatableNeuroLineage = data.frame(region=character(),
Neuro_lineage_discrete=character(),
Freq=double(),
slide=character(),
stringsAsFactors=FALSE)

for(var in names(datalist)){
    obj = seurat_obj_list[[var]]
    subset = obj@meta.data[c("Astro_lineage_discrete", "Neuro_lineage_discrete", "region")]
    tt = table(subset[, c(2,3)])/ colSums(table(subset[, c(3,2)])) * 100
    dataframeobj = as.data.frame(tt)
    dataframeobj$slide = var
    datatableNeuroLineage = rbind(datatableNeuroLineage, dataframeobj)
}
datatableNeuroLineage$condition = unlist(strsplit(datatableNeuroLineage$slide, "_"))[seq(2,length(datatableNeuroLineage$slide)*3,3)]


### astrocytes

datatableLineage = datatableAstroLineage
datatableLineageRegionYoung = datatableLineage[which(datatableLineage$region == region & datatableLineage$Astro_lineage_discrete == levels(unique(datatableLineage$Astro_lineage_discrete))[1]),]
datatableLineageRegionJuv = datatableLineage[which(datatableLineage$region == region & datatableLineage$Astro_lineage_discrete == levels(unique(datatableLineage$Astro_lineage_discrete))[2]),]
datatableLineageRegionAd = datatableLineage[which(datatableLineage$region == region & datatableLineage$Astro_lineage_discrete == levels(unique(datatableLineage$Astro_lineage_discrete))[3]),]

###

### plotting
region = "Granular_Layer"
region = "SGZ"
region = "Hilius"

### neurons
datatableLineage = datatableNeuroLineage
datatableLineageRegionYoung = datatableLineage[which(datatableLineage$region == region & datatableLineage$Neuro_lineage_discrete == levels(unique(datatableLineage$Neuro_lineage_discrete))[1]),]
datatableLineageRegionJuv = datatableLineage[which(datatableLineage$region == region & datatableLineage$Neuro_lineage_discrete == levels(unique(datatableLineage$Neuro_lineage_discrete))[2]),]
datatableLineageRegionAd = datatableLineage[which(datatableLineage$region == region & datatableLineage$Neuro_lineage_discrete == levels(unique(datatableLineage$Neuro_lineage_discrete))[3]),]



par(mfrow = c(2, 3))
tab = datatableLineageRegionYoung
p1 = ggplot(tab, aes(x=condition, y=Freq, fill = condition)) +ylim(0,100) + geom_boxplot(position=position_dodge(1))+ geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(1))+theme_classic() + stat_compare_means(method = "wilcox.test")+geom_text(aes(label=slide), hjust=0.4,vjust=-1.2)+ylab("% or total Nes+ cells detected in the region")
p2 = ggplot(tab, aes(x=condition, y=Freq, fill = condition)) +ylim(0,100) + geom_boxplot(position=position_dodge(1))+ geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(1))+theme_classic() + stat_compare_means(method = "wilcox.test")+ylab("% or total Nes+ cells detected in the region")
tab = datatableLineageRegionJuv
p3 = ggplot(tab, aes(x=condition, y=Freq, fill = condition)) +ylim(0,100) + geom_boxplot(position=position_dodge(1))+ geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(1))+theme_classic() + stat_compare_means(method = "wilcox.test")+geom_text(aes(label=slide), hjust=0.4,vjust=-1.2)+ylab("% or total Nes+ cells detected in the region")
p4 = ggplot(tab, aes(x=condition, y=Freq, fill = condition)) +ylim(0,100) + geom_boxplot(position=position_dodge(1))+ geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(1))+theme_classic() + stat_compare_means(method = "wilcox.test")+ylab("% or total Nes+ cells detected in the region")
tab = datatableLineageRegionAd
p5 = ggplot(tab, aes(x=condition, y=Freq, fill = condition)) +ylim(0,100) + geom_boxplot(position=position_dodge(1))+ geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(1))+theme_classic() + stat_compare_means(method = "wilcox.test")+geom_text(aes(label=slide), hjust=0.4,vjust=-1.2)+ylab("% or total Nes+ cells detected in the region")
p6 = ggplot(tab, aes(x=condition, y=Freq, fill = condition)) +ylim(0,100) + geom_boxplot(position=position_dodge(1))+ geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(1))+theme_classic() + stat_compare_means(method = "wilcox.test")+ylab("% or total Nes+ cells detected in the region")

grid.arrange( p2, p4, p6, p1, p3, p5,nrow = 2)




### end of barplots

mergeclusterssubset = SubsetData(mergeclusters, ident.remove = c(0,1, 14, 18,9,8,20,15,12, 19))
mergeclusterssubset = SubsetData(mergeclusters, ident.remove = c(12,7,11,0,1,4))
mergeclusterssubset = SubsetData(mergeclusters, ident.remove = c(6,1,2,12,0,13))
mergeclusterssubset = SubsetData(mergeclusters, ident.remove = c(8,7))
mergeclusterssubset = SubsetData(mergeclusters, ident.remove = c(10,11,8,1,12,2,0))

data_read_subset = data_read[, intersect(colnames(mergeclusterssubset), colnames(data_read))]


objSpatial <- CreateSeuratObject(data_read_subset, min.cells = 10, min.features = 2)
l=length(colnames(objSpatial))

objSpatial@meta.data$orig.ident = unlist(strsplit(colnames(objSpatial), "zz"))[seq(1,l*2,2)]

objSpatial@meta.data$slide = unlist(strsplit(colnames(data_read), "_"))[seq(1,l*3,3)]
objSpatial@meta.data$condition = unlist(strsplit(colnames(data_read), "_"))[seq(2,l*3,3)]

objSpatial <- FindVariableFeatures(objSpatial, selection.method = "vst", nfeatures = 5000)
objSpatial <- NormalizeData(objSpatial)
objSpatial <- ScaleData(objSpatial, features = unique(c( rownames(objSpatial))) )
objSpatial = RunPCA(objSpatial ,do.print=FALSE, pc.genes = unique(VariableFeatures(objSpatial) ) )
PCnumber = 8
res = 0.8
objSpatial <- FindNeighbors(objSpatial, dims = 1:PCnumber)
objSpatial <- FindClusters(objSpatial, resolution = res, print.output = 0, save.SNN = T )
objSpatial <- RunUMAP(objSpatial, dims = 1:PCnumber)
# obj <- RunTSNE(obj, dims = 1:PCnumber)
DimPlot(objSpatial, label = TRUE)


for(var in names(seurat_obj_list))
{
  
  obj <- seurat_obj_list[[var]]
  
  obj@meta.data$clustername = objSpatial@meta.data[colnames(obj),]$seurat_clusters
  obj@meta.data$clustername = objSpatial@meta.data[colnames(obj),]$seurat_clusters
  obj@meta.data$clustername = objSpatial@meta.data[colnames(obj),]$seurat_clusters
  obj@meta.data$clustername = objSpatial@meta.data[colnames(obj),]$seurat_clusters
  obj@meta.data$clustername = objSpatial@meta.data[colnames(obj),]$seurat_clusters

#  obj@meta.data[is.na(obj@meta.data)] <- 0
  seurat_obj_list[[var]] <- obj
}

DimPlot(seurat_obj_list$slide000_TBI_A2, group.by = "clustername" , reduction = "coord" , cols = c("gray", "gray", "red", "blue", "green", "orange"))

####

library(slingshot)
sce <- runSlingshot(mergeclusters, "UMAP", "7", c( "5" ) )
plot(reducedDims(sce)$UMAP, col = "gray", pch=16, asp = 1)
lines(SlingshotDataSet(sce), lwd=2, col='black', type = 'lineages')

mergeclusters@meta.data$slingPseudotime_1 = sce$slingPseudotime_1
mergeclusters@meta.data$slingPseudotime_2 = sce$slingPseudotime_2
mergeclusters@meta.data$slingPseudotime_3 = sce$slingPseudotime_3
mergeclusters@meta.data$slingPseudotime_4 = sce$slingPseudotime_4
mergeclusters@meta.data$slingPseudotime_5 = sce$slingPseudotime_5

####

for(var in names(seurat_obj_list))
{
  
  obj <- seurat_obj_list[[var]]
  
  obj@meta.data$Pseudotime_1 = mergeclusters@meta.data[colnames(obj),]$slingPseudotime_1
  obj@meta.data$Pseudotime_2 = mergeclusters@meta.data[colnames(obj),]$slingPseudotime_2
  obj@meta.data$Pseudotime_3 = mergeclusters@meta.data[colnames(obj),]$slingPseudotime_3
  obj@meta.data$Pseudotime_4 = mergeclusters@meta.data[colnames(obj),]$slingPseudotime_4
  obj@meta.data$Pseudotime_5 = mergeclusters@meta.data[colnames(obj),]$slingPseudotime_5
#  obj@meta.data[is.na(obj@meta.data)] <- 0
  seurat_obj_list[[var]] <- obj
}

FeaturePlot(seurat_obj_list$slide000_Control_A2, c("Pseudotime_1", "Pseudotime_2", "Pseudotime_3", "Pseudotime_4", "Pseudotime_5"), reduction = "coord" , cols = c("yellow", "red") , pt.size = 1.8)
