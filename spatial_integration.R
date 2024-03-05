
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
#filenames = filenames[grepl("Original", filenames)]
#filenames = filenames[grepl("/Preview", filenames)]

filenames = setdiff( filenames[grepl("[O|o]riginal", filenames)] ,filenames[grepl("[O|o]ld", filenames)])

datalist = lapply(filenames,
                    function(x){read.csv(file = x,
                                         row.names = 1,
                                         stringsAsFactors = FALSE)})

names(datalist) <- c("slide000_Control_A1", "slide000_Control_A2", "slide000_TBI_A1", "slide000_TBI_A2",  "slide011_Control_D2", "slide011_TBI_D1", "slide016_TBI_C1", "slide010_Control_D2", "slide010_TBI_D1")

#names(datalist) <- c("slide000_Control_A1", "slide000_Control_A2", "slide000_TBI_A1", "slide000_TBI_A2")

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






#data_read <- Read10X(data.dir = "ResolveDataToAnalyze/sc_S5/"  )
#colnames(data_read) = paste0("sc_Control_S5zz", colnames(data_read))

#datalist[["sc_Control_S5"]] = data_read

#data_read <- Read10X(data.dir = "ResolveDataToAnalyze/sc_S6/"  )
#colnames(data_read) = paste0("sc_Control_S6zz", colnames(data_read))

#datalist[["sc_Control_S6"]] = data_read


commongenes =  Reduce(intersect, lapply(datalist, rownames))
datalistcommon = datalist
for(var in names(datalistcommon) ){
    datalistcommon[[var]] = datalistcommon[[var]][commongenes,]
}



#spatialdata = Reduce(function(x,y) {merge(x, y, all = TRUE)}, datalist)

#folders = paste0("ResolveDataToAnalyze/", sampleslide, "/*measurements.csv")

###### integration

data_read = datalistcommon [[names(datalistcommon)[1]]]
for(var in names(datalistcommon)[2:length(names(datalistcommon))] ){
    data_read = cbind( data_read, datalistcommon[[var]])
}

###selection
#selectedSections = names(datalistcommon)[c(1, 6, 9)]
#selectedSections = names(datalistcommon)[c(1:4)]

#selectedSections = names(datalistcommon)[c(1:5,9)]

#data_read = datalistcommon [[selectedSections[1]]]
#for(var in selectedSections[2:length(selectedSections)] ){
#    data_read = cbind( data_read, datalistcommon[[var]])
#}
### end of the selection

#l=length(colnames(data_read))


#data_sep <- CreateSeuratObject(counts = data_read)
#data_sep@meta.data$orig.ident = unlist(strsplit(colnames(data_read), "zz"))[seq(1,l*2,2)]
#l=length(colnames(data_sep))

#data_sep@meta.data$slide = unlist(strsplit(colnames(data_read), "_"))[seq(1,l*3,3)]
#data_sep@meta.data$condition = unlist(strsplit(colnames(data_read), "_"))[seq(2,l*3,3)]
#data_sep@meta.data$section = unlist(strsplit(colnames(data_read), "_"))[seq(3,l*3,3)]

### old way works better
## data_read = data_read[, setdiff(colnames(data_read), "slide000_Control_A2zzPathCellObject.1604")]
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

source("~/Dropbox/araks/codeperproject/SeuratRun.R")

mergeclusters <-SeuratV3MergeSpatial( obj10x@assays$RNA@counts, objSpatial@assays$RNA@counts,  "10x", "Genexyz", number.cc = 10, dimensions.align = 10)
DimPlot(mergeclusters, group.by = "sample")
DimPlot(mergeclusters, label= TRUE)

FeaturePlot(mergeclusters, c("Ascl1", "Aldoc", "Neurod1", "Frzb"))
mergeclusters@meta.data$previousID = "spatial"
mergeclusters@meta.data[intersect(rownames(mergeclusters@meta.data), rownames(AstroNeurontest@meta.data)),]$previousID = AstroNeurontest@meta.data[intersect(rownames(mergeclusters@meta.data), rownames(AstroNeurontest@meta.data)),]$previousIdent
DimPlot(mergeclusters, group.by = "previousID", label = TRUE)
DimPlot(mergeclusters, group.by = "previousID", split.by = "sample", label=TRUE)

###subsetdata
#objSpatial1 = SubsetData(objSpatial, cells = intersect( colnames(objSpatial), WhichCells(mergeclusters, idents = setdiff(c(1:12), c(0,1,2,8) ) )))
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
#Astro <- SubsetData(mergeclusters, ident.use =  c(3,15,18,7,6,1,20,21,9,12))
#Neuro<- SubsetData(mergeclusters, ident.use =  c( 3,15,19,14,13,10,8,4,5,11))

#Astro <- SubsetData(mergeclusters, ident.use =  c(10,7,6,0,9,8))
#Neuro<- SubsetData(mergeclusters, ident.use =  c( 10,7,2,3,5,4,11))

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
#  obj@meta.data[is.na(obj@meta.data)] <- 0
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

pdf("hist_nFeature_RNA.pdf")
hist(data_sep@meta.data$nFeature_RNA)
dev.off()

pdf("hist_percent.mt.pdf")
hist(data_sep@meta.data$percent.mt)
dev.off()


## Set lower threshold quite low because of low seq depth:
data_sep <- subset(data_sep, subset = nFeature_RNA > 0 & nFeature_RNA < 6000)
dim(data_sep)
## [1] 33537 110314


## Split object by sample:
data_list = SplitObject(data_sep, split.by = "orig.ident")

## Run sctransform for each sample:
for(i in names(data_list))
{
  cat("Sample ", i, "....\n")
  data_list[[i]] <- SCTransform(data_list[[i]], 
                                verbose = FALSE)
}

## Select features for downstream integration, and run PrepSCTIntegration, which ensures that all necessary Pearson residuals have been calculated:
data_features <- SelectIntegrationFeatures(object.list = data_list, 
                                           nfeatures = 500)
data_features = setdiff(data_features, c("Pcdha11", "Pcdhgb2", "Pcdhga4", "Pcdha2", "Serpina3n") )


data_list <- PrepSCTIntegration(object.list = data_list, 
                                anchor.features = data_features)
## Run PCA so I can run the FindIntegrationAchors with reduction = "rpca":
data_list <- lapply(X = data_list, FUN = RunPCA, verbose = FALSE, features = data_features)
## Identify anchors and integrate the datasets:
## Find integration vectors:
data_anchors <- FindIntegrationAnchors(object.list = data_list, 
                                       normalization.method = "SCT", 
                                       reduction = "cca", ## normally this would be "cca", but it's recommended to switch to reciprocal pca for very big datasets. But to run "rpca" I need to run PCA beforehand (see above).
                                       anchor.features = data_features,
                                       reference = 1) ## chose which dataset should be reference for integration (that speeds up the integration, also not doing it gives me an error (see below)
## If I don't set a reference dataset, I get this error:
#Error in validityMethod(as(object, superClass)) :
#  long vectors not supported yet: ../../src/include/Rinlinedfuns.h:522 "")

## On default, all non-anchor genes are being thrown away during data integration. 
## Create list of common genes to keep:
genes_to_integrate = Reduce(intersect, lapply(data_list, function(x) rownames(x@assays$SCT@scale.data)))
length(genes_to_integrate)

#genes_to_integrate  = c("Ascl1","Aurkb", "Ube2c","Neurog2", "Frzb", "Lpar1", "Eomes", "Fabp7", "Aldoc", "Fgfr3", "Ogt", "Fam107a", "Hapln1", "Neat1", "Sparc", "Sned1", "Sox4", "Neurod1", "Calb2", "Snap25", "Unc5d", "Vim", "Mpped1", "Nell2", "Kitl", "Pcp4","Ogt", "Syt1")


## 6000
## save workspace in case integratedata doesn't work:
save.image(file = "seurat_preprocess_image.RData")

## This step takes a lot of memory if done on all genes:
## For some reason with Seurat v.3.1.1 it doesn't work anymore on genes_to_integrate. Figure this out. For now I remove it.
data <- IntegrateData(anchorset = data_anchors, 
                      normalization.method = "SCT",
                      features.to.integrate = genes_to_integrate)
## Warning: Adding a command log without an assay associated with it

## Proceed with downstream analysis:
data <- RunPCA(data)
## Run UMAP on first couple PCs:

data <- JackStraw(data, num.replicate = 100)
data <- ScoreJackStraw(data, dims = 1:40)
ElbowPlot(data)

data <- RunUMAP(data, dims = 1:15) ## for full study also try dims = 1:30 (that's being used in the more recent vignettes)
# data = RunTSNE(data)


## Cluster the cells based on PCs:
data = FindNeighbors(data, dims = 1:15)
data = FindClusters(data, resolution = 0.6)

DimPlot(data, pt.size = 1.0, label.size = 5, label = TRUE)
DefaultAssay(data) = "SCT"
FeaturePlot(data, c("Aldoc", "Ascl1", "Neurod1", "Snap25"))

metadata =  data@meta.data
tablecounts = table(metadata$slide, metadata$seurat_clusters)

barplot(tablecounts,
        main = "representation in each category",
        xlab = "Class",
        col = c("red2","green2", "blue2", "orange", "green4")
)
legend("topright",
       rownames(tablecounts),
       fill = c("red2","green2", "blue2", "orange", "green4")
)
save(data, file = "seurat_data.RData")

pdf("umap.pdf",
    width = 7.5, 
    height = 6.5)
DimPlot(object = data, 
               reduction = "umap",
               label = FALSE)
dev.off()

pdf("umap_cluster.pdf",
    width = 10, 
    height = 7)
DimPlot(object = data, 
        reduction = "umap",
        label = TRUE)
dev.off()


pdf("umap_sample.pdf",
    width = 7.5, 
    height = 6.5)
DimPlot(object = data, 
        group.by = "orig.ident",
        reduction = "umap",
        label = FALSE)

dev.off()

pdf("umap_nFeature_RNA.pdf",
    width = 7.5, 
    height = 6.5)
FeaturePlot(object = data, 
        reduction = "umap",
        features = "nFeature_RNA")
dev.off()

pdf("umap_nCount_RNA.pdf",
    width = 7.5, 
    height = 6.5)
FeaturePlot(object = data, 
            reduction = "umap",
            features = "nCount_RNA")
dev.off()


data@assays$RNA@counts["Prox1", ] = 0

for(var in names(data_list) ) {
    data@assays$RNA@counts["Prox1",intersect(colnames(data@assays$RNA@counts), colnames(datalist[[var]] ))] = as.numeric(datalist[[var]]["Prox1", intersect(colnames(data@assays$RNA@counts), colnames(datalist[[var]]))] )

}

#data@assays$RNA@counts["Prox1",intersect(colnames(data@assays$RNA@counts), colnames(datalist$sc_Control_S6))] = as.numeric(datalist$sc_Control_S6["Prox1", intersect(colnames(data@assays$RNA@counts), colnames(datalist$sc_Control_S6))] )

#data@meta.data[intersect(rownames(data@meta.data), colnames(datalist$sc_Control_S5)),]$Prox1   = as.numeric(datalist$sc_Control_S5["Prox1", intersect(rownames(data@meta.data), colnames(datalist$sc_Control_S5))] )


prox1 = data@assays$RNA@counts["Prox1",] >5
cellsProx1 = names(prox1[which(prox1)])
DimPlot(data, cells = setdiff( colnames(data), cellsProx1), label = TRUE)

ll = length(unique(data@active.ident)) - 1
nspc = subset( data, cells = WhichCells(data, idents = setdiff(c(0:ll), c(0) )))
DefaultAssay(nspc) = "SCT"
VariableFeatures(nspc) =  genes_to_integrate
nspc <- RunPCA(nspc)
nspc <- JackStraw(nspc, num.replicate = 100)
nspc <- ScoreJackStraw(nspc, dims = 1:20)
ElbowPlot(nspc)
nspc <- RunUMAP(nspc, dims = 1:20)
nspc <- FindNeighbors(nspc, dims = 1:20)
nspc <- FindClusters(nspc, resolution = 0.6)
DimPlot(nspc, label = TRUE, label.size = 5, pt.size = 0.8)
DimPlot(nspc,  label.size = 5, group.by = "slide")

FeaturePlot(nspc, c("Aldoc", "Ascl1", "Neurod1", "Snap25"))

metadata =  data@meta.data
table(metadata$slide, metadata$seurat_clusters)



nspc2 = subset( nspc, cells = WhichCells(nspc, idents = setdiff(c(1:29), c(2,14, 29, 27, 25, 10, 4, 24) )))
nspc2 = subset( nspc2, features = rownames(nspc)[rownames(nspc) != "Prox1" ]  )

nspc2 <- RunPCA(nspc2)
nspc2 <- RunUMAP(nspc2, dims = 1:10)
nspc2 = FindNeighbors(nspc2, dims = 1:10)
nspc2 = FindClusters(nspc2, resolution = 0.6)
DimPlot(nspc2,  label = TRUE, label.size = 5)
DimPlot(nspc2,  label.size = 5, group.by = "slide")
FeaturePlot(nspc2, c("Aldoc", "Ascl1", "Neurod1", "Snap25"))


nspc3 = subset( nspc2, cells = WhichCells(nspc, idents = setdiff(c(1:16), c(2,14, 29, 27, 25, 10, 4, 24) )))
nspc3 = subset( nspc, features = rownames(nspc)[rownames(nspc) != "Prox1" ]  )

nspc2 <- RunPCA(nspc2)
nspc2 <- RunUMAP(nspc2, dims = 1:10)
nspc2 = FindNeighbors(nspc2, dims = 1:10)
nspc2 = FindClusters(nspc2, resolution = 0.6)
DimPlot(nspc2,  label = TRUE, label.size = 5)
FeaturePlot(nspc2, c("Aldoc", "Ascl1", "Neurod1"))




### Average comparisons
AllSpatialAverage = Reduce(rbind, lapply(datalist, rowMeans))
rownames(AllSpatialAverage) =names(datalist)
AllSpatialAverage = as.data.frame(AllSpatialAverage)
AllSpatialAverage = t.data.frame(AllSpatialAverage)
AllSpatialAverage = as.data.frame(AllSpatialAverage)

## c(1,2,5,7) c(3,4,6,8,9))
par(mfrow=c(4, 4))
for(var in names(datalist) [c(3,4,6,8,9)]){
    for(var2 in names(datalist) [c(3,4,6,8,9)] ){
            # Creating the plot
            y = AllSpatialAverage[[var]]
            x = AllSpatialAverage[[var2]]
            plot(x, y, pch = 19, col = "lightblue3", xlab=var, ylab=var2) + abline(a = 0, b = 1, col = "blue", lwd = 1, xlab="", ylab="")+abline(lm(y ~ x), col = "red", lwd = 1) + title(paste("Correlation:", round(cor(x, y), 2)))
    }
}



obj@meta.data$slide = unlist(strsplit(colnames(AllSpatialAverage), "_"))[seq(1,27,3)]
obj@meta.data$condition = unlist(strsplit(colnames(AllSpatialAverage), "_"))[seq(2,27,3)]
obj@meta.data$section = unlist(strsplit(colnames(AllSpatialAverage), "_"))[seq(3,27,3)]
