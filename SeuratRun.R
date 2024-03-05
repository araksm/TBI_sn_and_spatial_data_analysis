#
#  SeuratRun.R
#  
#
#  Created by Araks Martirosyan on 09/11/18.
#
#

SeuratV3Cluster <- function (FilteredDatabase, PCnumber = 10, res = 0.8, bias = c(), nFeature_RNAmin = 300)
{
    obj <- CreateSeuratObject(FilteredDatabase, min.cells = 5, min.features = 1)
    obj <- subset(obj, subset = (nFeature_RNA > nFeature_RNAmin & nFeature_RNA < 5000 )  ) #500 #10
    obj <- NormalizeData(obj, verbose = FALSE, normalization.method = "LogNormalize", scale.factor = 10000)
    obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 5000)
    VariableFeatures(obj) <- unique(VariableFeatures(obj), bias)
    print("Shh" %in% VariableFeatures(obj) )
    obj <- ScaleData(obj, features = unique(c(bias, rownames(obj))) )
    obj = RunPCA(obj,do.print=FALSE, pc.genes = unique(VariableFeatures(obj), bias ) )

#  obj=JackStraw(obj,num.replicate = 200)
    #  PCElbowPlot(obj)
    
    obj <- FindNeighbors(obj, dims = 1:PCnumber)
    
    obj <- FindClusters(obj, resolution = res, print.output = 0, save.SNN = T )
    obj <- RunUMAP(obj, dims = 1:PCnumber)
    # obj <- RunTSNE(obj, dims = 1:PCnumber)
    
    #DimPlot(obj)
    
    return(obj)
    
}



SeuratV3Merge <- function(data1, data2, S1label="S5", S2label="S6", number.cc = 20, dimensions.align = 20, res = 0.8, minCells = 0, minGenes = 0, bias = c(), npcs.n = 30) {
    
    # Set up the first object
    obj1 <- CreateSeuratObject(counts = data1, project = S1label, min.cells = 5)
    obj1$sample <- S1label
    obj1 <- subset(obj1, subset = (nFeature_RNA > 300 & nFeature_RNA < 5000 )  ) #500 #10
    obj1 <- NormalizeData(obj1, verbose = FALSE)
    obj1 <- FindVariableFeatures(obj1, selection.method = "vst", nfeatures = 5000)
    VariableFeatures(obj1) <- unique(c(bias,  VariableFeatures(obj1) ))
    
    # Set up the second object
    obj2 <- CreateSeuratObject(counts = data2, project = S2label, min.cells = 5)
    obj2$sample <- S2label
    obj2 <- subset(obj2, subset = (nFeature_RNA > 200 & nFeature_RNA < 4500) ) # 500 #10
    obj2 <- NormalizeData(obj2, verbose = FALSE)
    obj2 <- FindVariableFeatures(obj2, selection.method = "vst", nfeatures = 5000)
    VariableFeatures(obj2) <- unique(c(bias,  VariableFeatures(obj2) ))
    print(c(obj1, obj2))
    # Gene selection for input to CCA
    anchors <- FindIntegrationAnchors(object.list = list(obj1, obj2), dims = 1:number.cc)
    combined <- IntegrateData(anchorset = anchors, dims = 1:number.cc)
    
    DefaultAssay(combined) <- "integrated"
    VariableFeatures(combined) <- unique(c(bias,  VariableFeatures(combined) ))
    
    combined <- ScaleData(combined, verbose = FALSE, features = unique(rownames(combined), bias) )
    combined <- RunPCA(combined, npcs = npcs.n, verbose = FALSE)
    
    # t-SNE and Clustering
    combined <- RunUMAP(combined, reduction = "pca", dims = 1:dimensions.align)
    combined <- FindNeighbors(combined, reduction = "pca", dims = 1:dimensions.align)
    combined <- FindClusters(combined, resolution = res)

    return (combined)
    
}


read10xData <- function(path){
    
    library(Matrix)
    matrix_dir = path
    barcode.path <- paste0(matrix_dir, "barcodes.tsv.gz")
    features.path <- paste0(matrix_dir, "features.tsv.gz")
    matrix.path <- paste0(matrix_dir, "matrix.mtx.gz")
    mat <- readMM(file = matrix.path)
    feature.names = read.delim(features.path, header = FALSE, stringsAsFactors = FALSE)
    barcode.names = read.delim(barcode.path, header = FALSE, stringsAsFactors = FALSE)
    colnames(mat) = barcode.names$V1
    rownames(mat) = feature.names$V1
    
    return(mat)
}


runSlingshot <- function (SeuratObj, dimMethod, startClusterName, endClusterName) {

    return( slingshot(as.SingleCellExperiment(SeuratObj), clusterLabels = 'ident', reducedDim = dimMethod,  start.clus = startClusterName, end.clus = endClusterName) )

}



#get seurat colors

getColors <- function ( seuratplot ) {
    
    pbuild <- ggplot2::ggplot_build(seuratplot) # Use ggplot_build to deconstruct the ggplot object
    pdata <- pbuild$data[[1]] # Pull the data used for the plot
    pdata <-  pdata[order(pdata$group), ] # Order the plot data by group
    ucols <- unique(pdata$colour) # Get a vector of unique colors
    names(ucols) <- unique(pdata$group) # Add the groups to the vector of colors as names
    
    return (ucols)
}



### runslingshot

runSlingshotLineage <- function ( sce, lineageN ) {

    require(gam)
    LineageN = 5
    t <- sce$slingPseudotime_5
    Y <- assays(sce)$logcounts
    var1K <- unique( as.vector( unlist(sapply(as.integer(sce@int_metadata$slingshot@lineages[LineageN]$Lineage5), function(x) top10Astro[which(top10Astro$cluster == x),]$gene))))


    Y <- Y[var1K,]
    gam.pval <- apply(Y,1,function(z){
        d <- data.frame(z=z, t=t)
        tmp <- gam(z ~ lo(t), data=d)
        p <- summary(tmp)[4][[1]][1,5]
        p
    })

    require(clusterExperiment)
    topgenes <- var1K
    heatdata <- assays(sim)$logcounts[rownames(assays(sim)$logcounts) %in% topgenes, order(t, na.last = NA)]

    heatclus <- sce$ident[order(t, na.last = NA)]
    ce <- ClusterExperiment(as.matrix(heatdata), heatclus)
    ce@clusterLegend[[1]][,"color"] <- as.character(ucols[ce@clusterLegend[[1]][,"name"] ])
    
    plotHeatmap(ce, clusterFeatures = FALSE, clusterSamplesData = "orderSamplesValue", clusterFeaturesData = var1K )
    
    
    png(filename = paste0("~/Documents/Pascal/lineage",LineageN,".png"), width = 580, height = 780, units = "px")
    plotHeatmap(ce, clusterFeatures = FALSE, clusterSamplesData = "orderSamplesValue", clusterFeaturesData = var1K )
    dev.off()
    
    
    pdf( paste0("~/Documents/Pascal/lineage", LineageN, "-barplot.pdf"), width=2,height=6)
    datasubtype <- Fdata[which(Fdata$ident %in% as.integer(sce@int_metadata$slingshot@lineages[LineageN]$Lineage5) ),]
    barplot(table(datasubtype$ident), col = ucols, horiz = TRUE,  face = "bold", family = "Times", las = 1)
    dev.off()
    
    
    write(var1K, file = paste0("~/Documents/Pascal/genelistlineage", LineageN))
    
    top10TF <- AstroMarkersPosTFs %>% group_by(cluster) %>% top_n(15, avg_logFC)
     
    DoHeatmap(mergeCutAstros, top10TF$gene[1:15]) + theme(axis.text.y  = element_text(size=12, face = 4, family = "ArialMT"),axis.text.x  = element_text(size=12, face = "bold", family = "ArialMT"), axis.title.y = element_text(size=12, face = "bold", family = "ArialMT"), axis.title.x = element_text(size=12, face = "bold", family = "ArialMT"), plot.title = element_text(hjust = 0.5, size=12, face = "bold", family = "ArialMT"), legend.position = "right", panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.line = element_blank()
    )


}



vulcanoPlot <- function (res, str1, str2, th1, th2, genestr, xlim, ylim) {
    res[,str2] <- as.numeric(res[,str2])
    with(res, plot(res[,str2], -log10(res[, str1]), pch=20, main="Volcano plot", xlim = xlim, ylim = ylim, xlab = str2, ylab = str1,  col = res[, "color"]))
    # Add colored points: red if padj<th1, orange of log2FC>th2, green if both)
    #  with(subset(res, res[,str1] < th1 ), points(subset(res, res[,str1] < th1 )[,str2], -log10(subset(res, res[,str1] < th1 )[,str1]), pch=20, col="red3"))
    # with(subset(res, abs(res[,str2]) > th2), points(subset(res, abs(res[,str2]) > th2)[,str2], -log10(subset(res, abs(res[,str2]) > th2)[,str1]), pch=20, col="orange3"))
    #  with(subset(res, res[,str1] < th1 & abs(res[,str2])>th2), points(subset(res, res[,str1] < th1 & abs(res[,str2])>th2)[,str2], -log10(subset(res, res[,str1] < th1 & abs(res[,str2])>th2)[,str1]), pch=20, col=ucols[as.character(astroSubMarkers$cl)]))

    # Label points with the textxy function from the calibrate plot
    # library(calibrate)
      with(subset(res, res[,str1] < th1 & abs(res[,str2])>th2), textxy(subset(res, res[,str1] < th1 & abs(res[,str2])>th2)[,str2], -log10(subset(res, res[,str1] < th1 & abs(res[,str2])>th2)[,str1]), labs=subset(res, res[,str1] < th1 & abs(res[,str2])>th2)[, genestr], cex=.8, offset = sample(1:9, 3, replace = TRUE)/10, srt=30 ))
    # return (subset(res, res[,str1] < th1 & abs(res[,str2])>th2 & res[, "color"] != "gray"))
    
}




stripeChart <- function (res, str1, str2, th1, th2, genestr, xlim, ylim, pch = 20, cex = 1) {
    res[,str2] <- as.numeric(res[,str2])
    with(res, plot(res[,str2], res[, str1], pch=pch, cex=cex, face = "plain", family = "ArialMT", main="DE genes between TBI and Sham", xlim = xlim, ylim = ylim, xlab = "ln fold change", ylab = "",  col = res[, "color"],  horiz = TRUE, yaxt='n') ) # "↑TBI (green)    ↓TBI (orange),
    # with(res, axis(2, at=c(1:27), labels= FALSE ))
     with(res, axis(2, at=c(1:length(unique(res$Lineage))),labels= unique(res$lineage_to_plot), cex.lab=1, las = 2, family = "ArialMT" ))
    
    
    #dataframe <- data.frame (y = c(1:27), x = rep(0, 27) )
    # geom_text(aes ( c(1:27), rep(0, 27) ),  as.character(unique(obj@active.ident)))
    
    #  cap <- grid.cap()
    #  grid.newpage()
    #  grid.raster(cap, vp=viewport(angle=-90))
    
    # Add colored points: red if padj<th1, orange of log2FC>th2, green if both)
    # with(subset(res, res[,str1] > th1 ), points(subset(res, res[,str1] > th1 )[,str2], subset(res, res[,str1] > th1 )[,str1], pch=20, col="gray"))
    #   with(subset(res, abs(res[,str2]) < th2), points(subset(res, abs(res[,str2]) < th2)[,str2], subset(res, abs(res[,str2]) < th2)[,str1], pch=20, col="gray"))
     #  with(subset(res, res[,str1] > th1 & abs(res[,str2])>th2), points(subset(res, res[,str1] > th1 & abs(res[,str2])>th2)[,str2], subset(res, res[,str1] > th1 & abs(res[,str2])>th2)[,str1], pch=20, col="green4"))
    
    # Label points with the textxy function from the calibrate plot
    #   library(calibrate)
    #  with(subset(res, res[,str1] > th1 & abs(res[,str2])>th2 & res[, "color"] != "gray"), textxy(subset(res, res[,str1] > th1 & abs(res[,str2])>th2 & res[, "color"] != "gray")[,str2], subset(res, res[,str1] > th1 & abs(res[,str2])>th2 & res[, "color"] != "gray")[,str1], labs=subset(res, res[,str1] > th1 & abs(res[,str2])>th2 & res[, "color"] != "gray")[, genestr], cex=.5, offset = sample(1:9, 10, replace = TRUE)/10, srt=30, xpd = TRUE))
    
}

findTBIgenesByReclustering <- function  (obj, clusterN) {
    genelist <- diseasedgenesAstroNeuro[which(diseasedgenesAstroNeuro$clusterName == clusterN),]$gene

    
    cells <- WhichCells(obj, idents = as.character(clusterN)  )
    S5namesANtest <- intersect(cells, paste0(colnames(S5), "_1") )
    S6namesANtest <- intersect(cells, paste0(colnames(S6), "_2") )
    
    
    S5cutANtest <- S5 [ , matrix(unlist(strsplit(S5namesANtest, "_")), ncol = 2, byrow=TRUE)[, 1] ]
    S6cutANtest <- S6 [ , matrix(unlist(strsplit(S6namesANtest, "_")), ncol = 2, byrow=TRUE)[, 1] ]

    objreturn <- SeuratV3Merge(S5cutANtest, S6cutANtest, "Sham", "TBI", bias = genelist,  number.cc = 5, dimensions.align = 5,  npcs.n = 5)


#   markers <- FindMarkers(astroSub, ident.1 = "TBI", ident.2 = "Sham", test.use = "MAST", logfc.threshold = 0.00002)
#   genelist <- rownames(markers)
    
    return(objreturn)
    
}

MASTtest <- function  (obj, clusterN, threshold = 0.00002, gene = "Ppp1r14b") {
    
   # astroSub <- subset(object = obj, ident.use = clusterN)
    astroSub <- subset(AstroNeurontest, cells = WhichCells(AstroNeurontest, idents = clusterN))
    #  print(astroSub)
    astroSub <- SetIdent(astroSub, value = astroSub@meta.data$sample)
    #astroSub <- RenameIdents(astroSub,  "S5" = "Sham", "S6" = "TBI" )
    markers <- FindMarkers(astroSub, ident.1 = "TBI", ident.2 = "Sham", test.use = "MAST", logfc.threshold = threshold)
    markers$p.value = markers$p_val_adj
    return( markers[gene,]$p_val_adj)
    
}


findTBIgenes <- function  (obj, clusterN, threshold = 0.2) {
    genelist <- c()
    
   # astroSub <- subset(object = obj, ident.use = clusterN)
    astroSub <- subset(AstroNeurontest, cells = WhichCells(AstroNeurontest, idents = clusterN))
    #  print(astroSub)
    astroSub <- SetIdent(astroSub, value = astroSub@meta.data$sample)
    #astroSub <- RenameIdents(astroSub,  "S5" = "Sham", "S6" = "TBI" )
    markers <- FindMarkers(astroSub, ident.1 = "TBI", ident.2 = "Sham", test.use = "MAST", logfc.threshold = threshold)
    genelist <- rownames(markers)
    return(markers)
    
}



VlnPlotPerCluster <- function  (obj, clusterN, genelist) {

    astroSub <- SubsetData(object = Neurons, ident.use = clusterN)
    print(astroSub)
    astroSub <- SetIdent(astroSub, value = astroSub@meta.data$sample)
    astroSub <- RenameIdents(astroSub,  "S5" = "Sham", "S6" = "TBI" )
    
    genelist <- rownames(FindMarkers(astroSub, ident.1 = "TBI", ident.2 = "Sham", test.use = "MAST"))
    
    print(genelist)
    
    VlnPlot(astroSub, genelist , ncol = 4)+ theme(axis.text.y  = element_text(size=12, face = 4, family = "ArialMT"),axis.text.x  = element_text(size=12, face = "bold", family = "ArialMT"), axis.title.y = element_text(size=12, face = "bold", family = "ArialMT"), axis.title.x = element_text(size=12, face = "bold", family = "ArialMT"), plot.title = element_text(hjust = 0.5, size=24, face = "bold", family = "ArialMT"), legend.position = "right", panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank() )


    
}


PlotAverages <- function(obj, geneNames){

    cluster.averages <- AverageExpression(obj, return.seurat = TRUE)
    DoHeatmap(cluster.averages,  features = geneNames, draw.lines = FALSE, slot = "data") + scale_fill_gradientn(colors = c("black", "darkseagreen4", "gray" ,"orange3")) + theme(axis.text.y  = element_text(size=12, face = 4, family = "ArialMT"),axis.text.x  = element_text(size=12, face = "bold", family = "ArialMT"), axis.title.y = element_text(size=12, face = "bold", family = "ArialMT"), axis.title.x = element_text(size=12, face = "bold", family = "ArialMT"), plot.title = element_text(hjust = 0.5, size=12, face = "bold", family = "ArialMT"), legend.position = "right", panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.line = element_blank() )
}



RunGO <- function(genelist){

    funcmatrix <- matrix( unlist(strsplit ( unlist(lapply (1:length(genelist), function(x) paste0( genelist[x], "&" , paste0("NNN", as.character(BP$Term[grepl(toupper(genelist[x]), BP$Genes)]) ), "&", paste0("NNN", as.character(BP$P.value[grepl(toupper(genelist[x]), BP$Genes)]) ) )  ) ), "&")), ncol = 3, byrow = TRUE )


    funcmatrix <- as.data.frame(funcmatrix)
    funcmatrix <- funcmatrix[which(funcmatrix$V2 != "NNN"),]
    genelist <- unique(funcmatrix$V1)

    funcmatrix <- matrix( unlist(strsplit ( unlist(lapply (1:length(genelist), function(x) paste0( genelist[x], "&" ,  as.character(BP$Term[grepl(toupper(genelist[x]), BP$Genes)]) , "&",  as.character(BP$P.value[grepl(toupper(genelist[x]), BP$Genes)])  )  ) ), "&")), ncol = 3, byrow = TRUE )

    funcmatrix <- as.data.frame(funcmatrix)



# funcmatrix <- funcmatrix[which(as.double(as.character(funcmatrix$V3) ) < 0.05),]


    funcmatrix$cl <- unlist(lapply(1: dim(funcmatrix)[1] , function(x) clustergene[which(clustergene$gene== as.character(funcmatrix$V1[x] )),]$cluster [1]  ))





    astrofunctions <- unique(unlist(lapply ( c( "neuro", "blood", "vassel", "Ca2+", "calcium", "channel", "glu", "synap", "neurogenesis", "gli",   "axon", "gaba", "lactate", "tight junction", "pruning", "aqua"), function(x) funcmatrix$V2[grepl( x , funcmatrix$V2 )] )))
    astrofunctions <- astrofunctions[-1]
    
    funcmatrix <- funcmatrix[funcmatrix$V2 %in% astrofunctions,]


    destination <- funcmatrix$V2
    origin <- funcmatrix$cl


    data <- data.frame(origin, destination)

    # Transform input data in a adjacency matrix
    adjacencyData <- with(data, table(origin, destination))

    # Charge the circlize library
    library(circlize)

    # Make the circular plot
    chordDiagram(adjacencyData, order = unique(union( origin, destination)), annotationTrack = "grid", transparency = 0.5, preAllocateTracks=list(link.sort = TRUE,link.decreasing = TRUE), link.decreasing = TRUE,  )+circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
        xlim = get.cell.meta.data("xlim")
        ylim = get.cell.meta.data("ylim")
        sector.name = get.cell.meta.data("sector.index")
        circos.text( mean(xlim), ylim[1] + .1, sector.name, family="ArialMT", facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.1))
        circos.axis(h = "top", labels.cex = 0.001, major.tick.percentage = 0, sector.index = sector.name, track.index = 2, minor.ticks = 0)
    }, bg.border = NA)



}



firstup <- function(x) {
    substr(x, 1, 1) <- toupper(substr(x, 1, 1))
    x
}

firstdown <- function(x) {
    substr(x, 1, 1) <- tolower(substr(x, 1, 1))
    x
}


Binomiltest <- function (obj, comp = threshold) {
    metadata <- FetchData(obj, c("ident", "sample"))

    aggdata <- aggregate(metadata, by = list(metadata$sample, metadata$ident), function(x) sum(x != 0))
    aggdata$ident <- NULL
    colnames(aggdata) <- c("Sample", "Ident", "NumberOfCells")
    aggdataTBI = aggdata[which(aggdata$Sample=="TBI"),]
    aggdataSham = aggdata[which(aggdata$Sample=="Sham"),]
    colnames(aggdataTBI) <- c("Sample", "Ident", "NumberOfTBI")
    aggdataTBI$NumberOfCells = aggdataTBI$NumberOfTBI + aggdataSham$NumberOfCells
    
    listPvalue = c()
    for(i in c(1:dim(aggdataTBI)[1])){
        x= aggdataTBI[i,]
        listPvalue = c(listPvalue, binom.test(as.integer(x$NumberOfTBI),x$NumberOfCells, p=as.numeric(threshold)/100 )$p.value)
    }
    aggdataTBI$binomialtest = listPvalue
    
    return (aggdataTBI)
    
}


plotOriginSampleBarplots <- function(obj, group1 = "TBI", group2 = "Sham", fontsize=15){


    metadata_tmp <- FetchData(obj, c("ident", "sample"))

    aggdata <- aggregate(metadata_tmp, by = list(metadata_tmp$sample, metadata_tmp$ident), function(x) sum(x != 0))
    #aggdata[24,] <- aggdata[23,]
    #aggdata[23,] <- c("S5", 11, 0)

    aggdata$ident <- NULL
    colnames(aggdata) <- c("Sample", "Ident", "NumberOfCells")
    ggplot( aggdata, aes(y = NumberOfCells, x = Ident,  fill = Sample)) +  geom_bar(stat="identity", position=position_dodge())

    listL <- c()
    for(i in seq(1, length(aggdata$NumberOfCells), 2)) {
        
        listL <- c( listL, aggdata$NumberOfCells[i]/ (aggdata$NumberOfCells[i] + aggdata$NumberOfCells[i+1]), aggdata$NumberOfCells[i+1]/ (aggdata$NumberOfCells[i] + aggdata$NumberOfCells[i+1]))
        
    }
    aggdata$Per <- listL*100

    threshold = dim(metadata_tmp[which(metadata_tmp$sample == group1),])[1] / (dim(metadata_tmp[which(metadata_tmp$sample == group1),])[1] + dim(metadata_tmp[which(metadata_tmp$sample == group2 ),])[1]) * 100
    compP = round(threshold/100, 2)
    aggdata$pval = rep(Binomiltest(obj, comp = compP )$binomialtest, each = 2)
    aggdata$Ident = ifelse(aggdata$pval < 0.05, paste0(aggdata$Ident, "*"), paste0(as.character(aggdata$Ident), "  " ) )
    aggdata$Ident = as.character(aggdata$Ident)
    aggdata$Ident = factor(aggdata$Ident, levels = unique(as.character(aggdata$Ident)))
    

    ggplot( aggdata, aes(y = Per, x = Ident,  fill = Sample)) +  geom_bar(stat="identity")+  labs(title = "", x = "", y =  "Origin %") +  geom_hline (yintercept = threshold, linetype="dashed", size = 1)+scale_y_continuous(breaks = c(0, round(threshold, 2), 100) ) +
    theme(axis.text.y  = element_text(size=fontsize, face = "bold", family = "ArialMT"), axis.text.x  = element_text(size=fontsize, face = "bold", family = "ArialMT", angle = 90, hjust = 1), axis.title.y = element_text(size=fontsize, face = "bold", family = "ArialMT"), axis.title.x = element_text(size=fontsize, face = "bold", family = "ArialMT"),plot.title = element_text(hjust = 0.5, size=fontsize, face = "bold", family = "ArialMT"), legend.title = element_blank(), legend.text = element_text(size=fontsize, face = "bold", family = "ArialMT"),legend.position = "top", panel.background = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line = element_blank(),panel.border = element_rect(colour = "black", fill=NA, size=1, linetype="solid")
    )

    
 #   return(aggdata)


}


barplotOfPercentages <- function(metadata, percentageIn =  "previousID", ids = "seurat_clusters" ){
    metadata_tmp = as.data.frame.matrix ( table(metadata[, c( percentageIn, ids)]) )
    metadata_tmp = metadata_tmp/rowSums(metadata_tmp)*100
    metadata_tmp = t.data.frame(as.data.frame(metadata_tmp))
    metadata_tmp = as.data.frame(metadata_tmp)
    metadata_tmp$spatial <- NULL
    
    library(gplots)
    heatmap.2(t(as.matrix(metadata_tmp)), Rowv=FALSE, dendrogram='none',Colv=FALSE, density.info = 'none', key = TRUE, keysize = 1, col="viridis", trace='none', margins = c(11,13), key.xlab = "% of representation", key.title = "", family = "Times", cexRow = 1.05, cexCol = 1.05,  key.par = list(cex=1))
}


#### spatial data


SeuratV3MergeSpatial <- function(data1, data2, S1label, S2label, number.cc = 20, dimensions.align = 20, res = 0.8, minCells = 0, minGenes = 0, bias = c(), npcs.n = 30) {
  
  # Set up the first object
  obj1 <- CreateSeuratObject(counts = data1, project = S1label, min.cells = 1)
  obj1$sample <- S1label
  obj1 <- subset(obj1, subset = (nFeature_RNA > 1 & nFeature_RNA < 5000 )  ) #500 #10
  obj1 <- NormalizeData(obj1, verbose = FALSE)
  obj1 <- FindVariableFeatures(obj1, selection.method = "vst", nfeatures = 5000)
  VariableFeatures(obj1) <- unique(c(bias,  VariableFeatures(obj1) ))
  
  # Set up the second object
  obj2 <- CreateSeuratObject(counts = data2, project = S2label, min.cells = 1)
  obj2$sample <- S2label
  obj2 <- subset(obj2, subset = (nFeature_RNA > 1 & nFeature_RNA < 4500) ) # 500 #10
  obj2 <- NormalizeData(obj2, verbose = FALSE)
  obj2 <- FindVariableFeatures(obj2, selection.method = "vst", nfeatures = 5000)
  VariableFeatures(obj2) <- unique(c(bias,  VariableFeatures(obj2) ))
  print(c(obj1, obj2))
  # Gene selection for input to CCA
  anchors <- FindIntegrationAnchors(object.list = list(obj1, obj2), dims = 1:number.cc)
  combined <- IntegrateData(anchorset = anchors, dims = 1:number.cc)
  
  DefaultAssay(combined) <- "integrated"
  VariableFeatures(combined) <- unique(c(bias,  VariableFeatures(combined) ))
  
  combined <- ScaleData(combined, verbose = FALSE, features = unique(rownames(combined), bias) )
  combined <- RunPCA(combined, npcs = npcs.n, verbose = FALSE)
  
  # t-SNE and Clustering
  combined <- RunUMAP(combined, reduction = "pca", dims = 1:dimensions.align)
  combined <- FindNeighbors(combined, reduction = "pca", dims = 1:dimensions.align)
  combined <- FindClusters(combined, resolution = res)
  
  return (combined)
  
}

####

SeuratV3MergeSpatial_SCT <- function(data1, data2, S1label, S2label, number.cc = 20, dimensions.align = 20, res = 0.8, minCells = 0, minGenes = 0, bias = c(), npcs.n = 30) {
  
  # Set up the first object
  obj1 <- CreateSeuratObject(counts = data1, project = S1label, min.cells = 1)
  obj1$sample <- S1label
  obj1 <- subset(obj1, subset = (nFeature_RNA > 1 & nFeature_RNA < 5000 )  ) #500 #10
  #obj1 <- NormalizeData(obj1, verbose = FALSE)
  obj1 <- SCTransform(obj1, ncells = 3000, verbose = FALSE)
  obj1 <- FindVariableFeatures(obj1, selection.method = "vst", nfeatures = 5000)
  VariableFeatures(obj1) <- unique(c(bias,  VariableFeatures(obj1) ))
  
  # Set up the second object
  obj2 <- CreateSeuratObject(counts = data2, project = S2label, min.cells = 1)
  obj2$sample <- S2label
  obj2 <- subset(obj2, subset = (nFeature_RNA > 1 & nFeature_RNA < 4500) ) # 500 #10
  obj2 <- SCTransform(obj2, ncells = 3000, verbose = FALSE)

  #obj2 <- NormalizeData(obj2, verbose = FALSE)
  obj2 <- FindVariableFeatures(obj2, selection.method = "vst", nfeatures = 5000)
  VariableFeatures(obj2) <- unique(c(bias,  VariableFeatures(obj2) ))
  print(c(obj1, obj2))
  # Gene selection for input to CCA
  
  objlist =c(obj1, obj2)
  features <- SelectIntegrationFeatures(objlist)

  objlist <- PrepSCTIntegration(
    objlist,
    anchor.features = features
  )
  
  anchors <- FindIntegrationAnchors(object.list = objlist, dims = 1:number.cc, normalization.method = "SCT", anchor.features=features )
  combined <- IntegrateData(anchorset = anchors, dims = 1:number.cc)
  
  DefaultAssay(combined) <- "SCT"
  combined = FindVariableFeatures(combined)
  VariableFeatures(combined) <- unique(c(bias,  VariableFeatures(combined) ))
  
 # combined <- ScaleData(combined, verbose = FALSE, features = unique(rownames(combined), bias) )
  combined <- RunPCA(combined, npcs = npcs.n, verbose = FALSE)
  
  # t-SNE and Clustering
  combined <- RunUMAP(combined, reduction = "pca", dims = 1:dimensions.align)
  combined <- FindNeighbors(combined, reduction = "pca", dims = 1:dimensions.align)
  combined <- FindClusters(combined, resolution = res)
  
  return (combined)
  
}



WilcoxonTest  <- function( datatable ###columns correspond to variable
            ){
                require(reshape2)
                require(ggplot2)
                
                meltData = melt(datatable)
                
                SIGN = FALSE
                test <- wilcox.test(meltData$value ~ meltData$variable)
                SIGN = ifelse(test$p.value < 0.05, TRUE, FALSE)

                ggplot(meltData, aes(x = variable, y = value)) + geom_violin(scale="width",adjust = 1,width = 0.5,fill = c( "red4", "blue4"))+  geom_jitter(shape = 16, position = position_jitter(0.07))+ theme_classic() + geom_signif(comparisons = list(unique(meltData$variable)), map_signif_level=SIGN)
                              
            }


#### clusterProfiler

PathwayAnalysis_GO_KEGG <- function(genelist, correction = "BH"){
    
    library(clusterProfiler)
    library(org.Hs.eg.db)

    
    gene.df <- bitr(genelist, toType = "ENTREZID",
    fromType = c("SYMBOL"),
    OrgDb = org.Hs.eg.db)


    ego2 <- enrichGO(gene = gene.df$SYMBOL,
    OrgDb         = "org.Hs.eg.db",
    keyType       = 'SYMBOL',
    ont           = "BP",
    pAdjustMethod = correction,
    pvalueCutoff  = 0.1,
    qvalueCutoff  = 0.05)


    result_up = ego2@result

    result_up$ONT = "BP"
    result_up1 = result_up[which(result_up$p.adjust < 0.1),]



    ego2 <- enrichGO(gene         = gene.df$SYMBOL,
    OrgDb         = "org.Hs.eg.db",
    keyType       = 'SYMBOL',
    ont           = "CC",
    pAdjustMethod = correction,
    pvalueCutoff  = 0.1,
    qvalueCutoff  = 0.05)


    result_up = ego2@result

    result_up$ONT = "CC"
    result_up2 = result_up[which(result_up$p.adjust < 0.1),]



    ego2 <- enrichGO(gene         = gene.df$SYMBOL,
    OrgDb         = "org.Hs.eg.db",
    keyType       = 'SYMBOL',
    ont           = "MF",
    pAdjustMethod = correction,
    pvalueCutoff  = 0.1,
    qvalueCutoff  = 0.05)


    result_up = ego2@result

    result_up$ONT = "MF"
    result_up3 = result_up[which(result_up$p.adjust < 0.1),]

    
   # kk <- enrichKEGG(gene = gene.df$ENTREZID,
   # organism     = 'human',
   # pvalueCutoff = 0.1, pAdjustMethod = correction)
    
   # result_up4 = kk@result
    
   # entID = lapply( strsplit(result_up4$geneID, "/"), FUN = function(x) as.numeric(x))
   # l = lapply(entID, function(x) gene.df[which(gene.df$ENTREZID %in% x),]$SYMBOL  )
   # result_up4$geneID = lapply(l, function(x) paste0(x, collapse = "/") )
    
   # result_up4$ONT = "KEGG"
   # result = rbind(result_up1,result_up2,result_up3, result_up4)
    result = rbind(result_up1,result_up2,result_up3)

    rm(result_up1,result_up2,result_up3, result_up4)

    return(result)
    
    
}
