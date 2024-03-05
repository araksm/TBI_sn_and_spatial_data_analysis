load("~/Documents/Pascal/test.RData")
source("~//SeuratRun.R")

Combinedtest <- SeuratV3Merge(S5[rownames(S5)[rownames(S5) != "Malat1"],], S6[rownames(S6)[rownames(S6) != "Malat1"],], "S5", "S6",  number.cc = 20, dimensions.align = 20, npcs.n = 30, res = 2)
ucols <- c("coral", "goldenrod1", "yellowgreen", "seagreen3", "darkolivegreen3", "darkolivegreen4", "mediumseagreen", "lightskyblue2", "lightskyblue3", "skyblue4", "lightslateblue", "turquoise4", "khaki3", "palevioletred", "tan", "tan3", "orange2")


Combinedtest <- CombinedtestV4
Combinedtest <- SetIdent( Combinedtest, WhichCells(object = Combinedtest, idents = c(10,9,8,4)), value = "Oligo" )
Combinedtest <- SetIdent( Combinedtest, WhichCells(object = Combinedtest, idents = c(0,16)), value = "Endo" )
Combinedtest <- SetIdent( Combinedtest, WhichCells(object = Combinedtest, idents = c(20)), value = "Mgl1" )
Combinedtest <- SetIdent( Combinedtest, WhichCells(object = Combinedtest, idents = c(23)), value = "Mgl2" )
Combinedtest <- SetIdent( Combinedtest, WhichCells(object = Combinedtest, idents = c(17,22,12,7)), value = "Mural" )
Combinedtest <- SetIdent( Combinedtest, WhichCells(object = Combinedtest, idents = c(13)), value = "NSC" )
Combinedtest <- SetIdent( Combinedtest, WhichCells(object = Combinedtest, idents = c(5)), value = "RG-like" )
Combinedtest <- SetIdent( Combinedtest, WhichCells(object = Combinedtest, idents = c(19)), value = "A1" )
Combinedtest <- SetIdent( Combinedtest, WhichCells(object = Combinedtest, idents = c(11)), value = "A2" )
Combinedtest <- SetIdent( Combinedtest, WhichCells(object = Combinedtest, idents = c(6)), value = "A3" )
Combinedtest <- SetIdent( Combinedtest, WhichCells(object = Combinedtest, idents = c(3)), value = "A4" )
Combinedtest <- SetIdent( Combinedtest, WhichCells(object = Combinedtest, idents = c(21)), value = "OPC" )
Combinedtest <- SetIdent( Combinedtest, WhichCells(object = Combinedtest, idents = c(14)), value = "U" )
Combinedtest <- SetIdent( Combinedtest, WhichCells(object = Combinedtest, idents = c(2)), value = "N1" )
Combinedtest <- SetIdent( Combinedtest, WhichCells(object = Combinedtest, idents = c(1)), value = "N2" )
Combinedtest <- SetIdent( Combinedtest, WhichCells(object = Combinedtest, idents = c(15)), value = "N3" )
Combinedtest <- SetIdent( Combinedtest, WhichCells(object = Combinedtest, idents = c(18)), value = "N4" )
Combinedtest@active.ident <- factor(Combinedtest@active.ident, levels = c(    "NSC"  ,   "RG-like", "A1" ,     "A2"    ,  "A3" ,     "A4" ,  "OPC"  , "U",  "N1"  ,    "N2" ,     "N3"  ,    "N4"  , "Oligo" ,  "Endo"   ,
"Mgl1"  , "Mgl2"   ,  "Mural") , ordered = T)
DimPlot(Combinedtest, label = TRUE, cols = ucols) + theme(axis.text.y  = element_text(size=15, face = "bold", family = "Times New Roman"),
axis.text.x  = element_text(size=15, face = "bold", family = "Times New Roman"),
axis.title.y = element_text(size=15, face = "bold", family = "Times New Roman"),
axis.title.x = element_text(size=15, face = "bold", family = "Times New Roman"),
plot.title = element_text(hjust = 0.5, size=15, face = "bold", family = "Times New Roman"),
legend.position = "right",
panel.background = element_blank(),
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
axis.line = element_blank(),
panel.border = element_rect(colour = "black", fill=NA, size=1, linetype="solid")
)

Combinedtest@meta.data [which(Combinedtest@meta.data$sample == "S6"),]$sample <- "TBI"
Combinedtest@meta.data [which(Combinedtest@meta.data$sample == "S5"),]$sample <- "Sham"
Combinedtest@meta.data <- na.omit(Combinedtest@meta.data)
source("~/Documents/AstrocyteDevelopment/SeuratRun.R")
plotOriginSampleBarplots(Combinedtest)

testAllMarkers <- FindAllMarkers(Combinedtest, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)


load("~/Documents/Pascal/TBIwork/nscv3/CombinedtestV4.RData")


highlevelgenestest <- as.character(read.table("~/Documents/Pascal/HighLevelGenestest copy.txt")$V1)
DoHeatmap(Combinedtest, intersect( testAllMarkers$gene, highlevelgenestest), lines.width = 10) + theme(axis.text.y  = element_text(size=12, face = 3, family = "Times New Roman"),axis.text.x  = element_text(size=12, face = "plain", family = "Times New Roman"), axis.title.y = element_text(size=12, face = "plain", family = "Times New Roman"), axis.title.x = element_text(size=12, face = "plain", family = "Times New Roman"), plot.title = element_text(hjust = 0.5, size=12, face = "plain", family = "Times New Roman"), legend.position = "right", panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.line = element_blank() )
PlotAverages(Combinedtest, intersect( testAllMarkers$gene, highlevelgenestest))

save(Combinedtest, ucols , testAllMarkers, file = "~/Documents/Pascal/TBIwork/nscv3/CombinedtestV4.RData", version = 2 )

##### Astrocytes and Neurons


##### Astrocytes and Neurons
load("~/Documents/Pascal/test.RData")
load("~/Documents/Pascal/TBIwork/nscv3/CombinedtestV4.RData")

AstroNeurontestV3 <- AstroNeurontest
cells <- WhichCells(Combinedtest, idents = levels(Combinedtest)[ ! (levels(Combinedtest)  %in% c("OPC" , "U", "Oligo" ,  "Endo"   ,
                                                                                                     "Mgl1"  , "Mgl2"   ,  "Mural"
))
]  )

S5namesAtest <- intersect(cells, paste0(colnames(S5), "_1") )
S6namesAtest <- intersect(cells, paste0(colnames(S6), "_2") )

source("~/Documents/AstrocyteDevelopment/SeuratRun.R")
set.seed(700)
S5cutAtest <- S5 [ rownames(S5)[rownames(S5) != "Malat1" ], matrix(unlist(strsplit(S5namesAtest, "_")), ncol = 2, byrow=TRUE)[, 1] ]
S6cutAtest <- S6 [ rownames(S5)[rownames(S5) != "Malat1" ], matrix(unlist(strsplit(S6namesAtest, "_")), ncol = 2, byrow=TRUE)[, 1] ]
AstroNeurontest <- SeuratV3Merge(S5cutAtest, S6cutAtest, "S5Cut", "S6Cut", res = 1.5, number.cc = 30, dimensions.align = 16, npcs.n = 30)




DimPlot(AstroNeurontest, label = TRUE)
DimPlot(AstroNeurontest, label = TRUE, split.by = "sample")
save(AstroNeurontest, file = "~/Documents/Pascal/TBIwork/nscv3/AstroNeurontestV5.RData", version = 2)


#AstroNeurontest <- JackStraw(AstroNeurontest, num.replicate = 100)
#AstroNeurontest <- ScoreJackStraw(AstroNeurontest, dims = 1:20)
#ElbowPlot(AstroNeurontest)
###16PC is good for astros
### 
mergeCutAstros@meta.data$active.ident <- mergeCutAstros@active.ident
AstroNeurontest@meta.data$previosIdent <- "100"
#AstroNeurontest@meta.data[colnames(mergeCutAstros),]$previosIdent <- as.character(mergeCutAstros@meta.data[intersect(colnames(mergeCutAstros), matrix(unlist(strsplit(colnames(AstroNeurontest), "_")), ncol = 2, byrow=TRUE)[, 1]),]$active.ident)
S5namesAst <- intersect(matrix(unlist(strsplit(colnames(AstroNeurontest)[which(AstroNeurontest@meta.data$sample == "S5")], "_")), ncol = 2, byrow=TRUE)[, 1], colnames(mergeCutAstros)[which(mergeCutAstros@meta.data$sample == "S5")] )
S6namesAst <- intersect(matrix(unlist(strsplit(colnames(AstroNeurontest)[which(AstroNeurontest@meta.data$sample == "S6")], "_")), ncol = 2, byrow=TRUE)[, 1], colnames(mergeCutAstros)[which(mergeCutAstros@meta.data$sample == "S6")] )
AstroNeurontest@meta.data[paste0(S5namesAst,  "_1"),]$previosIdent <- as.character(mergeCutAstros@meta.data[S5namesAst,]$active.ident)
AstroNeurontest@meta.data[paste0(S6namesAst,  "_2"),]$previosIdent <- as.character(mergeCutAstros@meta.data[S6namesAst,]$active.ident)
#AstroNeurontest@meta.data <- na.omit(AstroNeurontest@meta.data)
DimPlot(AstroNeurontest, group.by = "previosIdent",label = TRUE)

AstroNeurontest@meta.data$previousIdent <- as.character(Combinedtest[,colnames(AstroNeurontest)]@active.ident)
DimPlot(AstroNeurontest, group.by = "previousIdent", label = TRUE)


DimPlot(AstroNeurontest, label = TRUE)
FeaturePlot(AstroNeurontest, c("Sparc","Sned1","Cadps", "Neat1", "Hapln1", "Cdh2", "Gfap"))
FeaturePlot(AstroNeurontest, c("Aurkb", "Pclaf", "Add2", "Pou2f2", "Sstr2", "Mycn", "Neurog2") )
RidgePlot(AstroNeurontest, c(  "Neurod1" , "Sparc","Sned1","Cadps", "Neat1", "Hapln1"), pt.size = -1, group.by = "previosIdent"  )
testAstroMarkers <- FindAllMarkers(AstroNeurontest, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)

library(slingshot)
sce <- runSlingshot(AstroNeurontest, "UMAP", 8, c( 6 ) )
plot(reducedDims(sce)$UMAP, col = "gray", pch=16, asp = 1)
lines(SlingshotDataSet(sce), lwd=2, col='black', type = 'lineages')

AstroNeurontest <- AstroNeurontestV5
AstroNeurontest <- SetIdent( AstroNeurontest, WhichCells(object = AstroNeurontest, idents = c(18)), value = "NSC_1" )
AstroNeurontest <- SetIdent( AstroNeurontest, WhichCells(object = AstroNeurontest, idents = c(10)), value = "NSC_2" )
AstroNeurontest <- SetIdent( AstroNeurontest, WhichCells(object = AstroNeurontest, idents = c(15,8)), value = "RG-like_1" )
AstroNeurontest <- SetIdent( AstroNeurontest, WhichCells(object = AstroNeurontest, idents = c(11)), value = "RG-like_2" )
AstroNeurontest <- SetIdent( AstroNeurontest, WhichCells(object = AstroNeurontest, idents = c(17)), value = "A1" )
AstroNeurontest <- SetIdent( AstroNeurontest, WhichCells(object = AstroNeurontest, idents = c(7)), value = "A2" )
AstroNeurontest <- SetIdent( AstroNeurontest, WhichCells(object = AstroNeurontest, idents = c(2)), value = "A3_1" )
AstroNeurontest <- SetIdent( AstroNeurontest, WhichCells(object = AstroNeurontest, idents = c(13)), value = "A3_2" )
AstroNeurontest <- SetIdent( AstroNeurontest, WhichCells(object = AstroNeurontest, idents = c(5)), value = "A4_1" )
AstroNeurontest <- SetIdent( AstroNeurontest, WhichCells(object = AstroNeurontest, idents = c(6)), value = "A4_2" )

AstroNeurontest <- SetIdent( AstroNeurontest, WhichCells(object = AstroNeurontest, idents = c(4, 20)), value = "N1_1" )
AstroNeurontest <- SetIdent( AstroNeurontest, WhichCells(object = AstroNeurontest, idents = c(3)), value = "N1_2" )
AstroNeurontest <- SetIdent( AstroNeurontest, WhichCells(object = AstroNeurontest, idents = c(0)), value = "N2_1" )
AstroNeurontest <- SetIdent( AstroNeurontest, WhichCells(object = AstroNeurontest, idents = c(1)), value = "N2_2" )
AstroNeurontest <- SetIdent( AstroNeurontest, WhichCells(object = AstroNeurontest, idents = c(16)), value = "N3_1" )
AstroNeurontest <- SetIdent( AstroNeurontest, WhichCells(object = AstroNeurontest, idents = c(9)), value = "N3_2" )
AstroNeurontest <- SetIdent( AstroNeurontest, WhichCells(object = AstroNeurontest, idents = c(14)), value = "N4" )
AstroNeurontest <- SetIdent( AstroNeurontest, WhichCells(object = AstroNeurontest, idents = c(19)), value = "U1" )
AstroNeurontest <- SetIdent( AstroNeurontest, WhichCells(object = AstroNeurontest, idents = c(12)), value = "U2" )

ucolsAstroNeuro <- getColors(DimPlot(AstroNeurontest))

AstroNeurontest@active.ident <- factor(AstroNeurontest@active.ident, levels = rev( unique(levels(AstroNeurontest@active.ident))), ordered = T)
c("gray", "gray", "gray" ,"red", rep("gray", 14))
DimPlot(AstroNeurontest, label = TRUE, cols = as.character(ucolsAstroNeuro) ) + theme(axis.text.y  = element_text(size=15, face = "bold", family = "Times New Roman"),
                                                                                     axis.text.x  = element_text(size=15, face = "bold", family = "Times New Roman"),
                                                                                     axis.title.y = element_text(size=15, face = "bold", family = "Times New Roman"),
                                                                                     axis.title.x = element_text(size=15, face = "bold", family = "Times New Roman"),
                                                                                     plot.title = element_text(hjust = 0.5, size=15, face = "bold", family = "Times New Roman"),
                                                                                     legend.position = "right",
                                                                                     panel.background = element_blank(),
                                                                                     panel.grid.major = element_blank(),
                                                                                     panel.grid.minor = element_blank(),
                                                                                     axis.line = element_blank(),
                                                                                     panel.border = element_rect(colour = "black", fill=NA, size=1, linetype="solid")
)

AstroNeurontest@meta.data [which(AstroNeurontest@meta.data$sample == "S6"),]$sample <- "TBI"
AstroNeurontest@meta.data [which(AstroNeurontest@meta.data$sample == "S5"),]$sample <- "Sham"
AstroNeurontest@meta.data <- na.omit(AstroNeurontest@meta.data)

testAstroneuroMarkers <- FindAllMarkers(AstroNeurontest, only.pos = FALSE, min.pct = 0.25, thresh.use = 0.25)
PlotAverages(obj = AstroNeurontest, geneNames = testAstroneuroMarkers$gene)
top10 <- testAstroneuroMarkers %>% group_by(cluster) %>% top_n(50, avg_logFC)
DoHeatmap(AstroNeurontest, top10$gene)
PlotAverages(obj = AstroNeurontest, geneNames = top10$gene)+ theme(axis.text.y  = element_text(size=12, face = 4, family = "Times New Roman"),axis.text.x  = element_text(size=12, face = "bold", family = "Times New Roman"), axis.title.y = element_text(size=12, face = "bold", family = "Times New Roman"), axis.title.x = element_text(size=12, face = "bold", family = "Times New Roman"), plot.title = element_text(hjust = 0.5, size=12, face = "bold", family = "Times New Roman"), legend.position = "right", panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.line = element_blank()
)

FeaturePlot(AstroNeurontest, c("Sparc","Sned1","Cadps", "Neat1", "Hapln1", "Cdh2", "Gfap"))

AstroNeurontest@active.ident <- ordered(AstroNeurontest@active.ident, levels = c(   "NSC_1",  "NSC_2",  "RG-like_1", "RG-like_2",  "A1" , "U1", "A2", "U2", "A3_1", "A3_2","A4_2", "A4_1", "N1_1", "N1_2", "N2_1", "N2_2",  "N3_1", "N3_2", "N4"
))
AstroNeurontest@active.ident <- ordered(AstroNeurontest@active.ident, levels = 
c("NSC_1","NSC_2","RG-like","A1_1","A1_2","A2","A3_1","A3_2","A4_1","A4_2","N1_1","N1_2","N2_1","N2_2","N3_1","N3_2","N4" ))

load("~/Documents/Pascal/TBIwork/nscv3/AstroNeurontestV5.RData")


testAstroneuroMarkers <- FindAllMarkers(AstroNeurontest, only.pos = TRUE, min.pct = 0.4, thresh.use = 0.4)
top10 <- testAstroneuroMarkers %>% group_by(cluster) %>% top_n(2, avg_logFC)
DoHeatmap(obj = AstroNeurontest, features  = top10$gene, draw.lines = FALSE)+ theme(axis.text.y  = element_text(size=12, face = 4, family = "Times New Roman"),axis.text.x  = element_text(size=12, face = "bold", family = "Times New Roman"), axis.title.y = element_text(size=12, face = "bold", family = "Times New Roman"), axis.title.x = element_text(size=12, face = "bold", family = "Times New Roman"), plot.title = element_text(hjust = 0.5, size=12, face = "bold", family = "Times New Roman"), legend.position = "right", panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.line = element_blank()
)

DoHeatmap(obj = AstroNeurontest, features  = (testAstroneuroMarkers[ which ( testAstroneuroMarkers$gene %in% firstup ( tolower(TFs$V1) ) ),]  %>% group_by(cluster) %>% top_n(4, avg_logFC))$gene , draw.lines = FALSE)+ theme(axis.text.y  = element_text(size=12, face = 4, family = "Times New Roman"),axis.text.x  = element_text(size=12, face = "bold", family = "Times New Roman"), axis.title.y = element_text(size=12, face = "bold", family = "Times New Roman"), axis.title.x = element_text(size=12, face = "bold", family = "Times New Roman"), plot.title = element_text(hjust = 0.5, size=12, face = "bold", family = "Times New Roman"), legend.position = "right", panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.line = element_blank()
)
#save(AstroNeurontest, testAstroneuroMarkers, ucolsAstroNeuro, file = "~/Documents/Pascal/TBIwork/nscv3/Astroneurotest.RData", version = 2 )

var <- c("Ascl1", "Aurkb", "Ube2c", "Kif4", "Ung", "Mxra7", "Wnt8b",  "Frzb", "Lpar1" , "Fabp7", "Slc1a3",  "Fgfr3", "Syne1" , "Aldoc",  "Glul", "Atp1b2",  "Ogt", "Fam107a",  "Gja1", "Mfge8" , "Aldh1l1", "Gfap", "Agt", "Unc13c", "Id3", "Gabrg1", "Gria1",
         "Pax6", "Snap25", "Sox11", "Sox4" ,"Nnat",  "Syt1", "Celf4", "L1cam" , "Neurod1", "Neurog2", "Tubb3", "Tubb5", "Rbfox3", "Myt1l", "Sparc","Sned1","Cadps", "Neat1", "Hapln1", "Cdh2", "Gfap")
DoHeatmap(obj = AstroNeurontest, features  = var , draw.lines = FALSE)+ theme(axis.text.y  = element_text(size=12, face = 4, family = "Times New Roman"),axis.text.x  = element_text(size=12, face = "bold", family = "Times New Roman"), axis.title.y = element_text(size=12, face = "bold", family = "Times New Roman"), axis.title.x = element_text(size=12, face = "bold", family = "Times New Roman"), plot.title = element_text(hjust = 0.5, size=12, face = "bold", family = "Times New Roman"), legend.position = "right", panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.line = element_blank()
                                                                              )
VlnPlot(AstroNeurontest, c("Ascl1", "Aurkb", "Ube2c", "Neurog2","Frzb", "Fabp7", "Lpar1"), cols = as.character(ucols), pt.size = 0.1, ncol = 1  )
VlnPlot(AstroNeurontest, c("Slc1a3", "Aldoc", "Fgfr3", "Ogt", "Fam107a","Sparc", "Sned1", "Neat1", "Hapln1"), cols = as.character(ucols), pt.size = 0.1, ncol = 1 )

VlnPlot(AstroNeurontest, c("Snap25", "Syt1", "L1cam", "Celf4","Neurod1", "Sox11", "Sox4", "Gria1"), cols = as.character(ucolsAstroNeuro), pt.size = 0.1, ncol = 1 )

top10 <- testAstroneuroMarkers[which(testAstroneuroMarkers$cluster %in% c("N1_2", "N2_1","N2_2") ),] %>% group_by(cluster) %>% top_n(10, avg_logFC)
DoHeatmap(obj = AstroNeurontest, features  = top10$gene, draw.lines = FALSE)+ theme(axis.text.y  = element_text(size=12, face = 4, family = "Times New Roman"),axis.text.x  = element_text(size=12, face = "bold", family = "Times New Roman"), axis.title.y = element_text(size=12, face = "bold", family = "Times New Roman"), axis.title.x = element_text(size=12, face = "bold", family = "Times New Roman"), plot.title = element_text(hjust = 0.5, size=12, face = "bold", family = "Times New Roman"), legend.position = "right", panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.line = element_blank()
)


top10 <- testAstroneuroMarkers[which(testAstroneuroMarkers$cluster %in% c("N1_2", "N2_1","N2_2") ),] %>% group_by(cluster) %>% top_n(10, avg_logFC)
DoHeatmap(obj = AstroNeurontest, features  = top10$gene, draw.lines = FALSE)+ theme(axis.text.y  = element_text(size=12, face = 4, family = "Times New Roman"),axis.text.x  = element_text(size=12, face = "bold", family = "Times New Roman"), axis.title.y = element_text(size=12, face = "bold", family = "Times New Roman"), axis.title.x = element_text(size=12, face = "bold", family = "Times New Roman"), plot.title = element_text(hjust = 0.5, size=12, face = "bold", family = "Times New Roman"), legend.position = "right", panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.line = element_blank()
)

VlnPlot(AstroNeurontest, c("Pcp4", "Dcx", "Nnat", "Tac2", "Calb2", "Sema3c", "Stmn4", "Sema5a"), cols = as.character(ucols), pt.size = 0.1, ncol = 1 )

VlnPlot(AstroNeurontest, c("Pcp4", "Tac2", "Calb2", "Tubb3", "Rfbox3"), cols = as.character(ucols), pt.size = 0.1, ncol = 1 )
VlnPlot(AstroNeurontest, c("Adamts18", "Camk2b", "Nell2", "Kitl", "Vim", "Unc5d", "Mpped1", "Ogt" ), cols = as.character(ucols), pt.size = 0.1, ncol = 2)

RidgePlot(AstroNeurontest, c("Ascl1", "Aurkb", "Ube2c", "Neurog2","Frzb", "Fabp7", "Lpar1"), pt.size = 1  )

VlnPlotFormatted <- function(obj, genename) {
  VlnPlot(obj, genename, cols = as.character(ucolsAstroNeuro), pt.size = -1, ncol = 1, combine = TRUE ) + 
    theme(axis.text.y  = element_text(size=10, face = "bold", family = "Times New Roman"),
          axis.text.x  = element_blank(), axis.title.x = element_blank(), 
          axis.title.y = element_blank() , plot.title = element_text(hjust = 0.5, size=12, face = "italic", family = "Times New Roman"), legend.position = "none", 
          panel.background = element_blank(), 
          panel.grid.major = element_blank(), panel.grid.minor = element_blank()
  )
  
}

genelist <- c("Ascl1", "Aurkb", "Ube2c", "Neurog2","Frzb", "Fabp7", "Lpar1", "Slc1a3", "Aldoc", "Fgfr3", "Ogt", "Fam107a","Sparc", "Sned1", "Neat1", "Hapln1")
genelist <- c("Snap25",   "Neurod1", "Calb2",  "Unc5d", "Vim", "Mpped1",  "Nell2", "Kitl", "Adamts18", "Sema5a", "Pcp4",  "Tac2", "Syt1", "Ogt", "Gria1", "Sox4"   )

genelist <- c("Neurod1",   "Neurod4", "Neurog2",  "Tubb3", "Myt1l", "Rbfox3",  "Syn1" , "Slc17a7", "Slc17a6", "Gfap" )

graphPlot <- list()
j = 1
for (var in genelist){
  graphPlot[[j]] <- VlnPlotFormatted(AstroNeurontest, var)
  j = j+1
}

multiplot(plotlist = graphPlot, cols = 1)


save(AstroNeurontest, testAstroneuroMarkers, file = "~/Documents/Pascal/TBIwork/nscv3/AstroNeurontestV4.RData", version = 2)


### plotting the %%%%

metadata <- FetchData(AstroNeurontest, c("ident", "sample"))

aggdata <- aggregate(metadata, by = list(metadata$sample, metadata$ident), function(x) sum(x != 0))
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

ggplot( aggdata, aes(y = Per, x = Ident,  fill = Sample)) +  geom_bar(stat="identity")+  labs(title = "", x = "cluster", y =  "%") +
  theme(axis.text.y  = element_text(size=20, face = "bold", family = "Times New Roman"),
        axis.text.x  = element_text(size=20, face = "bold", family = "Times New Roman", angle = 90),
        axis.title.y = element_text(size=20, face = "bold", family = "Times New Roman"),
        axis.title.x = element_text(size=20, face = "bold", family = "Times New Roman"),
        plot.title = element_text(hjust = 0.5, size=20, face = "bold", family = "Times New Roman"),
        legend.position = "top",
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1, linetype="solid")
  )+  geom_hline (yintercept = dim(metadata[which(metadata$sample == "TBI"),])[1] / (dim(metadata[which(metadata$sample == "TBI"),])[1] + dim(metadata[which(metadata$sample == "Sham"),])[1]) * 100, linetype="dashed", size = 1)+scale_y_continuous(breaks = c(0, 57.7, 62.9, 100) ) + geom_text(aes(label=paste0(NumberOfCells, "")), position = position_stack(vjust = 0.5), size = 8 )

###

data_list = SplitObject(AstroNeurontest, split.by = "ident")
AstroNeurontest <- PercentageFeatureSet(AstroNeurontest, pattern = "^MT-", col.name = "percent.mt")
AstroNeurontest <- SCTransform(AstroNeurontest, vars.to.regress = "percent.mt", verbose = FALSE)


##### DE genes AstroNeuro

source("~/Documents/AstrocyteDevelopment/SeuratRun.R")

obj <- AstroNeurontest

#testAllMarkers <- FindAllMarkers(obj, only.pos = FALSE, min.pct = 0.25, thresh.use = 0.0025)
DefaultAssay(AstroNeurontest) <- "SCT"
obj <-AstroNeurontest
i=1
tmp <- findTBIgenes(obj, levels(obj@active.ident)[1] )

tmp$gene <- rownames(tmp)
diseasedgenes <- cbind( tmp , cluster = i-0.25+sample(1:10, dim(tmp)[1], replace=T)/20, subtype = rep(levels(obj@active.ident)[i], dim(tmp)[1])  )
i = 2
for(var in levels(obj@active.ident)[-1] ){
  tmp <- findTBIgenes(obj, var)
  tmp$gene <- rownames(tmp)
  print(var)
  if( dim(tmp)[1] > 1)
  diseasedgenes <- rbind(diseasedgenes,  cbind(tmp, cluster = i-0.25+sample(1:10, dim(tmp)[1], replace=T)/20 , subtype = rep(var, dim(tmp)[1]))  )
  i = i+1
  rm(tmp)
}

### combined

i=1
tmp <- findTBIgenes(obj, levels(obj@active.ident)[1] )

tmp$gene <- rownames(tmp)
diseasedgenes <- cbind( tmp , cluster = i-0.25+sample(1:10, dim(tmp)[1], replace=T)/20, subtype = rep(levels(obj@active.ident)[i], dim(tmp)[1])  )
i = 2
for(var in levels(obj@active.ident)[-1] ){
  tmp <- findTBIgenes(obj, var)
  tmp$gene <- rownames(tmp)
  print(var)
  if( dim(tmp)[1] > 1)
  diseasedgenes <- rbind(diseasedgenes,  cbind(tmp, cluster = i-0.25+sample(1:10, dim(tmp)[1], replace=T)/20 , subtype = rep(var, dim(tmp)[1]))  )
  i = i+1
  rm(tmp)
}





###### with Logistic regression LR
source("~/Documents/AstrocyteDevelopment/SeuratRun.R")

DefaultAssay(AstroNeurontest) <- "SCT"
obj = AstroNeurontest
i=1
tmp <- findTBIgenes(obj, levels(obj@active.ident)[1] )

tmp$gene <- rownames(tmp)
diseasedgenesbyLR <- cbind( tmp , cluster = i-0.25+sample(1:10, dim(tmp)[1], replace=T)/20 )
i = 2
for(var in levels(obj@active.ident)[-1] ){
  tmp <- findTBIgenes(obj, var)
  tmp$gene <- rownames(tmp)
  diseasedgenesbyLR <- rbind(diseasedgenesbyLR,  cbind(tmp, cluster = i-0.25+sample(1:10, dim(tmp)[1], replace=T)/20 ) )
  i = i+1
  rm(tmp)
}


diseasedgenes = diseasedgenesbyLR

diseasedgenes$cl <- round(diseasedgenes$cluster)
#diseasedgenes$gene <- rownames(diseasedgenes)

diseasedgenes$color <- "gray"
avg_logFCth = 0.25
diseasedgenes[which(diseasedgenes$avg_logFC < -avg_logFCth & diseasedgenes$p_val_adj < 0.05),]$color <- "orange1"
diseasedgenes[which(diseasedgenes$avg_logFC > avg_logFCth & diseasedgenes$p_val_adj < 0.05),]$color <- "darkseagreen4"
par(mar = c(5, 8, 5, 5))
stripeChart(diseasedgenes, "cluster", "avg_logFC", 0, avg_logFCth, "gene", c(-1.2, 2), c(1,17) )

diseasedgenes$marker <- 0
diseasedgenes$marker <- ifelse ( ( diseasedgenes$p_val_adj < 0.05 & diseasedgenes$avg_logFC < 0 ), -1, 0)
diseasedgenes$marker <- ifelse ( ( diseasedgenes$p_val_adj < 0.05 & diseasedgenes$avg_logFC > 0 ), 1, 0)

### down and up regulated genes
upgenes <- sapply( 1:17, function(x) dim(diseasedgenes[which(diseasedgenes$cl == x & diseasedgenes$p_val_adj < 0.05 & diseasedgenes$avg_logFC < 0.25 ),])[1] )
downgenes <- sapply( 1:17, function(x) dim(diseasedgenes[which(diseasedgenes$cl == x & diseasedgenes$p_val_adj < 0.05 & diseasedgenes$avg_logFC > -0.25 ),])[1] )

par(mfrow=c(1,2), mar = c(3,2,3,3))
barplot(-downgenes, horiz = TRUE, col = ucolsAstroNeuro,  , las = 2, main = "downregulated in TBI", xlim = c(-400, 0))
barplot(upgenes, horiz = TRUE, col = ucolsAstroNeuro, names.arg = levels(AstroNeurontest), las = 2, main = "upregulated in TBI",xlim = c(0,400) )

# diseasedgenesAstroNeuro
diseasedgenesAstroNeuro <- diseasedgenes[which(diseasedgenes$color != "gray" &  diseasedgenes$p_val_adj < 0.05),]
diseasedgenesAstroNeuro$clusterName <- levels(AstroNeurontest@active.ident)[diseasedgenesAstroNeuro$cl]

save(diseasedgenesAstroNeuro, file = "~/Documents/Pascal/TBIwork/nscv3/diseasedgenesAstroNeuroV4.Rdata")

#### checking the robustness
cluster = c("A2", "A3_1", "A3_2", "A4_1" )
cluster = c("N1_2", "N2_1", "N2_2")
cluster = c("NSC_1", "NSC_2")
cluster = c("RG-like", "A1_1", "A1_2", "N1_1")
cluster = c("NSC_1", "NSC_2","RG-like", "N1_1")

cluster = c("NSC_1", "NSC_2","RG-like", "N1_1")
cluster = c("NSC_1", "NSC_2","RG-like", "A1_1", "A1_2")

cluster = c("RG-like", "A1_1", "A1_2")
cluster = c("A1_1", "A1_2")
cluster = c("RG-like", "N1_1")
cluster = c("A1_2","A4_1", "A4_2")
cluster = c("A2")

cluster = c("A1_1","A2", "A3_1", "A3_2")
cluster = c("A2", "A3_1")


cluster = c("N1_1", "N1_2", "N2_1", "N2_2", "N3_1", "N3_2", "N4")

cluster = c("N3_1", "N3_2")


cluster = c( "N3_2")

clustervector = list(c("NSC_1", "NSC_2","RG-like"), c() )

i = 1+i
cluster = levels(AstroNeurontest)[i]

DefaultAssay(AstroNeurontest) <- "integrated"
tmp3 = findTBIgenes(AstroNeurontest, cluster)
tmp3 <- tmp2[which(tmp3$p_val_adj < 0.05),]

DefaultAssay(AstroNeurontest) <- "SCT"
tmp1 = findTBIgenes(AstroNeurontest, cluster)
tmp1 <- tmp1[which(tmp1$p_val_adj < 0.05),]

tmp1$Lineage = paste0(cluster, collapse = "_")
tmp1$Gene = rownames(tmp1)


lineageGenes = rbind(tmp1, lineageGenes)



DefaultAssay(AstroNeurontest) <- "RNA"
tmp2 = findTBIgenes(AstroNeurontest, clusterN = cluster)
tmp2 <- tmp2[which(tmp2$p_val_adj < 0.05),]

intersect(intersect(rownames(tmp1), rownames(tmp2)), rownames(tmp2))

LinegaeGenes = intersect(intersect(rownames(tmp1), rownames(tmp2)), rownames(tmp2))

View(tmp2[LinegaeGenes,])

VlnPlot(AstroNeurontest, "Gfap", idents = cluster, split.by = "sample")

clusters = c("NSC_1", "NSC_2","RG-like")
DotPlot( SubsetData(object = AstroNeurontest, ident.use = clusters), features = c("Psmc3", "Psmc4", "Psmd14", "Psmd6", "Psmd12") ,  group.by = "sample") + theme(axis.text.x = element_text(angle = 90, hjust=1))

### for the final figures lineage genes

cluster = list(levels(AstroNeurontest@active.ident)[1:3], levels(AstroNeurontest@active.ident)[4:9], levels(AstroNeurontest@active.ident)[11:14], levels(AstroNeurontest@active.ident)[15:17])

lineageGenes = data.frame("p_val"=double(),
"avg_logFC"=double(),
"pct.1"=double(),
"pct.1"=double(),
"p_val_adj" =double(),
"Lineage"=character(),
"Lineage"=character(),
stringsAsFactors=FALSE)

for(var in cluster) {
    DefaultAssay(AstroNeurontest) <- "SCT"
    tmp1 = findTBIgenes(AstroNeurontest, var, threshold = 0)
    #tmp1 <- tmp1[which(tmp1$p_val_adj < 0.05),]

    tmp1$Lineage = paste0(var, collapse = "_")
    tmp1$Gene = rownames(tmp1)


    lineageGenes = rbind(tmp1, lineageGenes)
}



lineageGenes$color <- "gray"
avg_logFCth = 0.2
lineageGenes[which(lineageGenes$avg_log2FC < -avg_logFCth & lineageGenes$p_val_adj < 0.05),]$color <- "orange1"
lineageGenes[which(lineageGenes$avg_log2FC > avg_logFCth & lineageGenes$p_val_adj < 0.05),]$color <- "darkseagreen4"

lineageGenes$cl = case_when(
    lineageGenes$Lineage == unique(lineageGenes$Lineage)[1] ~ 1,
    lineageGenes$Lineage == unique(lineageGenes$Lineage)[2] ~ 2,
    lineageGenes$Lineage == unique(lineageGenes$Lineage)[3] ~ 3,
    lineageGenes$Lineage == unique(lineageGenes$Lineage)[4] ~ 4
)

lineageGenes$lineage_to_plot = case_when(
    lineageGenes$Lineage == unique(lineageGenes$Lineage)[4] ~ "NSC and RG-like cells",
    lineageGenes$Lineage == unique(lineageGenes$Lineage)[3] ~ "Neuronal progenitors",
    lineageGenes$Lineage == unique(lineageGenes$Lineage)[2] ~ "Astrocytic progenitors",
    lineageGenes$Lineage == unique(lineageGenes$Lineage)[1] ~ "Remaining astrocytes"
)

lineageGenes$cluster = lineageGenes$cl -0.1 + sample(1:100, length(lineageGenes$cl), replace = TRUE)/1000
par(mar = c(5, 15, 5, 5))
stripeChart(lineageGenes, "cluster", "avg_logFC", 0, avg_logFCth, "p_val_adj", c(-0.5,0.5), c(0.5,4.5) )

### looking into literature 
diseasedgenesAstroNeuro <- diseasedgenesAstroNeuro [which(diseasedgenesAstroNeuro$color != "gray" & diseasedgenesAstroNeuro$p_val_adj < 0.05),] 
TBIgenesliterature <- read.table("~/Documents/Pascal/TBIvsSHAM_literature10.1016slashj.ebiom.2017.01.txt", sep = "\t")
diseasedgenesAstroNeuro[ which ( diseasedgenesAstroNeuro$gene %in% intersect(diseasedgenesAstroNeuro$gene, TBIgenesliterature$V2)),]
diseasedgenesAstroNeuro[ which ( diseasedgenesAstroNeuro$gene %in% intersect(diseasedgenesAstroNeuro$gene, TBIgenesliterature$V2)),]

TFs <- read.table("~/Documents/Pascal/TF_names_v_1.01.txt")
mouseTFs <- firstup(tolower(TFs$V1))


### heatmapPlot
TFDiseaseRelevant <- diseasedgenesAstroNeuro[ which(diseasedgenesAstroNeuro$gene %in% c( mouseTFs, "Tubb4b", "Tubb3", "Tubb5", "Tubb2b", "Tubb6", "Tubb4a", "Fos", "Jun", "Junb", "Dusp", "Ascl1", "Neurod1", "Neurog2", "Sox9", "Pax6", "Myt1l", "Slc1a3", "Tubb3", "Sox4", "Sox11" ) ),]# read.table("~/Documents/Pascal/diseasedTFsAstroNeuro.txt", header = TRUE)

TFDiseaseRelevant <- diseasedgenesAstroNeuro[ which(diseasedgenesAstroNeuro$gene %in% c( "Tubb4b", "Tubb3", "Tubb5", "Tubb2b", "Tubb6", "Tubb4a", "Fos", "Jun", "Junb", "Dusp", "Ascl1", "Neurod1", "Neurog2", "Sox9", "Pax6", "Myt1l", "Slc1a3", "Tubb3", "Sox4", "Sox11" ) ),]# read.table("~/Documents/Pascal/diseasedTFsAstroNeuro.txt", header = TRUE)
TFDiseaseRelevant$gene <- as.character(TFDiseaseRelevant$gene)



TFDiseaseRelevant <- TFDiseaseRelevant[order( TFDiseaseRelevant$marker , TFDiseaseRelevant$cl, TFDiseaseRelevant$avg_logFC  ),]
rownames(TFDiseaseRelevant) <- paste0 (TFDiseaseRelevant$gene, "." , TFDiseaseRelevant$cl)

library(gplots)
datacluster <- as.data.frame(matrix(0, ncol = length(TFDiseaseRelevant$gene), nrow = length(unique(TFDiseaseRelevant$cl)) ))
colnames(datacluster) <- rownames(TFDiseaseRelevant)
rownames(datacluster) <- unique(TFDiseaseRelevant$cl)
for( var in colnames(datacluster) ) {
  print(var)
  datacluster[ as.character(TFDiseaseRelevant[var,]$cl), var] <- TFDiseaseRelevant[var,]$avg_logFC
  # datacluster[var,TFDiseaseRelevant [which(TFDiseaseRelevant$cl == var),]$gene ] <- TFDiseaseRelevant [which(TFDiseaseRelevant$cl == var),]$avg_logFC
}
#datacluster$cl <- NULL
#datacluster$cl <- as.integer(rownames(datacluster))
#datacluster <- datacluster[order( datacluster$cl ),]
#datacluster$cl <- NULL
#datacluster <- datacluster[,colSums(datacluster[, -1])>0 ]


heatmap.2(t(as.matrix(datacluster)), dendrogram='none', Rowv=FALSE, Colv=FALSE,trace='none', 
          density.info = 'none', key = TRUE, keysize = 0.8, xlab = NULL,col=colorRampPalette(c( "orange1", "gray89", "darkseagreen4"))(n = 30),
          margins = c(7,7), family = "Times", cexRow = 1, labRow=as.expression(lapply(TFDiseaseRelevant$gene, function(a) bquote(italic(.(a))))),
          labCol =  as.character(unique(TFDiseaseRelevant$clusterName)) )



#### monocle
DefaultAssay(spatialfinal) <- "integrated"
Neuro  <- SubsetData(spatialfinal, ident.use =  c("S0", "S1" ,"S10" ,"S9"  ,"S8" , "S7" ,"S6" , "S5" , "S4" , "S3" , "S2"  ))
Astro <- SubsetData(spatialfinal, ident.remove =  c( "S10" ,"S9"  ,"S8" , "S7" ,"S6" , "S5" , "S4" , "S3" , "S2"  ))
#Astro <- SubsetData(AstroNeurontest, ident.use =  levels(AstroNeurontest)[c(12,11,10)] ) for ast4 pseudotime

Astro.cds <- as.cell_data_set(Astro)
Astro.cds <- cluster_cells(cds = Astro.cds, reduction_method = "UMAP")
Astro.cds <- learn_graph(Astro.cds, use_partition = TRUE)

Neuro.cds <- as.cell_data_set(Neuro)
Neuro.cds <- cluster_cells(cds = Neuro.cds, reduction_method = "UMAP")
Neuro.cds <- learn_graph(Neuro.cds, use_partition = TRUE)

# order cells
Astro.cds <- order_cells(Astro.cds, reduction_method = "UMAP", root_cells = WhichCells(spatialfinal, idents = c("S0", "S1") ) )
Neuro.cds <- order_cells(Neuro.cds, reduction_method = "UMAP", root_cells = WhichCells(spatialfinal, idents = c("S0", "S1") ) )


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


spatialfinal <- AddMetaData(
  object = spatialfinal,
  metadata = Astro.cds@principal_graph_aux@listData$UMAP$pseudotime,
  col.name = "AstroLineage"
)


spatialfinal <- AddMetaData(
  object = spatialfinal,
  metadata = Neuro.cds@principal_graph_aux@listData$UMAP$pseudotime,
  col.name = "NeuroLineage"
)



FeaturePlot(spatialfinal, c("AstroLineage", "NeuroLineage"), pt.size = 0.1) & scale_color_viridis_c()


Fetchdata = FetchData(spatialfinal, c("NeuroLineage", "AstroLineage") )

seurat_obj_list$A1control@meta.data$NeuroLineage = 0
seurat_obj_list$A1control@meta.data[intersect(colnames(Astro.cds),colnames(seurat_obj_list$A1control)),]$NeuroLineage = as.character(Fetchdata[intersect(colnames(Astro.cds),colnames(seurat_obj_list$A1control)),]$NeuroLineage)
FeaturePlot(seurat_obj_list$A1control,  c("NeuroLineage"))

seurat_obj_list$A1control <- AddMetaData(
  object = seurat_obj_list$A1control,
  metadata = Astro.cds@principal_graph_aux@listData$UMAP$pseudotime[intersect(colnames(Astro.cds),colnames(seurat_obj_list$A1control))],
  col.name = "AstroLineage"
)

seurat_obj_list$A1control <- AddMetaData(
  object = seurat_obj_list$A1control,
  metadata = Neuro.cds@principal_graph_aux@listData$UMAP$pseudotime[intersect(colnames(Neuro.cds),colnames(seurat_obj_list$A1control))],
  col.name = "NeuroLineage"
)
FeaturePlot(seurat_obj_list$A1control,  c("AstroLineage", "NeuroLineage"), reduction = "coord", pt.size = 2.9) &  scale_colour_gradient( low = "blue", high = "khaki1",
                                                                                                                   space = "Lab", na.value = t_col("gray95"), guide = "colourbar",
                                                                                                                   aesthetics = "colour")


seurat_obj_list$A2control <- AddMetaData(
  object = seurat_obj_list$A2control,
  metadata = Astro.cds@principal_graph_aux@listData$UMAP$pseudotime[intersect(colnames(Astro.cds),colnames(seurat_obj_list$A2control))],
  col.name = "AstroLineage"
)

seurat_obj_list$A2control <- AddMetaData(
  object = seurat_obj_list$A2control,
  metadata = Neuro.cds@principal_graph_aux@listData$UMAP$pseudotime[intersect(colnames(Neuro.cds),colnames(seurat_obj_list$A2control))],
  col.name = "NeuroLineage"
)
FeaturePlot(seurat_obj_list$A2control,  c("AstroLineage", "NeuroLineage"), reduction = "coord", pt.size = 2.9) &  scale_colour_gradient( low = "blue", high = "khaki1",space = "Lab", na.value = t_col("gray95"), guide = "colourbar",
                                                                                                                                         aesthetics = "colour")


seurat_obj_list$A1tbi <- AddMetaData(
  object = seurat_obj_list$A1tbi,
  metadata = Astro.cds@principal_graph_aux@listData$UMAP$pseudotime[intersect(colnames(Astro.cds),colnames(seurat_obj_list$A1tbi))],
  col.name = "AstroLineage"
)

seurat_obj_list$A1tbi <- AddMetaData(
  object = seurat_obj_list$A1tbi,
  metadata = Neuro.cds@principal_graph_aux@listData$UMAP$pseudotime[intersect(colnames(Neuro.cds),colnames(seurat_obj_list$A1tbi))],
  col.name = "NeuroLineage"
)
FeaturePlot(seurat_obj_list$A1tbi,  c("AstroLineage", "NeuroLineage"), reduction = "coord", pt.size = 2.9) &  scale_colour_gradient( low = "blue", high = "khaki1",space = "Lab", na.value = t_col("gray95"), guide = "colourbar",
                                                                                                                                         aesthetics = "colour")


seurat_obj_list$A2tbi <- AddMetaData(
  object = seurat_obj_list$A2tbi,
  metadata = Astro.cds@principal_graph_aux@listData$UMAP$pseudotime[intersect(colnames(Astro.cds),colnames(seurat_obj_list$A2tbi))],
  col.name = "AstroLineage"
)

seurat_obj_list$A2tbi <- AddMetaData(
  object = seurat_obj_list$A2tbi,
  metadata = Neuro.cds@principal_graph_aux@listData$UMAP$pseudotime[intersect(colnames(Neuro.cds),colnames(seurat_obj_list$A2tbi))],
  col.name = "NeuroLineage"
)
FeaturePlot(seurat_obj_list$A2tbi,  c("AstroLineage", "NeuroLineage"), reduction = "coord", pt.size = 2.9) &  scale_colour_gradient( low = "blue", high = "khaki1",space = "Lab", na.value = t_col("gray95"), guide = "colourbar",
                                                                                                                                     aesthetics = "colour")

                                                                                                                                         



Neuro  <- SubsetData(spatialfinal, ident.use =  c("S0", "S1" ,"S10" ,"S9"  ,"S8" , "S7" ,"S6" , "S5" , "S4" , "S3" , "S2"  ))

fneurodata = FetchData(Neuro, c("Mpped1", "NeuroLineage" ))
fneurodata$condition = substr(as.character( unlist( strsplit ( rownames(fneurodata), "=") )) [seq(1, 2*length(rownames(fneurodata)), 2) ], 1, 3)
## develop a model
fneurodata$pred1 <- predict(lm(Mpped1 ~ poly(NeuroLineage, 3), data=fneurodata))

## look at the data

ggplot(fneurodata, aes(x = NeuroLineage, y=Mpped1, color=condition)) +
  geom_point( ) +
  geom_hline(aes(yintercept=0))


p1 <- ggplot(fneurodata, aes(x = NeuroLineage, y=Snap25, color=condition)) +
  geom_point( ) +
  geom_hline(aes(yintercept=0))+scale_color_gradient( low = "blue", high = "khaki1", space ="Lab" , na.value = t_col("gray95"), guide = "colourbar",) +  theme(axis.text.y  = element_text(size=15, face = "bold", family = "Times New Roman"),
                                                                                                                                                                axis.text.x  = element_text(size=15, face = "bold", family = "Times New Roman"),
                                                                                                                                                                axis.title.y = element_text(size=15, face = "bold", family = "Times New Roman"),
                                                                                                                                                                axis.title.x = element_text(size=15, face = "bold", family = "Times New Roman"),
                                                                                                                                                                plot.title = element_text(hjust = 0.5, size=15, face = "bold", family = "Times New Roman"),
                                                                                                                                                                legend.position = "right",
                                                                                                                                                                panel.background = element_blank(),
                                                                                                                                                                panel.grid.major = element_blank(),
                                                                                                                                                                panel.grid.minor = element_blank(),
                                                                                                                                                                axis.line = element_blank(),
                                                                                                                                                                panel.border = element_rect(colour = "black", fill=NA, size=1, linetype="solid")
  )+
  geom_line(aes(y = pred1), color="red")

print(p1)

## check the model








##############  monocle




#### monocle
DefaultAssay(AstroNeurontest) <- "integrated"
Neurons <- SubsetData(AstroNeurontest, ident.use = c("NSC_1", "NSC_2", "N1_1"  ,  "N1_2"  ,  "N2_1"  ,  "N2_2"   , "N3_1"    ,"N3_2"   ) )
Astrocytes <- SubsetData(AstroNeurontest, ident.use = c("NSC_1", "NSC_2", "RG-like","A1_1"  ,    "A2" ,     "A3_1"  ,  "A3_2", "A1_2"  , "A4_1"  ,  "A4_2"    ) )

Astrocytes.cds <- as.cell_data_set(Astrocytes)
Astrocytes.cds <- preprocess_cds(Astrocytes.cds,method = "PCA" )
Astrocytes.cds <- reduce_dimension(Astrocytes.cds, reduction_method = "UMAP", preprocess_method = 'PCA')
Astrocytes.cds <- cluster_cells(cds = Astrocytes.cds, reduction_method = "UMAP")
Astrocytes.cds <- learn_graph(Astrocytes.cds, use_partition = TRUE)

Neurons.cds <- preprocess_cds(Neurons.cds,method = "PCA" )
Neurons.cds <- reduce_dimension(Neurons.cds, reduction_method = "UMAP")

Neurons.cds <- as.cell_data_set(Neurons)
Neurons.cds <- cluster_cells(cds = Neurons.cds, reduction_method = "UMAP")
Neurons.cds <- learn_graph(Neurons.cds, use_partition = TRUE)

AstroNeurontest.cds <- as.cell_data_set(AstroNeurontest)

AstroNeurontest.cds <- preprocess_cds(AstroNeurontest.cds,method = "PCA" )
AstroNeurontest.cds <- reduce_dimension(AstroNeurontest.cds, reduction_method = "UMAP")
AstroNeurontest.cds <- as.cell_data_set(AstroNeurontest)
AstroNeurontest.cds <- cluster_cells(cds = AstroNeurontest.cds, reduction_method = "UMAP")
AstroNeurontest.cds <- learn_graph(AstroNeurontest.cds, use_partition = TRUE)


# order cells
Astrocytes.cds <- order_cells(Astrocytes.cds, reduction_method = "UMAP", root_cells = WhichCells(AstroNeurontest, idents = c("NSC_1", "NSC_2") ) )
Neurons.cds <- order_cells(Neurons.cds, reduction_method = "UMAP", root_cells = WhichCells(AstroNeurontest, idents = c("N1_1") ) )
AstroNeurontest.cds <- order_cells(AstroNeurontest.cds, reduction_method = "UMAP", root_cells = WhichCells(AstroNeurontest, idents = c("NSC_1", "NSC_2") ) )


plot_cells(
  cds = Astrocytes.cds,
  color_cells_by = "pseudotime",
  label_cell_groups=FALSE,
  label_leaves=TRUE,
  label_branch_points=TRUE,
  show_trajectory_graph = TRUE
)


AstroNeurontest <- AddMetaData(
  object = AstroNeurontest,
  metadata = Astrocytes.cds@principal_graph_aux@listData$UMAP$pseudotime,
  col.name = "AstroLineage"
)


AstroNeurontest <- AddMetaData(
  object = AstroNeurontest,
  metadata = Neurons.cds@principal_graph_aux@listData$UMAP$pseudotime,
  col.name = "NeuroLineage"
)

AstroNeurontest <- AddMetaData(
  object = AstroNeurontest,
  metadata = AstroNeurontest.cds@principal_graph_aux@listData$UMAP$pseudotime,
  col.name = "Lineages"
)






FeaturePlot(AstroNeurontest, c( "AstroLineage"), pt.size = 0.1) +scale_color_gradient( low = "blue", high = "khaki1", space ="Lab" , na.value = t_col("gray77"), guide = "colourbar",)  +  theme(axis.text.y  = element_text(size=15, face = "bold", family = "Times New Roman"),
         axis.text.x  = element_text(size=15, face = "bold", family = "Times New Roman"),
         axis.title.y = element_text(size=15, face = "bold", family = "Times New Roman"),
         axis.title.x = element_text(size=15, face = "bold", family = "Times New Roman"),
         plot.title = element_text(hjust = 0.5, size=15, face = "bold", family = "Times New Roman"),
         legend.position = "right",
         panel.background = element_blank(),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         axis.line = element_blank(),
         panel.border = element_rect(colour = "black", fill=NA, size=1, linetype="solid")
)


Neurons <- SubsetData(AstroNeurontest, ident.use = c("NSC_1", "NSC_2", "RG-like","N1_1"  ,  "N1_2"  ,  "N2_1"  ,  "N2_2"   , "N3_1"    ,"N3_2" , "N4"  ) )
Astrocytes <- SubsetData(AstroNeurontest, ident.remove = c("A4_1", "A4_2", "N1_1"  ,  "N1_2"  ,  "N2_1"  ,  "N2_2"   , "N3_1"    ,"N3_2" , "N4"  ) )

gene = "Neurog2"
fneurodata = FetchData(Neurons, c(gene, "NeuroLineage", "sample" ))
fneurodata$condition = fneurodata$sample
fneurodata$pred1 <- predict(lm(fneurodata[,gene] ~ poly(NeuroLineage, 3), data=fneurodata))

## look at the data

ggplot(fneurodata, aes(x = NeuroLineage, y=fneurodata[,gene], color=condition)) +
  geom_point( ) + ylab(gene)+
  geom_hline(aes(yintercept=0))+theme_classic()


 ggplot(fneurodata, aes(x = NeuroLineage, y=fneurodata[,gene], colour = NeuroLineage)) +
  geom_point( ) + ylab(gene)+
  geom_hline(aes(yintercept=0))+scale_color_gradient( low = "blue", high = "khaki1", space ="Lab" , na.value = t_col("gray95"), guide = "colourbar",) +  theme(axis.text.y  = element_text(size=15, face = "bold", family = "Times New Roman"),
                                                                                                                                                               axis.text.x  = element_text(size=15, face = "bold", family = "Times New Roman"),
                                                                                                                                                               axis.title.y = element_text(size=15, face = "bold", family = "Times New Roman"),
                                                                                                                                                               axis.title.x = element_text(size=15, face = "bold", family = "Times New Roman"),
                                                                                                                                                               plot.title = element_text(hjust = 0.5, size=15, face = "bold", family = "Times New Roman"),
                                                                                                                                                               legend.position = "right",
                                                                                                                                                               panel.background = element_blank(),
                                                                                                                                                               panel.grid.major = element_blank(),
                                                                                                                                                               panel.grid.minor = element_blank(),
                                                                                                                                                               axis.line = element_blank(),
                                                                                                                                                               panel.border = element_rect(colour = "black", fill=NA, size=1, linetype="solid")
  )+
  geom_line(aes(y = pred1), color="red")


 
 
 
 gene = "Fgfr3"
 fastrodata = FetchData(Astrocytes, c(gene, "AstroLineage", "sample" ))
 fastrodata$condition = fastrodata$sample
 fastrodata <- fastrodata[!is.infinite(fastrodata$AstroLineage),]
 
 fastrodata$pred1 <- predict(lm(fastrodata[,gene] ~ poly(AstroLineage,4), data=fastrodata))
 
 ## look at the data
 
 ggplot(fastrodata, aes(x = AstroLineage, y=fastrodata[,gene], color=condition)) +
   geom_point( ) + ylab(gene)+
 geom_hline(aes(yintercept=0))
 
 
 ggplot(fastrodata, aes(x = AstroLineage, y=fastrodata[,gene], colour = AstroLineage)) +
   geom_point( ) + ylab(gene)+
   geom_hline(aes(yintercept=0))+scale_color_gradient( low = "blue", high = "khaki1", space ="Lab" , na.value = t_col("gray95"), guide = "colourbar",) +  theme(axis.text.y  = element_text(size=15, face = "bold", family = "Times New Roman"),
                                                                                                                                                                axis.text.x  = element_text(size=15, face = "bold", family = "Times New Roman"),
                                                                                                                                                                axis.title.y = element_text(size=15, face = "bold", family = "Times New Roman"),
                                                                                                                                                                axis.title.x = element_text(size=15, face = "bold", family = "Times New Roman"),
                                                                                                                                                                plot.title = element_text(hjust = 0.5, size=15, face = "bold", family = "Times New Roman"),
                                                                                                                                                                legend.position = "right",
                                                                                                                                                                panel.background = element_blank(),
                                                                                                                                                                panel.grid.major = element_blank(),
                                                                                                                                                                panel.grid.minor = element_blank(),
                                                                                                                                                                axis.line = element_blank(),
                                                                                                                                                                panel.border = element_rect(colour = "black", fill=NA, size=1, linetype="solid")
   )+ geom_line(aes(y = pred1), color="red")
 
 
vargenes = c( "Ascl1", "Aurkb", "Ube2c", "Neurog2", "Frzb", "Fabp7", "Lpar1" ,
"Aldoc", "Ogt", "Hapln1", "Neat1", "Sparc", "Sned1", "Fgfr3" ,
"Sox4", "Neurod1", "Calb2", "Unc5d", "Mpped1", "Nell2", "Syt1")
 


### GO analysis 

clusterN = "NSC_1"
markerN = 0

diseasedgenesAstroNeuroSubset = diseasedgenesAstroNeuro[which(diseasedgenesAstroNeuro$clusterName == clusterN & diseasedgenesAstroNeuro$marker  == markerN),]

genelist <- queryMany(diseasedgenesAstroNeuroSubset$gene, scopes="symbol", fields="entrezgene", species="mouse")

go <- goana(genelist@listData$entrezgene, species="Mm")
topGO(go, n=15)

Ntop = 100
termsOfInterest = "[S|s]ynapse|[A|a]xon|[N|n]euron|[N|n]ervous|[G|g]lia|[A|a]strocyt|[M|m]icroglia|[P|p]arkinson|[A|a]lzheimer|[E|e]pilepcy|TNF|[N|n]euro|[A|a]ge"

topGO(go, n = Ntop)[grepl(termsOfInterest, topGO(go, n = Ntop)$Term) & topGO(go, n = Ntop)$P.DE <0.05,]


GOgenes <- list()
for(var in rownames(topGO(go, n = Ntop)[grepl(termsOfInterest, topGO(go, n = Ntop)$Term) & topGO(go, n = Ntop)$P.DE <0.05,]) )
{
  GOTerms <- var 
  allegs = get(GOTerms, org.Mm.egGO2ALLEGS)
  genes = unlist(mget(allegs,org.Mm.egSYMBOL))
  GOgenes[var] <- list( var = intersect(  diseasedgenesAstroNeuroSubset$gene,  as.character(genes) ) )
}
cbind(topGO(go, n = Ntop)[grepl(termsOfInterest, topGO(go, n = Ntop)$Term) & topGO(go, n = Ntop)$P.DE <0.05,], as.character(GOgenes) )





 
### TGF-smad genes
TGF_beta_genes = c("Chrd", "Nog", "Nbl1", "Thbs1", "Dcn", "Lefty", "Fst", "Bmp2", "Bmp4", "Bmp5", "Bmp6", "Bmp7", "Bmp8", "Gdf5", "Gdf6", "Gdf7", "Amh", "Ltbp1", "Tgfb1", "Tgfb2", "Tgfb3", "Inhbb", "Inhbc", "Inhbe", "Nodal", "Bmpr2", "Amhr2", "Tgfbr2", "Acvr2a", "Acvr2b", "Bmpr1a", "Bmpr1b", "Acvr1", "Tgfbr1", "Acvr1b", "Acvr1c", "Bambi", "Smad1", "Smad5", "Smad9", "Smad2", "Smad4", "Smad6", "Smad7", "Smurf", "Madhip", "Id1", "Id2", "Id3", "Id4", "Emc", "Rbl1", "E2f4", "Tfdp1", "Ep300", "Sp1", "Tgif1", "Tgif2", "Myc", "Cdkn2b", "Pitx2", "Rbx1", "Cul1", "Skp1", "Erk", "Ifng", "Tnf", "Rhoa", "Rock1", "Ppp2r1", "Ppp2c", "Rps6kb")

clusters = c("NSC_1", "NSC_2", "RG-like","A1_1"  ,  "A1_2"  ,  "A2" ,     "A3_1"  ,  "A3_2", "A4_1"  ,  "A4_2"    )
clusters = c("NSC_1", "NSC_2", "N1_1"  ,  "N1_2"  ,  "N2_1"  ,  "N2_2"   , "N3_1"    ,"N3_2"  , "N4" )


bioasedfile = read.csv("/Volumes/Transcend/Pascal/TBI_Biased.xlsx - Sheet 1 - TBI_BiasedNew2.csv", header = TRUE, stringsAsFactors = FALSE)

genes =  intersect( markers$gene, unlist(strsplit( bioasedfile[grepl( "neuroblast proliferation", bioasedfile$Description),]$geneID, ", ")))

genes =  intersect( markers[which(markers$cluster %in% clusters),]$gene, unlist(strsplit( bioasedfile[grepl( "astrocyte", bioasedfile$Description),]$geneID, ", ")))

DotPlot(AstroNeurontest, features = unique(genes), idents = clusters)+ theme(axis.text.x = element_text(angle = 90, hjust=1))





### violin plots
adjpval =  list()
for( var in levels(Idents(AstroNeurontest))) {
    
    adjpval[[var]] =  MASTtest(AstroNeurontest, var)
}

dataToPlot = FetchData(AstroNeurontest, c("Ppp1r14b", "idents", "sample"), slot = "data")
ylim = max(dataToPlot$Ppp1r14b)+0.01
var = names(adjpval)[1]
ggplot(dataToPlot[which(dataToPlot$idents == var),], aes(x = sample, y = Ppp1r14b)) +
geom_violin(aes(color=sample, fill=sample), trim = TRUE) +
geom_jitter(height = 0, width = 0.1) +
theme_classic() + geom_segment(aes(x=1,xend=2,y=ylim,yend=ylim)) + geom_text(aes(x=1.5,y=ylim+0.5), label =adjpval[var], )
