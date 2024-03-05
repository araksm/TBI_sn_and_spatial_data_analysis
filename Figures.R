
## Run with Seurat v3.2.0

#setwd("/Users/u0114327/Dropbox/")
setwd("/Volumes/Transcend/Pascal/Figures/")
load("RevisionOfAnnotations.RData")

Combinedtest = SetIdent(Combinedtest, cells = WhichCells(Combinedtest, idents = c("N1", "N2", "N3", "N4")), value = "Neuronal lineage")
Combinedtest = SetIdent(Combinedtest, cells = WhichCells(Combinedtest, idents =  c("A1", "A2", "A3", "A4")), value = "Astrocyte lineage")

AstroNeurontest = SetIdent(AstroNeurontest, cells = WhichCells(AstroNeurontest, idents = c("A1_1")), value = "A-stage1")
AstroNeurontest = SetIdent(AstroNeurontest, cells = WhichCells(AstroNeurontest, idents = c("A2")), value = "A-stage2")
AstroNeurontest = SetIdent(AstroNeurontest, cells = WhichCells(AstroNeurontest, idents = c("A3_1")), value = "A-stage3")
AstroNeurontest = SetIdent(AstroNeurontest, cells = WhichCells(AstroNeurontest, idents = c("A3_2")), value = "A-stage4")
AstroNeurontest = SetIdent(AstroNeurontest, cells = WhichCells(AstroNeurontest, idents = c("A1_2")), value = "A-stage5")
AstroNeurontest = SetIdent(AstroNeurontest, cells = WhichCells(AstroNeurontest, idents = c("A4_2")), value = "A-stage6")
AstroNeurontest = SetIdent(AstroNeurontest, cells = WhichCells(AstroNeurontest, idents = c("A4_1")), value = "A-stage7")
AstroNeurontest = SetIdent(AstroNeurontest, cells = WhichCells(AstroNeurontest, idents = c("N1_1")), value = "N-stage1")
AstroNeurontest = SetIdent(AstroNeurontest, cells = WhichCells(AstroNeurontest, idents = c("N1_2")), value = "N-stage2")
AstroNeurontest = SetIdent(AstroNeurontest, cells = WhichCells(AstroNeurontest, idents = c("N2_1")), value = "N-stage3")
AstroNeurontest = SetIdent(AstroNeurontest, cells = WhichCells(AstroNeurontest, idents = c("N2_2")), value = "N-stage4")
AstroNeurontest = SetIdent(AstroNeurontest, cells = WhichCells(AstroNeurontest, idents = c("N3_1")), value = "N-stage5")
AstroNeurontest = SetIdent(AstroNeurontest, cells = WhichCells(AstroNeurontest, idents = c("N3_2")), value = "N-stage6")
AstroNeurontest = SetIdent(AstroNeurontest, cells = WhichCells(AstroNeurontest, idents = c("N4")), value = "N-stage7")
AstroNeurontest = SetIdent(AstroNeurontest, cells = WhichCells(AstroNeurontest, idents = c("NSC_1")), value = "NSC-stage1")
AstroNeurontest = SetIdent(AstroNeurontest, cells = WhichCells(AstroNeurontest, idents = c("NSC_2")), value = "NSC-stage2")


#genes.viz = c(genes.viz[1:19], genes.viz[22:25],genes.viz[c(20,21,36)], genes.viz[26:35],genes.viz[37:63] )

AstroNeurontest@active.ident = factor(x=AstroNeurontest@active.ident, levels = c("NSC-stage1", "NSC-stage2", "RG-like", "N-stage1","N-stage2", "N-stage3", "N-stage4", "N-stage5", "N-stage6", "N-stage7", "A-stage1","A-stage2","A-stage3","A-stage4","A-stage5","A-stage6","A-stage7") )
Neurons = SubsetData( AstroNeurontest, ident.use = levels(AstroNeurontest@active.ident)[1:10])
Astros = SubsetData( AstroNeurontest, ident.use = levels(AstroNeurontest@active.ident)[c(1,2,3, 11:17)])


Combinedtest = SetIdent(Combinedtest, cells = WhichCells(Combinedtest, idents = c("Neuronal lineage")), value = "Neuron-like cells")
Combinedtest = SetIdent(Combinedtest, cells = WhichCells(Combinedtest, idents =  c("Astrocyte lineage")), value = "Astrocyte-like cells")


library(ggplot2)
library(Seurat)
fontsize = 15

#emf(file = "majorCellTypeClusters.emf",width = 7, height = 7)
pdf(file = "majorCellTypeClusters.pdf",   # The directory you want to save the file in
    width = 580, # The width of the plot in inches
    height = 351) # The height of the plot in inches

DimPlot(Combinedtest) + theme(axis.text.y  = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
axis.text.x  = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
axis.title.y = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
axis.title.x = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
plot.title = element_text(hjust = 0.5, size=fontsize, face = "bold", family = "Times New Roman"),
legend.position = "right",
panel.background = element_blank(),
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
axis.line = element_blank(),
panel.border = element_rect(colour = "black", fill=NA, size=1, linetype="solid")
)

dev.off()
par(mar = c(5, 5, 5, 10))

pdf(file = "majorCellTypeHeatmap.pdf",   # The directory you want to save the file in
    width = 1527#1492, # The width of the plot in inches
    height = 793) #830) # The height of the plot in inches
DoHeatmap(Combinedtest, features = genes.viz)+ theme(axis.text.y  = element_text(size=12, face = 4, family = "Times New Roman"),axis.text.x  = element_text(size=12, face = "bold", family = "Times New Roman"), axis.title.y = element_text(size=12, face = "bold", family = "Times New Roman"), axis.title.x = element_text(size=12, face = "bold", family = "Times New Roman"), plot.title = element_text(hjust = 0.5, size=12, face = "bold", family = "Times New Roman"), legend.position = "right", panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.line = element_blank())

dev.off()


pdf(file = "AstroNeuroCellTypeClusters.pdf",   # The directory you want to save the file in
    width = 642, # The width of the plot in inches
    height = 425) # The height of the plot in inches

fontsize = 15
DimPlot(AstroNeurontest, label = TRUE) + theme(axis.text.y  = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
axis.text.x  = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
axis.title.y = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
axis.title.x = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
plot.title = element_text(hjust = 0.5, size=fontsize, face = "bold", family = "Times New Roman"),
legend.position = "right",
panel.background = element_blank(),
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
axis.line = element_blank(),
panel.border = element_rect(colour = "black", fill=NA, size=1, linetype="solid")
)

dev.off()


pdf(file = "AstroPseudotime.pdf",   # The directory you want to save the file in
    width = 543, # The width of the plot in inches
    height = 436) # The height of the plot in inches

FeaturePlot(AstroNeurontest, features = c("AstroLineage"), cols = c("yellow2", "red2")) + theme(axis.text.y  = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
axis.text.x  = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
axis.title.y = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
axis.title.x = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
plot.title = element_text(hjust = 0.5, size=fontsize, face = "bold", family = "Times New Roman"),
legend.position = "right",
panel.background = element_blank(),
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
axis.line = element_blank(),
panel.border = element_rect(colour = "black", fill=NA, size=1, linetype="solid")
)

dev.off()



pdf(file = "NeuroPseudotime.pdf",   # The directory you want to save the file in
    width = 543, # The width of the plot in inches
    height = 436) # The height of the plot in inches

FeaturePlot(AstroNeurontest, features = c("NeuroLineage"), cols = c("yellow2", "red2")) + theme(axis.text.y  = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
axis.text.x  = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
axis.title.y = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
axis.title.x = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
plot.title = element_text(hjust = 0.5, size=fontsize, face = "bold", family = "Times New Roman"),
legend.position = "right",
panel.background = element_blank(),
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
axis.line = element_blank(),
panel.border = element_rect(colour = "black", fill=NA, size=1, linetype="solid")
)

dev.off()


pdf(file = "A4Pseudotime.pdf",   # The directory you want to save the file in
    width = 543, # The width of the plot in inches
    height = 436) # The height of the plot in inches

FeaturePlot(AstroNeurontest, features = c("A4Lineage"), cols = c("yellow2", "red2")) + theme(axis.text.y  = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
axis.text.x  = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
axis.title.y = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
axis.title.x = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
plot.title = element_text(hjust = 0.5, size=fontsize, face = "bold", family = "Times New Roman"),
legend.position = "right",
panel.background = element_blank(),
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
axis.line = element_blank(),
panel.border = element_rect(colour = "black", fill=NA, size=1, linetype="solid")
)

dev.off()


term = "glial cell development"
genes_glial_cell_development = strsplit("Id4, Kcnj10, Ntrk2, Gfap, Cspg5, Phgdh, Lamb2, Clu, Mt3, Ascl1, Vim, Hes5, Dag1, Smo, Mt3, Vim, Phgdh, Ntrk2, Clu, Lamb2, Id4, Gfap, Ascl1, Fgfr3, Dmd, Cspg5, Ntrk2, Lgi4, Plp1, Lrp1, Ldlr, Wasf3, Tcf7l2, Agt, Fgfr3, Omg, Cspg5, Clu, Ntrk2, Dmd, Mt3, Hes5, Lgi4, Dag1, Lamb2, Olig1, Phgdh, Id4, Egfr, Nr1d1, Kcnj10, Lrp1, Tcf7l2, Ldlr, Wasf3, Pou3f2, Omg, Cspg5, Fgfr3, Clu, Dmd, Pmp22, Id4, Mt3, Lgi4, Gfap, Ntrk2, Wasf3, Phgdh, Agt, Dag1, Ldlr, Hes5, Pou3f2, Akt2, Id2, Lrp1, Clu, Mt3, Cspg5, Ntrk2, Id4, Gfap, Phgdh, Fgfr3, Plp1, Kcnj10, Mt3, Id4, Phgdh, Gfap, Clu, Id2, Nr1d1, Akt2", ", " )[[1]]


pdf(file = "DotPlot_glial_cell_development.pdf",   # The directory you want to save the file in
    width = 915, # The width of the plot in inches
    height = 420) # The height of the plot in inches

fontsize = 15
DotPlot(Astros, features = unique(intersect(AstroNeuroMarkersSCT[which(AstroNeuroMarkersSCT$pct.1 > 0),]$gene, genes_glial_cell_development))) + theme(axis.text.x  = element_text(size=fontsize, face = 4, family = "Times", angle = 90, hjust = 1),axis.text.y  = element_text(size=fontsize, face = "bold", family = "Times"), axis.title.y = element_text(size=fontsize, face = "bold", family = "Times"), axis.title.x = element_text(size=fontsize, face = "bold", family = "Times"), plot.title = element_text(hjust = 0.5, size=fontsize, face = "bold", family = "Times"), legend.position = "right", panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

dev.off()


term = "positive regulation of neuron differentiation"
genes_positive_regulation_of_neuron_differentiation = strsplit("Nbl1, Ascl1, Neurog2, Ranbp1, Mif, Ptbp1, Tgif2, Fxn, Elavl4, Rpl4, D16Ertd472e, Nfe2l2, Cdh4, Nme1, Tcf3, Itgb1, Hmg20b, Eef1a1, Ezh2, Mdk, Fezf2, Ncoa1, Zeb1, Mdk, Elavl4, Sox11, Pcp4, Robo2, Neurog2, Eef1a1, Epha4, Dcx, Flna, Rpl4, Arhgef2, Bcl11a, Itgb1, Sez6, Dpysl3, Ptbp1, Ehd1, Enc1, Cbfa2t2, Llph, Nedd4l, Rnd2, Foxg1, Hspa5, Dcx, Pcp4, Sox11, Neurod1, Tubb2b, Dpysl3, Epha4, Rnd2, Sez6, Enc1, Elavl4, Map1b, Tcf4, Nedd4l, Cd24a, Dixdc1, Ehd1, Stmn2, Limk1, Rtn4, Ptprd, Plxna3, Dbn1, Mdk, Arhgef2, Palm, Bcl11a, Parp6, Gsk3b, Eef1a1, Trak1, Rufy3, Kdm1a, Foxg1, Tubb2b, Stmn2, Neurod2, Neurod1, Plxna4, Cd24a, Dcx, Dpysl3, Ptprd, Islr2, Map1b, Sox11, Epha4, Tcf4, Cnr1, Prox1, Nrp1, Gdpd5, Adra2c, Zeb2, Socs2, Enc1, Snap91, Mapt, Dbn1, Ehd1, L1cam, Palm, Rtn4, Pacsin1, Plxna2, Trim32, Dixdc1, Eef1a1, Apbb1, Kif3c, App, Gsk3b, Nedd4l, Trak1, Fyn, Itsn1, Lrp8, Foxg1, Sema5a, Camk2b, Prox1, Plxna4, Neurod2, Cnr1, Islr2, Syt2, Adra2c, Socs2, Cd24a, Plxna2, Neurod1, Sh3gl3, Lrrc7, Stmn2, Mapt, Eef1a1, L1cam, Map1b, Shank1, Tubb2b, Tcf4, Pak3, Zeb2, Dcx, Nin, Kif3c, Trim67, Apbb1, Sox11, Ntrk3, Ptprd, Camk1, Epha4, Enc1, Fgfr1, Rpl4, Fyn, Trim32, Skil, Ehd1, Nrp1, Dbn1, Tcf4, Sema5a, Plxna2, Zeb2, Camk2b, Prox1, Fgfr1, Map1b, Epha4, Ahi1, Ptprd, Dhx36, Rtn4, Caprin1, Islr2, Eef1a1, Tubb2b, Rpl4, Stmn2, Dpysl3, Eif4g2, Fez1, Neurod1, Pcp4, Arf1, Mif, Ranbp1, Map1b, Mapt, Ndnf, Reln, Cxcl12, Syt4, Gdf5, Rit2, Syt1, Plxnd1, Scn1b, Ache, Nrg1, Baiap2, Robo2, App, Map1b, Bend6, Cdh4, Ndrg4, Bcl11a, Impact, Pcp4, Epha3, Stmn2, L1cam, Pacsin1, Islr2, Ahi1, Dab2ip, Nap1l2, Ptprf, Timp2, Dixdc1, Rufy3, Elavl4, Magi2, Actr3, Ncoa1, Negr1, Rtn4", ", " )[[1]]



pdf(file = "DotPlot_positive_regulation_of_neuron_differentiation.pdf",   # The directory you want to save the file in
    width = 1616, # The width of the plot in inches
    height = 420) # The height of the plot in inches

fontsize = 15
DotPlot(Neurons, features = unique(intersect(AstroNeuroMarkersSCT[which(AstroNeuroMarkersSCT$pct.1 > 0.65),]$gene, genes_positive_regulation_of_neuron_differentiation))) + theme(axis.text.x  = element_text(size=fontsize, face = 4, family = "Times", angle = 90, hjust = 1),axis.text.y  = element_text(size=fontsize, face = "bold", family = "Times"), axis.title.y = element_text(size=fontsize, face = "bold", family = "Times"), axis.title.x = element_text(size=fontsize, face = "bold", family = "Times"), plot.title = element_text(hjust = 0.5, size=fontsize, face = "bold", family = "Times"), legend.position = "right", panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

dev.off()



#### Barplots of origins


fontsize = 15

emf(file = "Origin_barplot_astrocytes_and_neurons.emf",width = 7, height = 7
#pdf("Origin_barplot_astrocytes_and_neurons.pdf",)
plotOriginSampleBarplots(AstroNeurontest) ## in SeuratRun.R

dev.off()


pdf("Origin_barplot_high_level_cell_types.pdf",
    width = 512,
    height = 431)
plotOriginSampleBarplots(Combinedtest) ## in SeuratRun.R

dev.off()


### stipecharts
lineageGenes$pch = 1
lineageGenes[which(lineageGenes$color != "gray"),]$pch = 16
lineageGenes$cex = 0.3
lineageGenes[which(lineageGenes$color != "gray"),]$cex = 1

par(mar = c(5, 10, 5, 5))
pdf("Stripechart_DEGs.pdf",
    width = 437,
    height = 418)
#stripeChart(lineageGenes[which(lineageGenes$p_val_adj < 0.05),],  "cluster", "avg_logFC", 0, avg_logFCth, "p_val_adj", c(-0.5,0.5), c(0.5,4.5) )
#stripeChart(lineageGenes[which(lineageGenes$p_val_adj < 0.05 & ( lineageGenes$avg_logFC < -0.2 | lineageGenes$avg_logFC > 0.2)),],  "cluster", "avg_logFC", 0, avg_logFCth, "p_val_adj", c(-0.5,0.5), c(0.5,4.5) )
#stripeChart(lineageGenes,  "cluster", "avg_logFC", 0, avg_logFCth, "p_val_adj", c(-0.5,0.5), c(0.5,4.5) )
stripeChart(lineageGenes,  "cluster", "avg_logFC", 0, avg_logFCth, "p_val_adj", c(-0.5,0.5), c(0.5,4.5), pch = lineageGenes$pch, cex = lineageGenes$cex )

dev.off()



### heatmap
TFDiseaseRelevant = lineageGenes[which( lineageGenes$color != "gray"),]
TFDiseaseRelevantall = TFDiseaseRelevant
TFDiseaseRelevant =  rbind(TFDiseaseRelevantall[which(TFDiseaseRelevantall$cl ==1 & TFDiseaseRelevantall$pct.1 > 0.7 & abs(TFDiseaseRelevantall$avg_logFC)> 0.3 ),] , TFDiseaseRelevantall[which(TFDiseaseRelevantall$cl > 1 ),])
TFDiseaseRelevant =  TFDiseaseRelevantall[which(TFDiseaseRelevantall$cl > 1 ),]

TFDiseaseRelevant = TFDiseaseRelevant[which( !grepl( "Rpl" , TFDiseaseRelevant$Gene) ), ]
TFDiseaseRelevant = TFDiseaseRelevant[which( !grepl( "Rps" , TFDiseaseRelevant$Gene) ), ]
TFDiseaseRelevant = TFDiseaseRelevant[which( !grepl( "mt" , TFDiseaseRelevant$Gene) ), ]
TFDiseaseRelevant = TFDiseaseRelevant[which( !grepl( "Mir" , TFDiseaseRelevant$Gene) ), ]
#TFDiseaseRelevant = TFDiseaseRelevant[which( !grepl( "Gm" , TFDiseaseRelevant$Gene) ), ]
TFDiseaseRelevant = TFDiseaseRelevant[which( !grepl( "AC" , TFDiseaseRelevant$Gene) ), ]
TFDiseaseRelevant = TFDiseaseRelevant[which( !grepl( "CT" , TFDiseaseRelevant$Gene) ), ]

TFDiseaseRelevant <- TFDiseaseRelevant[order( TFDiseaseRelevant$cl, TFDiseaseRelevant$avg_logFC  ),]

rownames(TFDiseaseRelevant) <- paste0 (TFDiseaseRelevant$Gene, "." , TFDiseaseRelevant$cl)
datacluster <- as.data.frame(matrix(0, ncol = length(TFDiseaseRelevant$Gene), nrow = length(unique(TFDiseaseRelevant$cl)) ))
colnames(datacluster) <- rownames(TFDiseaseRelevant)
rownames(datacluster) <- unique(TFDiseaseRelevant$cl)
for( var in colnames(datacluster) ) {
  print(var)
  datacluster[ as.character(TFDiseaseRelevant[var,]$cl), var] <- TFDiseaseRelevant[var,]$avg_logFC
}

pdf("Heatmap_DEGs.pdf",
width = 463,
height = 793)

heatmap.2(t(as.matrix(datacluster)), dendrogram='none', Rowv=FALSE, Colv=FALSE,trace='none',
          density.info = 'none', key = TRUE, keysize = 1, xlab = NULL,col=colorRampPalette(c("orange1", "gray89", "darkseagreen4"))(n = 30),
          margins = c(11,13), family = "Times", cexRow = 1.05, cexCol = 1.05, labRow=as.expression(lapply(TFDiseaseRelevant$Gene, function(a) bquote(italic(.(a))))), key.xlab = "ln fold change", key.title = "",  lhei=c(1,4.5), lwid=c(2,3.5), key.par = list(cex=1),
          labCol =  as.character(unique(TFDiseaseRelevant$lineage_to_plot)) )
dev.off()




#### spatial data

### without ML and devided in 3 intervals
load("/Volumes/Transcend/Pascal/Figures/seurat_spatial_final_withoutML_final.RData")
load("/Volumes/Transcend/Pascal/Figures/AstroNeurontest_newAnnotations.RData")


var = names(seurat_obj_list)[3]


pdf(paste0("spatial_segmentation_", var, ".pdf"),
width = 278,
height = 358)

DimPlot(seurat_obj_list[[var]], reduction = "coord", group.by = "region", cols = c("red3", "red4",  "orange3", "blue4") ) +  ggtitle(var)+ theme(axis.text.y  = element_blank(),
axis.text.x  = element_blank(),
axis.title.y = element_blank(),
axis.title.x = element_blank(),
plot.title = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
panel.background = element_blank(),
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
axis.line = element_blank(),
panel.border = element_blank(),
axis.ticks.x=element_blank(),
axis.ticks.y=element_blank()
) + scale_color_discrete( type = c( "indianred1", "hotpink4", "orange3", "blue4"), labels=c("Granular_Layer_1" = "GL1", "Granular_Layer_2" = "GL2","Hilius" = "Hilus", "outer" = "ML", "SGZ" = "SGZ" ))

dev.off()



### merge clusters

pdf("Merging_cluster_round1.pdf"),
width = 648,
height = 466)

DimPlot(mergeclusters_round1, group.by = "sample") + theme(axis.text.y  = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
axis.text.x  = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
axis.title.y = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
axis.title.x = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
plot.title = element_text(hjust = 0.5, size=fontsize, face = "bold", family = "Times New Roman"),
legend.position = "right",
panel.background = element_blank(),
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
axis.line = element_blank(),
panel.border = element_rect(colour = "black", fill=NA, size=1, linetype="solid")
)

dev.off()


pdf("Merging_cluster_round1_labeled.pdf"),
width = 648,
height = 466)

DimPlot(mergeclusters_round1,label=TRUE) + theme(axis.text.y  = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
axis.text.x  = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
axis.title.y = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
axis.title.x = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
plot.title = element_text(hjust = 0.5, size=fontsize, face = "bold", family = "Times New Roman"),
legend.position = "right",
panel.background = element_blank(),
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
axis.line = element_blank(),
panel.border = element_rect(colour = "black", fill=NA, size=1, linetype="solid")
)

dev.off()



pdf("Merging_clusters_grouped.pdf"),
width = 476,
height = 357)

DimPlot(mergeclusters, group.by = "sample") + theme(axis.text.y  = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
axis.text.x  = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
axis.title.y = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
axis.title.x = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
plot.title = element_text(hjust = 0.5, size=fontsize, face = "bold", family = "Times New Roman"),
legend.position = "right",
panel.background = element_blank(),
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
axis.line = element_blank(),
panel.border = element_rect(colour = "black", fill=NA, size=1, linetype="solid")
)

dev.off()


pdf("Merging_cluster_labeled.pdf"),
width = 476,
height = 357)

DimPlot(mergeclusters,label=TRUE) + theme(axis.text.y  = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
axis.text.x  = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
axis.title.y = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
axis.title.x = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
plot.title = element_text(hjust = 0.5, size=fontsize, face = "bold", family = "Times New Roman"),
legend.position = "right",
panel.background = element_blank(),
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
axis.line = element_blank(),
panel.border = element_rect(colour = "black", fill=NA, size=1, linetype="solid")
)

dev.off()

IdentsToPlot = setdiff(Idents(AstroNeurontest), c("A-stage7", "A-stage6", "A-stage5", "N-stage7"))
CellsToPlot = rownames(mergeclusters@meta.data[mergeclusters@meta.data$previousID %in% c(IdentsToPlot, "spatial"),])
mergeclusters_subset = subset(mergeclusters, cells = CellsToPlot)

IdentsToPlot = setdiff(Idents(mergeclusters), c(11,5,0))
CellsToPlot = WhichCells(mergeclusters_subset, idents =  IdentsToPlot)
mergeclusters_subset = subset(mergeclusters_subset, cells = CellsToPlot)

pdf("Merging_cluster_labeled_previousID.pdf"),
width = 720,
height = 509)

DimPlot(mergeclusters_subset, label=TRUE, group.by = "previousID", label.size = 5, pt.size = 1.5)  + theme(axis.text.y  = element_text(size=fontsize, face = "bold", family = "Times New Roman"), axis.text.x  = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
axis.title.y = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
axis.title.x = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
plot.title = element_text(hjust = 0.5, size=fontsize, face = "bold", family = "Times New Roman"),
legend.position = "right",
panel.background = element_blank(),
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
axis.line = element_blank(),
panel.border = element_rect(colour = "black", fill=NA, size=1, linetype="solid")
)

dev.off()


pdf("Merging_cluster_labeled_split_by_tech.pdf"),
width = 1168,
height = 508)

DimPlot(mergeclusters_subset, group.by = "previousID", split.by = "sample", label=TRUE, label.size = 4.5, pt.size = 1.7)  + theme(axis.text.y  = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
axis.text.x  = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
axis.title.y = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
axis.title.x = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
plot.title = element_text(hjust = 0.5, size=fontsize, face = "bold", family = "Times New Roman"),
legend.position = "right",
panel.background = element_blank(),
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
axis.line = element_blank(),
panel.border = element_rect(colour = "black", fill=NA, size=1, linetype="solid")
)


dev.off()


pdf("Merging_cluster_Neuro_Lineage.pdf"),
width = 510,
height = 413)


FeaturePlot(mergeclusters, c( "NeuroLineage"), pt.size = 0.7, cols =  c("yellow2", "red2")) + theme(axis.text.y  = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
axis.text.x  = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
axis.title.y = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
axis.title.x = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
plot.title = element_text(hjust = 0.5, size=fontsize, face = "bold", family = "Times New Roman"),
legend.position = "right",
panel.background = element_blank(),
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
axis.line = element_blank(),
panel.border = element_rect(colour = "black", fill=NA, size=1, linetype="solid")
)


dev.off()




pdf("Merging_cluster_Astro_Lineage.pdf"),
width = 510,
height = 413)


FeaturePlot(mergeclusters, c( "AstroLineage"), pt.size = 0.7, cols =  c("yellow2", "red2")) + theme(axis.text.y  = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
axis.text.x  = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
axis.title.y = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
axis.title.x = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
plot.title = element_text(hjust = 0.5, size=fontsize, face = "bold", family = "Times New Roman"),
legend.position = "right",
panel.background = element_blank(),
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
axis.line = element_blank(),
panel.border = element_rect(colour = "black", fill=NA, size=1, linetype="solid")
)


dev.off()


### deviding in 4
mergeclusters@meta.data$Neuro_lineage_discrete = cut(mergeclusters@meta.data$NeuroLineage, breaks = c(0, max(NSCRG@meta.data$NeuroLineage) , max(NSCRG@meta.data$NeuroLineage)+ (NeuroLineageMax - max(NSCRG@meta.data$NeuroLineage))/3, max(NSCRG@meta.data$NeuroLineage)+ 2*(NeuroLineageMax - max(NSCRG@meta.data$NeuroLineage))/3, NeuroLineageMax ) )
mergeclusters@meta.data$Astro_lineage_discrete = cut(mergeclusters@meta.data$AstroLineage, breaks = c(0, max(NSCRG@meta.data$AstroLineage) , max(NSCRG@meta.data$AstroLineage)+ (AstroLineageMax - max(NSCRG@meta.data$AstroLineage))/3 , max(NSCRG@meta.data$AstroLineage)+2*(AstroLineageMax - max(NSCRG@meta.data$AstroLineage))/3, AstroLineageMax) )


pdf("Merging_cluster_Neuro_Lineage_Discrete_4parts.pdf"),
width = 582,
height = 413)

DimPlot(mergeclusters, group.by = "Neuro_lineage_discrete", pt.size = 0.7) + theme(axis.text.y  = element_text(size=fontsize, face = "bold", family = "Times New Roman"), axis.text.x  = element_text(size=fontsize, face = "bold", family = "Times New Roman"), axis.title.y = element_text(size=fontsize, face = "bold", family = "Times New Roman"), axis.title.x = element_text(size=fontsize, face = "bold", family = "Times New Roman"), plot.title = element_text(hjust = 0.5, size=fontsize, face = "bold", family = "Times New Roman"), legend.position = "right",panel.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1, linetype="solid") )+ggtitle("Neuro_lineage_discrete")

dev.off()



pdf("Merging_cluster_Astro_Lineage_Discrete_4parts.pdf"),
width = 582,
height = 413)

DimPlot(mergeclusters, group.by = "Astro_lineage_discrete", pt.size = 0.7) + theme(axis.text.y  = element_text(size=fontsize, face = "bold", family = "Times New Roman"), axis.text.x  = element_text(size=fontsize, face = "bold", family = "Times New Roman"), axis.title.y = element_text(size=fontsize, face = "bold", family = "Times New Roman"), axis.title.x = element_text(size=fontsize, face = "bold", family = "Times New Roman"), plot.title = element_text(hjust = 0.5, size=fontsize, face = "bold", family = "Times New Roman"), legend.position = "right",panel.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1, linetype="solid") )+ggtitle("Astro_lineage_discrete_4")

dev.off()

i=1
pdf(paste0("Neuro_Lineage_part",i,"_dist.pdf")),
width = 775,
height = 296)

ggplot(datatableNeuroLineage[which(datatableNeuroLineage$Neuro_lineage_discrete ==  levels(unique(datatableNeuroLineage$Neuro_lineage_discrete))[i]),], aes(x=region, y=Freq, fill = Neuro_lineage_discrete, group = condition, color = condition)) + stat_compare_means(aes(group = condition), label = "p.format", label.y = 100) + scale_x_discrete(labels=c("Granular_Layer_1" = "GL1", "Granular_Layer_2" = "GL2","Hilius" = "Hilus", "outer" = "ML", "SGZ" = "SGZ" )) + theme(axis.text.y  = element_text(size=15, face = "bold", family = "Times New Roman"), axis.text.x  = element_text(size=15, face = "bold", family = "Times New Roman"), axis.title.y = element_text(size=15, face = "bold", family = "Times New Roman"),
axis.title.x = element_text(size=15, face = "bold", family = "Times New Roman"),
plot.title = element_text(hjust = 0.5, size=15, face = "bold", family = "Times New Roman"),
legend.position = "right",
panel.background = element_blank(),
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
axis.line = element_blank(),
legend.title = element_text(size=20, color = "black"),
legend.text = element_text(size=20),
legend.key=element_rect(fill='white', size = 20),
panel.border = element_rect(colour = "black", fill=NA, size=1, linetype="solid")) + guides(colour = guide_legend(override.aes = list(size=0, stroke=4)))  + ylim(-20,120) + xlab("") + ylab("% of cells")+
    geom_point(position = position_jitterdodge(jitter.width = .2, dodge.width = .7), alpha = .5) +
    stat_summary(fun = mean, na.rm = TRUE, geom = "point", shape = "diamond",
                 size = 4, color = "black", position = position_dodge(width = .7)) +
    stat_summary(fun.data = mean_cl_normal, na.rm = TRUE, geom = "errorbar", width = .2, color = "black",
                 position = position_dodge(width = .7)) +
    scale_color_brewer(palette = "Set1")

dev.off()



i=1
pdf(paste0("Astro_Lineage_part",i,"_dist.pdf")),
width = 775,
height = 296)

ggplot(datatableAstroLineage[which(datatableAstroLineage$Astro_lineage_discrete ==  levels(unique(datatableAstroLineage$Astro_lineage_discrete))[1]),], aes(x=region, y=Freq, fill = Astro_lineage_discrete, group = condition, color = condition)) + stat_compare_means(aes(group = condition), label = "p.format", label.y = 110) + scale_x_discrete(labels=c("Granular_Layer_1" = "GL1", "Granular_Layer_2" = "GL2","Hilius" = "Hilus", "outer" = "ML", "SGZ" = "SGZ" )) + theme(axis.text.y  = element_text(size=15, face = "bold", family = "Times New Roman"), axis.text.x  = element_text(size=15, face = "bold", family = "Times New Roman"), axis.title.y = element_text(size=15, face = "bold", family = "Times New Roman"),
axis.title.x = element_text(size=15, face = "bold", family = "Times New Roman"),
plot.title = element_text(hjust = 0.5, size=15, face = "bold", family = "Times New Roman"),
legend.position = "right",
panel.background = element_blank(),
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
axis.line = element_blank(),
legend.title = element_text(size=20, color = "black"),
legend.text = element_text(size=20),
legend.key=element_rect(fill='white', size = 20),
panel.border = element_rect(colour = "black", fill=NA, size=1, linetype="solid")) + guides(colour = guide_legend(override.aes = list(size=0, stroke=4)))  + ylim(-20,120) + xlab("") + ylab("% of cells")+
    geom_point(position = position_jitterdodge(jitter.width = .2, dodge.width = .7), alpha = .5) +
    stat_summary(fun = mean, na.rm = TRUE, geom = "point", shape = "diamond",
                 size = 4, color = "black", position = position_dodge(width = .7)) +
    stat_summary(fun.data = mean_cl_normal, na.rm = TRUE, geom = "errorbar", width = .2, color = "black",
                 position = position_dodge(width = .7)) +
    scale_color_brewer(palette = "Set1")


dev.off()





#### deviding in 3
NeuroLineageMax = max(na.omit(mergeclusters@meta.data$NeuroLineage))
AstroLineageMax = max(na.omit(mergeclusters@meta.data$AstroLineage))
NSCRG = SubsetData(mergeclusters, ident.use = c(8))

mergeclusters@meta.data$Neuro_lineage_discrete = cut(mergeclusters@meta.data$NeuroLineage, breaks = c(0, max(NSCRG@meta.data$NeuroLineage),max(NSCRG@meta.data$NeuroLineage)+ (NeuroLineageMax - max(NSCRG@meta.data$NeuroLineage))/2, NeuroLineageMax ) )
mergeclusters@meta.data$Astro_lineage_discrete = cut(mergeclusters@meta.data$AstroLineage, breaks = c(0, max(NSCRG@meta.data$AstroLineage) , max(NSCRG@meta.data$AstroLineage)+ (AstroLineageMax - max(NSCRG@meta.data$AstroLineage))/2 , AstroLineageMax) )


# With ML and devided in 3 intervals

load("/Volumes/Transcend/Pascal/Figures/seurat_spatial_final_withML.RData")
var = names(seurat_obj_list)[3]


pdf(paste0("spatial_segmentation_", var, ".pdf"),
width = 278,
height = 358)

DimPlot(seurat_obj_list[[var]], reduction = "coord", group.by = "region", cols = c("red3", "red4",  "orange3", "blue4") ) +  ggtitle(var) + theme(axis.text.y  = element_blank(),axis.text.x  = element_blank(),axis.title.y = element_blank(),axis.title.x = element_blank(),plot.title = element_text(size=fontsize, face = "bold", family = "Times New Roman"),panel.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line = element_blank(),panel.border = element_blank(), axis.ticks.x=element_blank(),axis.ticks.y=element_blank()
) + scale_color_discrete( type = c( "indianred1", "hotpink4", "orange3", "darkcyan" , "blue4"), labels=c("Granular_Layer_1" = "GL1", "Granular_Layer_2" = "GL2","Hilius" = "Hilus", "outer" = "ML", "SGZ" = "SGZ" ))

dev.off()


fwrite(NumberOfCellsInHighLevel, file = "Table1_NumberOfCellsInHighLevelClustering.csv")
fwrite(NumberOfCellsInAstrocyticAndNeuronalPopulations, file = "Table2_NumberOfCellsInAstrocyticAndNeuronalPopulations.csv")
fwrite(AstroNeuroMarkersSCT, file = "Table3_AstroNeuroMarkersSCTnormalization.csv")
fwrite(TFDiseaseRelevantall, file = "Table6_disease_related_DGE_list.csv")

datatableNeuroLineage = data.frame(region=character(),
                                   Neuro_lineage_discrete=character(),
                                   Freq=double(),
                                   slide=character(),
                                   stringsAsFactors=FALSE)

for(var in names(datalist)){
    obj = seurat_obj_list[[var]]
    subset = obj@meta.data[c("Astro_lineage_discrete", "Neuro_lineage_discrete", "region")]
    tt = table(subset[, c(2,3)])
    dataframeobj = as.data.frame(tt)
    dataframeobj$slide = var
    datatableNeuroLineage = rbind(datatableNeuroLineage, dataframeobj)
}
datatableNeuroLineage$condition = unlist(strsplit(datatableNeuroLineage$slide, "_"))[seq(2,length(datatableNeuroLineage$slide)*3,3)]

fwrite(datatableNeuroLineage, file = "Table8_The_absolute_number_of_neuronal_progenitors_per_pseudotime_location_tissue.csv")


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
    tt = table(subset[, c(1,3)])
    dataframeobj = as.data.frame(tt)
    dataframeobj$slide = var
    datatableAstroLineage = rbind(datatableAstroLineage, dataframeobj)
}
datatableAstroLineage$condition = unlist(strsplit(datatableAstroLineage$slide, "_"))[seq(2,length(datatableAstroLineage$slide)*3,3)]


fwrite(datatableAstroLineage, file = "Table9_The_absolute_of_astrocytic_progenitors_per_pseudotime_location_tissue.csv")



### additional Figures for the supplementary information

Shamcells = rownames(Combinedtest@meta.data[which(Combinedtest@meta.data$sample == "Sham"),])
TBIcells = rownames(Combinedtest@meta.data[which(Combinedtest@meta.data$sample == "TBI"),])

TBICombined_test = SubsetData(Combinedtest, cells = TBIcells)
ShamCombined_test = SubsetData(Combinedtest, cells = Shamcells)

histData = data.frame("TBI"=tbihist )

ggplot(histData, aes(x = TBI)) +
    geom_density(aes(y = ..density..))  + theme(axis.text.y  = element_text(size=15, face = "bold", family = "Times New Roman"), axis.text.x  = element_text(size=15, face = "bold", family = "Times New Roman"), axis.title.y = element_text(size=15, face = "bold", family = "Times New Roman"),
                            axis.title.x = element_text(size=15, face = "bold", family = "Times New Roman"),
                            plot.title = element_text(hjust = 0.5, size=15, face = "bold", family = "Times New Roman"),
                            legend.position = "right",
                            panel.background = element_blank(),
                            panel.grid.major = element_blank(),
                            panel.grid.minor = element_blank(),
                            axis.line = element_blank(),
                            legend.title = element_text(size=20, color = "black"),
                            legend.text = element_text(size=20),
                            legend.key=element_rect(fill='white', size = 20),
                            panel.border = element_rect(colour = "black", fill=NA, size=1, linetype="solid")) + guides(colour = guide_legend(override.aes = list(size=0, stroke=4)))  + xlab("Transcript count per cells in TBI animals") + ylab("Density")+geom_vline(xintercept = 200, col = "red" )



### correlation between samples

pdf(paste0("correlation_between_spatial_TBI_samples", var, ".pdf"), width = 1397, height = 605)

par(mfrow=c(2,5))
listToconsider = names(seurat_obj_list)[grepl("TBI", names(seurat_obj_list))]
l = length(listToconsider)
for( i in 1:(l-1)) {
    
    var1 = listToconsider[i]
    for(var2  in listToconsider [(i+1):l]) {
        
        sample1 =  as.numeric(rowMeans( seurat_obj_list[[var1]]@assays$RNA@counts[ggenelist,] ))
        sample2 =  as.numeric(rowMeans( seurat_obj_list[[var2]]@assays$RNA@counts[ggenelist,] ))
        
        plot(sample1 ~ sample2, main = paste0( "Pearson corr: ", round(cor(sample1, sample2), 3) ), xlab = var2, ylab= var1, pch = 16) + abline(lm(sample1 ~  sample2), col="red" )
        
    }
}

dev.off()


pdf(paste0("correlation_between_spatial_control_samples", var, ".pdf"), width = 838, height = 358)

par(mfrow=c(2,3))
ggenelist = rownames(objSpatial)
listToconsider = names(seurat_obj_list)[grepl("ontrol", names(seurat_obj_list))]
l = length(listToconsider)
for( i in 1:(l-1)) {
    
    var1 = listToconsider[i]
    for(var2  in listToconsider [(i+1):l]) {
        
        sample1 =  as.numeric(rowMeans( seurat_obj_list[[var1]]@assays$RNA@counts[ggenelist,] ))
        sample2 =  as.numeric(rowMeans( seurat_obj_list[[var2]]@assays$RNA@counts[ggenelist,] ))
        
        plot(sample1 ~ sample2, main = paste0( "Pearson corr: ", round(cor(sample1, sample2), 3) ), xlab = var2, ylab= var1, pch = 16) + abline(lm(sample1 ~  sample2), col="red" )
        
    }
}

dev.off()


par(mfrow=c(2,3))
ggenelist = rownames(objSpatial)
listToconsider = names(seurat_obj_list)[grepl("ontrol", names(seurat_obj_list))]
l = length(listToconsider)
for( i in 1:l) {
    
        var1 = listToconsider[i]
        sample1 =  as.numeric(rowMeans( seurat_obj_list[[var1]]@assays$RNA@counts[ggenelist,] ))
        sample2 =  as.numeric(rowMeans( AstroNeurontest@assays$RNA@counts[ggenelist,] ))
        
        plot(sample1 ~ sample2, main = paste0( "Pearson corr: ", round(cor(sample1, sample2), 3) ), xlab = var2, ylab= "SCRNA-seq", pch = 16) + abline(lm(sample1 ~  sample2), col="red" )
        

}




#### Figure of markers



pdf("FeaturePlot_HighLevelMarkers.pdf"),
width = 1310,
height = 767)

FeaturePlot(Combinedtest, c("Mki67","Ascl1", "Snap25", "Neurod1","Aldoc", "Slc1a3", "Pdgfra", "Mog", "Pecam1", "Cd68", "Cx3cr1", "Des")) & theme( plot.title = element_text( face = "italic") )

dev.off()



pdf("FeaturePlot_RNA_scopeMarkers.pdf"),
width = 983,
height = 767)

FeaturePlot(AstroNeurotest, c("Hapln1","Frzb", "Ascl1", "Sparc","Sned1", "Neat1", "Slc1a3"))& theme( plot.title = element_text( face = "italic") )

dev.off()



### checking GFAP expression


var = names(seurat_obj_list)[1]


pdf(paste0("spatial_Gfap_expression_", var, ".pdf"),
width = 278,
height = 358)

FeaturePlot(seurat_obj_list[[var]], "Gfap", reduction = "coord") +  ggtitle(paste0(var, "_Gfap"))+ theme(axis.text.y  = element_blank(),
axis.text.x  = element_blank(),
axis.title.y = element_blank(),
axis.title.x = element_blank(),
plot.title = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
panel.background = element_blank(),
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
axis.line = element_blank(),
panel.border = element_blank(),
axis.ticks.x=element_blank(),
axis.ticks.y=element_blank()
)

dev.off()

### Lineage plots

var = names(seurat_obj_list)[1]

pdf(paste0("spatial_AstroLineage_", var, ".pdf"),
width = 278,
height = 358)

seurat_obj_list[[var]]@meta.data$pt.size.astro = ifelse(is.na(seurat_obj_list[[var]]@meta.data$AstroLineage), 0.5, 1)

FeaturePlot(seurat_obj_list[[var]], "AstroLineage", reduction = "coord", cols = c("yellow2", "red2"), pt.size = seurat_obj_list[[var]]@meta.data$pt.size.astro )  +  ggtitle(paste0(var, "_Astro"))+ theme(axis.text.y  = element_blank(),
axis.text.x  = element_blank(),
axis.title.y = element_blank(),
axis.title.x = element_blank(),
plot.title = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
panel.background = element_blank(),
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
axis.line = element_blank(),
panel.border = element_blank(),
axis.ticks.x=element_blank(),
axis.ticks.y=element_blank()
)

dev.off()


var = names(seurat_obj_list)[1]

pdf(paste0("spatial_NeuroLineage_", var, ".pdf"),
width = 278,
height = 358)

seurat_obj_list[[var]]@meta.data$pt.size.neuro = ifelse(is.na(seurat_obj_list[[var]]@meta.data$NeuroLineage), 0.5, 1)

FeaturePlot(seurat_obj_list[[var]], "NeuroLineage", reduction = "coord", cols = c("yellow2", "red2"), pt.size = seurat_obj_list[[var]]@meta.data$pt.size.neuro )  +  ggtitle(paste0(var, "_Neuro"))+ theme(axis.text.y  = element_blank(),
axis.text.x  = element_blank(),
axis.title.y = element_blank(),
axis.title.x = element_blank(),
plot.title = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
panel.background = element_blank(),
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
axis.line = element_blank(),
panel.border = element_blank(),
axis.ticks.x=element_blank(),
axis.ticks.y=element_blank()
)

dev.off()


### cell displacement in Neuronal Lineage

var = names(seurat_obj_list)[1]
lineageRange = 3

pdf(paste0("spatial_Neuro_Lineage_", lineageRange, "_", var, ".pdf"),
width = 278,
height = 358)

cellsToHighlights = rownames(seurat_obj_list[[var]]@meta.data[which(seurat_obj_list[[var]]@meta.data$Neuro_lineage_discrete == levels(seurat_obj_list[[var]]@meta.data$Neuro_lineage_discrete)[lineageRange] ),])

cols.final = as.character(getColors(DimPlot(seurat_obj_list[[var]], group.by = "Neuro_lineage_discrete")))[lineageRange]

DimPlot(seurat_obj_list[[var]], group.by = "Neuro_lineage_discrete", reduction = "coord", cells.highlight = cellsToHighlights,  sizes.highlight = 2) + scale_color_manual(labels  = c("Unselected", levels(seurat_obj_list[[var]]@meta.data$Neuro_lineage_discrete)[lineageRange]), values = c("grey81", cols.final) ) +  ggtitle(paste0(var, "_Neuro")) + theme(axis.text.y  = element_blank(),
axis.text.x  = element_blank(),
axis.title.y = element_blank(),
axis.title.x = element_blank(),
plot.title = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
panel.background = element_blank(),
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
axis.line = element_blank(),
panel.border = element_blank(),
axis.ticks.x=element_blank(),
axis.ticks.y=element_blank()
)

dev.off()


### cell displacement in Astrocytic Lineage

var = names(seurat_obj_list)[1]
lineageRange = 1

pdf(paste0("spatial_Astro_Lineage_", lineageRange, "_", var, ".pdf"),
width = 278,
height = 358)

cellsToHighlights = rownames(seurat_obj_list[[var]]@meta.data[which(seurat_obj_list[[var]]@meta.data$Astro_lineage_discrete == levels(seurat_obj_list[[var]]@meta.data$Astro_lineage_discrete)[lineageRange] ),])

cols.final = as.character(getColors(DimPlot(seurat_obj_list[[var]], group.by = "Astro_lineage_discrete")))[lineageRange]

DimPlot(seurat_obj_list[[var]], group.by = "Astro_lineage_discrete", reduction = "coord", cells.highlight = cellsToHighlights,  sizes.highlight = 2) + scale_color_manual(labels  = c("Unselected", levels(seurat_obj_list[[var]]@meta.data$Astro_lineage_discrete)[lineageRange]), values = c("grey81", cols.final) ) +  ggtitle(paste0(var, "_Astro")) + theme(axis.text.y  = element_blank(),
axis.text.x  = element_blank(),
axis.title.y = element_blank(),
axis.title.x = element_blank(),
plot.title = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
panel.background = element_blank(),
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
axis.line = element_blank(),
panel.border = element_blank(),
axis.ticks.x=element_blank(),
axis.ticks.y=element_blank()
)

dev.off()



### heatmap of the mapping data

pdf("Heatmap_spatial_mapping.pdf",
width = 917,
height = 742)

barplotOfPercentages(metadata = metadataCombined)

dev.off()

### save the table

metadata_tmp = as.data.frame.matrix (table(mergeclusters@meta.data[, c( "seurat_clusters","previousID")]) )
write.csv(metadata_tmp, file = "/Volumes/Transcend/Pascal/Figures/Table12_spatial_to_10x_mapping.csv")


#### spatial mapping
mergeclusters@meta.data$ID10x = "Unmapped"
mergeclusters@meta.data[WhichCells(mergeclusters, idents = c(8)), ]$ID10x = "8, NSC-stages1,2"
mergeclusters@meta.data[WhichCells(mergeclusters, idents = c(7)), ]$ID10x = "7, RG-like"
mergeclusters@meta.data[WhichCells(mergeclusters, idents = c(1)), ]$ID10x = "1, N-stages1,2"
mergeclusters@meta.data[WhichCells(mergeclusters, idents = c(2)), ]$ID10x = "2, N-stage3"
mergeclusters@meta.data[WhichCells(mergeclusters, idents = c(6)), ]$ID10x = "6, N-stages4,5"
mergeclusters@meta.data[WhichCells(mergeclusters, idents = c(10)), ]$ID10x = "10, N-stages5,6"
mergeclusters@meta.data[WhichCells(mergeclusters, idents = c(11)), ]$ID10x = "11, N-stage7"
mergeclusters@meta.data[WhichCells(mergeclusters, idents = c(9)), ]$ID10x= "9, A-stages1,5"
mergeclusters@meta.data[WhichCells(mergeclusters, idents = c(3)), ]$ID10x = "3, A-stage2"
mergeclusters@meta.data[WhichCells(mergeclusters, idents = c(4)), ]$ID10x = "4, A-stages3,4"
mergeclusters@meta.data[WhichCells(mergeclusters, idents = c(0)), ]$ID10x = "0, A-stages6,7,5"
mergeclusters@meta.data[WhichCells(mergeclusters, idents = c(5)), ]$ID10x = "5, Unmapped"

cols.mergedclusters = getColors(DimPlot(mergeclusters))


for(var in names(seurat_obj_list)){
    obj = seurat_obj_list[[var]]
    seurat_obj_list[[var]]@meta.data$ID10x = "Unmapped"
    seurat_obj_list[[var]]@meta.data[intersect(colnames(seurat_obj_list[[var]]), colnames(mergeclusters)),]$ID10x = mergeclusters@meta.data[intersect(colnames(seurat_obj_list[[var]]), colnames(mergeclusters)),]$ID10x
    seurat_obj_list[[var]]@meta.data$ID10x = factor(seurat_obj_list[[var]]@meta.data$ID10x, levels = c("8, NSC-stages1,2", "7, RG-like",  "1, N-stages1,2" , "2, N-stage3", "6, N-stages4,5",  "10, N-stages5,6", "11, N-stage7", "9, A-stages1,5", "3, A-stage2", "4, A-stages3,4", "0, A-stages6,7,5", "5, Unmapped", "Unmapped"))

}


var = names(seurat_obj_list)[1]

pdf(paste0("spatial_location_of_NSC_like_on_", var, ".pdf"),
width = 492,
height = 519)

DimPlot(seurat_obj_list[[var]], group.by = "ID10x", reduction = "coord", cols = c(cols.mergedclusters [[9]], cols.mergedclusters [[8]], "gray","gray","gray","gray","gray","gray","gray", "gray","gray","gray"), pt.size = 1.5) + ggtitle(var)+ theme(axis.text.y  = element_blank(),
axis.text.x  = element_blank(),
axis.title.y = element_blank(),
axis.title.x = element_blank(),
plot.title = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
panel.background = element_blank(),
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
axis.line = element_blank(),
panel.border = element_blank(),
axis.ticks.x=element_blank(),
axis.ticks.y=element_blank()
)

dev.off()




pdf(paste0("spatial_location_of_Neuronal_lineage_on_", var, ".pdf"),
width = 492,
height = 519)

DimPlot(seurat_obj_list[[var]], group.by = "ID10x", reduction = "coord", cols = c( "gray","gray", cols.mergedclusters [[2]], cols.mergedclusters [[3]],cols.mergedclusters [[7]],cols.mergedclusters [[11]], "gray","gray","gray", "gray","gray","gray"), pt.size = 1.5) + ggtitle(var)+ theme(axis.text.y  = element_blank(),
axis.text.x  = element_blank(),
axis.title.y = element_blank(),
axis.title.x = element_blank(),
plot.title = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
panel.background = element_blank(),
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
axis.line = element_blank(),
panel.border = element_blank(),
axis.ticks.x=element_blank(),
axis.ticks.y=element_blank()
)

dev.off()


pdf(paste0("spatial_location_of_Astrocytic_lineage_on_", var, ".pdf"),
width = 492,
height = 519)

DimPlot(seurat_obj_list[[var]], group.by = "ID10x", reduction = "coord", cols = c( "gray","gray",  "gray","gray","gray", "gray",cols.mergedclusters [[10]], cols.mergedclusters [[4]],cols.mergedclusters [[5]],cols.mergedclusters [[1]],"gray","gray"), pt.size = 1.5) + ggtitle(var)+ theme(axis.text.y  = element_blank(),
axis.text.x  = element_blank(),
axis.title.y = element_blank(),
axis.title.x = element_blank(),
plot.title = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
panel.background = element_blank(),
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
axis.line = element_blank(),
panel.border = element_blank(),
axis.ticks.x=element_blank(),
axis.ticks.y=element_blank()
)

dev.off()


pdf("Marker_genes_in_spatial_data.pdf",
width = 1422,
height = 678)

spatialcells = rownames(meta_data_combined[which(meta_data_combined$sample == "Genexyz"),])
onlyspatialmergeddata = SubsetData(mergeclusters, cells = spatialcells)



FeaturePlot(onlyspatialmergeddata, features = c("Ascl1", "Aurkb", "Ube2c", "Mki67", "Pax6",
                                                "Sox4", "Neurog2", "Unc5d", "Calb2", "Syt1",
                                                "Aldoc", "Sned1", "Fgfr3", "Ogt", "Fam107a"
), ncol = 5, ) & theme( plot.title = element_text( face = "italic") )
dev.off()


#### Elbow plots

pdf("AstroNeurons_ElbowPlot.pdf",
width = 488,
height = 353)

ElbowPlot(AstroNeurontest) + theme(axis.text.y  = element_text(size=fontsize, face = "bold", family = "Times New Roman"), axis.text.x  = element_text(size=fontsize, face = "bold", family = "Times New Roman"), axis.title.y = element_text(size=fontsize, face = "bold", family = "Times New Roman"), axis.title.x = element_text(size=fontsize, face = "bold", family = "Times New Roman"), plot.title = element_text(hjust = 0.5, size=fontsize, face = "bold", family = "Times New Roman"), legend.position = "right",panel.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1, linetype="solid") )+ggtitle("")
dev.off()

pdf("HigherOrderCell_ElbowPlot.pdf",
width = 488,
height = 353)

ElbowPlot(Combinedtest) + theme(axis.text.y  = element_text(size=fontsize, face = "bold", family = "Times New Roman"), axis.text.x  = element_text(size=fontsize, face = "bold", family = "Times New Roman"), axis.title.y = element_text(size=fontsize, face = "bold", family = "Times New Roman"), axis.title.x = element_text(size=fontsize, face = "bold", family = "Times New Roman"), plot.title = element_text(hjust = 0.5, size=fontsize, face = "bold", family = "Times New Roman"), legend.position = "right",panel.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1, linetype="solid") )+ggtitle("")
dev.off()

#### Feature number filters

S5data  = Read10X ("~/Dropbox/TBI_DATA/10xdata/S5_10X_05042019/filtered_feature_bc_matrix/")
S6data  = Read10X ("~/Dropbox/TBI_DATA/10xdata/S6_10X_05042019/filtered_feature_bc_matrix/")


pdf("ControlData_S5_TranscriptCount.pdf",
width = 488,
height = 353)

data10x = list("S5" = colSums(S5data!=0))
data10x = as.data.frame(data10x)
ggplot(data10x, aes(x=S5))+geom_density( color="blue", size = 2)+ stat_function(fun = dnorm, colour = "blue", args = list(mean = 0.3, sd = 1))+ scale_x_continuous(breaks = c(300,2500,5000,6000)) + xlab ("Transcript count per nuclei" ) + ylab ("Density" ) + geom_vline(xintercept = 300, col = "black",linetype = "dashed")+geom_vline(xintercept = 5000, col = "black", linetype = "dashed")+ theme(axis.text.y  = element_text(size=fontsize, face = "bold", family = "Times New Roman"), axis.text.x  = element_text(size=fontsize, face = "bold", family = "Times New Roman"), axis.title.y = element_text(size=fontsize, face = "bold", family = "Times New Roman"), axis.title.x = element_text(size=fontsize, face = "bold", family = "Times New Roman"), plot.title = element_text(hjust = 0.5, size=fontsize, face = "bold", family = "Times New Roman"), legend.position = "right",panel.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1, linetype="solid") ) + ggtitle("Transcript count distribution in control")

dev.off()



pdf("ControlData_S6_TranscriptCount.pdf",
width = 488,
height = 353)

data10x = list("S6" = colSums(S6data!=0))
data10x = as.data.frame(data10x)
ggplot(data10x, aes(x=S6))+geom_density( color="blue", size = 2)+ stat_function(fun = dnorm, colour = "blue", args = list(mean = 0.3, sd = 1))+ scale_x_continuous(breaks = c(200,2500,5000,6000)) + xlab ("Transcript count per nuclei" ) + ylab ("Density" ) + geom_vline(xintercept = 200, col = "black",linetype = "dashed")+geom_vline(xintercept = 5000, col = "black", linetype = "dashed")+ theme(axis.text.y  = element_text(size=fontsize, face = "bold", family = "Times New Roman"), axis.text.x  = element_text(size=fontsize, face = "bold", family = "Times New Roman"), axis.title.y = element_text(size=fontsize, face = "bold", family = "Times New Roman"), axis.title.x = element_text(size=fontsize, face = "bold", family = "Times New Roman"), plot.title = element_text(hjust = 0.5, size=fontsize, face = "bold", family = "Times New Roman"), legend.position = "right",panel.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1, linetype="solid") ) + ggtitle("Transcript count distribution in TBI sample")

dev.off()

dataToPlot = FetchData(AstroNeurontest, c("Ppp1r14b", "idents", "sample"), slot = "data")

pdf("Ppp1r14bEXPinSub.pdf",
width = 401,
height = 349)

var = names(adjpval)[2]
ylim = max([which([which(dataToPlot$idents == var),]$idents == var),]$Ppp1r14b)
ggplot(dataToPlot[which(dataToPlot$idents == var),], aes(x = sample, y = Ppp1r14b)) +
geom_violin(aes(color=sample, fill=sample), trim = TRUE) +xlab(var)+ylab("")+ggtitle("Ppp1r14b")+
geom_jitter(height = 0, width = 0.1)  + geom_segment(aes(x=1,xend=2,y=ylim,yend=ylim)) + geom_text(aes(x=1.5,y=ylim+0.2), label = ifelse(adjpval[var]<0.001, "***", ifelse(adjpval[var]<0.01, "**", ifelse(adjpval[var]<0.05, "*", "ns") ))) + theme(axis.text.y  = element_text(size=fontsize, face = "bold", family = "Times New Roman"), axis.text.x  = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
axis.title.y = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
axis.title.x = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
plot.title = element_text(hjust = 0.5, size=fontsize, face = "bold.italic", family = "Times New Roman"),
legend.position = "right",
panel.border = element_rect(colour = "black", fill=NA, size=1, linetype="solid"),
legend.title = element_text( size=fontsize, face = "bold", family = "Times New Roman"),
legend.text = element_text( size=fontsize, face = "bold", family = "Times New Roman"),
panel.grid.major = element_blank(),panel.grid.minor = element_blank(), panel.background = element_blank()
)

dev.off()




### FIRST VERSION


## Run with Seurat v3.2.0

#setwd("/Users/u0114327/Dropbox/")
setwd("/Volumes/Transcend/Pascal/Figures/")
load("RevisionOfAnnotations.RData")

Combinedtest = SetIdent(Combinedtest, cells = WhichCells(Combinedtest, idents = c("N1", "N2", "N3", "N4")), value = "Neuronal lineage")
Combinedtest = SetIdent(Combinedtest, cells = WhichCells(Combinedtest, idents =  c("A1", "A2", "A3", "A4")), value = "Astrocyte lineage")

AstroNeurontest = SetIdent(AstroNeurontest, cells = WhichCells(AstroNeurontest, idents = c("A1_1")), value = "A-stage1")
AstroNeurontest = SetIdent(AstroNeurontest, cells = WhichCells(AstroNeurontest, idents = c("A2")), value = "A-stage2")
AstroNeurontest = SetIdent(AstroNeurontest, cells = WhichCells(AstroNeurontest, idents = c("A3_1")), value = "A-stage3")
AstroNeurontest = SetIdent(AstroNeurontest, cells = WhichCells(AstroNeurontest, idents = c("A3_2")), value = "A-stage4")
AstroNeurontest = SetIdent(AstroNeurontest, cells = WhichCells(AstroNeurontest, idents = c("A1_2")), value = "A-stage5")
AstroNeurontest = SetIdent(AstroNeurontest, cells = WhichCells(AstroNeurontest, idents = c("A4_2")), value = "A-stage6")
AstroNeurontest = SetIdent(AstroNeurontest, cells = WhichCells(AstroNeurontest, idents = c("A4_1")), value = "A-stage7")
AstroNeurontest = SetIdent(AstroNeurontest, cells = WhichCells(AstroNeurontest, idents = c("N1_1")), value = "N-stage1")
AstroNeurontest = SetIdent(AstroNeurontest, cells = WhichCells(AstroNeurontest, idents = c("N1_2")), value = "N-stage2")
AstroNeurontest = SetIdent(AstroNeurontest, cells = WhichCells(AstroNeurontest, idents = c("N2_1")), value = "N-stage3")
AstroNeurontest = SetIdent(AstroNeurontest, cells = WhichCells(AstroNeurontest, idents = c("N2_2")), value = "N-stage4")
AstroNeurontest = SetIdent(AstroNeurontest, cells = WhichCells(AstroNeurontest, idents = c("N3_1")), value = "N-stage5")
AstroNeurontest = SetIdent(AstroNeurontest, cells = WhichCells(AstroNeurontest, idents = c("N3_2")), value = "N-stage6")
AstroNeurontest = SetIdent(AstroNeurontest, cells = WhichCells(AstroNeurontest, idents = c("N4")), value = "N-stage7")
AstroNeurontest = SetIdent(AstroNeurontest, cells = WhichCells(AstroNeurontest, idents = c("NSC_1")), value = "NSC-stage1")
AstroNeurontest = SetIdent(AstroNeurontest, cells = WhichCells(AstroNeurontest, idents = c("NSC_2")), value = "NSC-stage2")


#genes.viz = c(genes.viz[1:19], genes.viz[22:25],genes.viz[c(20,21,36)], genes.viz[26:35],genes.viz[37:63] )

AstroNeurontest@active.ident = factor(x=AstroNeurontest@active.ident, levels = c("NSC-stage1", "NSC-stage2", "RG-like", "N-stage1","N-stage2", "N-stage3", "N-stage4", "N-stage5", "N-stage6", "N-stage7", "A-stage1","A-stage2","A-stage3","A-stage4","A-stage5","A-stage6","A-stage7") )
Neurons = SubsetData( AstroNeurontest, ident.use = levels(AstroNeurontest@active.ident)[1:10])
Astros = SubsetData( AstroNeurontest, ident.use = levels(AstroNeurontest@active.ident)[c(1,2,3, 11:17)])


Combinedtest = SetIdent(Combinedtest, cells = WhichCells(Combinedtest, idents = c("Neuronal lineage")), value = "Neuron-like cells")
Combinedtest = SetIdent(Combinedtest, cells = WhichCells(Combinedtest, idents =  c("Astrocyte lineage")), value = "Astrocyte-like cells")


library(ggplot2)
library(Seurat)
fontsize = 15

#emf(file = "majorCellTypeClusters.emf",width = 7, height = 7)
pdf(file = "majorCellTypeClusters.pdf",   # The directory you want to save the file in
    width = 580, # The width of the plot in inches
    height = 351) # The height of the plot in inches

DimPlot(Combinedtest) + theme(axis.text.y  = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
axis.text.x  = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
axis.title.y = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
axis.title.x = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
plot.title = element_text(hjust = 0.5, size=fontsize, face = "bold", family = "Times New Roman"),
legend.position = "right",
panel.background = element_blank(),
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
axis.line = element_blank(),
panel.border = element_rect(colour = "black", fill=NA, size=1, linetype="solid")
)

dev.off()


pdf(file = "majorCellTypeClusters_grouped_by_sample.pdf",   # The directory you want to save the file in
    width = 500, # The width of the plot in inches
    height = 351) # The height of the plot in inches

DimPlot(Combinedtest, group.by = "sample") + theme(axis.text.y  = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
axis.text.x  = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
axis.title.y = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
axis.title.x = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
plot.title = element_text(hjust = 0.5, size=fontsize, face = "bold", family = "Times New Roman"),
legend.position = "right",
panel.background = element_blank(),
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
axis.line = element_blank(),
panel.border = element_rect(colour = "black", fill=NA, size=1, linetype="solid")
)

dev.off()


par(mar = c(5, 5, 5, 10))

pdf(file = "majorCellTypeHeatmap.pdf",   # The directory you want to save the file in
    width = 1527#1492, # The width of the plot in inches
    height = 793) #830) # The height of the plot in inches
DoHeatmap(Combinedtest, features = genes.viz)+ theme(axis.text.y  = element_text(size=12, face = 4, family = "Times New Roman"),axis.text.x  = element_text(size=12, face = "bold", family = "Times New Roman"), axis.title.y = element_text(size=12, face = "bold", family = "Times New Roman"), axis.title.x = element_text(size=12, face = "bold", family = "Times New Roman"), plot.title = element_text(hjust = 0.5, size=12, face = "bold", family = "Times New Roman"), legend.position = "right", panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.line = element_blank())

dev.off()


pdf(file = "AstroNeuroCellTypeClusters.pdf",   # The directory you want to save the file in
    width = 642, # The width of the plot in inches
    height = 425) # The height of the plot in inches

fontsize = 15
DimPlot(AstroNeurontest, label = TRUE) + theme(axis.text.y  = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
axis.text.x  = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
axis.title.y = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
axis.title.x = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
plot.title = element_text(hjust = 0.5, size=fontsize, face = "bold", family = "Times New Roman"),
legend.position = "right",
panel.background = element_blank(),
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
axis.line = element_blank(),
panel.border = element_rect(colour = "black", fill=NA, size=1, linetype="solid")
)

dev.off()




pdf(file = "AstroNeuroCellTypeClusters_grouped_by_sample.pdf",   # The directory you want to save the file in
    width = 600, # The width of the plot in inches
    height = 425) # The height of the plot in inches

fontsize = 15
DimPlot(AstroNeurontest, group.by="sample") + theme(axis.text.y  = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
axis.text.x  = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
axis.title.y = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
axis.title.x = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
plot.title = element_text(hjust = 0.5, size=fontsize, face = "bold", family = "Times New Roman"),
legend.position = "right",
panel.background = element_blank(),
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
axis.line = element_blank(),
panel.border = element_rect(colour = "black", fill=NA, size=1, linetype="solid")
)

dev.off()

load("/Volumes/Transcend/Pascal/Figures/AstroNeurontest_newAnnotations.RData")

IdentsToPlot = setdiff(Idents(AstroNeurontest), c("A-stage7", "A-stage6", "A-stage5", "N-stage7"))
CellsToPlot = WhichCells(AstroNeurontest, idents =  IdentsToPlot)
pdf(file = "AstroPseudotime.pdf",   # The directory you want to save the file in
    width = 543, # The width of the plot in inches
    height = 436) # The height of the plot in inches

FeaturePlot(AstroNeurontest, cells = CellsToPlot,  features = c("AstroLineage"), cols = c("yellow2", "red2")) + theme(axis.text.y  = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
axis.text.x  = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
axis.title.y = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
axis.title.x = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
plot.title = element_text(hjust = 0.5, size=fontsize, face = "bold", family = "Times New Roman"),
legend.position = "right",
panel.background = element_blank(),
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
axis.line = element_blank(),
panel.border = element_rect(colour = "black", fill=NA, size=1, linetype="solid")
)

dev.off()



pdf(file = "NeuroPseudotime.pdf",   # The directory you want to save the file in
    width = 543, # The width of the plot in inches
    height = 436) # The height of the plot in inches

FeaturePlot(AstroNeurontest,  cells = CellsToPlot, features = c("NeuroLineage"), cols = c("yellow2", "red2")) + theme(axis.text.y  = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
axis.text.x  = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
axis.title.y = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
axis.title.x = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
plot.title = element_text(hjust = 0.5, size=fontsize, face = "bold", family = "Times New Roman"),
legend.position = "right",
panel.background = element_blank(),
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
axis.line = element_blank(),
panel.border = element_rect(colour = "black", fill=NA, size=1, linetype="solid")
)

dev.off()


pdf(file = "A4Pseudotime.pdf",   # The directory you want to save the file in
    width = 543, # The width of the plot in inches
    height = 436) # The height of the plot in inches

FeaturePlot(AstroNeurontest, features = c("A4Lineage"), cols = c("yellow2", "red2")) + theme(axis.text.y  = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
axis.text.x  = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
axis.title.y = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
axis.title.x = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
plot.title = element_text(hjust = 0.5, size=fontsize, face = "bold", family = "Times New Roman"),
legend.position = "right",
panel.background = element_blank(),
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
axis.line = element_blank(),
panel.border = element_rect(colour = "black", fill=NA, size=1, linetype="solid")
)

dev.off()


term = "glial cell development"
genes_glial_cell_development = strsplit("Mki67, Id4, Kcnj10, Ntrk2, Gfap, Cspg5, Phgdh, Lamb2, Clu, Mt3, Ascl1, Vim, Hes5, Dag1, Smo, Mt3, Vim, Phgdh, Ntrk2, Clu, Lamb2, Id4, Gfap, Ascl1, Fgfr3, Dmd, Cspg5, Ntrk2, Lgi4, Plp1, Lrp1, Ldlr, Wasf3, Tcf7l2, Agt, Fgfr3, Omg, Cspg5, Clu, Ntrk2, Dmd, Mt3, Hes5, Lgi4, Dag1, Lamb2, Olig1, Phgdh, Id4, Egfr, Nr1d1, Kcnj10, Lrp1, Tcf7l2, Ldlr, Wasf3, Pou3f2, Omg, Cspg5, Fgfr3, Clu, Dmd, Pmp22, Id4, Mt3, Lgi4, Gfap, Ntrk2, Wasf3, Phgdh, Agt, Dag1, Ldlr, Hes5, Pou3f2, Akt2, Id2, Lrp1, Clu, Mt3, Cspg5, Ntrk2, Id4, Gfap, Phgdh, Fgfr3, Plp1, Kcnj10, Mt3, Id4, Phgdh, Gfap, Clu, Id2, Nr1d1, Akt2", ", " )[[1]]


pdf(file = "DotPlot_glial_cell_development.pdf",   # The directory you want to save the file in
    width = 915, # The width of the plot in inches
    height = 420) # The height of the plot in inches

fontsize = 15
DotPlot(Astros,  idents = IdentsToPlot, features = unique(intersect(AstroNeuroMarkersSCT[which(AstroNeuroMarkersSCT$pct.1 > 0),]$gene, genes_glial_cell_development))) + theme(axis.text.x  = element_text(size=fontsize, face = 4, family = "Times New Roman", angle = 90, hjust = 1),axis.text.y  = element_text(size=fontsize, face = "bold", family = "Times New Roman"), axis.title.y = element_text(size=fontsize, face = "bold", family = "Times New Roman"), axis.title.x = element_text(size=fontsize, face = "bold", family = "Times New Roman"), plot.title = element_text(hjust = 0.5, size=fontsize, face = "bold", family = "Times New Roman"), legend.position = "right", panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

dev.off()

pdf(file = "DotPlot_glial_cell_development_flipped.pdf",   # The directory you want to save the file in
    width = 501, # The width of the plot in inches
    height = 791) # The height of the plot in inches

DotPlot(Astros,  idents = IdentsToPlot, features = unique(intersect(AstroNeuroMarkersSCT[which(AstroNeuroMarkersSCT$pct.1 > 0),]$gene, genes_glial_cell_development))) + theme(axis.text.x  = element_text(size=fontsize, face = "bold", family = "Times New Roman", angle = 90, hjust = 1),axis.text.y  = element_text(size=fontsize, face = 4, family = "Times New Roman"), axis.title.y = element_text(size=fontsize, face = "bold", family = "Times New Roman"), axis.title.x = element_text(size=fontsize, face = "bold", family = "Times New Roman"), plot.title = element_text(hjust = 0.5, size=fontsize, face = "bold", family = "Times New Roman"), legend.position = "right", panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())+coord_flip()

dev.off()


term = "positive regulation of neuron differentiation"
genes_positive_regulation_of_neuron_differentiation = strsplit("Mki67, Nbl1, Ascl1, Neurog2, Ranbp1, Mif, Ptbp1, Tgif2, Fxn, Elavl4, Rpl4, D16Ertd472e, Nfe2l2, Cdh4, Nme1, Tcf3, Itgb1, Hmg20b, Eef1a1, Ezh2, Mdk, Fezf2, Ncoa1, Zeb1, Mdk, Elavl4, Sox11, Pcp4, Robo2, Neurog2, Eef1a1, Epha4, Dcx, Flna, Rpl4, Arhgef2, Bcl11a, Itgb1, Sez6, Dpysl3, Ptbp1, Ehd1, Enc1, Cbfa2t2, Llph, Nedd4l, Rnd2, Foxg1, Hspa5, Dcx, Pcp4, Sox11, Neurod1, Tubb2b, Dpysl3, Epha4, Rnd2, Sez6, Enc1, Elavl4, Map1b, Tcf4, Nedd4l, Cd24a, Dixdc1, Ehd1, Stmn2, Limk1, Rtn4, Ptprd, Plxna3, Dbn1, Mdk, Arhgef2, Palm, Bcl11a, Parp6, Gsk3b, Eef1a1, Trak1, Rufy3, Kdm1a, Foxg1, Tubb2b, Stmn2, Neurod2, Neurod1, Plxna4, Cd24a, Dcx, Dpysl3, Ptprd, Islr2, Map1b, Sox11, Epha4, Tcf4, Cnr1, Prox1, Nrp1, Gdpd5, Adra2c, Zeb2, Socs2, Enc1, Snap91, Mapt, Dbn1, Ehd1, L1cam, Palm, Rtn4, Pacsin1, Plxna2, Trim32, Dixdc1, Eef1a1, Apbb1, Kif3c, App, Gsk3b, Nedd4l, Trak1, Fyn, Itsn1, Lrp8, Foxg1, Sema5a, Camk2b, Prox1, Plxna4, Neurod2, Cnr1, Islr2, Syt2, Adra2c, Socs2, Cd24a, Plxna2, Neurod1, Sh3gl3, Lrrc7, Stmn2, Mapt, Eef1a1, L1cam, Map1b, Shank1, Tubb2b, Tcf4, Pak3, Zeb2, Dcx, Nin, Kif3c, Trim67, Apbb1, Sox11, Ntrk3, Ptprd, Camk1, Epha4, Enc1, Fgfr1, Rpl4, Fyn, Trim32, Skil, Ehd1, Nrp1, Dbn1, Tcf4, Sema5a, Plxna2, Zeb2, Camk2b, Prox1, Fgfr1, Map1b, Epha4, Ahi1, Ptprd, Dhx36, Rtn4, Caprin1, Islr2, Eef1a1, Tubb2b, Rpl4, Stmn2, Dpysl3, Eif4g2, Fez1, Neurod1, Pcp4, Arf1, Mif, Ranbp1, Map1b, Mapt, Ndnf, Reln, Cxcl12, Syt4, Gdf5, Rit2, Syt1, Plxnd1, Scn1b, Ache, Nrg1, Baiap2, Robo2, App, Map1b, Bend6, Cdh4, Ndrg4, Bcl11a, Impact, Pcp4, Epha3, Stmn2, L1cam, Pacsin1, Islr2, Ahi1, Dab2ip, Nap1l2, Ptprf, Timp2, Dixdc1, Rufy3, Elavl4, Magi2, Actr3, Ncoa1, Negr1, Rtn4", ", " )[[1]]



pdf(file = "DotPlot_positive_regulation_of_neuron_differentiation.pdf",   # The directory you want to save the file in
    width = 1616, # The width of the plot in inches ## 1655
    height = 420) # The height of the plot in inches ##395

fontsize = 15
DotPlot(Neurons, features = unique(intersect(AstroNeuroMarkersSCT[which(AstroNeuroMarkersSCT$pct.1 > 0.65),]$gene, genes_positive_regulation_of_neuron_differentiation))) + theme(axis.text.x  = element_text(size=fontsize, face = 4, family = "Times", angle = 90, hjust = 1),axis.text.y  = element_text(size=fontsize, face = "bold", family = "Times"), axis.title.y = element_text(size=fontsize, face = "bold", family = "Times"), axis.title.x = element_text(size=fontsize, face = "bold", family = "Times"), plot.title = element_text(hjust = 0.5, size=fontsize, face = "bold", family = "Times"), legend.position = "right", panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

dev.off()

genes_positive_regulation_of_neuron_differentiation = strsplit("Mki67, Nbl1, Ascl1, Neurog2, Ranbp1, Mif, Ptbp1, Tgif2, Fxn, Elavl4, Rpl4, D16Ertd472e, Nfe2l2, Cdh4, Nme1, Tcf3, Itgb1, Hmg20b, Eef1a1, Ezh2, Mdk, Fezf2, Ncoa1, Zeb1, Mdk, Elavl4, Sox11, Pcp4, Neurog2, Eef1a1, Epha4, Dcx, Flna, Rpl4, Arhgef2, Bcl11a, Itgb1, Sez6, Dpysl3, Ptbp1, Ehd1, Enc1, Cbfa2t2, Llph, Nedd4l, Rnd2, Foxg1, Hspa5, Dcx, Pcp4, Sox11, Neurod1, Tubb2b, Dpysl3, Epha4, Rnd2, Sez6, Enc1, Elavl4, Map1b, Tcf4, Nedd4l, Cd24a, Dixdc1, Ehd1, Stmn2, Limk1, Rtn4, Ptprd, Plxna3, Dbn1, Mdk, Arhgef2, Palm, Bcl11a, Parp6, Gsk3b, Eef1a1, Trak1, Rufy3, Kdm1a, Foxg1, Tubb2b, Stmn2, Neurod2, Neurod1, Plxna4, Cd24a, Dcx, Dpysl3, Ptprd, Islr2, Map1b, Sox11, Epha4, Tcf4, Cnr1, Prox1, Nrp1, Gdpd5, Adra2c, Zeb2, Socs2, Enc1, Snap91, Mapt, Dbn1, Ehd1, L1cam, Palm, Rtn4, Pacsin1, Plxna2, Trim32, Dixdc1, Eef1a1, Apbb1, Kif3c, App, Gsk3b, Nedd4l, Trak1, Fyn, Itsn1, Lrp8, Foxg1, Sema5a, Camk2b, Prox1, Plxna4, Neurod2, Cnr1, Islr2, Syt2, Adra2c, Socs2, Cd24a, Plxna2, Neurod1, Sh3gl3, Lrrc7, Stmn2, Mapt, Eef1a1, L1cam, Map1b, Shank1, Tubb2b, Tcf4, Pak3, Zeb2, Dcx, Nin, Kif3c, Trim67, Apbb1, Sox11, Ntrk3, Ptprd, Camk1, Epha4, Enc1, Fgfr1, Rpl4, Fyn, Trim32, Skil, Ehd1, Nrp1, Dbn1, Tcf4, Sema5a, Plxna2, Zeb2, Camk2b, Prox1, Fgfr1, Map1b, Epha4, Ahi1, Ptprd, Dhx36, Rtn4, Caprin1, Islr2, Eef1a1, Tubb2b, Rpl4, Stmn2, Dpysl3, Eif4g2, Fez1, Neurod1, Pcp4, Arf1, Mif, Ranbp1, Map1b, Mapt, Syt4, Gdf5, Rit2, Plxnd1, Scn1b, Ache, Nrg1, Baiap2,  App, Map1b, Bend6, Cdh4, Ndrg4, Bcl11a, Impact, Pcp4, Epha3, Stmn2, L1cam, Pacsin1, Islr2, Ahi1, Dab2ip, Nap1l2, Ptprf, Timp2, Dixdc1, Rufy3, Elavl4, Magi2, Actr3, Ncoa1, Negr1, Rtn4", ", " )[[1]]


pdf(file = "DotPlot_positive_regulation_of_neuron_differentiation_flipped.pdf",   # The directory you want to save the file in
    width = 501, # The width of the plot in inches
    height = 1394) # The height of the plot in inches

DotPlot(Neurons,  idents = IdentsToPlot, features = unique(intersect(AstroNeuroMarkersSCT[which(AstroNeuroMarkersSCT$pct.1 > 0.65),]$gene, genes_positive_regulation_of_neuron_differentiation)))+ theme(axis.text.x  = element_text(size=fontsize, face = "bold", family = "Times New Roman", angle = 90, hjust = 1),axis.text.y  = element_text(size=fontsize, face = 4, family = "Times New Roman"), axis.title.y = element_text(size=fontsize, face = "bold", family = "Times New Roman"), axis.title.x = element_text(size=fontsize, face = "bold", family = "Times New Roman"), plot.title = element_text(hjust = 0.5, size=fontsize, face = "bold", family = "Times New Roman"), legend.position = "right", panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())+coord_flip()

dev.off()


#### Barplots of origins


fontsize = 15

emf(file = "Origin_barplot_astrocytes_and_neurons.emf",width = 7, height = 7
#pdf("Origin_barplot_astrocytes_and_neurons.pdf",)
plotOriginSampleBarplots(AstroNeurontest) ## in SeuratRun.R

dev.off()


pdf("Origin_barplot_high_level_cell_types.pdf",
    width = 512,
    height = 431)
plotOriginSampleBarplots(Combinedtest) ## in SeuratRun.R

dev.off()


### stipecharts
lineageGenes$pch = 1
lineageGenes[which(lineageGenes$color != "gray"),]$pch = 16
lineageGenes$cex = 0.3
lineageGenes[which(lineageGenes$color != "gray"),]$cex = 1

par(mar = c(5, 10, 5, 5))
pdf("Stripechart_DEGs.pdf",
    width = 437,
    height = 418)
#stripeChart(lineageGenes[which(lineageGenes$p_val_adj < 0.05),],  "cluster", "avg_logFC", 0, avg_logFCth, "p_val_adj", c(-0.5,0.5), c(0.5,4.5) )
#stripeChart(lineageGenes[which(lineageGenes$p_val_adj < 0.05 & ( lineageGenes$avg_logFC < -0.2 | lineageGenes$avg_logFC > 0.2)),],  "cluster", "avg_logFC", 0, avg_logFCth, "p_val_adj", c(-0.5,0.5), c(0.5,4.5) )
#stripeChart(lineageGenes,  "cluster", "avg_logFC", 0, avg_logFCth, "p_val_adj", c(-0.5,0.5), c(0.5,4.5) )
stripeChart(lineageGenes,  "cluster", "avg_logFC", 0, avg_logFCth, "p_val_adj", c(-0.5,0.5), c(0.5,4.5), pch = lineageGenes$pch, cex = lineageGenes$cex )

dev.off()

###
TFDiseaseRelevant$cl = 0
for (i in c(1:13)){
    TFDiseaseRelevant[which(TFDiseaseRelevant$Lineage == intersect( levels(AstroNeurontest), unique(TFDiseaseRelevant$Lineage))[i]),]$cl = i
}

### heatmap
TFDiseaseRelevant = lineageGenes[which( lineageGenes$color != "gray"),]
TFDiseaseRelevantall = TFDiseaseRelevant
TFDiseaseRelevant =  rbind(TFDiseaseRelevantall[which(TFDiseaseRelevantall$cl ==1 & TFDiseaseRelevantall$pct.1 > 0.7 & abs(TFDiseaseRelevantall$avg_logFC)> 0.3 ),] , TFDiseaseRelevantall[which(TFDiseaseRelevantall$cl > 1 ),])
TFDiseaseRelevant =  TFDiseaseRelevantall[which(TFDiseaseRelevantall$cl > 1 ),]

TFDiseaseRelevant = TFDiseaseRelevant[which( !grepl( "Rpl" , TFDiseaseRelevant$Gene) ), ]
TFDiseaseRelevant = TFDiseaseRelevant[which( !grepl( "Rps" , TFDiseaseRelevant$Gene) ), ]
TFDiseaseRelevant = TFDiseaseRelevant[which( !grepl( "mt" , TFDiseaseRelevant$Gene) ), ]
TFDiseaseRelevant = TFDiseaseRelevant[which( !grepl( "Mir" , TFDiseaseRelevant$Gene) ), ]
TFDiseaseRelevant = TFDiseaseRelevant[which( !grepl( "Gm" , TFDiseaseRelevant$Gene) ), ]
TFDiseaseRelevant = TFDiseaseRelevant[which( !grepl( "AC" , TFDiseaseRelevant$Gene) ), ]
TFDiseaseRelevant = TFDiseaseRelevant[which( !grepl( "CT" , TFDiseaseRelevant$Gene) ), ]
TFDiseaseRelevant = TFDiseaseRelevant[which( !grepl( "ik" , TFDiseaseRelevant$Gene) ), ]

TFDiseaseRelevant <- TFDiseaseRelevant[order( TFDiseaseRelevant$cl, TFDiseaseRelevant$avg_logFC  ),]
#TFDiseaseRelevant <- TFDiseaseRelevant[order( TFDiseaseRelevant$lineage, TFDiseaseRelevant$avg_logFC  ),]

rownames(TFDiseaseRelevant) <- paste0 (TFDiseaseRelevant$Gene, "." , TFDiseaseRelevant$cl)
datacluster <- as.data.frame(matrix(0, ncol = length(TFDiseaseRelevant$Gene), nrow = length(unique(TFDiseaseRelevant$cl)) ))
colnames(datacluster) <- rownames(TFDiseaseRelevant)
rownames(datacluster) <- unique(TFDiseaseRelevant$cl)
for( var in colnames(datacluster) ) {
  print(var)
  datacluster[ as.character(TFDiseaseRelevant[var,]$cl), var] <- TFDiseaseRelevant[var,]$avg_logFC
}

pdf("Heatmap_DEGs.pdf",
width = 463,
height = 793)

heatmap.2(t(as.matrix(datacluster)), dendrogram='none', Rowv=FALSE, Colv=FALSE,trace='none',
          density.info = 'none', key = TRUE, keysize = 1, xlab = NULL,col=colorRampPalette(c("orange1", "gray89", "darkseagreen4"))(n = 30),
          margins = c(11,13), family = "Times", cexRow = 1.05, cexCol = 1.05, labRow=as.expression(lapply(TFDiseaseRelevant$Gene, function(a) bquote(italic(.(a))))), key.xlab = "ln fold change", key.title = "",  lhei=c(1,4.5), lwid=c(2,3.5), key.par = list(cex=1),
          labCol =  as.character(unique(TFDiseaseRelevant$lineage_to_plot)) )
dev.off()


### new version

upregulated = TFDiseaseRelevant[which(TFDiseaseRelevant$avg_log2FC > 0),]
upregulated = upregulated %>% group_by(Lineage) %>% top_n(10, avg_log2FC)

## upregulated
datacluster <- as.data.frame(matrix(0, ncol = length( unique(upregulated$Gene)), nrow = length(unique(upregulated$cl)) ))
colnames(datacluster) <- unique(upregulated$Gene)
rownames(datacluster) <- intersect( levels(AstroNeurontest), unique(upregulated$Lineage))


for(cluster in rownames(datacluster)){
    for( var in colnames(datacluster) ) {
        print(var)
        varrr = upregulated[which(upregulated$Gene ==var & upregulated$Lineage == cluster),]
        ifelse( length(varrr$avg_log2FC) > 0,
        datacluster[ cluster, var] <- upregulated[which(upregulated$Gene ==var & upregulated$Lineage == cluster),]$avg_log2FC, 0)
    }
}


pdf("SubFigure13_Heatmap_DEGs_upregulated.pdf",
width = 601,
height = 821)
heatmap.2(t(as.matrix(datacluster)), dendrogram='none', Rowv=FALSE, Colv=FALSE,trace='none',
          density.info = 'none', key = TRUE, keysize = 1, xlab = NULL,col=colorRampPalette(c("orange1", "gray89", "darkseagreen4"))(n = 30),
          margins = c(11,13), family = "Times", cexRow = 1.05, cexCol = 1.05, labRow=as.expression(lapply(colnames(datacluster), function(a) bquote(italic(.(a))))), key.xlab = "ln fold change", key.title = "",  lhei=c(1,4.5), lwid=c(2,3.5), key.par = list(cex=1),
          labCol =  rownames(datacluster) )
dev.off()

## downregulated

downregulated = TFDiseaseRelevant[which(TFDiseaseRelevant$avg_log2FC < 0),]
downregulated = downregulated %>% group_by(Lineage) %>% top_n(20, avg_log2FC)

datacluster <- as.data.frame(matrix(0, ncol = length( unique(downregulated$Gene)), nrow = length(unique(downregulated$cl)) ))
colnames(datacluster) <- unique(downregulated$Gene)
rownames(datacluster) <- intersect( levels(AstroNeurontest), unique(downregulated$Lineage))



for(cluster in rownames(datacluster)){
    for( var in colnames(datacluster) ) {
        print(var)
        varrr = upregulated[which(downregulated$Gene ==var & downregulated$Lineage == cluster),]
        ifelse( length(varrr$avg_log2FC) > 0,
        datacluster[ cluster, var] <- downregulated[which(downregulated$Gene ==var & downregulated$Lineage == cluster),]$avg_log2FC, 0)
    }
}

pdf("SubFigure13_Heatmap_DEGs_downregulated.pdf",
width = 601,
height = 821)
heatmap.2(t(as.matrix(datacluster)), dendrogram='none', Rowv=FALSE, Colv=FALSE,trace='none',
          density.info = 'none', key = TRUE, keysize = 1, xlab = NULL,col=colorRampPalette(c("orange1", "gray89", "darkseagreen4"))(n = 30),
          margins = c(11,13), family = "Times", cexRow = 1.05, cexCol = 1.05, labRow=as.expression(lapply(colnames(datacluster), function(a) bquote(italic(.(a))))), key.xlab = "ln fold change", key.title = "",  lhei=c(1,4.5), lwid=c(2,3.5), key.par = list(cex=1),
          labCol =  rownames(datacluster) )
dev.off()




pdf("Heatmap_DEGs.pdf",
width = 463,
height = 793)

heatmap.2(t(as.matrix(datacluster)), dendrogram='none', Rowv=FALSE, Colv=FALSE,trace='none',
          density.info = 'none', key = TRUE, keysize = 1, xlab = NULL,col=colorRampPalette(c("orange1", "gray89", "darkseagreen4"))(n = 30),
          margins = c(11,13), family = "Times", cexRow = 1.05, cexCol = 1.05, labRow=as.expression(lapply(TFDiseaseRelevant$Gene, function(a) bquote(italic(.(a))))), key.xlab = "ln fold change", key.title = "",  lhei=c(1,4.5), lwid=c(2,3.5), key.par = list(cex=1),
          labCol =  as.character(unique(TFDiseaseRelevant$lineage_to_plot)) )
dev.off()



#### spatial data

### without ML and devided in 3 intervals
load("/Volumes/Transcend/Pascal/Figures/seurat_spatial_final_withoutML_final.RData")
load("/Volumes/Transcend/Pascal/Figures/AstroNeurontest_newAnnotations.RData")


var = names(seurat_obj_list)[3]


pdf(paste0("spatial_segmentation_", var, ".pdf"),
width = 278,
height = 358)

DimPlot(seurat_obj_list[[var]], reduction = "coord", group.by = "region", cols = c("red3", "red4",  "orange3", "blue4") ) +  ggtitle(var)+ theme(axis.text.y  = element_blank(),
axis.text.x  = element_blank(),
axis.title.y = element_blank(),
axis.title.x = element_blank(),
plot.title = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
panel.background = element_blank(),
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
axis.line = element_blank(),
panel.border = element_blank(),
axis.ticks.x=element_blank(),
axis.ticks.y=element_blank()
) + scale_color_discrete( type = c( "indianred1", "hotpink4", "orange3", "blue4"), labels=c("Granular_Layer_1" = "GL1", "Granular_Layer_2" = "GL2","Hilius" = "Hilus", "outer" = "ML", "SGZ" = "SGZ" ))

dev.off()



### merge clusters

pdf("Merging_cluster_round1.pdf"),
width = 648,
height = 466)

DimPlot(mergeclusters_round1, group.by = "sample") + theme(axis.text.y  = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
axis.text.x  = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
axis.title.y = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
axis.title.x = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
plot.title = element_text(hjust = 0.5, size=fontsize, face = "bold", family = "Times New Roman"),
legend.position = "right",
panel.background = element_blank(),
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
axis.line = element_blank(),
panel.border = element_rect(colour = "black", fill=NA, size=1, linetype="solid")
)

dev.off()


pdf("Merging_cluster_round1_labeled.pdf"),
width = 648,
height = 466)

DimPlot(mergeclusters_round1,label=TRUE) + theme(axis.text.y  = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
axis.text.x  = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
axis.title.y = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
axis.title.x = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
plot.title = element_text(hjust = 0.5, size=fontsize, face = "bold", family = "Times New Roman"),
legend.position = "right",
panel.background = element_blank(),
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
axis.line = element_blank(),
panel.border = element_rect(colour = "black", fill=NA, size=1, linetype="solid")
)

dev.off()



pdf("Merging_clusters_grouped.pdf"),
width = 476,
height = 357)

DimPlot(mergeclusters, group.by = "sample") + theme(axis.text.y  = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
axis.text.x  = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
axis.title.y = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
axis.title.x = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
plot.title = element_text(hjust = 0.5, size=fontsize, face = "bold", family = "Times New Roman"),
legend.position = "right",
panel.background = element_blank(),
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
axis.line = element_blank(),
panel.border = element_rect(colour = "black", fill=NA, size=1, linetype="solid")
)

dev.off()


pdf("Merging_cluster_labeled.pdf"),
width = 476,
height = 357)

DimPlot(mergeclusters,label=TRUE) + theme(axis.text.y  = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
axis.text.x  = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
axis.title.y = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
axis.title.x = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
plot.title = element_text(hjust = 0.5, size=fontsize, face = "bold", family = "Times New Roman"),
legend.position = "right",
panel.background = element_blank(),
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
axis.line = element_blank(),
panel.border = element_rect(colour = "black", fill=NA, size=1, linetype="solid")
)

dev.off()


pdf("Merging_cluster_labeled_previousID.pdf"),
width = 720,
height = 509)

DimPlot(mergeclusters, label=TRUE, group.by = "previousID", label.size = 5, pt.size = 1.5)  + theme(axis.text.y  = element_text(size=fontsize, face = "bold", family = "Times New Roman"), axis.text.x  = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
axis.title.y = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
axis.title.x = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
plot.title = element_text(hjust = 0.5, size=fontsize, face = "bold", family = "Times New Roman"),
legend.position = "right",
panel.background = element_blank(),
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
axis.line = element_blank(),
panel.border = element_rect(colour = "black", fill=NA, size=1, linetype="solid")
)

dev.off()


pdf("Merging_cluster_labeled_split_by_tech.pdf"),
width = 1168,
height = 508)

DimPlot(mergeclusters, group.by = "previousID", split.by = "sample", label=TRUE, label.size = 4.5, pt.size = 1.7)  + theme(axis.text.y  = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
axis.text.x  = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
axis.title.y = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
axis.title.x = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
plot.title = element_text(hjust = 0.5, size=fontsize, face = "bold", family = "Times New Roman"),
legend.position = "right",
panel.background = element_blank(),
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
axis.line = element_blank(),
panel.border = element_rect(colour = "black", fill=NA, size=1, linetype="solid")
)


dev.off()


pdf("Merging_cluster_Neuro_Lineage.pdf"),
width = 510,
height = 413)


FeaturePlot(mergeclusters, c( "NeuroLineage"), pt.size = 0.7, cols =  c("yellow2", "red2")) + theme(axis.text.y  = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
axis.text.x  = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
axis.title.y = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
axis.title.x = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
plot.title = element_text(hjust = 0.5, size=fontsize, face = "bold", family = "Times New Roman"),
legend.position = "right",
panel.background = element_blank(),
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
axis.line = element_blank(),
panel.border = element_rect(colour = "black", fill=NA, size=1, linetype="solid")
)


dev.off()




pdf("Merging_cluster_Astro_Lineage.pdf"),
width = 510,
height = 413)


FeaturePlot(mergeclusters, c( "AstroLineage"), pt.size = 0.7, cols =  c("yellow2", "red2")) + theme(axis.text.y  = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
axis.text.x  = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
axis.title.y = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
axis.title.x = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
plot.title = element_text(hjust = 0.5, size=fontsize, face = "bold", family = "Times New Roman"),
legend.position = "right",
panel.background = element_blank(),
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
axis.line = element_blank(),
panel.border = element_rect(colour = "black", fill=NA, size=1, linetype="solid")
)


dev.off()


### deviding in 4
mergeclusters@meta.data$Neuro_lineage_discrete = cut(mergeclusters@meta.data$NeuroLineage, breaks = c(0, max(NSCRG@meta.data$NeuroLineage) , max(NSCRG@meta.data$NeuroLineage)+ (NeuroLineageMax - max(NSCRG@meta.data$NeuroLineage))/3, max(NSCRG@meta.data$NeuroLineage)+ 2*(NeuroLineageMax - max(NSCRG@meta.data$NeuroLineage))/3, NeuroLineageMax ) )
mergeclusters@meta.data$Astro_lineage_discrete = cut(mergeclusters@meta.data$AstroLineage, breaks = c(0, max(NSCRG@meta.data$AstroLineage) , max(NSCRG@meta.data$AstroLineage)+ (AstroLineageMax - max(NSCRG@meta.data$AstroLineage))/3 , max(NSCRG@meta.data$AstroLineage)+2*(AstroLineageMax - max(NSCRG@meta.data$AstroLineage))/3, AstroLineageMax) )


pdf("Merging_cluster_Neuro_Lineage_Discrete_4parts.pdf"),
width = 582,
height = 413)

DimPlot(mergeclusters, group.by = "Neuro_lineage_discrete", pt.size = 0.7) + theme(axis.text.y  = element_text(size=fontsize, face = "bold", family = "Times New Roman"), axis.text.x  = element_text(size=fontsize, face = "bold", family = "Times New Roman"), axis.title.y = element_text(size=fontsize, face = "bold", family = "Times New Roman"), axis.title.x = element_text(size=fontsize, face = "bold", family = "Times New Roman"), plot.title = element_text(hjust = 0.5, size=fontsize, face = "bold", family = "Times New Roman"), legend.position = "right",panel.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1, linetype="solid") )+ggtitle("Neuro_lineage_discrete")

dev.off()



pdf("Merging_cluster_Astro_Lineage_Discrete_4parts.pdf"),
width = 582,
height = 413)

DimPlot(mergeclusters, group.by = "Astro_lineage_discrete", pt.size = 0.7) + theme(axis.text.y  = element_text(size=fontsize, face = "bold", family = "Times New Roman"), axis.text.x  = element_text(size=fontsize, face = "bold", family = "Times New Roman"), axis.title.y = element_text(size=fontsize, face = "bold", family = "Times New Roman"), axis.title.x = element_text(size=fontsize, face = "bold", family = "Times New Roman"), plot.title = element_text(hjust = 0.5, size=fontsize, face = "bold", family = "Times New Roman"), legend.position = "right",panel.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1, linetype="solid") )+ggtitle("Astro_lineage_discrete_4")

dev.off()

i=1
pdf(paste0("Neuro_Lineage_part",i,"_dist.pdf")),
width = 775,
height = 296)

ggplot(datatableNeuroLineage[which(datatableNeuroLineage$Neuro_lineage_discrete ==  levels(unique(datatableNeuroLineage$Neuro_lineage_discrete))[i]),], aes(x=region, y=Freq, fill = Neuro_lineage_discrete, group = condition, color = condition)) + stat_compare_means(aes(group = condition), label = "p.format", label.y = 100) + scale_x_discrete(labels=c("Granular_Layer_1" = "GL1", "Granular_Layer_2" = "GL2","Hilius" = "Hilus", "outer" = "ML", "SGZ" = "SGZ" )) + theme(axis.text.y  = element_text(size=15, face = "bold", family = "Times New Roman"), axis.text.x  = element_text(size=15, face = "bold", family = "Times New Roman"), axis.title.y = element_text(size=15, face = "bold", family = "Times New Roman"),
axis.title.x = element_text(size=15, face = "bold", family = "Times New Roman"),
plot.title = element_text(hjust = 0.5, size=15, face = "bold", family = "Times New Roman"),
legend.position = "right",
panel.background = element_blank(),
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
axis.line = element_blank(),
legend.title = element_text(size=20, color = "black"),
legend.text = element_text(size=20),
legend.key=element_rect(fill='white', size = 20),
panel.border = element_rect(colour = "black", fill=NA, size=1, linetype="solid")) + guides(colour = guide_legend(override.aes = list(size=0, stroke=4)))  + ylim(-20,120) + xlab("") + ylab("% of cells")+
    geom_point(position = position_jitterdodge(jitter.width = .2, dodge.width = .7), alpha = .5) +
    stat_summary(fun = mean, na.rm = TRUE, geom = "point", shape = "diamond",
                 size = 4, color = "black", position = position_dodge(width = .7)) +
    stat_summary(fun.data = mean_cl_normal, na.rm = TRUE, geom = "errorbar", width = .2, color = "black",
                 position = position_dodge(width = .7)) +
    scale_color_brewer(palette = "Set1")

dev.off()



i=1
pdf(paste0("Astro_Lineage_part",i,"_dist.pdf")),
width = 775,
height = 296)

ggplot(datatableAstroLineage[which(datatableAstroLineage$Astro_lineage_discrete ==  levels(unique(datatableAstroLineage$Astro_lineage_discrete))[1]),], aes(x=region, y=Freq, fill = Astro_lineage_discrete, group = condition, color = condition)) + stat_compare_means(aes(group = condition), label = "p.format", label.y = 110) + scale_x_discrete(labels=c("Granular_Layer_1" = "GL1", "Granular_Layer_2" = "GL2","Hilius" = "Hilus", "outer" = "ML", "SGZ" = "SGZ" )) + theme(axis.text.y  = element_text(size=15, face = "bold", family = "Times New Roman"), axis.text.x  = element_text(size=15, face = "bold", family = "Times New Roman"), axis.title.y = element_text(size=15, face = "bold", family = "Times New Roman"),
axis.title.x = element_text(size=15, face = "bold", family = "Times New Roman"),
plot.title = element_text(hjust = 0.5, size=15, face = "bold", family = "Times New Roman"),
legend.position = "right",
panel.background = element_blank(),
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
axis.line = element_blank(),
legend.title = element_text(size=20, color = "black"),
legend.text = element_text(size=20),
legend.key=element_rect(fill='white', size = 20),
panel.border = element_rect(colour = "black", fill=NA, size=1, linetype="solid")) + guides(colour = guide_legend(override.aes = list(size=0, stroke=4)))  + ylim(-20,120) + xlab("") + ylab("% of cells")+
    geom_point(position = position_jitterdodge(jitter.width = .2, dodge.width = .7), alpha = .5) +
    stat_summary(fun = mean, na.rm = TRUE, geom = "point", shape = "diamond",
                 size = 4, color = "black", position = position_dodge(width = .7)) +
    stat_summary(fun.data = mean_cl_normal, na.rm = TRUE, geom = "errorbar", width = .2, color = "black",
                 position = position_dodge(width = .7)) +
    scale_color_brewer(palette = "Set1")


dev.off()





#### deviding in 3
NeuroLineageMax = max(na.omit(mergeclusters@meta.data$NeuroLineage))
AstroLineageMax = max(na.omit(mergeclusters@meta.data$AstroLineage))
NSCRG = SubsetData(mergeclusters, ident.use = c(8))

mergeclusters@meta.data$Neuro_lineage_discrete = cut(mergeclusters@meta.data$NeuroLineage, breaks = c(0, max(NSCRG@meta.data$NeuroLineage),max(NSCRG@meta.data$NeuroLineage)+ (NeuroLineageMax - max(NSCRG@meta.data$NeuroLineage))/2, NeuroLineageMax ) )
mergeclusters@meta.data$Astro_lineage_discrete = cut(mergeclusters@meta.data$AstroLineage, breaks = c(0, max(NSCRG@meta.data$AstroLineage) , max(NSCRG@meta.data$AstroLineage)+ (AstroLineageMax - max(NSCRG@meta.data$AstroLineage))/2 , AstroLineageMax) )


# With ML and devided in 3 intervals

load("/Volumes/Transcend/Pascal/Figures/seurat_spatial_final_withML.RData")
var = names(seurat_obj_list)[3]


pdf(paste0("spatial_segmentation_", var, ".pdf"),
width = 278,
height = 358)

DimPlot(seurat_obj_list[[var]], reduction = "coord", group.by = "region", cols = c("red3", "red4",  "orange3", "blue4") ) +  ggtitle(var)+ theme(axis.text.y  = element_blank(),
                                                                                                                                                 axis.text.x  = element_blank(),
                                                                                                                                                 axis.title.y = element_blank(),
                                                                                                                                                 axis.title.x = element_blank(),
                                                                                                                                                 plot.title = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
                                                                                                                                                 panel.background = element_blank(),
                                                                                                                                                 panel.grid.major = element_blank(),
                                                                                                                                                 panel.grid.minor = element_blank(),
                                                                                                                                                 axis.line = element_blank(),
                                                                                                                                                 panel.border = element_blank(),
                                                                                                                                                 axis.ticks.x=element_blank(),
                                                                                                                                                 axis.ticks.y=element_blank()
) + scale_color_discrete( type = c( "indianred1", "hotpink4", "orange3", "darkcyan" , "blue4"), labels=c("Granular_Layer_1" = "GL1", "Granular_Layer_2" = "GL2","Hilius" = "Hilus", "outer" = "ML", "SGZ" = "SGZ" ))

dev.off()


fwrite(NumberOfCellsInHighLevel, file = "Table1_NumberOfCellsInHighLevelClustering.csv")
fwrite(NumberOfCellsInAstrocyticAndNeuronalPopulations, file = "Table2_NumberOfCellsInAstrocyticAndNeuronalPopulations.csv")
fwrite(AstroNeuroMarkersSCT, file = "Table3_AstroNeuroMarkersSCTnormalization.csv")
fwrite(TFDiseaseRelevantall, file = "Table6_disease_related_DGE_list.csv")

datatableNeuroLineage = data.frame(region=character(),
                                   Neuro_lineage_discrete=character(),
                                   Freq=double(),
                                   slide=character(),
                                   stringsAsFactors=FALSE)

for(var in names(datalist)){
    obj = seurat_obj_list[[var]]
    subset = obj@meta.data[c("Astro_lineage_discrete", "Neuro_lineage_discrete", "region")]
    tt = table(subset[, c(2,3)])
    dataframeobj = as.data.frame(tt)
    dataframeobj$slide = var
    datatableNeuroLineage = rbind(datatableNeuroLineage, dataframeobj)
}
datatableNeuroLineage$condition = unlist(strsplit(datatableNeuroLineage$slide, "_"))[seq(2,length(datatableNeuroLineage$slide)*3,3)]

fwrite(datatableNeuroLineage, file = "Table8_The_absolute_number_of_neuronal_progenitors_per_pseudotime_location_tissue.csv")


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
    tt = table(subset[, c(1,3)])
    dataframeobj = as.data.frame(tt)
    dataframeobj$slide = var
    datatableAstroLineage = rbind(datatableAstroLineage, dataframeobj)
}
datatableAstroLineage$condition = unlist(strsplit(datatableAstroLineage$slide, "_"))[seq(2,length(datatableAstroLineage$slide)*3,3)]


fwrite(datatableAstroLineage, file = "Table9_The_absolute_of_astrocytic_progenitors_per_pseudotime_location_tissue.csv")



### additional Figures for the supplementary information

Shamcells = rownames(Combinedtest@meta.data[which(Combinedtest@meta.data$sample == "Sham"),])
TBIcells = rownames(Combinedtest@meta.data[which(Combinedtest@meta.data$sample == "TBI"),])

TBICombined_test = SubsetData(Combinedtest, cells = TBIcells)
ShamCombined_test = SubsetData(Combinedtest, cells = Shamcells)

histData = data.frame("TBI"=tbihist )

ggplot(histData, aes(x = TBI)) +
    geom_density(aes(y = ..density..))  + theme(axis.text.y  = element_text(size=15, face = "bold", family = "Times New Roman"), axis.text.x  = element_text(size=15, face = "bold", family = "Times New Roman"), axis.title.y = element_text(size=15, face = "bold", family = "Times New Roman"),
                            axis.title.x = element_text(size=15, face = "bold", family = "Times New Roman"),
                            plot.title = element_text(hjust = 0.5, size=15, face = "bold", family = "Times New Roman"),
                            legend.position = "right",
                            panel.background = element_blank(),
                            panel.grid.major = element_blank(),
                            panel.grid.minor = element_blank(),
                            axis.line = element_blank(),
                            legend.title = element_text(size=20, color = "black"),
                            legend.text = element_text(size=20),
                            legend.key=element_rect(fill='white', size = 20),
                            panel.border = element_rect(colour = "black", fill=NA, size=1, linetype="solid")) + guides(colour = guide_legend(override.aes = list(size=0, stroke=4)))  + xlab("Transcript count per cells in TBI animals") + ylab("Density")+geom_vline(xintercept = 200, col = "red" )



#### Figure of markers



pdf("FeaturePlot_HighLevelMarkers.pdf"),
width = 1310,
height = 767)

FeaturePlot(Combinedtest, c("Mki67","Ascl1", "Snap25", "Neurod1","Aldoc", "Slc1a3", "Pdgfra", "Mog", "Pecam1", "Cd68", "Cx3cr1", "Des")) & theme( plot.title = element_text( face = "italic") )

dev.off()



pdf("FeaturePlot_RNA_scopeMarkers.pdf"),
width = 983,
height = 767)

FeaturePlot(AstroNeurotest, c("Hapln1","Frzb", "Ascl1", "Sparc","Sned1", "Neat1", "Slc1a3"))& theme( plot.title = element_text( face = "italic") )

dev.off()



### checking GFAP expression


var = names(seurat_obj_list)[1]


pdf(paste0("spatial_Gfap_expression_", var, ".pdf"),
width = 278,
height = 358)

FeaturePlot(seurat_obj_list[[var]], "Gfap", reduction = "coord") +  ggtitle(paste0(var, "_Gfap"))+ theme(axis.text.y  = element_blank(),
axis.text.x  = element_blank(),
axis.title.y = element_blank(),
axis.title.x = element_blank(),
plot.title = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
panel.background = element_blank(),
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
axis.line = element_blank(),
panel.border = element_blank(),
axis.ticks.x=element_blank(),
axis.ticks.y=element_blank()
)

dev.off()

### Lineage plots

var = names(seurat_obj_list)[1]

pdf(paste0("spatial_AstroLineage_", var, ".pdf"),
width = 278,
height = 358)

seurat_obj_list[[var]]@meta.data$pt.size.astro = ifelse(is.na(seurat_obj_list[[var]]@meta.data$AstroLineage), 0.5, 1)

FeaturePlot(seurat_obj_list[[var]], "AstroLineage", reduction = "coord", cols = c("yellow2", "red2"), pt.size = seurat_obj_list[[var]]@meta.data$pt.size.astro )  +  ggtitle(paste0(var, "_Astro"))+ theme(axis.text.y  = element_blank(),
axis.text.x  = element_blank(),
axis.title.y = element_blank(),
axis.title.x = element_blank(),
plot.title = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
panel.background = element_blank(),
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
axis.line = element_blank(),
panel.border = element_blank(),
axis.ticks.x=element_blank(),
axis.ticks.y=element_blank()
)

dev.off()


var = names(seurat_obj_list)[1]

pdf(paste0("spatial_NeuroLineage_", var, ".pdf"),
width = 278,
height = 358)

seurat_obj_list[[var]]@meta.data$pt.size.neuro = ifelse(is.na(seurat_obj_list[[var]]@meta.data$NeuroLineage), 0.5, 1)

FeaturePlot(seurat_obj_list[[var]], "NeuroLineage", reduction = "coord", cols = c("yellow2", "red2"), pt.size = seurat_obj_list[[var]]@meta.data$pt.size.neuro )  +  ggtitle(paste0(var, "_Neuro"))+ theme(axis.text.y  = element_blank(),
axis.text.x  = element_blank(),
axis.title.y = element_blank(),
axis.title.x = element_blank(),
plot.title = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
panel.background = element_blank(),
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
axis.line = element_blank(),
panel.border = element_blank(),
axis.ticks.x=element_blank(),
axis.ticks.y=element_blank()
)

dev.off()


### cell displacement in Neuronal Lineage

var = names(seurat_obj_list)[1]
lineageRange = 3

pdf(paste0("spatial_Neuro_Lineage_", lineageRange, "_", var, ".pdf"),
width = 278,
height = 358)

cellsToHighlights = rownames(seurat_obj_list[[var]]@meta.data[which(seurat_obj_list[[var]]@meta.data$Neuro_lineage_discrete == levels(seurat_obj_list[[var]]@meta.data$Neuro_lineage_discrete)[lineageRange] ),])

cols.final = as.character(getColors(DimPlot(seurat_obj_list[[var]], group.by = "Neuro_lineage_discrete")))[lineageRange]

DimPlot(seurat_obj_list[[var]], group.by = "Neuro_lineage_discrete", reduction = "coord", cells.highlight = cellsToHighlights,  sizes.highlight = 2) + scale_color_manual(labels  = c("Unselected", levels(seurat_obj_list[[var]]@meta.data$Neuro_lineage_discrete)[lineageRange]), values = c("grey81", cols.final) ) +  ggtitle(paste0(var, "_Neuro")) + theme(axis.text.y  = element_blank(),
axis.text.x  = element_blank(),
axis.title.y = element_blank(),
axis.title.x = element_blank(),
plot.title = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
panel.background = element_blank(),
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
axis.line = element_blank(),
panel.border = element_blank(),
axis.ticks.x=element_blank(),
axis.ticks.y=element_blank()
)

dev.off()


### cell displacement in Astrocytic Lineage

var = names(seurat_obj_list)[1]
lineageRange = 1

pdf(paste0("spatial_Astro_Lineage_", lineageRange, "_", var, ".pdf"),
width = 278,
height = 358)

cellsToHighlights = rownames(seurat_obj_list[[var]]@meta.data[which(seurat_obj_list[[var]]@meta.data$Astro_lineage_discrete == levels(seurat_obj_list[[var]]@meta.data$Astro_lineage_discrete)[lineageRange] ),])

cols.final = as.character(getColors(DimPlot(seurat_obj_list[[var]], group.by = "Astro_lineage_discrete")))[lineageRange]

DimPlot(seurat_obj_list[[var]], group.by = "Astro_lineage_discrete", reduction = "coord", cells.highlight = cellsToHighlights,  sizes.highlight = 2) + scale_color_manual(labels  = c("Unselected", levels(seurat_obj_list[[var]]@meta.data$Astro_lineage_discrete)[lineageRange]), values = c("grey81", cols.final) ) +  ggtitle(paste0(var, "_Astro")) + theme(axis.text.y  = element_blank(),
axis.text.x  = element_blank(),
axis.title.y = element_blank(),
axis.title.x = element_blank(),
plot.title = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
panel.background = element_blank(),
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
axis.line = element_blank(),
panel.border = element_blank(),
axis.ticks.x=element_blank(),
axis.ticks.y=element_blank()
)

dev.off()



### heatmap of the mapping data

pdf("Heatmap_spatial_mapping.pdf",
width = 917,
height = 742)

barplotOfPercentages(metadata = metadataCombined)

dev.off()

### save the table

metadata_tmp = as.data.frame.matrix (table(mergeclusters@meta.data[, c( "seurat_clusters","previousID")]) )
write.csv(metadata_tmp, file = "/Volumes/Transcend/Pascal/Figures/Table12_spatial_to_10x_mapping.csv")


#### spatial mapping
mergeclusters@meta.data$ID10x = "Unmapped"
mergeclusters@meta.data[WhichCells(mergeclusters, idents = c(8)), ]$ID10x = "8, NSC-stages1,2"
mergeclusters@meta.data[WhichCells(mergeclusters, idents = c(7)), ]$ID10x = "7, RG-like"
mergeclusters@meta.data[WhichCells(mergeclusters, idents = c(1)), ]$ID10x = "1, N-stages1,2"
mergeclusters@meta.data[WhichCells(mergeclusters, idents = c(2)), ]$ID10x = "2, N-stage3"
mergeclusters@meta.data[WhichCells(mergeclusters, idents = c(6)), ]$ID10x = "6, N-stages4,5"
mergeclusters@meta.data[WhichCells(mergeclusters, idents = c(10)), ]$ID10x = "10, N-stages5,6"
mergeclusters@meta.data[WhichCells(mergeclusters, idents = c(11)), ]$ID10x = "11, N-stage7"
mergeclusters@meta.data[WhichCells(mergeclusters, idents = c(9)), ]$ID10x= "9, A-stages1,5"
mergeclusters@meta.data[WhichCells(mergeclusters, idents = c(3)), ]$ID10x = "3, A-stage2"
mergeclusters@meta.data[WhichCells(mergeclusters, idents = c(4)), ]$ID10x = "4, A-stages3,4"
mergeclusters@meta.data[WhichCells(mergeclusters, idents = c(0)), ]$ID10x = "0, A-stages6,7,5"
mergeclusters@meta.data[WhichCells(mergeclusters, idents = c(5)), ]$ID10x = "5, unknown"

cols.mergedclusters = getColors(DimPlot(mergeclusters))

for(var in names(seurat_obj_list)){
    obj = seurat_obj_list[[var]]
    seurat_obj_list[[var]]@meta.data$ID10x = "Unmapped"
    seurat_obj_list[[var]]@meta.data[intersect(colnames(seurat_obj_list[[var]]), colnames(mergeclusters)),]$ID10x = mergeclusters@meta.data[intersect(colnames(seurat_obj_list[[var]]), colnames(mergeclusters)),]$ID10x
    seurat_obj_list[[var]]@meta.data$ID10x = factor(seurat_obj_list[[var]]@meta.data$ID10x, levels = c("8, NSC-stages1,2", "7, RG-like",  "1, N-stages1,2" , "2, N-stage3", "6, N-stages4,5",  "10, N-stages5,6", "11, N-stage7", "9, A-stages1,5", "3, A-stage2", "4, A-stages3,4", "0, A-stages6,7,5", "5, unknown", "Unmapped"))

}


var = names(seurat_obj_list)[2]

identstokeep = c("8, NSC-stages1,2", "7, RG-like", "1, N-stages1,2", "2, N-stage3", "6, N-stages4,5", "10, N-stages5,6", "3, A-stage2", "4, A-stages3,4", "5, unknown", "Unmapped")
#identstokeep = c("8, NSC-stages1,2", "7, RG-like")
#identstokeep = c("1, N-stages1,2", "2, N-stage3", "6, N-stages4,5", "10, N-stages5,6")
#identstokeep = c("3, A-stage2", "4, A-stages3,4",)

cellstokeep = rownames(mergeclusters@meta.data[ which(mergeclusters@meta.data$ID10x %in% identstokeep),])


pdf(paste0("spatial_location_of_NSC_like_on_", var, ".pdf"),
width = 492,
height = 519)

DimPlot(seurat_obj_list[[var]], group.by = "ID10x", reduction = "coord", cols = c(cols.mergedclusters [[9]], cols.mergedclusters [[8]], "gray","gray","gray","gray","gray","gray","gray", "gray","gray","gray"), pt.size = 1.5, cells = intersect(colnames(seurat_obj_list[[var]]), union( rownames(seurat_obj_list[[var]]@meta.data[which(seurat_obj_list[[var]]@meta.data$ID10x == "Unmapped"),]), cellstokeep))) + ggtitle(var)+ theme(axis.text.y  = element_blank(),
axis.text.x  = element_blank(),
axis.title.y = element_blank(),
axis.title.x = element_blank(),
plot.title = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
panel.background = element_blank(),
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
axis.line = element_blank(),
panel.border = element_blank(),
axis.ticks.x=element_blank(),
axis.ticks.y=element_blank()
)

dev.off()




pdf(paste0("spatial_location_of_Neuronal_lineage_on_", var, ".pdf"),
width = 492,
height = 519)

DimPlot(seurat_obj_list[[var]], group.by = "ID10x", reduction = "coord", cols = c( "gray","gray", cols.mergedclusters [[2]], cols.mergedclusters [[3]],cols.mergedclusters [[7]],cols.mergedclusters [[11]], "gray","gray","gray", "gray","gray","gray"), pt.size = 1.5, cells = intersect(colnames(seurat_obj_list[[var]]), union( rownames(seurat_obj_list[[var]]@meta.data[which(seurat_obj_list[[var]]@meta.data$ID10x == "Unmapped"),]), cellstokeep))) + ggtitle(var)+ theme(axis.text.y  = element_blank(),
axis.text.x  = element_blank(),
axis.title.y = element_blank(),
axis.title.x = element_blank(),
plot.title = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
panel.background = element_blank(),
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
axis.line = element_blank(),
panel.border = element_blank(),
axis.ticks.x=element_blank(),
axis.ticks.y=element_blank()
)

dev.off()


pdf(paste0("spatial_location_of_Astrocytic_lineage_on_", var, ".pdf"),
width = 492,
height = 519)

DimPlot(seurat_obj_list[[var]], group.by = "ID10x", reduction = "coord", cols = c( "gray","gray",  "gray","gray","gray", "gray",cols.mergedclusters [[10]], cols.mergedclusters [[4]],cols.mergedclusters [[5]],cols.mergedclusters [[1]],"gray","gray"), pt.size = 1.5, cells = intersect(colnames(seurat_obj_list[[var]]), union( rownames(seurat_obj_list[[var]]@meta.data[which(seurat_obj_list[[var]]@meta.data$ID10x == "Unmapped"),]), cellstokeep))) + ggtitle(var)+ theme(axis.text.y  = element_blank(),
axis.text.x  = element_blank(),
axis.title.y = element_blank(),
axis.title.x = element_blank(),
plot.title = element_text(size=fontsize, face = "bold", family = "Times New Roman"),
panel.background = element_blank(),
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
axis.line = element_blank(),
panel.border = element_blank(),
axis.ticks.x=element_blank(),
axis.ticks.y=element_blank()
)

dev.off()


pdf("Marker_genes_in_spatial_data.pdf",
width = 1422,
height = 678)

spatialcells = rownames(meta_data_combined[which(meta_data_combined$sample == "Genexyz"),])
onlyspatialmergeddata = SubsetData(mergeclusters, cells = spatialcells)



FeaturePlot(onlyspatialmergeddata, features = c("Ascl1", "Aurkb", "Ube2c", "Mki67", "Pax6",
                                                "Sox4", "Neurog2", "Unc5d", "Calb2", "Syt1",
                                                "Aldoc", "Sned1", "Fgfr3", "Ogt", "Fam107a"
), ncol = 5, ) & theme( plot.title = element_text( face = "italic") )
dev.off()




## spatial mapping

spatial_mapping = read.csv("/Volumes/Transcend/Pascal/Figures/Table12_spatial_to_10x_mapping.csv", row.names = 1)
norm_vector = colSums(spatial_mapping)

for (i in 1:length(colnames(spatial_mapping))){
    spatial_mapping[[i]] = spatial_mapping[[i]]/norm_vector[i]*100
}

spatial_mapping_to_plot = spatial_mapping[1:12, 0:17]


## Figure5_K_reordered_without_5.pdf

heatmap.2(as.matrix(t(spatial_mapping_to_plot[c(9,3,4,1,2,6,10,8,7),])), dendrogram='none', Rowv=F, Colv=F,trace='none',
          density.info = 'none', key = TRUE, keysize = 3.5, xlab = NULL,col=viridis(50),
          margins = c(5,7), family = "Times", cexRow = 1.05, cexCol = 1.05, key.xlab = "% of representation", key.title = "",  lhei=c(1,2.8), lwid=c(1,1.8), key.par = list(cex=1))

## Figure5_K_reordered_with_5.pdf
heatmap.2(as.matrix(t(spatial_mapping_to_plot[c(9,3,4,1,2,6,10,8,7, 5),])), dendrogram='none', Rowv=F, Colv=F,trace='none',
          density.info = 'none', key = TRUE, keysize = 3.5, xlab = NULL,col=viridis(50),
          margins = c(5,7), family = "Times", cexRow = 1.05, cexCol = 1.05, key.xlab = "% of representation", key.title = "",  lhei=c(1,2.8), lwid=c(1,1.8), key.par = list(cex=1))
