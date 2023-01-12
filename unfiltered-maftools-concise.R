## ----results='hide', message=FALSE---------------------------------------
require(maftools)
setwd(".")
unfiltered_mrdmm.maf = file.path(".", "unfilter_together.maf")
unfiltered_mrdmm = read.maf(maf = unfiltered_mrdmm.maf)#, clinicalData = unfiltered_mrdmm.clin)
## ------------------------------------------------------------------------
#Typing unfiltered_mrdmm shows basic summary of MAF file.
unfiltered_mrdmm
#Shows sample summry.
getSampleSummary(unfiltered_mrdmm)
#Shows gene summary.
getGeneSummary(unfiltered_mrdmm)
#Shows all fields in MAF
getFields(unfiltered_mrdmm)
#Writes maf summary to an output file with basename unfiltered_mrdmm.
write.mafSummary(maf = unfiltered_mrdmm, basename = 'unfiltered_mrdmm')
pdf(file="unfiltered_mrdmm_plot1.pdf",width=9.1,height=6.7)
plotmafSummary(maf = unfiltered_mrdmm, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
dev.off()
tiff(filename="unfiltered_mrdmm_plot1.tiff",width=9.1,height=6.7,units="in",res=300,compression="lzw")
plotmafSummary(maf = unfiltered_mrdmm, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
dev.off()
pdf(file="unfiltered_mrdmm_plot2.pdf",width=9.1,height=6.7)
oncoplot(maf = unfiltered_mrdmm, top = 10, fontSize = 12)
dev.off()
tiff(filename="unfiltered_mrdmm_plot2.tiff",width=9.1,height=6.7,units="in",res=300,compression="lzw")
oncoplot(maf = unfiltered_mrdmm, top = 10, fontSize = 12)
dev.off()
pdf(file="unfiltered_mrdmm_plot3.pdf",width=9.1,height=6.7)
oncostrip(maf = unfiltered_mrdmm, genes = c('EGFR','IGL','IRF4','LTB','TRAF','CSMD3','CYLD','FAT1','FGFR3','KDM6A','RB1','ATM','CCND1','DIS3','BRAF','TP53','FAM46C','NRAS','KRAS'))
dev.off()
tiff(filename="unfiltered_mrdmm_plot3.tiff",width=9.1,height=6.7,units="in",res=300,compression="lzw")
oncostrip(maf = unfiltered_mrdmm, genes = c('EGFR','IGL','IRF4','LTB','TRAF','CSMD3','CYLD','FAT1','FGFR3','KDM6A','RB1','ATM','CCND1','DIS3','BRAF','TP53','FAM46C','NRAS','KRAS'))
dev.off()
#tiff(filename="unfiltered_mrdmm_plot4.tiff",width=9.1,height=6.7,units="in",res=300,compression="lzw")
#unfiltered_mrdmm.mutsig <- file.path("/home/kevin/mrd_wes/processed_data/samples/28-04-2018_ensemble_vcf/filtered/mutsigCV", "2018-05-08.sig_genes.txt")#, package = "maftools")
#oncoplot(maf = unfiltered_mrdmm, colors = col, mutsig = unfiltered_mrdmm.mutsig, mutsigQval = 0.01, sortByAnnotation = FALSE)#, annotationColor = fabcolors)
#dev.off()
