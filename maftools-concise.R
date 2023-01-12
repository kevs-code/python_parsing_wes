## ----results='hide', message=FALSE---------------------------------------
require(maftools)
setwd(".")
mrdmm.maf = file.path(".", "together.maf")
mrdmm = read.maf(maf = mrdmm.maf)#, clinicalData = mrdmm.clin)
## ------------------------------------------------------------------------
#Typing mrdmm shows basic summary of MAF file.
mrdmm
#Shows sample summry.
getSampleSummary(mrdmm)
#Shows gene summary.
getGeneSummary(mrdmm)
#Shows all fields in MAF
getFields(mrdmm)
#Writes maf summary to an output file with basename mrdmm.
write.mafSummary(maf = mrdmm, basename = 'mrdmm')
pdf(file="mrdmm_plot1.pdf",width=9.1,height=6.7)
plotmafSummary(maf = mrdmm, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
dev.off()
tiff(filename="mrdmm_plot1.tiff",width=9.1,height=6.7,units="in",res=300,compression="lzw")
plotmafSummary(maf = mrdmm, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
dev.off()
pdf(file="mrdmm_plot2.pdf",width=9.1,height=6.7)
oncoplot(maf = mrdmm, top = 10, fontSize = 12)
dev.off()
tiff(filename="mrdmm_plot2.tiff",width=9.1,height=6.7,units="in",res=300,compression="lzw")
oncoplot(maf = mrdmm, top = 10, fontSize = 12)
dev.off()
pdf(file="mrdmm_plot3.pdf",width=9.1,height=6.7)
oncostrip(maf = mrdmm, genes = c('EGFR','IGL','IRF4','LTB','TRAF','CSMD3','CYLD','FAT1','FGFR3','KDM6A','RB1','ATM','CCND1','DIS3','BRAF','TP53','FAM46C','NRAS','KRAS'))
dev.off()
tiff(filename="mrdmm_plot3.tiff",width=9.1,height=6.7,units="in",res=300,compression="lzw")
oncostrip(maf = mrdmm, genes = c('EGFR','IGL','IRF4','LTB','TRAF','CSMD3','CYLD','FAT1','FGFR3','KDM6A','RB1','ATM','CCND1','DIS3','BRAF','TP53','FAM46C','NRAS','KRAS'))
dev.off()
#tiff(filename="mrdmm_plot4.tiff",width=9.1,height=6.7,units="in",res=300,compression="lzw")
#mrdmm.mutsig <- file.path("/home/kevin/mrd_wes/processed_data/samples/28-04-2018_ensemble_vcf/filtered/mutsigCV", "2018-05-08.sig_genes.txt")#, package = "maftools")
#oncoplot(maf = mrdmm, colors = col, mutsig = mrdmm.mutsig, mutsigQval = 0.01, sortByAnnotation = FALSE)#, annotationColor = fabcolors)
#dev.off()
