#!/usr/bin/env Rscript
################################################################################################
# Produces rainfall plots given a "simple text file" format as input. 
# This input file may have 1 or more samples included, and it also works on both SNVs and INDELs,
# however, you must separate the SNVs from the INDELs and place in separate files. The script 
# will ultimately produce a rainfall plot for each sample. 
#
#  
# Requires: regioneR  (https://bioconductor.org/packages/release/bioc/html/regioneR.html)
#			karyoploteR (https://bioconductor.org/packages/release/bioc/html/karyoploteR.html)
#
#
#
# This script was adpated from the following tutorial:
#	https://bernatgel.github.io/karyoploter_tutorial/Examples/Rainfall/Rainfall.html
################################################################################################

args <- commandArgs(trailingOnly = TRUE)

somatic.mutations <- read.table(args[1], header=FALSE, sep="\t", stringsAsFactors=FALSE)
somatic.mutations <- setNames(somatic.mutations, c("test", "sample", "ID", "Genome", "type", "chr", "start","end","ref","mut","subtype"))
library(regioneR)
somatic.mutations <- split(somatic.mutations, somatic.mutations$sample)



for (sample in ls(somatic.mutations)){
	#print(sample)
	sm <- somatic.mutations[[sample]]
	sm.gr <- toGRanges(sm[,c("chr", "start", "end", "type", "ref", "mut")])
	seqlevelsStyle(sm.gr) <- "UCSC"
	sm.gr

	library(karyoploteR)
	variant.colors <- getVariantsColors(sm.gr$ref, sm.gr$mut)
	pp <- getDefaultPlotParams(plot.type = 4)
	pp$data1inmargin <- 0
	pp$bottommargin <- 20
	
	file = paste(sample, "rainfall_plot.pdf", sep="_")
	pdf(file, width=15, height=10)
	
	kp <- plotKaryotype(plot.type=4, ideogram.plotter = NULL, labels.plotter = NULL, plot.params = pp)
	kpAddCytobandsAsLine(kp)
	kpAddChromosomeNames(kp, srt=45)
	kpAddMainTitle(kp, main=sample, cex=1.2)
	try(
		{kpPlotRainfall(kp, data = sm.gr, col=variant.colors, r0=0, r1=0.7)}, silent=T
	)
	kpAxis(kp, ymax = 7, tick.pos = 1:7, r0=0, r1=0.7)
	kpAddLabels(kp, labels = c("Distance between mutations (log10)"), srt=90, pos=1, label.margin = 0.04, r0=0, r1=0.7)
	kpPlotDensity(kp, data = sm.gr, r0=0.72, r1=1)
	kpAddLabels(kp, labels = c("Density"), srt=90, pos=1, label.margin = 0.04, r0=0.71, r1=1)


	dev.off()

}
