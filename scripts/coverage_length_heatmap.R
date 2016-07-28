#!/usr/bin/env Rscript

# ----------------------------------------------------------------------------------------
# --- Plot heatmap of coverage and RADtag length
# ----------------------------------------------------------------------------------------

library(ggplot2)
library(reshape)

cov.beds = list.files(path="results/", pattern="*.gt1.cov.bed")

cov.heatmap = function (this.bed) {

	cov.gt1 = read.table(paste0("results/", this.bed))
	
	names(cov.gt1) = c("chr", "start", "end", "info", "read.count")
	cov.gt1$len = cov.gt1$end - cov.gt1$start
	
	# Q&D plot
	plot(jitter(cov.gt1$read.count) ~ cov.gt1$len, xlim=c(0,1000))
	abline(v=c(275 - 65, 325 - 65), col='red')
	
	breaks = seq(from=0, to=max(cov.gt1$len)+10, by=50)
	bins = cut(cov.gt1$len, c(-Inf,breaks,Inf))
	cov.gt1$len.bin.l = breaks[as.numeric(bins) - 1]
	
	cov.gt1.t = table(cov.gt1$read.count, cov.gt1$len.bin.l)
	cov.gt1.m = melt(cov.gt1.t)
	names(cov.gt1.m) = c("read.count", "len.bin", "value")
	
	cov.gt1.m.sm = cov.gt1.m[cov.gt1.m$len.bin < 1001 & cov.gt1.m$read.count < 100,]
	
	p = ggplot(cov.gt1.m.sm, aes(len.bin, read.count)) + 
		geom_tile(aes(fill = value), color = "white") + 
		scale_fill_gradient(low = "steelblue", high = "red")
	
	ggsave(filename=paste0("reports/", gsub(".bed", ".pdf", this.bed)), plot=p)
}

lapply(cov.beds, cov.heatmap)
