#!/usr/bin/env Rscript

# ========================================================================================
# --- Explore how coverage varies by sample, sample type, mapped read count, etc.
# ========================================================================================

options(stringsAsFactors = FALSE)

# For testing, set to value other than -1.
to.read = -1

info.full = read.table("results/RADtags.info.bed",        nrows=to.read)
cvg.full  = read.table("results/RADtag.coverage.all.txt", nrows=to.read)

total.RADtags = nrow(info.full)

# ----------------------------------------------------------------------------------------
# --- Create dataset with entirely un-sequenced RADtags removed
# ----------------------------------------------------------------------------------------

no.hit.RADtags = which(rowSums(cvg.full[4:ncol(cvg.full)]) == 0)

info.full = info.full[-no.hit.RADtags,]
cvg.full  = cvg.full[-no.hit.RADtags,]

info = info.full
cvg  = cvg.full

# ----------------------------------------------------------------------------------------
# --- Bring in individual info
# ----------------------------------------------------------------------------------------

ind.info = read.csv("data/fecalRAD_individual_info.csv")

# Add sample names to individual coverage
# Warning: This assumes that the individuals are in the right order, 
# same as ls -v would produce

names(cvg) = c("chr", "start", "end", ind.info$NGS.ID)

# Make "chr" into a factor

cvg$chr = as.factor(cvg$chr)

# ----------------------------------------------------------------------------------------
# --- Parse RADtag info
# ----------------------------------------------------------------------------------------

rad.info = do.call(rbind, strsplit(info$V4, split='[;:]'))

rad.vars = rad.info[1,seq(from=1, to=ncol(rad.info), by=2)]
rad.info = data.frame(info[,1:3], rad.info[,seq(from=2, to=ncol(rad.info), by=2)])
names(rad.info) = c("chr", "start", "end", rad.vars)

rad.info[,4:ncol(rad.info)] = apply(rad.info[,4:ncol(rad.info)], 2, 
									function(x) as.numeric(x))

# Make "chr" into a factor
rad.info$chr = as.factor(rad.info$chr)

# ----------------------------------------------------------------------------------------
# --- For each individual, find number of RADtags with at least N reads
# ----------------------------------------------------------------------------------------

get.num.tags.gt.N = function (N) {
	return (colSums(cvg[,4:ncol(cvg)] > N))
}

hits.by.cov = do.call(rbind, lapply(1:50, get.num.tags.gt.N))
hits.by.cov.na = hits.by.cov
hits.by.cov.na[hits.by.cov.na == 0] = NA

pdf(file="reports/RADtag_count_gtNreads.pdf")
	matplot(hits.by.cov.na, type='l', 
		col=c('red', 'darkgoldenrod')[factor(ind.info$Sample.type)], 
		lty=as.numeric(factor(ind.info$Sample.type)), 
		ylab="Number of RADtags with more than N reads", 
		xlab="Number of reads (N)")
dev.off()

# ----------------------------------------------------------------------------------------
# --- Do more coverage curves, this time plotting against number of mapped reads
# ----------------------------------------------------------------------------------------

reads.mapped2 = ind.info$Reads.mapped ^ 2

pdf(file="reports/RADtag_count_curves_by_total_mapped.pdf")

	plot(hits.by.cov[1,] ~ ind.info$Reads.mapped, ylim = c(0,max(hits.by.cov[1,])),
		type='n', xlab="Reads Mapped", ylab="Number of RADtags with more than N reads")
	
	cols = vector("character", 20)
	cols[1]  = rgb(red=1, blue=0, green=0, alpha=0.5)
	cols[5]  = rgb(red=0, blue=1, green=0, alpha=0.5)
	cols[10] = rgb(red=0, blue=0, green=1, alpha=0.5)
	cols[20] = rgb(red=1, blue=1, green=0, alpha=0.5)
	
	add.cov.curve = function(N) {
	
		# Add black outline around blood
		points(hits.by.cov[N,] ~ ind.info$Reads.mapped,
			pch=1,
			col=c('black', NA)[factor(ind.info$Sample.type)])
	
		points(hits.by.cov[N,] ~ ind.info$Reads.mapped,
			pch=16, col=cols[N])
		fit = lm(hits.by.cov[N,] ~ ind.info$Reads.mapped + reads.mapped2)
		fit.l = cbind(ind.info$Reads.mapped, predict(fit))
		fit.l = fit.l[order(fit.l[,1]),]
	
		lines(fit.l, col = cols[N], lwd = 2)
	}
	
	lapply(c(1,5,10,20), add.cov.curve)

dev.off()

# ----------------------------------------------------------------------------------------
# --- Plot residuals by various variables
# ----------------------------------------------------------------------------------------

# Work in progress...

fit = lm(hits.by.cov[1,] ~ ind.info$Reads.mapped + reads.mapped2)
pairs(cbind(resid(fit), ind.info[,c(11:15,18)]), lower.panel=NULL)
