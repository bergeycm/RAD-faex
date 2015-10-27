#!/usr/bin/env Rscript

# ========================================================================================
# --- Model RADtag coverage (to explain inter-RADtag variation)
# ========================================================================================

options(stringsAsFactors = FALSE)

if(!require(reshape)){
	install.packages("reshape", dependencies=TRUE, repos='http://cran.rstudio.com/')
	library(reshape)
}
if(!require(lme4)){
	install.packages("lme4", dependencies=TRUE, repos='http://cran.rstudio.com/')
	library(lme4)
}
if(!require(car)){
	install.packages("car", dependencies=TRUE, repos='http://cran.rstudio.com/')
	library(car)
}
if(!require(R2admb)){
	install.packages("R2admb", dependencies=TRUE, repos='http://cran.rstudio.com/')
	library(R2admb)
}
if(!require(coda)){
	install.packages("coda", dependencies=TRUE, repos='http://cran.rstudio.com/')
	library(coda)
}
if(!require(glmmADMB)){
	install.packages("glmmADMB", repos="http://glmmadmb.r-forge.r-project.org/repos")
	library(glmmADMB)
}

# Get subset size from user
args = commandArgs(trailingOnly = TRUE)
subset.size = as.numeric(args[1])
write(paste("Processing subset of", subset.size, "lines") , stderr())

info.full = read.table("results/RADtags.info.bed")
cvg.full  = read.table("results/RADtag.coverage.all.txt")

total.RADtags = nrow(info.full)

# ----------------------------------------------------------------------------------------
# --- Create dataset with entirely un-sequenced RADtags removed
# ----------------------------------------------------------------------------------------

no.hit.RADtags = which(rowSums(cvg.full[4:ncol(cvg.full)]) == 0)

info.full = info.full[-no.hit.RADtags,]
cvg.full  = cvg.full[-no.hit.RADtags,]

# ----------------------------------------------------------------------------------------
# --- Randomly subset RADtags
# ----------------------------------------------------------------------------------------

subset.ids = sample(1:nrow(info.full), subset.size, replace=FALSE)

info = info.full[subset.ids,]
cvg  = cvg.full[subset.ids,]

rm(list=c("info.full", "cvg.full"))

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
# --- Normalize coverage of RADtags for each individual
# ----------------------------------------------------------------------------------------

ind.mean.cov = colMeans(cvg[,4:ncol(cvg)])

# Normalize coverage by dividing RADtag's coverage by individual's mean coverage
normalized.cov = data.frame(cbind(	cvg[,1:3],
									sweep(cvg[,4:ncol(cvg)], 2, ind.mean.cov, `/`)))

# ----------------------------------------------------------------------------------------
# --- Compute stats for each RADtag
# ----------------------------------------------------------------------------------------

# --- First overall, without regard to sample type

rad.mean.cov = rowMeans(cvg[,4:ncol(cvg)])
rad.med.cov  = apply(cvg[,4:ncol(cvg)], 1, median)

rad.norm.mean.cov = rowMeans(normalized.cov[,4:ncol(normalized.cov)])
rad.norm.med.cov  = apply(normalized.cov[,4:ncol(normalized.cov)], 1, median)

# --- Split by sample type to model them separately

# --- First blood.
blood.ind.idx = c(1:3, (3 + which(ind.info$Sample.type == "blood")))
cvg.blood = cvg[,blood.ind.idx]
normalized.cov.blood = normalized.cov[,blood.ind.idx]

rad.mean.cov.blood = rowMeans(cvg[,4:ncol(cvg.blood)])
rad.med.cov.blood  = apply(cvg.blood[,4:ncol(cvg.blood)], 1, median)

rad.norm.mean.cov.blood = rowMeans(normalized.cov.blood[,4:ncol(normalized.cov.blood)])
rad.norm.med.cov.blood  = apply(normalized.cov.blood[,4:ncol(normalized.cov.blood)], 1, 
								median)

# --- Now poo.
feces.ind.idx = c(1:3, (3 + which(ind.info$Sample.type == "feces")))
cvg.feces = cvg[,feces.ind.idx]
normalized.cov.feces = normalized.cov[,feces.ind.idx]

rad.mean.cov.feces = rowMeans(cvg[,4:ncol(cvg.feces)])
rad.med.cov.feces  = apply(cvg.feces[,4:ncol(cvg.feces)], 1, median)

rad.norm.mean.cov.feces = rowMeans(normalized.cov.feces[,4:ncol(normalized.cov.feces)])
rad.norm.med.cov.feces  = apply(normalized.cov.feces[,4:ncol(normalized.cov.feces)], 1, 
								median)

# ----------------------------------------------------------------------------------------
# --- Combine RADtag info and coverage summary stats
# ----------------------------------------------------------------------------------------

rad.info.all = cbind(	rad.info, 
						rad.mean.cov, rad.med.cov, 
						rad.norm.mean.cov, rad.norm.med.cov)

rad.info.blood = cbind(	rad.info, 
						rad.mean.cov.blood, rad.med.cov.blood, 
						rad.norm.mean.cov.blood, rad.norm.med.cov.blood)

rad.info.feces = cbind(	rad.info, 
						rad.mean.cov.feces, rad.med.cov.feces, 
						rad.norm.mean.cov.feces, rad.norm.med.cov.feces)

names(rad.info.blood) = gsub("\\.blood", "", names(rad.info.blood))
names(rad.info.feces) = gsub("\\.feces", "", names(rad.info.feces))

# ----------------------------------------------------------------------------------------
# --- Improve length variable
# ----------------------------------------------------------------------------------------

# Final median, modal length for RADtags that have coverage

hit.lens = rad.info.all[rad.info.all$rad.mean.cov > 7,]$length
hit.lens.uniq = unique(hit.lens)
mode.len = hit.lens.uniq[which.max(tabulate(match(hit.lens, hit.lens.uniq)))]
med.len  = median(rad.info.all[rad.info.all$rad.mean.cov > 7,]$length)

pdf("reports/RADtag_length_distribution.pdf")
par(mfrow=c(3,1))
hist(rad.info.all[rad.info.all$rad.mean.cov > 7,]$length, breaks=500, xlim=c(0,500), 
	main="All Samples", xlab="Length of RADtags with average read count > 7")
abline(v=med.len, lwd=2, col='red')
hist(rad.info.blood[rad.info.blood$rad.mean.cov > 7,]$length, breaks=500, xlim=c(0,500),
	main="Blood Samples", xlab="Length of RADtags with average read count > 7")
abline(v=med.len, lwd=2, col='red')
hist(rad.info.feces[rad.info.feces$rad.mean.cov > 7,]$length, breaks=500, xlim=c(0,500),
	main="Fecal Samples", xlab="Length of RADtags with average read count > 7")
abline(v=med.len, lwd=2, col='red')
par(mfrow=c(1,1))
dev.off()

# NOTE: How many bp do the adapters actually add to the reaction?
#abline(v=300 - 65, lwd=2, col='blue', lty=2)

# Compute deviation from median
rad.info.all$len.deviation = abs(rad.info.all$length - med.len)
rad.info.blood$len.deviation = rad.info.all$len.deviation
rad.info.feces$len.deviation = rad.info.all$len.deviation

# Rearrange so it's near its friends
len.col = which(names(rad.info.all) == "length")
new.col.ord = c(1:len.col, ncol(rad.info.all), (len.col+1):(ncol(rad.info.all)-1))

rad.info.all = rad.info.all[,new.col.ord]
rad.info.blood = rad.info.blood[,new.col.ord]
rad.info.feces = rad.info.feces[,new.col.ord]

# Recompute dnorm based on median
rad.info.all$len_dnorm = rad.info.blood$len_dnorm = rad.info.feces$len_dnorm = 
	dnorm(rad.info.all$length, mean = med.len, sd = 25)

# ----------------------------------------------------------------------------------------
# --- Visualize relationships
# ----------------------------------------------------------------------------------------

# Currently skipped
# pairs(rad.info.all[(rad.info.all$rad.mean.cov > 7) & 
#						(rad.info.all$rad.mean.cov < 50),
#						7:15])
# pairs(rad.info.blood[(rad.info.blood$rad.mean.cov > 7) & 
#						(rad.info.blood$rad.mean.cov < 50),
#						7:15])
# pairs(rad.info.all[(rad.info.all$rad.mean.cov > 7) & 
#						(rad.info.all$rad.mean.cov < 50),
#						7:15])

# ----------------------------------------------------------------------------------------
# --- Do simple multiple regression with everying thrown in
# ----------------------------------------------------------------------------------------

lm.all = lm(rad.mean.cov ~ length + len.deviation + len_dnorm + gc_perc + N_count + 
				gc_perc_5000 + N_count_5000 + 
				CpG_dist + CpG_ct + CpG_is_dist + CpG_5000, data=rad.info.all)

lm.blood = lm(rad.mean.cov ~ length + len.deviation + len_dnorm + gc_perc + N_count + 
				gc_perc_5000 + N_count_5000 + 
				CpG_dist + CpG_ct + CpG_is_dist + CpG_5000, data=rad.info.blood)

lm.feces = lm(rad.mean.cov ~ length + len.deviation + len_dnorm + gc_perc + N_count + 
				gc_perc_5000 + N_count_5000 + 
				CpG_dist + CpG_ct + CpG_is_dist + CpG_5000, data=rad.info.feces)

sink("reports/RADtag.lm.blood.summary.txt")
	summary(lm.blood)
sink()
sink("reports/RADtag.lm.feces.summary.txt")
	summary(lm.feces)
sink()

# ----------------------------------------------------------------------------------------
# --- Restructure to include individual-level info
# ----------------------------------------------------------------------------------------

# Copy in length deviation and len_dnorm
rad.info.lens = rad.info
rad.info.lens$len.deviation = rad.info.all$len.deviation
rad.info.lens$len_dnorm     = rad.info.all$len_dnorm

# Merge info and coverage into one massive data.frame
info.cvg = merge(rad.info.lens, cvg, by=c("chr", "start", "end"))

info.cvg.m = melt(info.cvg, id=names(rad.info.lens))

# Fix names
names(info.cvg.m)[which(names(info.cvg.m) == "value")]    = "num.reads"
names(info.cvg.m)[which(names(info.cvg.m) == "variable")] = "NGS.ID"

# ----------------------------------------------------------------------------------------
# --- Bring in relevant info about individuals
# ----------------------------------------------------------------------------------------

info.cvg.m.ind = merge(info.cvg.m, ind.info, by="NGS.ID")

# Factorify relevant columns

info.cvg.m.ind = within(info.cvg.m.ind,{
	NGS.ID = factor(NGS.ID)
	Pool.ID = factor(Pool.ID)
	Sample.ID = factor(Sample.ID)
	Individual.ID = factor(Individual.ID)
	Sample.type = factor(Sample.type)
})

# ----------------------------------------------------------------------------------------
# --- Test for Poisson fit (No, it does not fit)
# ----------------------------------------------------------------------------------------

# See:
# http://ase.tufts.edu/gsc/gradresources/guidetomixedmodelsinr/mixed%20model%20guide.html

# poisson.fit = fitdistr(info.cvg.m.ind$num.reads, "Poisson")
# qqp(info.cvg.m.ind$num.reads, "pois", poisson.fit$estimate)

num.reads = info.cvg.m.ind$num.reads + 0.000001

#	nbinom.fit = fitdistr(num.reads, "Negative Binomial")
#	pdf(file="results/qqp_nbinom_test.pdf")
#		qqp(num.reads, "nbinom", 
#			size = nbinom.fit$estimate[[1]], 
#			mu = nbinom.fit$estimate[[2]])
#	dev.off()

# ----------------------------------------------------------------------------------------
# --- Do multiple regression with individual-level info
# ----------------------------------------------------------------------------------------

info.cvg.m.ind.subset = info.cvg.m.ind[,c(	'num.reads','len_dnorm','gc_perc',
											'CpG_5000','Sample.type','NGS.ID')]

# Center predictors
info.cvg.m.ind.subset = within(info.cvg.m.ind.subset,{
	len_dnorm = len_dnorm - mean(len_dnorm)
#	len.deviation = len.deviation - mean(len.deviation)
	gc_perc = gc_perc - mean(gc_perc)
#	CpG_is_dist = CpG_is_dist - mean(CpG_is_dist)
	CpG_5000 = CpG_5000 - mean(CpG_5000)
})

# Scale predictors
info.cvg.m.ind.subset = within(info.cvg.m.ind.subset,{
	len_dnorm = len_dnorm * 1000			## Arbitrary rescaling
#	len.deviation = len.deviation / 1000	## kilobases
	gc_perc = gc_perc * 100					## Percentages instead of fraction
#	CpG_is_dist = CpG_is_dist / 1000 		## kilobases
	CpG_5000 = CpG_5000 /100				## CpG / 100 bases in the greater +- 5kb area
})

# Create temporary directory to hold files
tmp.dir = tempfile(tmpdir = "./tmp", pattern=paste0("subset_",subset.size,"_"))

ptm = proc.time()
lm.ind = glmmadmb(num.reads ~ len_dnorm + (gc_perc + CpG_5000) * Sample.type + (1|NGS.ID),
					data=info.cvg.m.ind.subset, 
					extra.args="-ndi 200000",
					admb.opts=admbControl(	shess=FALSE, noinit=FALSE, impSamp=200,
											maxfn=1000, imaxfn=500, maxph=5),
					save.dir=tmp.dir,
					family="nbinom1", 
					zeroInflation=TRUE) ## Try also with zeroInflation=FALSE
elapsed = (proc.time() - ptm)['elapsed']
write(paste(elapsed,'seconds passed\n'), stderr())

sink("reports/RADtag.lm.indiv.summary.txt")
	summary(lm.ind)
sink()

# Save glmmadmb results
# For now, includes the number of rows read in filename
save(list="lm.ind", file=paste0("results/lm.ind.indiv.", subset.size, ".Rdata"))

# Remove temporary directory
#unlink(tmp.dir)
