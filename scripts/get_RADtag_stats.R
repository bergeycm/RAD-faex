options(stringsAsFactors = FALSE)

library(reshape2)

info.full = read.table("results/RADtags.info.bed")
cvg.full  = read.table("results/RADtag.coverage.all.txt")

total.RADtags = nrow(info.full)


# ----------------------------------------------------------------------------------------
# --- Create dataset with entirely un-sequenced RADtags removed
# ----------------------------------------------------------------------------------------

no.hit.RADtags = which(rowSums(cvg.full[4:ncol(cvg.full)]) == 0)

info.full = info.full[-no.hit.RADtags,]
cvg.full  = cvg.full[-no.hit.RADtags,]

info = info.full
cvg = cvg.full

# ----------------------------------------------------------------------------------------
# --- Bring in individual info
# ----------------------------------------------------------------------------------------

ind.info = read.csv("data/individual_info.csv")


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


# Compute deviation from median
rad.info.all$len.deviation = abs(rad.info.all$length - med.len)
rad.info.blood$len.deviation = rad.info.all$len.deviation
rad.info.feces$len.deviation = rad.info.all$len.deviation

# Final median, modal length for RADtags that have coverage

hit.lens = rad.info.all[rad.info.all$rad.mean.cov > 7,]$length
hit.lens.uniq = unique(hit.lens)
mode.len = hit.lens.uniq[which.max(tabulate(match(hit.lens, hit.lens.uniq)))]
med.len  = median(rad.info.all[rad.info.all$rad.mean.cov > 7,]$length)

# Add a "type" column for ggplot2

rad.info.blood$type = 'blood'
rad.info.feces$type = 'feces'

rad.info.both = rbind(rad.info.blood,rad.info.feces)

# Rearrange so it's near its friends
len.col = which(names(rad.info.all) == "length")
new.col.ord = c(1:len.col, ncol(rad.info.all), (len.col+1):(ncol(rad.info.all)-1))

rad.info.all = rad.info.all[,new.col.ord]
rad.info.blood = rad.info.blood[,new.col.ord]
rad.info.feces = rad.info.feces[,new.col.ord]


# Recompute dnorm based on median
rad.info.all$len_dnorm = rad.info.blood$len_dnorm = rad.info.feces$len_dnorm = 
	dnorm(rad.info.all$length, mean = med.len, sd = 25)

rad.info.lens = rad.info
rad.info.lens$len.deviation = rad.info.all$len.deviation
rad.info.lens$len_dnorm     = rad.info.all$len_dnorm

# Merge info and coverage into one massive data.frame
info.cvg = merge(rad.info.lens, cvg, by=c("chr", "start", "end"))

info.cvg.m = melt(info.cvg, id=names(rad.info.lens))

# Fix names
names(info.cvg.m)[which(names(info.cvg.m) == "value")]    = "num.reads"
names(info.cvg.m)[which(names(info.cvg.m) == "variable")] = "NGS.ID"

info.cvg.m.ind = merge(info.cvg.m, ind.info, by="NGS.ID")

# Factorify relevant columns

info.cvg.m.ind = within(info.cvg.m.ind,{
	NGS.ID = factor(NGS.ID)
	Pool.ID = factor(Pool.ID)
	Sample.ID = factor(Sample.ID)
	Individual.ID = factor(Individual.ID)
	Sample.type = factor(Sample.type)
})

# Save as R Data file because the file is pretty big

# write.csv(rad.info.both,file='results/radtag_stats.csv',row.names=FALSE,quote=FALSE)
save(list=c('rad.info.both','info.cvg.m.ind'),file='results/radtag_stats.RData')






