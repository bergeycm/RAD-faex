#!/usr/bin/Rscript

# ========================================================================================
# --- Model RADtag coverage (to explain inter-RADtag variation)
# ========================================================================================

options(stringsAsFactors = FALSE)

info = read.table("results/RADtags.info.bed")
cvg  = read.table("results/RADtag.coverage.all.txt")

# ----------------------------------------------------------------------------------------
# --- Bring in individual info
# ----------------------------------------------------------------------------------------

ind.info = read.csv("data/fecalRAD_individual_info.csv")

# Add sample names to individual coverage
# Warning: This assumes that the individuals are in the right order, 
# same as ls -v would produce

names(cvg) = c("chr", "start", "end", ind.info$NGS.ID)

# ----------------------------------------------------------------------------------------
# --- Parse RADtag info
# ----------------------------------------------------------------------------------------

rad.info = do.call(rbind, strsplit(info$V4, split='[;:]'))

rad.vars = rad.info[1,seq(from=1, to=ncol(rad.info), by=2)]
rad.info = data.frame(info[,1:3], rad.info[,seq(from=2, to=ncol(rad.info), by=2)])
names(rad.info) = c("chr", "start", "end", rad.vars)

rad.info[,4:ncol(rad.info)] = apply(rad.info[,4:ncol(rad.info)], 2, 
									function(x) as.numeric(x))

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

pairs(rad.info.all[(rad.info.all$rad.mean.cov > 7) & 
						(rad.info.all$rad.mean.cov < 50),
						7:15])
pairs(rad.info.blood[(rad.info.blood$rad.mean.cov > 7) & 
						(rad.info.blood$rad.mean.cov < 50),
						7:15])
pairs(rad.info.all[(rad.info.all$rad.mean.cov > 7) & 
						(rad.info.all$rad.mean.cov < 50),
						7:15])

# ----------------------------------------------------------------------------------------
# --- Do multiple regression
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
