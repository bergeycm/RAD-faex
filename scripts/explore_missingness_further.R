#!/usr/bin/env Rscript

ind.info = read.csv("data/fecalRAD_individual_info.csv")

# -------------------------------------------------------------------------------------- #
# --- Make MDS plot of IBS and IBM
# -------------------------------------------------------------------------------------- #

# === Creating MDS plot of IBS and IBM ================================================

mds = read.table("results/all.cleaned.LDpruned.mds", header=TRUE)
ibm = read.table("results/all.cleaned.LDpruned.missing.mds", header=TRUE)

mds$FID = gsub(".PE", "", mds$FID)

names(ibm) = paste0(names(ibm), ".missing")
mds = cbind(mds, ibm$C1.missing, ibm$C2.missing)
names(mds) = gsub("ibm\\$", "", names(mds))

mds.ind = merge(ind.info, mds, by.x="NGS.ID", by.y="FID")

pal = brewer.pal(length(unique(mds.ind$Individual.ID)), "Set3")

pdf("results/mds_ibs.pdf")
	plot(mds.ind$C1 ~ mds.ind$C2, 
		pch=c(16,1)[factor(mds.ind$Sample.type)], 
		col=pal[factor(mds.ind$Individual.ID)], 
		cex=1.5)
	#identify(mds.ind$C1 ~ mds.ind$C2, label=mds.ind$Sample.ID, cex=0.5, pos=4)
dev.off()

pdf("results/mds_ibm.pdf")
	plot(mds.ind$C1.missing ~ mds.ind$C2.missing,
		pch=c(16,1)[factor(mds.ind$Sample.type)], 
		col=pal[factor(mds.ind$Individual.ID)],
		cex=1.5)
	#identify(mds.ind$C1.missing ~ mds.ind$C2.missing, label=mds.ind$Sample.ID, cex=0.5, pos=4)
dev.off()

# ----------------------------------------------------------------------------------------
# --- Parse association test results
# ----------------------------------------------------------------------------------------

# LD and missingness pruned dataset
assoc = read.table("results/all.cleaned.LDpruned.assoc", header=TRUE)
table(assoc[which(assoc$P < 0.05),]$A2)
nrow(assoc[which(assoc$P < 0.05),])
#ggplot(assoc, aes(x=P, color=factor(A2), group=factor(A2))) + geom_density()

# All SNPs
assoc = read.table("results/all.assoc", header=TRUE)
table(assoc[which(assoc$P < 0.05),]$A2)
nrow(assoc[which(assoc$P < 0.05),])
#ggplot(assoc, aes(x=P, color=factor(A2), group=factor(A2))) + geom_density()

# ----------------------------------------------------------------------------------------
# --- Parse missingness chi-sq results
# ----------------------------------------------------------------------------------------

# All SNPs
miss = read.table("results/all.missing", header=TRUE)
sig.miss = miss[miss$P < 0.05 & miss$F_MISS_A > miss$F_MISS_U,]
nrow(sig.miss)
# Output BED file of locations where allele missingness is significantly different
sig.miss.bed = do.call(rbind, strsplit(sig.miss$SNP, split=":"))
sig.miss.bed = cbind(sig.miss.bed, as.numeric(sig.miss.bed[,2]) + 1)
write.table(sig.miss.bed, "results/significantly.differently.missing.all.bed", 
	quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

# LD and missingness pruned dataset
miss = read.table("results/all.cleaned.LDpruned.missing", header=TRUE)
sig.miss = miss[miss$P < 0.05 & miss$F_MISS_A > miss$F_MISS_U,]
nrow(sig.miss)
# Output BED file of locations where allele missingness is significantly different
sig.miss.bed = do.call(rbind, strsplit(sig.miss$SNP, split=":"))
sig.miss.bed = cbind(sig.miss.bed, as.numeric(sig.miss.bed[,2]) + 1)
write.table(sig.miss.bed, "results/significantly.differently.missing.LDpruned.bed", 
	quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

# ----------------------------------------------------------------------------------------
# --- Analyze distance between samples
# ----------------------------------------------------------------------------------------

# All SNPs
dist = read.table("results/all.dist")
inds = read.table("results/all.nosex")

# or 

# LD and missingness pruned dataset
dist = read.table("results/all.cleaned.LDpruned.dist")
inds = read.table("results/all.cleaned.LDpruned.nosex")

# --- 

inds$V1 = gsub(".PE", "", inds$V1)
inds = merge(inds, ind.info, by.x="V1", by.y="NGS.ID")

names(dist) = inds$V1
inds.dist = cbind(inds, dist)

inds.dist.m = melt(inds.dist, id=names(inds.dist)[1:ncol(inds)])

tx.inds = grep("Tx", unique(inds$Individual.ID), value=TRUE)
tx.inds = tx.inds[tx.inds != "Tx07"]

feces.to.skip = c("Tx05H2", "Tx02HH")
feces.to.skip.ids = unique(inds.dist.m[inds.dist.m$Sample.ID %in% feces.to.skip,1])

pull.out.dists = function (this.tx.ind) {

	this.inds.ids = inds[inds$Individual.ID == this.tx.ind,]$V1

	this.inds.dist.m = inds.dist.m[which(inds.dist.m$V1 %in% this.inds.ids & inds.dist.m$variable %in% this.inds.ids),]

	blood.ids = unique(this.inds.dist.m[this.inds.dist.m$Sample.type == "blood",]$V1)
	feces.ids = unique(this.inds.dist.m[this.inds.dist.m$Sample.type == "feces",]$V1)
	
	# Remove double enrichments from comparison
	feces.ids = feces.ids[!feces.ids %in% feces.to.skip.ids]	

	blood.to.poop = this.inds.dist.m[which(this.inds.dist.m$V1 %in% blood.ids & this.inds.dist.m$variable %in% feces.ids),]
	poop.to.poop  = this.inds.dist.m[which(this.inds.dist.m$V1 %in% feces.ids & this.inds.dist.m$variable %in% feces.ids),]
	
	blood.to.poop = blood.to.poop[blood.to.poop$value != 0,]
	
	poop.to.poop = poop.to.poop[poop.to.poop$value != 0,]
	poop.to.poop = poop.to.poop[order(poop.to.poop$value),]
	poop.to.poop = poop.to.poop[seq(from=1, to=nrow(poop.to.poop), by=2),]

	return(list(blood.to.poop$value, poop.to.poop$value))
}

dist.list = lapply(tx.inds, pull.out.dists)

blood.poop.dists = do.call(c, do.call(rbind, lapply(tx.inds, pull.out.dists))[,1])
poop.poop.dists  = do.call(c, do.call(rbind, lapply(tx.inds, pull.out.dists))[,2])

boxplot(c(blood.poop.dists, poop.poop.dists) ~ 
	c(rep("Blood to feces distance", length(blood.poop.dists)), 
	rep("Feces to feces distance", length(poop.poop.dists))))

wilcox.test(blood.poop.dists, poop.poop.dists)

