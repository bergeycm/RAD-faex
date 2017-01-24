#!/usr/bin/env Rscript

library(ggplot2)
library(reshape2)

ind.info = read.csv("../data/individual_info.csv")

sink("results/missingness_info.txt")

# -------------------------------------------------------------------------------------- #
# --- Make MDS plot of IBS and IBM
# -------------------------------------------------------------------------------------- #

mds = read.table("results/multi.LDpruned.mds", header=TRUE)
ibm = read.table("results/multi.LDpruned.missing.mds", header=TRUE)

mds$FID = gsub(".PE", "", mds$FID)

names(ibm) = paste0(names(ibm), ".missing")
mds = cbind(mds, ibm$C1.missing, ibm$C2.missing)
names(mds) = gsub("ibm\\$", "", names(mds))

mds.ind = merge(ind.info, mds, by.x="NGS.ID", by.y="FID")

inds = levels(factor(mds.ind$Individual.ID))
new.tx.names = c("SNPRC-13245", "SNPRC-14068", "SNPRC-25567",
                 "SNPRC-27278", "SNPRC-27958", "SNPRC-28064")
inds.fixed = inds
inds.fixed[5:10] = new.tx.names

p = ggplot(mds.ind, aes(C1, C2, color=Individual.ID, pch=Sample.type)) +
        geom_point(size=3) +
        xlab("Dimension 1") + ylab("Dimension 2") +
        coord_fixed() +
        scale_color_discrete(name="Individual",
            breaks=inds,
            labels=inds.fixed) +
        scale_shape_manual(name="Sample Type", values = c(16, 17, 1))



ggsave("results/mds_ibs.pdf", plot=p)

q = ggplot(mds.ind, aes(C1.missing, C2.missing, color=Individual.ID, pch=Sample.type)) +
        geom_point(size=3) +
        xlab("Dimension 1") + ylab("Dimension 2") +
        coord_fixed() +
        scale_color_discrete(name="Individual",
            breaks=inds,
            labels=inds.fixed) +
        scale_shape_manual(name="Sample Type", values = c(16, 17, 1))

ggsave("results/mds_ibm.pdf", plot=q)

# ----------------------------------------------------------------------------------------
# --- Parse association test results
# ----------------------------------------------------------------------------------------

# LD and missingness pruned dataset
assoc = read.table("results/multi.LDpruned.assoc", header=TRUE)
table(assoc[which(assoc$P < 0.05),]$A2)
sig.assoc = nrow(assoc[which(assoc$P < 0.05),])
cat(paste(sig.assoc, "SNPs out of", nrow(assoc),
    "significantly associated with sample type. (Filtered dataset)"))

# All SNPs
assoc = read.table("results/all.assoc", header=TRUE)
table(assoc[which(assoc$P < 0.05),]$A2)
sig.assoc = nrow(assoc[which(assoc$P < 0.05),])
cat(paste(sig.assoc, "SNPs out of", nrow(assoc),
    "significantly associated with sample type. (Full dataset)"))

# ----------------------------------------------------------------------------------------
# --- Parse missingness chi-sq results
# ----------------------------------------------------------------------------------------

options(stringsAsFactors = FALSE)

# All SNPs
miss = read.table("results/all.missing", header=TRUE)
sig.miss = miss[miss$P < 0.05 & miss$F_MISS_A > miss$F_MISS_U,]
cat(paste(nrow(sig.miss), "SNPs out of", nrow(miss),
    "are significant for missingness chi-sq. (Full dataset)"))
# Output BED file of locations where allele missingness is significantly different
sig.miss.bed = do.call(rbind, strsplit(sig.miss$SNP, split=":"))
sig.miss.bed = cbind(sig.miss.bed, as.numeric(sig.miss.bed[,2]) + 1)
write.table(sig.miss.bed, "results/significantly.differently.missing.all.bed",
    quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

# LD and missingness pruned dataset
miss = read.table("results/multi.LDpruned.missing", header=TRUE)
sig.miss = miss[miss$P < 0.05 & miss$F_MISS_A > miss$F_MISS_U,]
cat(paste(nrow(sig.miss), "SNPs out of", nrow(miss),
    "are significant for missingness chi-sq. (Filtered dataset)"))
# Output BED file of locations where allele missingness is significantly different
sig.miss.bed = do.call(rbind, strsplit(sig.miss$SNP, split=":"))
sig.miss.bed = cbind(sig.miss.bed, as.numeric(sig.miss.bed[,2]) + 1)
write.table(sig.miss.bed, "results/significantly.differently.missing.LDpruned.bed",
    quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

# ----------------------------------------------------------------------------------------
# --- Analyze distance between samples
# ----------------------------------------------------------------------------------------

# LD and missingness pruned dataset
dist = read.table("results/multi.LDpruned.dist")
inds = read.table("results/multi.LDpruned.nosex")

# ---

# Remove downsampled samples
ds = grep("samp", inds$V1)
inds = inds[-ds,]
dist = dist[-ds,-ds]

inds$V1 = gsub(".PE", "", inds$V1)
inds = merge(inds, ind.info, by.x="V1", by.y="NGS.ID")

names(dist) = inds$V1
inds.dist = cbind(inds, dist)

inds.dist.m = melt(inds.dist, id=names(inds.dist)[1:ncol(inds)])

tx.inds = grep("^T", unique(inds$Individual.ID), value=TRUE)

# Skip Tx07
tx.inds = tx.inds[tx.inds != "Tx07"]

# Skip double enrichments
feces.to.skip = c("Tx05H2", "Tx02HH")
feces.to.skip.ids = unique(inds.dist.m[inds.dist.m$Sample.ID %in% feces.to.skip,1])

pull.out.dists = function (this.tx.ind) {

    this.inds.ids = inds[inds$Individual.ID == this.tx.ind,]$V1

    this.inds.dist.m = inds.dist.m[which(inds.dist.m$V1 %in% this.inds.ids &
                                         inds.dist.m$variable %in% this.inds.ids),]

    blood.ids = unique(this.inds.dist.m[this.inds.dist.m$Sample.type == "blood",]$V1)
    feces.ids = unique(this.inds.dist.m[this.inds.dist.m$Sample.type == "feces",]$V1)

    # Remove double enrichments from comparison
    feces.ids = feces.ids[!feces.ids %in% feces.to.skip.ids]

    blood.to.poop = this.inds.dist.m[which(this.inds.dist.m$V1 %in% blood.ids &
                                           this.inds.dist.m$variable %in% feces.ids),]
    poop.to.poop  = this.inds.dist.m[which(this.inds.dist.m$V1 %in% feces.ids &
                                           this.inds.dist.m$variable %in% feces.ids),]

    blood.to.poop = blood.to.poop[blood.to.poop$value != 0,]

    poop.to.poop = poop.to.poop[poop.to.poop$value != 0,]
    poop.to.poop = poop.to.poop[order(poop.to.poop$value),]
    poop.to.poop = poop.to.poop[seq(from=1, to=nrow(poop.to.poop), by=2),]

    return(list(blood.to.poop$value, poop.to.poop$value))
}

dist.list = lapply(tx.inds, pull.out.dists)

blood.poop.dists = do.call(c, do.call(rbind, lapply(tx.inds, pull.out.dists))[,1])
poop.poop.dists  = do.call(c, do.call(rbind, lapply(tx.inds, pull.out.dists))[,2])

dists = data.frame(
            distance = c(blood.poop.dists, poop.poop.dists),
            category = c(rep("Blood to feces", length(blood.poop.dists)),
                         rep("Feces to feces", length(poop.poop.dists))))

p = ggplot(dists, aes(category, distance)) +
        geom_boxplot() +
        geom_point() +
        ylab("Distance") + xlab("")

ggsave("results/intra.ind.distance.comparison.pdf", plot=p)

wilcox.test(blood.poop.dists, poop.poop.dists)
