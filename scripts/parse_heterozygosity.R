#!/usr/bin/env Rscript

# ----------------------------------------------------------------------------------------
# --- Parse heterozygosity info
# ----------------------------------------------------------------------------------------

library(stringr)
library(ggplot2)

# Get heterozygosity file (vcftools output) from user
args = commandArgs(trailingOnly = TRUE)
het.file = args[1]	# e.g. "results/baboon.pass.snp.het"

out.pdf.f       = gsub('results','reports',paste0(het.file, ".pdf"))
out.pdf.hom_obs = gsub('results','reports',paste0(het.file, ".obs_hom.pdf"))
out.txt = gsub('results','reports',paste0(het.file, ".wilcox.txt"))
out.table = gsub('results','reports',paste0(het.file, ".table.csv"))

het = read.table(het.file, header=TRUE)

ind.info = read.csv("data/individual_info.csv")

het$NGS.ID = het$INDV
het$NGS.ID = gsub("_samp.*", "", het$NGS.ID)
het$NGS.ID = gsub(".PE", "", het$NGS.ID)

het.info = merge(ind.info, het, by="NGS.ID")

# Skipping. Very coarse boxplot of inbreeding coefficient, separated by sample type
# ggplot(het.info, aes(Sample.type, F)) + geom_boxplot()

het.info$ds = c("Raw", "Downsampled")[grepl("samp", het.info$INDV) + 1]
het.info$is.special = grepl("(H2|HH)", het.info$Sample.ID)

het.info$ds.pair = str_to_upper(gsub(".*_samp-", "", het.info$INDV))
het.info[which(het.info$ds == "Raw"),]$ds.pair = NA
het.info$is.ds.pair = het.info$INDV %in% het.info$ds.pair

het.info.pairs = subset(het.info, (het.info$ds == "Downsampled" | het.info$is.ds.pair) & 
	het.info$is.special == FALSE)

het.info.pairs = droplevels(het.info.pairs)

p = ggplot(het.info.pairs, aes(Sample.type, O.HOM.)) + 
	geom_boxplot() + 
	geom_point(aes(color=Individual.ID, size=N_SITES), pch=1) +
	xlab("Sample type") + 
	ylab("Observed Homozygous Sites") + 
	scale_color_discrete(name="Individual",
		breaks=c("Tx01", "Tx02", "Tx05", "Tx09", "Tx10"),
		labels=c('SNPRC-13245','SNPRC-14068','SNPRC-25567','SNPRC-27958','SNPRC-28064')) + 
	scale_size_continuous(name="Num. sites")

q = ggplot(het.info.pairs, aes(Sample.type, F)) + 
	geom_boxplot() + 
	geom_point(aes(color=Individual.ID, size=N_SITES), pch=1) +
	xlab("Sample type") + 
	ylab("F") + 
	scale_color_discrete(name="Individual",
		breaks=c("Tx01", "Tx02", "Tx05", "Tx09", "Tx10"),
		labels=c('SNPRC-13245','SNPRC-14068','SNPRC-25567','SNPRC-27958','SNPRC-28064')) + 
	scale_size_continuous(name="Num. sites")

ggsave(out.pdf.hom_obs, plot=p)
ggsave(out.pdf.f,       plot=q)

write.csv(het.info.pairs,file=out.table,row.names=FALSE)

# Test to see if F differs
sink(out.txt)
	wilcox.test(het.info.pairs[het.info.pairs$Sample.type == "feces",]$F, 
		het.info.pairs[het.info.pairs$Sample.type == "blood",]$F)
sink()
