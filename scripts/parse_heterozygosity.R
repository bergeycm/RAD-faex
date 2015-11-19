#!/usr/bin/env Rscript

# ----------------------------------------------------------------------------------------
# --- Parse heterozygosity info
# ----------------------------------------------------------------------------------------

library(stringr)

het = read.table("results/baboon.INDIV_DIPLOID.pass.snp.het", header=TRUE)

ind.info = read.csv("data/fecalRAD_individual_info.csv")

het$NGS.ID = het$INDV
het$NGS.ID = gsub("_samp.*", "", het$NGS.ID)
het$NGS.ID = gsub(".PE", "", het$NGS.ID)

ggplot(het.info, aes(Sample.type, F)) + geom_boxplot()

het.info$ds = c("Raw", "Downsampled")[grepl("samp", het.info$INDV) + 1]
het.info$is.special = grepl("(H2|HH)", het.info$Sample.ID)

het.info$ds.pair = paste0("fecalRAD-", 
							str_to_upper(gsub(".*_samp-", "", het.info$INDV)), 
							".PE")
het.info[which(het.info$ds == "Raw"),]$ds.pair = NA
het.info$is.ds.pair = het.info$INDV %in% het.info$ds.pair

het.info.pairs = subset(het.info, (het.info$ds == "Downsampled" | het.info$is.ds.pair) & 
	het.info$is.special == FALSE)

p = ggplot(het.info.pairs, aes(Sample.type, F)) + 
	geom_boxplot() + 
	geom_point(aes(color=Individual.ID, size=N_SITES), pch=1)

ggsave("results/heterozygosity.pdf", plot=p)
