#!/usr/bin/env Rscript

# ----------------------------------------------------------------------------------------
# --- Look at discordance and symmetry in SNPs that are unique to one file
# ----------------------------------------------------------------------------------------

library(ggplot2)

# Get discordance summary file (output from parse_all_discordance_matrices.R) from user
args = commandArgs(trailingOnly = TRUE)
d.file = args[1]	# e.g. "results/discordance.multi.txt"

d = read.table(d.file)
names(d) = c("file1", "file2", "match", "total", "file1.uniq", "file2.uniq")

d = within(d,{
	file1 = gsub('-','_',as.character(file1))
	file2 = gsub('-','_',as.character(file2))
})

d$whoami = gsub('(Tx[0-9]{2}).+','\\1',as.character(d$file1))

aliases = c(Tx01='SNPRC #13245',Tx02='SNPRC #14068',Tx05='SNPRC #25567',Tx07='SNPRC #27278',Tx09='SNPRC #27958',Tx10='SNPRC #28064')
d$whoami = aliases[as.character(d$whoami)]

d$file1[grep('_samp',d$file1)] = paste0(d$file1[grep('_samp',d$file1)],'-',d$file2[grep('_samp',d$file1)])
d$file2[grep('_samp',d$file2)] = paste0(d$file2[grep('_samp',d$file2)],'-',d$file1[grep('_samp',d$file2)])

conversion = data.frame(new.name=c("TB01a","TB01a","TB02a","TB02a","TB02a","TB02a","TB10a","TF05a","TF05b","TF09a","TF09b","TF10a","TF01a","TF01b","TF02a","TF02c","TF02b","TF02d","TF10b","TB05a","TB05a","TB09a","TB09a","TB10a"),old.name=c("Tx01B","Tx01B","Tx02B","Tx02B","Tx02B","Tx02B","Tx10B","Tx05F","Tx05H","Tx09F","Tx09H","Tx10F","Tx01F","Tx01H","Tx02F","Tx02H","Tx02F","Tx02HH","Tx10H","Tx05B","Tx05B","Tx09B","Tx09B","Tx10B"))

for (i in 1:nrow(conversion)) {
	o = conversion$old.name[i]
	n = conversion$new.name[i]
	d = within(d,{
		file1 = gsub(paste0(o,'$'),n,file1)
		file2 = gsub(paste0(o,'$'),n,file2)
		file1 = gsub(paste0(o,'_'),paste0(n,'_'),file1)
		file2 = gsub(paste0(o,'_'),paste0(n,'_'),file2)
	})
}



d$match.perc = d$match / d$total

d$file1.uniq.perc = d$file1.uniq / d$total
d$file2.uniq.perc = d$file2.uniq / d$total

d$ds = c("Raw", "Downsampled")[(grepl("samp", d$file1) | grepl("samp", d$file2)) + 1]
d$ds.blood = grepl("samp", d$file1)
d$ds.feces = grepl("samp", d$file2)
d$is.special = grepl("(H2|HH)", d$file2)

d$comparison.id = 1:nrow(d)

# How downsampling affects the distance
boxplot(d$match.perc ~ d$ds)

# Unique SNPs in blood vs. unique SNPs in feces (without regard to downsampling)
boxplot(c(d$file1.uniq.perc, d$file2.uniq.perc) ~ 
	c(rep("blood", nrow(d)), rep("feces", nrow(d))))

# Split by downsampled vs. raw
d.ds = d[d$ds == "Downsampled" & d$is.special == FALSE,]
d.raw = d[d$ds == "Raw" & d$is.special == FALSE,]
boxplot(c(d.ds$file1.uniq.perc, d.ds$file2.uniq.perc) ~ 
	c(rep("blood", nrow(d.ds)), rep("feces", nrow(d.ds))))
boxplot(c(d.raw$file1.uniq.perc, d.raw$file2.uniq.perc) ~ 
	c(rep("blood", nrow(d.raw)), rep("feces", nrow(d.raw))))

d.ds.m = data.frame(uniq.perc = c(d.ds$file1.uniq.perc, d.ds$file2.uniq.perc),
	samp=c(rep("blood", nrow(d.ds)), rep("feces", nrow(d.ds))),
	file.id=c(as.character(d.ds$file1), as.character(d.ds$file2)),
	whoami=c(d.ds$whoami,d.ds$whoami),
	comparison.id=c(d.ds$comparison.id,d.ds$comparison.id))

message('p-value is ',wilcox.test(uniq.perc~samp,data=d.ds.m,paired=TRUE)$p.value)

write.csv(d.ds.m,file=paste0(d.file, ".table.csv"),row.names=FALSE)

p = ggplot(d.ds.m, aes(samp, uniq.perc)) + 
	geom_boxplot() + 
	geom_point(aes(color=whoami)) + 
	xlab("Sample type") + 
	ylab ("Proportion unique alleles")

ggsave(file=paste0(d.file, ".boxplot.pdf"))

sink(paste0(d.file, ".basic.stats.txt"))

	# Test if unique fecal allele counts are different than unique blood allele counts
	wilcox.test(d.ds$file1.uniq,  d.ds$file2.uniq, paired=TRUE)
	wilcox.test(d.raw$file1.uniq, d.raw$file2.uniq, paired=TRUE)
	wilcox.test(d.ds$file1.uniq.perc,  d.ds$file2.uniq.perc, paired=TRUE)	# For paper
	wilcox.test(d.raw$file1.uniq.perc, d.raw$file2.uniq.perc, paired=TRUE)

	summary(d.ds$file1.uniq.perc)	# For paper, blood
	summary(d.ds$file2.uniq.perc)	# For paper, poop

sink()

ggplot(d, aes(file1.uniq.perc, file2.uniq.perc, col=ds, size=total)) + 
	geom_point() + xlim(c(0,1)) + ylim(c(0,1)) + 
	geom_abline(intercept = 0, slope=1) +
	geom_text(aes(label=paste(file1, file2)), size=2, hjust=0, vjust=0)

# Explore how residuals change with proxy for read count
d$ratio = d$file1.uniq / d$file2.uniq
d$resids = 1 - d$ratio

ggplot(d[d$ds == "Downsampled" & d$is.special == FALSE,], 
		aes(total, ratio, size=match, col=ds.feces)) + 
	geom_point() + 
	geom_abline(slope=0, intercept=1)
