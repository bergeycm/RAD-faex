#!/usr/bin/env Rscript

# ----------------------------------------------------------------------------------------
# --- Look at discordance and symmetry in SNPs that are unique to one file
# ----------------------------------------------------------------------------------------

d = read.table("results/discordance.txt")
names(d) = c("file1", "file2", "match", "total", "file1.uniq", "file2.uniq")

d$match.perc = d$match / d$total

d$file1.uniq.perc = d$file1.uniq / d$total
d$file2.uniq.perc = d$file2.uniq / d$total

d$ds = c("Raw", "Downsampled")[(grepl("samp", d$file1) | grepl("samp", d$file2)) + 1]
d$ds.blood = grepl("samp", d$file1)
d$ds.feces = grepl("samp", d$file2)
d$is.special = grepl("(H2|HH)", d$file2)

# How downsampling affects the distance (increase slightly)
boxplot(d$match.perc ~ d$ds)

# Unique SNPs in blood vs. unique SNPs in feces
boxplot(c(d$file1.uniq.perc, d$file2.uniq.perc) ~ 
	c(rep("blood", nrow(d)), rep("feces", nrow(d))))

# Split by downsampled vs. raw
d.ds = d[d$ds == "Downsampled",]
d.raw = d[d$ds == "Raw",]
boxplot(c(d.ds$file1.uniq.perc, d.ds$file2.uniq.perc) ~ 
	c(rep("blood", nrow(d.ds)), rep("feces", nrow(d.ds))))
boxplot(c(d.raw$file1.uniq.perc, d.raw$file2.uniq.perc) ~ 
	c(rep("blood", nrow(d.raw)), rep("feces", nrow(d.raw))))

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
