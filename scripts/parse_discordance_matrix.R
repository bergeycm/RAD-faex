#!/usr/bin/env Rscript

# ----------------------------------------------------------------------------------------
# --- Parse discordance matrix output from vcf-tools
# ----------------------------------------------------------------------------------------

# Get subset size from user
args = commandArgs(trailingOnly = TRUE)
disc.matrix = args[1]

diff.m = read.table(disc.matrix, header=TRUE)
row.names(diff.m) = diff.m[,1]
diff.m = diff.m[2:4,3:5]

matching.sites = diff.m[1,1] + diff.m[2,2]
matching.N = 2 * diff.m[1,1] + diff.m[2,2]

num.N.poss  = matrix(c(2,2,2,2,1,1,2,1,NA), nrow=3)
num.N.match = matrix(c(2,1,0,1,1,0,0,0,NA), nrow=3)

total.N.poss  = sum(diff.m * num.N.poss, na.rm=TRUE)
total.N.match = sum(diff.m * num.N.match, na.rm=TRUE)

match.perc = total.N.match / total.N.poss

file1.m = matrix(c(0,0,0,1,0,0,2,1,0), nrow=3)
file2.m = matrix(c(0,1,2,0,0,1,0,0,0), nrow=3)

just.file1 = sum(diff.m * file1.m, na.rm=TRUE)
just.file2 = sum(diff.m * file2.m, na.rm=TRUE)

cat(total.N.match, total.N.poss, just.file1, just.file2, sep="\t")
cat("\n")
