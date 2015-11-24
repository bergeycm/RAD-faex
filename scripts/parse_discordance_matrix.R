#!/usr/bin/env Rscript

# ----------------------------------------------------------------------------------------
# --- Parse discordance matrix output from vcf-tools
# ----------------------------------------------------------------------------------------

# Get subset size from user
args = commandArgs(trailingOnly = TRUE)
disc.matrix = args[1]

diff.m = read.table(disc.matrix, header=TRUE)
row.names(diff.m) = diff.m[,1]

# Get rid of name column and make matrix square
diff.m = diff.m[,2:5]
dim(diff.m)

matching.sites = diff.m[1,1] + diff.m[2,2] + diff.m[3,3]
matching.N = diff.m[1,1] + (2 * diff.m[2,2]) + diff.m[3,3]

# Number of unique alleles possible
num.N.poss  = matrix(c(1,2,2,1, 2,2,2,2, 2,2,1,1, 1,2,1,NA), nrow=4)
# Number of unique alleles matching
num.N.match = matrix(c(1,1,0,0, 1,2,1,0, 0,1,1,0, 0,0,0,NA), nrow=4)

total.N.poss  = sum(diff.m * num.N.poss, na.rm=TRUE)
total.N.match = sum(diff.m * num.N.match, na.rm=TRUE)

match.perc = total.N.match / total.N.poss

file1.m = matrix(c(0,0,0,0, 1,0,0,0, 2,1,0,0, 0,0,0,0), nrow=4)
file2.m = t(file1.m)

just.file1 = sum(diff.m * file1.m, na.rm=TRUE)
just.file2 = sum(diff.m * file2.m, na.rm=TRUE)

cat(total.N.match, total.N.poss, just.file1, just.file2, sep="\t")
cat("\n")
