#!/usr/bin/env Rscript

source("http://bioconductor.org/biocLite.R")
library("ShortRead")
library("SimRAD")

# ========================================================================================
# --- Simulate ddRAD digestion -----------------------------------------------------------
# ========================================================================================

# ----------------------------------------------------------------------------------------
# --- Redefine adapt.select to return boolean of passing, rather than subset of tags
# ----------------------------------------------------------------------------------------

adapt.select.bool = function (sequences, type = "AB+BA", 
								cut_site_5prime1, cut_site_3prime1, 
								cut_site_5prime2, cut_site_3prime2) {
	if (type == "AB+BA") {
	
		# --- AB
		# Get those that match 3p1 ("AATT")
		dig1A.bool = as.vector(isMatchingStartingAt(cut_site_3prime1, sequences))
		
		# Reverse seqs and get those that match 5p2 ("GCATG", reversed)
		sequences.rev = reverseComplement(DNAStringSet(sequences))
		re2match.5p2 = reverseComplement(DNAStringSet(cut_site_5prime2))[[1]]
		dig1B.bool = as.vector(isMatchingStartingAt(re2match.5p2, sequences.rev))
		
		# Those that match both 3p1 and 5p2
		dig1 = dig1A.bool & dig1B.bool
		
		# --- BA
		# Get those that match 3p2 ("C", reversed)
		dig2A.bool = as.vector(isMatchingStartingAt(cut_site_3prime2, sequences))
		
		# Reverse seqs and get those that match 5p1 too ("", reversed)
		re2match.5p1 = reverseComplement(DNAStringSet(cut_site_5prime1))[[1]]
		dig2B.bool = as.vector(isMatchingStartingAt(re2match.5p1, sequences.rev))

		# Those that match both 3p2 and 5p1
		dig2 = dig2A.bool & dig2B.bool
		
		dig.bool = dig1 | dig2
		
		return(dig.bool)
		
	} else {
	
		# Call normal function if type is not AB+BA
		warn(paste("adapt.select.bool not implemented for types other than AB+BA.",
					"Calling original function."))
		return(adapt.select(sequences, type, 
								cut_site_5prime1, cut_site_3prime1, 
								cut_site_5prime2, cut_site_3prime2))
	}
}

# ----------------------------------------------------------------------------------------
# --- Read in sequence to digest
# ----------------------------------------------------------------------------------------

args = commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
	quit(save = "no", status = 1, runLast = FALSE)
}

for (seq.fa in args) {

	# Infer chromosome
	chr = gsub(".*(chr[0-9X]+).*", "\\1", seq.fa)

	# Read in sequence
	seq = ref.DNAseq(seq.fa, subselect.contigs = FALSE)

	# ------------------------------------------------------------------------------------
	# --- Set restriction enzymes
	# ------------------------------------------------------------------------------------
	
	# Restriction Enzyme 1
	# MluCI  |AATT
	cs_5p1 = ""
	cs_3p1 = "AATT"
	
	# Restriction Enzyme 2
	# SphI  GCATG|C
	cs_5p2 = "GCATG"
	cs_3p2 = "C"
	
	# ------------------------------------------------------------------------------------
	# --- Digest and remove AA and BB fragments
	# ------------------------------------------------------------------------------------
	
	seq.dig = insilico.digest(seq, cs_5p1, cs_3p1, cs_5p2, cs_3p2, verbose=TRUE)
	
	# Original way, no longer needed
	#seq.sel = adapt.select(seq.dig, type="AB+BA", cs_5p1, cs_3p1, cs_5p2, cs_3p2)
	
	# Find AB+BA fragments using function that returns boolean
	seq.sel.bool = adapt.select.bool(seq.dig, type="AB+BA", 
										cs_5p1, cs_3p1, cs_5p2, cs_3p2)
	seq.sel = seq.dig[seq.sel.bool]
	
	seq.dig.lens = nchar(seq.dig)
	seq.sel.lens = nchar(seq.sel)
	
	# ------------------------------------------------------------------------------------
	# --- Write BED file of RE locations
	# ------------------------------------------------------------------------------------
	
	RE.loc = cumsum(seq.dig.lens)
	RAD.tag.starts = c(0, RE.loc[-c(length(RE.loc))])
	RAD.tag.ends   = RE.loc
	
	# Also get height of normal probability density function that models size selection
	len.dnorm = dnorm(seq.sel.lens, mean = 300, sd = 25)
	
	RAD.bed = cbind(chr, RAD.tag.starts[seq.sel.bool], RAD.tag.ends[seq.sel.bool], 
					paste(	paste0("RAD_id:", which(seq.sel.bool)), 
							paste0("length:", seq.sel.lens),
							paste0("len_dnorm:", len.dnorm),		sep=';'))
	
	out.bed = paste0("results/", gsub(".+/([^/]+)\\.fa", "\\1", seq.fa), ".bed")
	write.table(RAD.bed, file=out.bed, sep="\t", 
				quote=FALSE, row.names=FALSE, col.names=FALSE)
}
