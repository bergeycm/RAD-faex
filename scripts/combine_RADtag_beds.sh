#!/bin/sh

# ----------------------------------------------------------------------------------------
# --- Combine BED files of RADtags
# ----------------------------------------------------------------------------------------

cat `ls -v results/chr*.bed` > results/RADtags.bed
