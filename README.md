# RAD-faex
This repository contains code used in the development of the **FecalSeq** method for enriching host DNA from feces, 
by Kenneth Chiou ([website](http://www.kennychiou.com)) and Christina Bergey ([website](http://www.christinabergey.com)).
A paper describing the method been submitted for publication and is [available as a preprint](http://biorxiv.org/content/early/2015/11/25/032870)
from [bioRxiv](http://biorxiv.org/about-biorxiv).

The most up to date accompanying labwork protocol can be found in [this repository](https://github.com/kchiou/fecalseq).

- Mapping and variant identification uses the commands 
from [`fecalRAD_analysis_cmds_through_variant_ID.sh`](fecalRAD_analysis_cmds_through_variant_ID.sh) and 
makes use of the mapping pipeline from the repo, [NGS-map](http://www.github.com/bergeycm/NGS-map).
- Subsequent analytical steps are summarized in [`RADfaex_cmds.sh`](RADfaex_cmds.sh).
