#!/bin/sh

# ----------------------------------------------------------------------------------------
# --- Download CpG island coordinates
# ----------------------------------------------------------------------------------------

CPG_URL=http://hgdownload.soe.ucsc.edu/goldenPath/papAnu2/
CPG_URL+=database/cpgIslandExtUnmasked.txt.gz

wget $CPG_URL -O data/cpgIslandExtUnmasked.txt.gz

# Table schema:
#  `bin` smallint(6) NOT NULL,
#  `chrom` varchar(255) NOT NULL,
#  `chromStart` int(10) unsigned NOT NULL,
#  `chromEnd` int(10) unsigned NOT NULL,
#  `name` varchar(255) NOT NULL,
#  `length` int(10) unsigned NOT NULL,
#  `cpgNum` int(10) unsigned NOT NULL,
#  `gcNum` int(10) unsigned NOT NULL,
#  `perCpg` float NOT NULL,
#  `perGc` float NOT NULL,
#  `obsExp` float NOT NULL,

gunzip -c data/cpgIslandExtUnmasked.txt.gz | \
	sed -ne '/^[0-9]*\tchr/p' \
	> data/cpgIslandExtUnmasked.chr.txt

# Also make BED format file
awk -v OFS='\t' '{ print $2,$3,$4,$5$6 }' data/cpgIslandExtUnmasked.chr.txt |
	> data/cpgIslandExtUnmasked.chr.bed

exit
