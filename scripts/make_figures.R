#!/usr/bin/env Rscript

options(stringsAsFactors=FALSE)

library(ggplot2)
library(reshape2)
library(scales)

short.names = structure(c("TF02a","TF05a","TF07a","TF09a","BZ06_051a","BZ06_053","TB02a","TB05a","TB07","TB09a","TF01a","TF02b","TF07b","TF10a","Chiou_14_030a","Chiou_14_050a","Chiou_14_065a","Chiou_14_069a","BZ06_066","BZ06_148","TB01a","TB10a","TF01b","TF02c","TF05b","TF02d","TF09b","TF10b","Chiou_14_050b","Chiou_14_050c","BZ06_051b","TF05c","Chiou_14_069b","BZ06_051c","BZ06_220","BZ07_001","BZ07_007","BZ07_029","BZ07_032","BZ07_034","BZ07_039","Chiou_14_057","Chiou_14_003","Chiou_14_044","Chiou_14_041","Chiou_14_042","BZ06_221","BZ06_227","BZ07_004","BZ07_005","BZ07_030","BZ07_035","BZ07_041","Chiou_14_001","Chiou_14_004","Chiou_14_065b","Chiou_14_030b","Chiou_14_039","BZ06_218","BZ06_224","BZ06_225","BZ07_042","BZ07_045","BZ07_047","BZ07_100","Chiou_14_036","Chiou_14_054","Chiou_14_056","Chiou_14_058","Chiou_14_059","Chiou_15_003","Chiou_15_004","Chiou_15_005","Chiou_14_005","TB01b","TB02b","TB05b","TB09b","TB10b"),.Names=c("A01F.T002","A02F.T005","A03F.T007","A04F.T009","A05F.B051","A06F.B053","A07B.T002","A08B.T005","A09B.T007","A10B.T009","B01F.T001","B02F.T002","B03F.T007","B04F.T010","B05F.C030","B06F.C050","B07F.C065","B08F.C069","B09F.B066","B10F.B148","B11B.T001","B12B.T010","C01F.T001","C02F.T002","C03F.T005","C04D.T002","C05F.T009","C06F.T010","C07F.C050","C08D.C050","C09F.B051","C10F.T005","C11D.C069","C12D.B051","D01D.B220","D02D.J001","D03D.J007","D04D.J029","D05D.J032","D06D.J034","D07D.J039","D08D.C057","D09D.C003","D10D.C044","D11D.C041","D12D.C042","D13D.B221","D14D.B227","D15D.J004","D16D.J005","D17D.J030","D18D.J035","D19D.J041","D20D.C001","D21D.C004","D22D.C065","D23D.C030","D24D.C039","D25D.B218","D26D.B224","D27D.B225","D28D.J042","D29D.J045","D30D.J047","D31D.J100","D32D.C036","D33D.C054","D34D.C056","D35D.C058","D36D.C059","D37D.H003","D38D.H004","D39D.H005","D40D.C005","E01B.T001","E02B.T002","E03B.T005","E04B.T009","E05B.T010"))

info = data.frame(id=c('A01F.T002','A02F.T005','A03F.T007','A04F.T009','A05F.B051','A06F.B053','A07B.T002','A08B.T005','A09B.T007','A10B.T009','B01F.T001','B02F.T002','B03F.T007','B04F.T010','B05F.C030','B06F.C050','B07F.C065','B08F.C069','B09F.B066','B10F.B148','B11B.T001','B12B.T010','C01F.T001','C02F.T002','C03F.T005','C04D.T002','C05F.T009','C06F.T010','C07F.C050','C08D.C050','C09F.B051','C10F.T005','C11D.C069','C12D.B051','D01D.B220','D02D.J001','D03D.J007','D04D.J029','D05D.J032','D06D.J034','D07D.J039','D08D.C057','D09D.C003','D10D.C044','D11D.C041','D12D.C042','D13D.B221','D14D.B227','D15D.J004','D16D.J005','D17D.J030','D18D.J035','D19D.J041','D20D.C001','D21D.C004','D22D.C065','D23D.C030','D24D.C039','D25D.B218','D26D.B224','D27D.B225','D28D.J042','D29D.J045','D30D.J047','D31D.J100','D32D.C036','D33D.C054','D34D.C056','D35D.C058','D36D.C059','D37D.H003','D38D.H004','D39D.H005','D40D.C005','E01B.T001','E02B.T002','E03B.T005','E04B.T009','E05B.T010'),ind.id=c('SNPRC #14068','SNPRC #25567','SNPRC #27278','SNPRC #27958','BZ06-051','BZ06-053','SNPRC #14068','SNPRC #25567','SNPRC #27278','SNPRC #27958','SNPRC #13245','SNPRC #14068','SNPRC #27278','SNPRC #28064','Chiou-14-030','Chiou-14-050','Chiou-14-065','Chiou-14-069','BZ06-066','BZ06-148','SNPRC #13245','SNPRC #28064','SNPRC #13245','SNPRC #14068','SNPRC #25567','SNPRC #14068','SNPRC #27958','SNPRC #28064','Chiou-14-050','Chiou-14-050','BZ06-051','SNPRC #25567','Chiou-14-069','BZ06-051','BZ06-220','BZ07-001','BZ07-007','BZ07-029','BZ07-032','BZ07-034','BZ07-039','Chiou-14-057','Chiou-14-003','Chiou-14-044','Chiou-14-041','Chiou-14-042','BZ06-221','BZ06-227','BZ07-004','BZ07-005','BZ07-030','BZ07-035','BZ07-041','Chiou-14-001','Chiou-14-004','Chiou-14-065','Chiou-14-030','Chiou-14-039','BZ06-218','BZ06-224','BZ06-225','BZ07-042','BZ07-045','BZ07-047','BZ07-100','Chiou-14-036','Chiou-14-054','Chiou-14-056','Chiou-14-058','Chiou-14-059','Chiou-15-003','Chiou-15-004','Chiou-15-005','Chiou-14-005','SNPRC #13245','SNPRC #14068','SNPRC #25567','SNPRC #27958','SNPRC #28064'),type=c('feces','feces','feces','feces','feces','feces','blood','blood','blood','blood','feces','feces','feces','feces','feces','feces','feces','feces','feces','feces','blood','blood','feces','feces','feces','feces','feces','feces','feces','feces','feces','feces','feces','feces','feces','feces','feces','feces','feces','feces','feces','feces','feces','feces','feces','feces','feces','feces','feces','feces','feces','feces','feces','feces','feces','feces','feces','feces','feces','feces','feces','feces','feces','feces','feces','feces','feces','feces','feces','feces','feces','feces','feces','feces','blood','blood','blood','blood','blood'),origin=c("captivity","captivity","captivity","captivity","wild","wild","captivity","captivity","captivity","captivity","captivity","captivity","captivity","captivity","wild","wild","wild","wild","wild","wild","captivity","captivity","captivity","captivity","captivity","captivity","captivity","captivity","wild","wild","wild","captivity","wild","wild","wild","wild","wild","wild","wild","wild","wild","wild","wild","wild","wild","wild","wild","wild","wild","wild","wild","wild","wild","wild","wild","wild","wild","wild","wild","wild","wild","wild","wild","wild","wild","wild","wild","wild","wild","wild","wild","wild","wild","wild","captivity","captivity","captivity","captivity","captivity"),protocol=structure(c(1,1,1,1,1,1,NA,NA,NA,NA,1,1,1,1,1,1,1,1,1,1,NA,NA,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,NA,NA,NA,NA,NA),.Label=c('manufacturer protocol','revised protocol'),class="factor"),ph=c(0.0015,0.044,0,0.1359,0.0159,0.0114,NA,NA,NA,NA,0.0355,0.0015,0,0.0987,0.0124,0.024,0.0076,0.0045,0.0485,0.0068,NA,NA,0.0433,0.0112,0.0838,1e-04,0.174,0.044,0.0299,0.0082,0.0333,0.1283,4.9e-05,0.0036,0.000289,0.000158,0.00122,3e-05,0.008793,0.004713,0.000629,0.002374,0.000339,0.000292,0.000271,0.00531,0.000894,0.000309,0.003858,0.00212,0.002748,0.00437,0.004441,7.9e-05,8.8e-05,0.00077,0.002962,0.0029,0.000581,0.000307,0.000573,0.001062,0.000262,0.000213,0.000793,0.003708,0.004415,0.001098,0.005733,0.008954,0.001017,0.031094,0.001113,0.02412,NA,NA,NA,NA,NA),dna=c(83,200,191,200,122,120,200,200,200,200,60.06,13.27,22.11,65.01,33.39,11.62,20.09,13.93,34.3,8.68,196.5,201.5,17.68,5.78,25.06,3.41,30.6,12.85,4.45,0.67,2.48,11.97,0.49,0.77,0.53,1.79,0.71,1.69,0.73,1.67,0.79,0.67,0.66,0.7,0.98,0.59,0.49,0.39,0.52,0.35,0.51,0.52,0.4,0.48,0.4,0.36,0.47,0.43,0.23,0.23,0.23,0.23,0.23,0.23,0.31,0.32,0.23,0.23,0.23,0.23,0.23,0.23,0.23,0.23,200,200,200,200,200),pid=c('A1','A1','A1','A2','A2','A2','A3','A3','A4','A4','B1','B1','B1','B1','B1','B1','B1','B1','B1','B1','B2','B2','C1','C2','C1','C3','C1','C1','C2','C4','C2','C1','C5','C6','D1','D1','D1','D1','D1','D1','D1','D1','D1','D1','D1','D1','D2','D2','D2','D2','D2','D2','D2','D2','D2','D2','D2','D2','D3','D3','D3','D3','D3','D3','D3','D3','D4','D4','D4','D4','D4','D4','D4','D4','E1','E1','E2','E1','E2'),pcr=c(24L,24L,24L,24L,24L,24L,24L,24L,24L,24L,24L,24L,24L,24L,24L,24L,24L,24L,24L,24L,24L,24L,20L,24L,20L,26L,20L,20L,24L,26L,24L,20L,26L,26L,22L,22L,22L,22L,22L,22L,22L,22L,22L,22L,22L,22L,22L,22L,22L,22L,22L,22L,22L,22L,22L,22L,22L,22L,26L,26L,26L,26L,26L,26L,26L,26L,22L,22L,22L,22L,22L,22L,22L,22L,12L,12L,12L,12L,12L))

mapping.results = lapply(short.names,function(x) {
	# First check if file exists and throw error if it doesn't
	if (!length(list.files('NGS-map/reports/',pattern=paste0(x,'.PE.bwa.baboon.aln_stats.txt')))) stop(paste('File',x,'not found!'))
	
	# Scan all the lines of the aln_stats file
	report.data = scan(file=paste0('NGS-map/reports/',x,'.PE.bwa.baboon.aln_stats.txt'),what='',sep='\n',quiet=TRUE)
	
	# Search for total reads and extract the number
	total.reads = report.data[grep('Total reads:',report.data)]
	total.reads = as.numeric(gsub('Total reads: +?([0-9]+)$','\\1',total.reads))
	
	# Search for mapped reads and extract the number
	mapped.reads = report.data[grep('Mapped reads:',report.data)]
	mapped.reads = as.numeric(gsub('^Mapped reads: +?([0-9]+).+','\\1',mapped.reads))
	
	# Return as data frame
	data.frame(total=total.reads,mapped=mapped.reads,pm=mapped.reads/total.reads)
})

# Bind all results into a data frame
mapping.results = do.call(rbind,mapping.results)

mapping.for.print = within(mapping.results,{
	total = prettyNum(total,big.mark=',')
	mapped = prettyNum(mapped,big.mark=',')
	pm = paste0(format(round(pm * 100,2),nsmall=2),'%')
})

# Now add experiment info to mapping info
info = data.frame(info,mapping.results)
rownames(info) = NULL

# Get library ID
info$lid = gsub('[0-9]+$','',info$pid)

# Plot qPCR by mapping percentage

p = ggplot(subset(info,type %in% 'feces'),aes(ph,pm,color=lid)) +
	geom_point() +
	geom_smooth(method='lm') +
	facet_wrap(~protocol,ncol=1) +
#	coord_fixed() +
	coord_cartesian(ylim=c(0,1)) +
	theme_classic() +
	scale_x_continuous(name='Detected by qPCR (pre-enrichment)',
		labels=scales::percent) +
	scale_y_continuous(name='Aligned to reference (post-enrichment)',
		labels=scales::percent) +
	scale_color_discrete(name='Library')
ggsave(p,file='reports/fecalseq_enrichment_regression.pdf',width=4,height=5,useDingbats=FALSE)


###################
# Use aln for now
###################

# Expand into two fractions
mapped = info
mapped$portion = 'host'

unmapped = info
unmapped$pm = 1 - unmapped$pm
unmapped$portion = 'nonhost'

info.all = rbind(mapped,unmapped)

# Refactor ID levels so that they sort according to mapping percentage
info.all$id = factor(info.all$id,levels=as.character(info$id[order(info$pm)]))
info.all$portion = factor(info.all$portion,levels=c('nonhost','host'))

p = ggplot(subset(info.all,type %in% 'feces'),aes(id,pm,fill=portion,alpha=origin)) +
	geom_bar(stat='identity') +
	geom_hline(aes(yintercept=mean(subset(info,type %in% 'blood')$pm)),linetype=2,color='red') +
	facet_wrap(~protocol,ncol=1,scales='free_x') +
	scale_alpha_manual(name='Origin',values=c(1/2,1)) +
	scale_fill_manual(name='Fraction',values=c('#b2df8a','#1f78b4')) +
	scale_y_continuous(labels=scales::percent) +
	theme_classic() +
	theme(axis.text.x=element_blank(),axis.ticks.x=element_blank()) +
	xlab('Sample') + ylab('Aligned to reference')
ggsave(p,file='reports/fecalseq_mapping_results.pdf',width=6,height=2.5,useDingbats=FALSE)

p = ggplot(subset(info.all,type %in% 'feces' & ph > 0 & portion %in% 'host'),aes(pm,log(pm/ph,base=10),color=protocol)) +
	geom_point(size=2) + stat_smooth(method=lm,fill='#cccccc') +
	scale_color_manual(name='Enrichment protocol',values=c('#E41A1C','#4DAF4A')) +
	scale_x_continuous(labels=scales::percent) +
	ylim(c(0.5,3.5)) +
	theme_classic() +
	theme(legend.title=element_text(size=10),legend.text=element_text(size=8),legend.margin=margin(-120,0,0,-88)) +
	xlab('Aligned to reference') + ylab(expression(log[10]('Fold enrichment')))
ggsave(p,file='reports/fecalseq_enrichment_magnitude_scatter.pdf',width=3,height=3,useDingbats=FALSE)


# colors: #fc8d62#8da0cb

p = ggplot(subset(info.all,type %in% 'feces' & ph > 0 & portion %in% 'host'),aes(protocol,log(pm/ph,base=10),color=protocol)) +
	geom_violin(draw_quantiles=0.5,size=1) + # geom_jitter(size=2) +
	scale_color_manual(name='Protocol',values=c('#E41A1C','#4DAF4A'),guide=FALSE) +
	ylim(c(0.5,3.5)) +
	theme_classic() +
	xlab('Enrichment protocol') + ylab(expression(log[10]('Fold enrichment')))
ggsave(p,file='reports/fecalseq_enrichment_magnitude_violin.pdf',width=3,height=3,useDingbats=FALSE)

mds = read.table("results/multi.LDpruned.mds", header=TRUE)
ibm = read.table("results/multi.LDpruned.missing.mds", header=TRUE)

mds$FID = gsub(".PE", "", mds$FID)

names(ibm) = paste0(names(ibm), ".missing")
mds = cbind(mds, ibm$C1.missing, ibm$C2.missing)
names(mds) = gsub("ibm\\$", "", names(mds))

info$FID = short.names[as.character(info$id)]

mds.ind = merge(info, mds, by="FID")


# p = ggplot(subset(mds.ind,origin %in% 'captivity'), aes(C1, C2, color=ind.id, shape=type)) +
# 	geom_point(size=3) +
# 	coord_fixed() +
# 	scale_color_manual(name="Individual",
# 		values=c('#7fc97f','#beaed4','#fdc086','#ffff99','#386cb0','#f0027f')) +
# 	scale_shape_manual(name="Sample Type", values = c(16, 17, 1),guide=FALSE) +
# 	scale_x_continuous(labels=function(x) format(x,nsmall=3)) +
# 	scale_y_continuous(labels=function(x) format(x,nsmall=3)) +
# 	xlab("Dimension 1") + ylab("Dimension 2") +
# 	theme_classic() +
# 	theme(legend.position='bottom',legend.margin=margin(0,0,0,-48))
# ggsave(p,file='reports/fecalseq_mds_snprc.pdf',width=4.5,height=4.5,useDingbats=FALSE)

# Keep the MDS faithful

p = ggplot(within(mds.ind,{ind.id=as.character(ind.id); ind.id[-grep('SNPRC',ind.id)] = 'Zambian (all)'}), aes(C1, C2, color=ind.id, shape=type)) +
	geom_point(size=2) +
	coord_fixed() +
	scale_color_manual(name="Individual",
		values=c('#7fc97f','#beaed4','#fdc086','#ffff99','#386cb0','#f0027f','#000000'),
		breaks=function(x) c(x[1:(length(x)-1)])) +
	scale_shape_manual(name="Sample Type", values = c(15, 17, 1),guide=FALSE) +
	xlab("Dimension 1") + ylab("Dimension 2") +
	theme_classic() +
	theme(legend.position='bottom',
		legend.key.width=unit(0.2,'in'),
		legend.key.height=unit(0.2,'in'),
		legend.title=element_text(size=10),
		legend.text=element_text(size=8))
ggsave(p,file='reports/fecalseq_mds_all.pdf',width=6,height=4.8,useDingbats=FALSE,scale=0.75)

p = ggplot(within(mds.ind,{ind.id=as.character(ind.id); ind.id[-grep('SNPRC',ind.id)] = 'Zambian (all)'}), aes(C1, C2, color=ind.id, shape=type)) +
	geom_point(size=2) +
	coord_fixed() +
	scale_color_manual(name="Individual",
		values=c('#7fc97f','#beaed4','#fdc086','#ffff99','#386cb0','#f0027f','#000000'),
		breaks=function(x) c(x[1:(length(x)-1)]),guide=FALSE) +
	scale_shape_manual(name="Sample Type", values = c(15, 17, 1)) +
	xlab("Dimension 1") + ylab("Dimension 2") +
	theme_classic() +
	theme(legend.key.width=unit(0.2,'in'),
		legend.key.height=unit(0.2,'in'),
		legend.title=element_text(size=10),
		legend.text=element_text(size=8))
ggsave(p,file='reports/fecalseq_mds_all2.pdf',width=6,height=4.8,useDingbats=FALSE,scale=0.75)

# Missingness

p = ggplot(subset(mds.ind,origin %in% 'captivity'), aes(C1.missing, C2.missing, color=ind.id, shape=type)) +
	geom_point(size=3) +
#	coord_fixed() +
	scale_color_manual(name="Individual",
		values=c('#7fc97f','#beaed4','#fdc086','#ffff99','#386cb0','#f0027f')) +
	scale_shape_manual(name="Sample Type", values = c(16, 17, 1), guide=FALSE) +
	xlab("Dimension 1") + ylab("Dimension 2") +
	theme_classic() +
	theme(legend.position='bottom')
ggsave(p,file='reports/fecalseq_ibm_snprc.pdf',width=6,height=6,useDingbats=FALSE)

# Discordance

uniq.alleles = read.csv('results/discordance.multi.txt.table.csv')
uniq.alleles$file.id = paste0(as.character(uniq.alleles$file.id),'.PE')

heterozygosity = read.csv('results/baboon.pass.snp.het.table.csv')
heterozygosity$ind.id = paste0('SNPRC #',heterozygosity$Individual.ID2)

het.pairs = data.frame(ind.id=heterozygosity$ind.id[1:(nrow(heterozygosity)/2)],down=heterozygosity$F[1:(nrow(heterozygosity)/2)],raw=heterozygosity$F[(nrow(heterozygosity)/2+1):nrow(heterozygosity)],down.type=heterozygosity$Sample.type[1:(nrow(heterozygosity)/2)],raw.type=heterozygosity$Sample.type[(nrow(heterozygosity)/2+1):nrow(heterozygosity)])

n.sites = unique(heterozygosity[,c('INDV','N_SITES')])
names(n.sites) = c('file.id','n.sites')

uniq.alleles = merge(uniq.alleles,n.sites,by='file.id',sort=FALSE,all.x=TRUE,all.y=FALSE)

p = ggplot(uniq.alleles, aes(samp, uniq.perc)) + 
	geom_boxplot() + 
	geom_point(aes(color=whoami,size=n.sites),shape=21) +
	scale_color_manual(name="Individual",
		values=c('#7fc97f','#beaed4','#fdc086','#386cb0','#f0027f')) +
	scale_size_area(name='Polymorphic sites',labels=scales::comma,max_size=4) +
	xlab("Sample type") + 
	ylab ("Proportion unique alleles") +
	guides(color = guide_legend(order=1),size = guide_legend(order=2)) +
#	guides(color = guide_legend(order=1,override.aes = list(size=3)), size = guide_legend(order=2)) +
	theme_classic() +
	theme(legend.title=element_text(size=10),legend.text=element_text(size=8),legend.margin=margin(0,0,0,-12))
ggsave(p,file='reports/fecalseq_multi_discordance.pdf',width=3,height=3,useDingbats=FALSE)


p = ggplot(heterozygosity, aes(Sample.type, F)) + 
	geom_boxplot() + 
	geom_point(aes(color=ind.id, size=N_SITES), shape=21) +
#	scale_color_manual(name="Individual",
#		values=c('#7fc97f','#beaed4','#fdc086','#386cb0','#f0027f'), guide=FALSE) +
#	geom_segment(aes(x=down.type,y=down,xend=raw.type,yend=raw,color=ind.id),data=het.pairs) +
	xlab("Sample type") + 
	ylab(expression(italic(F))) + 
	scale_color_manual(name="Individual",
		values=c('#7fc97f','#beaed4','#fdc086','#386cb0','#f0027f')) +
	scale_size_area(name='Polymorphic sites',labels=scales::comma,max_size=4) +
	scale_y_continuous(labels=function(x) format(x,nsmall=2)) +
	guides(color = guide_legend(order=1), size = guide_legend(order=2)) +
#	guides(color = guide_legend(order=1,override.aes = list(size=3)), size = guide_legend(order=2)) +
	theme_classic() +
	theme(legend.title=element_text(size=10),legend.text=element_text(size=8),legend.margin=margin(0,0,0,-12))
ggsave(p,file='reports/fecalseq_multi_inbreeding.pdf',width=3,height=3,useDingbats=FALSE)

wilcox.test(uniq.perc~samp,uniq.alleles,paired=TRUE)

wilcox.test(F~ds,heterozygosity,paired=TRUE)

# RAD stats

load('results/radtag_stats.RData')

rad.info = subset(rad.info.both,rad.mean.cov > 1)

p = ggplot(subset(rad.info, rad.mean.cov > 7),aes(length,fill=type)) +
	geom_histogram(bins=100,color='black') +
	facet_wrap(~type,ncol=1,scales='free_y') +
	xlim(0,500) +
	scale_fill_manual(values=c('red','darkgoldenrod'),guide=FALSE) +
	theme_classic() +
	xlab('Length (bp)') + ylab('Count')
ggsave(p,file='reports/fecalseq_radtag_length.pdf',width=6,height=3,useDingbats=FALSE)

info.cvg = info.cvg.m.ind

# Make a data frame of info.cvg row numbers and number of times to repeat them
reps = data.frame(id = 1:nrow(info.cvg),reads = info.cvg$num.reads)
subset.reps = apply(reps,1,function(x) {
	as.numeric(rep(x['id'],x['reads']))
})
reps.ind = do.call(c,subset.reps)

info.cvg = info.cvg[reps.ind,]

info.cvg$num.reads = NULL

p = ggplot(info.cvg,aes(gc_perc,fill=Sample.type)) +
	geom_histogram(bins=100,color='black') +
	facet_wrap(~Sample.type,ncol=1,scales='free_y') +
	scale_fill_manual(values=c('red','darkgoldenrod'),guide=FALSE) +
	scale_x_continuous(labels=scales::percent) +
	scale_y_continuous(labels=function(x) format(round(x/1e6,1),nsmall=1)) +
	theme_classic() +
	xlab('GC percentage') + ylab('Count (millions)')
ggsave(p,file='reports/fecalseq_radtag_gc.pdf',width=6,height=2.5,useDingbats=FALSE)

p = ggplot(info.cvg,aes(CpG_5000,fill=Sample.type)) +
	geom_histogram(bins=100,color='black') +
	facet_wrap(~Sample.type,ncol=1,scales='free_y') +
	scale_fill_manual(values=c('red','darkgoldenrod'),guide=FALSE) +
	scale_y_continuous(labels=function(x) format(round(x/1e6,1),nsmall=1)) +
	theme_classic() +
	xlab('CpGs') + ylab('Count (millions)')
ggsave(p,file='reports/fecalseq_radtag_cpg5000.pdf',width=6,height=2.5,useDingbats=FALSE)

p = ggplot(info.cvg,aes(length,fill=Sample.type)) +
	geom_histogram(bins=100,color='black') +
	facet_wrap(~Sample.type,ncol=1,scales='free_y') +
	xlim(0,500) +
	scale_fill_manual(values=c('red','darkgoldenrod'),guide=FALSE) +
	scale_y_continuous(labels=function(x) format(round(x/1e6,1),nsmall=1)) +
	theme_classic() +
	xlab('Length (bp)') + ylab('Count (millions)')
ggsave(p,file='reports/fecalseq_radtag_len.pdf',width=6,height=2.5,useDingbats=FALSE)



