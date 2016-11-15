#!/usr/bin/env Rscript

options(stringsAsFactors=FALSE)

library(ggplot2)
library(reshape2)

info = read.csv('data/mapping_info.csv')

###################
# Use aln for now
###################

# Calculate mapped and unmapped fractions for plotting

info$mapped = info$reads_mapped_aln / info$reads_aln

# Calculate starting endogenousity and exogenousity

info$start_end = as.numeric(gsub('%','',info$host_percentage)) / 100

# Calculate fold-increase in host fraction

info$fold = info$mapped / info$start_end

# Assign protocol

info$protocol = NA
info$protocol[info$library %in% c('A','B') & info$type %in% 'feces'] = 'manufacturer protocol'
info$protocol[info$library %in% c('C','E') & info$type %in% 'feces'] = 'revised protocol'

# Create a barplot dataset

host = info[info$type %in% 'feces',c('id','animal','origin','library','fold','protocol','mapped')]
host = host[order(host$mapped),]
names(host)[length(names(host))] = 'fraction'

host$number = 1:nrow(host)

nonhost = host
nonhost$fraction = 1 - nonhost$fraction
host$set = 'host'
nonhost$set = 'non-host'

barset = rbind(host,nonhost)

# Make stacked barplot showing the mapping percentage. Separate by captive/wild and new/old protocol.
ggplot(barset,aes(factor(number),fraction,fill=set,alpha=origin)) +
	geom_bar(stat='identity',position='stack') +
	facet_wrap(~protocol,ncol=1,scales='free_x') +
	theme(axis.ticks.x=element_blank(),axis.text.x=element_blank(),panel.grid=element_blank()) +
	scale_fill_discrete(name='Fraction') + 
	scale_alpha_manual(name='Origin',values=c(0.25,1)) +
	xlab('Sample') + ylab('Percentage')
ggsave(file='reports/mapping_results.pdf',useDingbats=FALSE,width=7,height=5)

# Make scatterplot showing the mapping success plotted against the fold enrichment

ggplot(info[info$type %in% 'feces' & !is.infinite(info$fold),],aes(mapped * 100,log(fold,base=10),color=protocol)) +
	geom_point() +
	theme(axis.ticks.x=element_blank(),panel.grid.minor=element_blank()) +
	scale_y_continuous(limits=c(0,4)) + scale_x_continuous(limits=c(0,100)) +
	scale_color_manual(name='Protocol',values=c('#e41a1c','#4daf4a')) +
	xlab('Mapping percentage') + ylab(expression(log[10]('Fold enrichment')))
ggsave(file='reports/enrichment_magnitude.pdf',useDingbats=FALSE,width=7,height=5)