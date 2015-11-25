#!/usr/bin/env Rscript

# ----------------------------------------------------------------------------------------
# --- Simulates RADtag recovery with differing levels of dropout
# ----------------------------------------------------------------------------------------

# Number of RADtags (defined in silico) falling between 150 and 450
total.radtags <- 100000
id <- 1:total.radtags
size <- round(runif(total.radtags,150,450))
radtags <- data.frame(id,size)
radtags$dropout <- FALSE

radtags$blood <- 0

# Simulate RAD reads
selection <- round(rnorm(100000,mean=300,sd=25))

# Unique RAD sizes in the size selection
selection.sizes <- sort(unique(selection))

# Number of RAD reads falling in each size
selection.number <- as.numeric(table(selection))

# Cycle through each unique RAD size
for (a in 1:length(selection.sizes)) {
	rad.size <- selection.sizes[a]
	rad.number <- selection.number[a]

	# Total number of defined RADtags that equal the current size
	num.radtags <- sum(radtags$size == rad.size)

	# Distribute the total reads matching this size among defined RADtags matching this size
	radtags.distributed <- as.numeric(table(sample(1:num.radtags,rad.number,replace=TRUE)))

	radtags$blood[radtags$size == rad.size][sample(1:num.radtags,length(radtags.distributed),replace=FALSE)] <- radtags.distributed
}

threshold <- round(max(radtags$blood) * 1.5)
sizes <- seq(1,threshold,1)

plot(sizes,sapply(sizes,function(x) sum(radtags$blood > x)),type='n',xlab='Number of reads (N)',ylab='Number of RADtags with more than N reads')
lines(sizes,sapply(sizes,function(x) sum(radtags$blood > x)),col='red',lwd=2)

# Vector of dropout percentages
percentages <- c(0.05,seq(0.1,0.50,0.1))

colors <- rainbow(percentages*1.8)

# Keep indices (j) to match colors vector
j <- 1
for (i in percentages) {

	dropout <- sample(1:total.radtags,total.radtags*i,replace=F)
	radtags$dropout <- FALSE
	radtags$dropout[dropout] <- TRUE

	radtags[[paste0('drop',i)]] <- 0
	
	for (a in 1:length(selection.sizes)) {
		rad.size <- selection.sizes[a]
		rad.number <- selection.number[a]

		# Total number of defined RADtags (excluding dropouts) that equal the current size
		num.radtags <- sum(radtags$size == rad.size & !radtags$dropout)

		# Distribute the number of size-selected RADtags among defined RADtags (excluding dropouts)
		radtags.distributed <- as.numeric(table(sample(1:num.radtags,rad.number,replace=TRUE)))

		radtags[[paste0('drop',i)]][radtags$size == rad.size & !radtags$dropout][sample(1:num.radtags,length(radtags.distributed),replace=FALSE)] <- radtags.distributed
	}

	lines(1:threshold,sapply(1:threshold,function(x) sum(radtags[[paste0('drop',i)]] > x)),col=colors[j],lty=2)
	j <- j + 1
}