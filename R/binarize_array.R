mskmeans <- function(data,k=2){

	classifications <- rep(1,length(data))
	centers <- mean(data,na.rm=TRUE)

	if(k>1){
	#remove NAs from data
	hasNAs <- sum(is.na(data))>0
	if(hasNAs){
		fulldata <- data
		data <- data[!is.na(data)]
	}
	# initialise cluster means
	sorteddata <- sort(data)
	avds <- sorteddata[-1]-sorteddata[-length(sorteddata)]
	
	maxkavds <- sort(avds,decreasing=TRUE,index.return=TRUE)$ix[1:(k-1)]
	
	cutoffs <- mapply(function(x,y)mean(c(x,y)),x=sorteddata[maxkavds[1:(k-1)]],y=sorteddata[maxkavds[1:(k-1)]+1])
	sortedcutoffs <- sort(cutoffs,decreasing=FALSE)
	
	for (i in 1:(length(sortedcutoffs))){
		classifications[data>=sortedcutoffs[i]] <- i+1
	}
	
	centroidmeans <- unlist(lapply(1:k,function(x)mean(data[classifications==x])))

	clus <- kmeans(data,centers=centroidmeans) # original
	
	# if(class(data) == "numeric"){
	# 	clus <- kmeans(data,centers=centroidmeans, algorithm=c("Lloyd")) # julien: avoid error "Error: number of cluster centres must lie between 1 and nrow(x)"
	# } else {
	# 	clus <- kmeans(data,centers=centroidmeans) # original
	# }
	
	classifications <- clus$cluster
	centers <- clus$centers
		
	if(hasNAs){
		fullassignments <- rep(NA,length(fulldata))
		fullassignments[!is.na(fulldata)] <- classifications
		classifications <- fullassignments
	}

	# if(min(table(classifications))<2) classifications <- rep(1,length(data))

	}
	list(cluster=classifications,centers=centers)
	
}

clusterDisc <- function(x,use.gap){

	outvals <- rep(NA,length(x))
	notNAs <- which(!is.na(x))
	if(length(notNAs)>1){
		x <- x[notNAs]
		clusters <- list()
		gap.stats <- list()
		if(use.gap){
			for (k in 1:2){
				clusters[[k]] <- mskmeans(x,k)
				gap.stats[[k]] <- gap(data=as.matrix(x),class=clusters[[k]]$cluster,B=100,cluster.func=mskmeans)
			}
			# use gap statistic to choose k=1 or k=2
			k <- ifelse(gap.stats[[1]][1]>=(gap.stats[[2]][1]-gap.stats[[2]][2]),1,2)
		}
		if(!use.gap){
			k <- 2
			clusters[[k]] <- mskmeans(x,k)
		}
		outvals[notNAs] <- clusters[[k]]$cluster-1
	}
	outvals

}

binarize.array <- function(x,min.filter=NA,var.filter=0,fc.filter=0,na.filter=FALSE,log.base=NA,use.gap=FALSE){

	filter <- c()
	if(!is.na(min.filter)){
		cat(paste("filtering all rows with no values greater than",min.filter,"\n"))
		filter <- c(filter,which(apply(x,MARGIN=1,max,na.rm=TRUE)>min.filter))
	}
	if (var.filter>0){
		cat(paste("filtering ",var.filter*100,"% of rows with lowest variation \n",sep=""))
		sds <- apply(x,MARGIN=1,sd,na.rm=TRUE)
		sd.order <- sort(sds,decreasing=FALSE,index.return=TRUE)$ix
		filter <- c(filter,sd.order[1:floor(var.filter*nrow(x))])
        }
	if(fc.filter>0){
		cat(paste("filtering all rows with no fold-change greater than",fc.filter,"\n"))
		if(is.na(log.base)){
			fcs <- apply(x,MARGIN=1,function(y)max(y,na.rm=TRUE)/min(y,na.rm=TRUE))
			filter <- c(filter,which(fcs<fc.filter))
		}
		if(!is.na(log.base)){
			fcs <- apply(x,MARGIN=1,function(y)max(y,na.rm=TRUE)-min(y,na.rm=TRUE))
			filter <- c(filter,which(fcs<log(fc.filter,base=log.base)))
		}
	}
	if(na.filter){
		cat("filtering out all rows with missing values \n")
		filter <- c(filter,which(apply(x,MARGIN=1,function(y)sum(is.na(y)))>0))
	}

	unfiltered <- setdiff(1:nrow(x),filter)
	
	output <- array(0,dim=dim(x))
	cat(paste("Julien edited: applying cluster-based binarization to",length(unfiltered),"rows of data. This may take some time... \n"))
	if(use.gap) cat("using gap-statistic to determine cluster number. if this takes too long, try setting 'use.gap=FALSE' \n")
	output[unfiltered,] <- t(apply(x[unfiltered,],MARGIN=1,clusterDisc,use.gap=use.gap))
	rownames(output) <- rownames(x)
	colnames(output) <- colnames(x)
	output
}
