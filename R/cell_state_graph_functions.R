
get_communities_from_graph <- function(
  graph,
  method                       = 'spin_glass_no_weight',
  communities_number           = NULL,
  checkCommunitiesConnectivity = TRUE,
  seed                         = 1234){

  set.seed(seed)

  if(is.null(graph))
    stop('Cell state graph must be generated first !')

  accepted.method <- c('spin_glass', 'spin_glass_no_weight', 'louvain_no_weight', 'louvain_weight', 'fast_greedy', 'fast_greedy_no_weight', 'newman_leading_eigenvector', 'newman_leading_eigenvector_no_weight', 'edge_betweenness', 'edge_betweenness_no_weight', 'info_map', 'info_map_no_weight')
  # accepted.method <- c('spin_glass', 'spin_glass_no_weight', 'louvain_no_weight', 'louvain_weight', 'fast_greedy', 'fast_greedy_no_weight', 'newman_leading_eigenvector', 'newman_leading_eigenvector_no_weight', 'edge_betweenness', 'edge_betweenness_no_weight', 'info_map', 'info_map_no_weight', 'all')

  if (! Reduce('&', method %in% accepted.method)) 
    stop("invalid communities detection method(s) (either '", paste0(accepted.method, collapse="', '"), "' must be used)")

  accepted.numbered.method = c('spin_glass', 'spin_glass_no_weight', 'fast_greedy', 'fast_greedy_no_weight')

  if(!is.null(communities_number) & (length(method) > 1 | !any(method %in% accepted.numbered.method)))
    stop("If a target number of communities is specified, only one of the following method is allowed: '", paste0(accepted.numbered.method, collapse="', '"), "'")

  # Spin glass
  # Finding communities in graphs based on statistical mechanics.
  # This function tries to find communities in graphs via a spin-glass model and simulated annealing.
  if(any(method %in% c('spin_glass', 'all'))){
    if(!is.null(communities_number)){
      communities_detection=cluster_spinglass(graph, weight=1/E(graph)$distance, spins=communities_number)
    } else {
      communities_detection=cluster_spinglass(graph, weight=1/E(graph)$distance)
    }
  }
  
  if(any(method %in% c('spin_glass_no_weight', 'all'))){
    if(!is.null(communities_number)){
      communities_detection = cluster_spinglass(graph, weight=NA, spins=communities_number)
    } else {
      communities_detection = cluster_spinglass(graph, weight=NA)
    }
  }
  
  
  # Louvain: based on the modularity measure and a hierarchial approach. Initially, each vertex is assigned to a community on its own. In every step, vertices are re-assigned to communities in a local, greedy way: each vertex is moved to the community with which it achieves the highest contribution to modularity. (igraph doc)
  if(any(method %in% c('louvain_no_weight', 'all')))
    communities_detection=cluster_louvain(graph, weights=NA)

  # weights influence ?
  if(any(method %in% c('louvain_weight', 'all')))
    communities_detection=cluster_louvain(graph, weights=1/E(graph)$distance)
  

  # Fast greedy: This algorithm merges individual nodes into communities in a way that greedily maximizes the modularity score of the graph. It can be proven that if no merge can increase the current modularity score, the algorithm can be stopped since no further increase can be achieved (igraph doc).
  # Modularity is the fraction of the edges that fall within the given groups minus the expected such fraction if edges were distributed at random (wikipedia).
  if(any(method %in% c('fast_greedy', 'all'))){
    communities_detection=cluster_fast_greedy(graph, weights=1/E(graph)$distance)
    if(!is.null(communities_number)){
      communities_detection=cutat(communities_detection, no=communities_number)
    }
  }


  if(any(method %in% c('fast_greedy_no_weight', 'all'))){
    communities_detection=cluster_fast_greedy(graph, weights=NULL)
    if(!is.null(communities_number)){
      communities_detection=cutat(communities_detection, no=communities_number)
    }
  }


  # Newman's leading eigenvector method for detecting community structure.
  # This function tries to find densely connected subgraphs in a graph 
  # by calculating the leading non-negative eigenvector 
  # of the modularity matrix of the graph.
  if(any(method %in% c('newman_leading_eigenvector', 'all')))
    communities_detection=cluster_leading_eigen(graph, weights=1/E(graph)$distance, options=list(maxiter=1000000, ncv=30))

  if(any(method %in% c('newman_leading_eigenvector_no_weight', 'all')))
    communities_detection=cluster_leading_eigen(graph, NULL, options=list(maxiter=1000000, ncv=30))

  
  # cluster_edge_betweenness
  if(any(method %in% c('edge_betweenness', 'all')))
    communities_detection=cluster_edge_betweenness(graph, weights=1/E(graph)$distance)
  
  if(any(method %in% c('edge_betweenness_no_weight', 'all')))
    communities_detection=cluster_edge_betweenness(graph, weights=NULL)
  

  # infomap ??
  if(any(method %in% c('info_map', 'all')))
    communities_detection=cluster_infomap(graph, e.weights=NULL, v.weights=1/E(graph)$distance, modularity = TRUE)
  
  if(any(method %in% c('info_map_no_weight', 'all')))
    communities_detection=cluster_infomap(graph, e.weights=NULL, v.weights=NULL, modularity = TRUE)
  
  # walktrap:  tries to find densely connected subgraphs, also called communities in a graph via random walks. The idea is that short random walks tend to stay in the same community.
  # steps <- The length of the random walks to perform.
  # communities_detection[['walktrap']]=cluster_walktrap(graph, weights=NA, steps=10)
  # # cluster_optimal: optimal community structure for a graph, in terms of maximal modularity score. (calculation takes too much time)
  # communities_detection[['optimal']]=cluster_optimal(graph, weights=NA)

  communities_with_names = as.numeric(membership(communities_detection))
  names(communities_with_names) =  V(graph)$name

  # Some community detection algorithms (eg cluster_spinglass) may define communities with multiple subcomponents. If "checkCommunitiesConnectivity" is TRUE, each community' subcomponents are forming their own communities.
  if(checkCommunitiesConnectivity){

    current_communities_number = length(unique(communities_with_names))

    for(c in sort(unique(communities_with_names))){

      cur_graph = igraph::delete_vertices(graph, v=names(communities_with_names[communities_with_names != c]))

      if(igraph::components(cur_graph)$no != 1){

        for(sc in seq(2, igraph::components(cur_graph)$no)) {

          current_communities_number <- current_communities_number + 1
          communities_with_names[names(which(igraph::components(cur_graph)$membership == sc))] <- current_communities_number

        }
      }
    }

  }

  return(communities_with_names)
}


# from http://blogs.msdn.com/b/gpalem/archive/2013/03/29/convert-igraph-r-objects-to-gexf-gephi-format.aspx
export_graph_as_gexf <- function(
  g,
  filepath = "converted_graph.gexf",
  data     = NA){

  require(igraph)
 
  # gexf nodes require two column data frame (id, label)
  # check if the input vertices has label already present
  # if not, just have the ids themselves as the label
  if(is.null(V(g)$label))
    V(g)$label <- as.character(V(g)$name)
  
  # similarily if edges does not have weight, add default 1 weight
  # if(is.null(E(g)$weight))
  E(g)$weight <- rep.int(1, ecount(g))
  
  nodes <- data.frame(cbind(V(g), V(g)$label))
  
  #edges <- t(Vectorize(get.edge, vectorize.args='id')(g, 1:ecount(g)))
  edges <- ends(g, E(g), names=FALSE)
    
  # combine all node attributes into a matrix (and take care of & for xml)
  vAttrNames <- setdiff(list.vertex.attributes(g), "label")
  
  nodesAtt <- data.frame(sapply(vAttrNames, function(attr) sub("&", "&",get.vertex.attribute(g, attr))),
                         stringsAsFactors = FALSE)
   
  # TODO: add popname in input data variable if needed
  # pop=as.matrix(getPopName(data), ncol=1)
  # colnames(pop)=c("Population")
  # nodesAtt=cbind(nodesAtt, pop)

  if(!identical(data, NA)){
    genelevels = as.data.frame(t(data))
    rownames(genelevels) = seq(1:dim(genelevels)[1])
    nodesAtt=cbind(nodesAtt, genelevels)
  }
  
  # combine all edge attributes into a matrix (and take care of & for xml)
  
  eAttrNames <- setdiff(list.edge.attributes(g), "weight")
  
  edgesAtt <- data.frame(sapply(eAttrNames, function(attr) sub("&", "&",get.edge.attribute(g, attr))),
                         
                         stringsAsFactors = FALSE)
  
  # generate the gexf object
  
  output <- rgexf::write.gexf(
    nodes,
    edges,
    edgesWeight = E(g)$weight,
    edgesAtt    = edgesAtt,
    nodesAtt    = nodesAtt)

  # suppressMessages(print(output, filepath, replace=T))
  fileConn <- file(filepath)
  writeLines(output$graph, fileConn)
  close(fileConn)
  
}

getGraphAverage <- function(graph, data_orig, rescale=.99, depth=1){

  if(depth==0){
    return(data_orig)
  }

  data_new = matrix(0, nrow=dim(data_orig)[1], ncol=dim(data_orig)[2])
  rownames(data_new) = rownames(data_orig)
  colnames(data_new) = colnames(data_orig)

  for(c in colnames(data_new)){
    neighb_cells_name = unlist(igraph::neighborhood(graph, order=depth, nodes=c, mode="all"))
    data_new[,c] <- rowMeans(data_orig[, neighb_cells_name, drop=FALSE])
  }

  if(rescale != 0){
    M99 = apply(data_orig, 1, quantile, probs=rescale)
    M100 = apply(data_orig, 1, max)
    indices = which(M99==0)
    M99[indices] = M100[indices]
    M99_new = apply(data_new, 1, quantile, probs=rescale)
    M100_new =  apply(data_new, 1, max)
    indices = which(M99_new==0)
    M99_new[indices] = M100_new[indices]
    M99_new[M99_new == 0] <- 1 # fix for null lines, avoid Nan at next line
    max_ratio = M99 / M99_new
    data_new = max_ratio * data_new
  }

  return(data_new)
}


getMagicAffinity <- function(
  weighted_adjacency_matrix,
  graph_autotune = 10,
  epsilon        = 1){

  W = weighted_adjacency_matrix

  # Autotuning distances: Normalize distances using the distance to the "knn_autotune" neighbor (here 10)
  for(j in rev(seq(nrow(W)))) {
    
    dj = W[j,]

    pos_idx = which(dj > 0)

    ordered_idx = pos_idx[order(dj[pos_idx])]
    
    ordered_d = dj[ordered_idx]

    lMaxTempIdxs = min(graph_autotune, length(ordered_d))
    
    if(lMaxTempIdxs == 0 | ordered_d[lMaxTempIdxs] == 0){
        W[j, ] <- 0
        print('Check whats up!')
    } else {
      W[j, ] = dj / ordered_d[lMaxTempIdxs]
    }
  }

  # Symmetrize W
  W <- W + t(W)

  # Affinity
  W = exp(- W / epsilon**2)
  W[W == 1] <- 0
  diag(W) <- 1 # "The distance between a cell and itself is zero, therefore its weight in the affinity matrix before normalization is 1"

  return(W)
}


getMagicAffinity_v2 <- function(
  weighted_adjacency_matrix,
  graph_autotune = 3,
  epsilon        = 1,
  kDiffmax       = 8){

  W <- weighted_adjacency_matrix

  cell_numNeighb <- colSums(W > 0)
  # we use rank to decrease the effect of outliers
  cell_numNeighb_rank = rank(cell_numNeighb, ties.method="min")/nrow(W)

  min_Neighb = min(colSums(W>0))

  kDiff = min_Neighb + as.integer((kDiffmax - min_Neighb) * cell_numNeighb_rank)

  W_discarded = matrix(0, nrow=nrow(W), ncol=nrow(W))

  # Autotuning distances: Normalize distances using the distance to the "knn_autotune" neighbor (here 10)
  for(j in rev(seq(nrow(W)))) {
    
    dj = W[j,]

    pos_idx = which(dj > 0)
    ordered_idx = pos_idx[order(dj[pos_idx])]
    ordered_d = dj[ordered_idx]

    curr_kDiff = kDiff[j]
    if(length(ordered_d) > curr_kDiff){
      print(paste0(j, " ", paste0(ordered_idx[(curr_kDiff+1):length(ordered_idx)], collapse=',')))
      W_discarded[j, ordered_idx[(curr_kDiff+1):length(ordered_idx)]] <- 1
    }

    lMaxTempIdxs = min(graph_autotune, length(ordered_d))
    
    if(lMaxTempIdxs == 0 | ordered_d[lMaxTempIdxs] == 0){
        W[j, ] <- 0
        print('Check whats up!')
    } else {
      W[j, ] = dj / ordered_d[lMaxTempIdxs]
    }
  }

  W[(W_discarded + t(W_discarded)) > 0] <- 0

  # Symmetrize W
  W <- W + t(W)

  # Affinity
  W = exp(- W / epsilon**2)
  W[W == 1] <- 0
  diag(W) <- 1 # "The distance between a cell and itself is zero, therefore its weight in the affinity matrix before normalization is 1"

  return(W)
}


# !!! Directed graph as input...
getMagicDiffusionData <- function(weighted_adjacency_matrix, data_orig, graph_autotune=10, epsilon=1, t=1, rescale=.99, kDiffmax=NA){

  if(t==0){
    return(data_orig)
  }

  # Affinity matrix
  if(identical(kDiffmax, NA)) {
    W = getMagicAffinity(weighted_adjacency_matrix, graph_autotune=graph_autotune, epsilon=epsilon)
  } else {
    W = getMagicAffinity_v2(weighted_adjacency_matrix, graph_autotune=graph_autotune, epsilon=epsilon, kDiffmax=kDiffmax)
  }
  
  #markov normalization
  D = colSums(W)
  D[D!=0] <- 1 / D[D!=0]

  L = D * W

  # L^t
  L_t = L %^% t

  # New data
  data_new = data_orig %*% t(L_t)

  if(rescale != 0){
    M99 = apply(data_orig, 1, quantile, probs=rescale)
    M100 = apply(data_orig, 1, max)
    indices = which(M99==0)
    M99[indices] = M100[indices]
    M99_new = apply(data_new, 1, quantile, probs=rescale)
    M100_new =  apply(data_new, 1, max)
    indices = which(M99_new==0)
    M99_new[indices] = M100_new[indices]
    M99_new[M99_new == 0] <- 1 # fix for null lines, avoid Nan at next line
    max_ratio = M99 / M99_new
    data_new = max_ratio * data_new
  }

  return(data_new)
}

plotColorGraph2 <- function(gf, lay, colors=NA, label=NA, shape="circle", vertex.size=8, rescale=TRUE, vertex.frame.color='black', edge.width=1, edge.color='grey', shuffle=FALSE, vertex.plot.order=NA){

  require(igraph)

  V(gf)$color=colors
 
  labels=NA
  if(!is.na(label[1])){
    labels=label
  }

  if(length(shape)==1){
    shape <- rep(shape, vcount(gf))
  }

  if(length(vertex.frame.color)==1){
    vertex.frame.color <- rep(vertex.frame.color, vcount(gf))
  }

  if(length(vertex.size)==1){
    vertex.size <- rep(vertex.size, vcount(gf))
  }

  # no rescaling for the graph, so the coordinates must belong to [-1, 1]**2
  lay <- t(apply(lay, 1, function(r){r - colMeans(lay)}))
  lay <- lay / max(abs(lay))

  # Change vertex z-order if needed
  if(!identical(vertex.plot.order, NA) | shuffle){

    if(shuffle)
      new_order = sample(vcount(gf))
    if(!identical(vertex.plot.order, NA))
      new_order = vertex.plot.order

    adj_mat = igraph::as_adjacency_matrix(gf, type="both", attr="weight", names=TRUE, sparse=TRUE)
    adj_mat <- adj_mat[new_order, new_order]
    gf <- igraph::graph_from_adjacency_matrix(adj_mat, mode="undirected", weighted="weight")

    lay <- lay[new_order, , drop=F]
    colors <- colors[new_order]
    vertex.size <- vertex.size[new_order]
    label <- label[new_order]
    shape <- shape[new_order]
    vertex.frame.color <- vertex.frame.color[new_order]    
  }

  igraph::plot.igraph(
                gf,
                layout=lay,
                vertex.color=colors,
                vertex.size=vertex.size,
                vertex.label=label,
                vertex.label.dist=.0, vertex.label.color='grey', vertex.label.cex=.2,
                vertex.shape = shape,
                vertex.frame.color=vertex.frame.color,
                edge.width=edge.width,
                edge.color=edge.color,
                rescale=rescale
                )
}


plotGeneGraph <- function(gf, lay, data, label=FALSE, title=TRUE, shape = "circle", palette = "redgreen", center_rescale = TRUE, vsize=8, colormap.truncate_outliers=F, gene.name=NULL, extra.title=NA){

  require(igraph)

  if(palette == "redgreen"){
    pal_grad <- colorRampPalette(c("red", "green"))(n = 101)
  } else if(palette == "blueyellow"){
    pal_grad <- colorRampPalette(c("blue", "yellow"))(n = 101)
  } else if(palette == "niceblueyellow"){
    pal_grad <- colorRampPalette(c("#0464DF", "#FFE800"))(n = 101)
  } else if(palette == "zscore"){
    pal_grad <- colorRampPalette(c("#191d73", "white", "#ed7901"))(n = 101)
  }

  gf$palette=pal_grad
  #currdata=norm01(data)

  if(min(data) >= 0){
    if(colormap.truncate_outliers){
      percentile.95 = quantile(x=data, .95)
      currdata=data / percentile.95
      currdata[which(currdata > 1)] <- 1

    } else {
      currdata=data / max(data)
    }
    currdata = as.integer(1+100*currdata)
  } else {

    # we assume this is z-scored data

    if(colormap.truncate_outliers){
      scaling_factor = quantile(x=abs(data), .95)
    } else {
      scaling_factor = max(abs(data))
    }
 
    currdata = data / scaling_factor
    currdata[which(currdata > 1)] <- 1
    currdata[which(currdata < -1)] <- -1

    # currdata = data / max(abs(range(data)))
    # values are in [-1, 1]
    currdata = as.integer(51+50*currdata)
  }
  
  labels=rep("", length(data))
  if(label==TRUE){
    labels=names(data)
  }

  plotColorGraph2(gf, lay, colors=gf$palette[currdata], label=labels, shape=shape, vertex.size=vsize, rescale=center_rescale, vertex.frame.color=gf$palette[currdata], shuffle=TRUE)

  # # no rescaling for the graph, so the coordinates must belong to [-1, 1]**2
  # lay <- t(apply(lay, 1, function(r){r - colMeans(lay)}))
  # lay <- lay / max(abs(lay))
  
  # igraph::plot.igraph(gf, layout=lay, vertex.label=labels, vertex.label.cex=.7, vertex.size=vsize, vertex.color=currdata, vertex.shape=shape, vertex.frame.color=currdata, rescale=center_rescale)

  if(title){
    tt = paste(gene.name, format(round(mean(as.matrix(data)), 2), nsmall = 2),format(round(max(data), 2), nsmall = 2), sep=" ")
    if(!identical(extra.title, NA))
      tt <- paste(tt, '(', extra.title, ')')
    title(tt, cex.main=.5)
  }
}

get_rMST_adjacency_matrix_edge_based <- function(distmat, RandomShotRatio=.2, numRand= NA, plot=TRUE, removedClass=2, allowed_adjacency=NA, numcores=12, debug_plot_basename=NULL){

  # require(classInt)

  if(is.na(numRand)){
      rMST_adjmatrix = aggregateRMSTs_LowMemory(distmat, conv=TRUE, RandomShotRatio=RandomShotRatio, plot=plot, allowed_adjacency=allowed_adjacency, numcores=numcores)
  } else {
    rMST_adjmatrix = aggregateRMSTs_LowMemory(distmat, conv=FALSE, Nrand=numRand, RandomShotRatio=RandomShotRatio, plot=plot, allowed_adjacency=allowed_adjacency, numcores=numcores)
  }

  if( (sum(rMST_adjmatrix[rMST_adjmatrix>0]) == length(rMST_adjmatrix[rMST_adjmatrix>0])) & (var(rMST_adjmatrix[rMST_adjmatrix>0]) == 0) ){
    print("A unique MST is generated.")
    return(rMST_adjmatrix)
  }

  # Cleaning -> Remove first class
  
  edge_norm_occurences = rMST_adjmatrix[rMST_adjmatrix>0]
  
  # Fisher clustering 
  # #################

  # cl <- classInt::classIntervals(edge_norm_occurences, style="fisher")
  # cutoff_val = cl$brks[[removedClass+1]]
  # # print(cutoff_val)
  # # print(length(intersect(which(rMST_adjmatrix>0), which(rMST_adjmatrix < cutoff_val))))
  # # print(length(which(rMST_adjmatrix>0)))
  # # cutoff_val = quantile(edge_norm_occurences, prob=.5)
  # # cutoff_val = .0
  
  # # Gaussian mixture
  # # ################
  # # Run EM algorithm multiple times as it performs quite inconsistently.
  # # select higher loglikelihood fit.
  # numfit = 10
  # k = 4

  # nm.all = lapply(seq(numfit), function(i){
  #   nm.res = mixtools::normalmixEM(edge_norm_occurences, k=k)
  #   return(list('loglik'=nm.res$loglik, 'posterior'=nm.res$posterior))
  # })

  # loglik.all = unlist(lapply(nm.all, function(x){x$loglik}))

  # kept = which.max(loglik.all)

  # for(kid in kept){

  #   posterior = nm.all[[kid]]$posterior

  #   higher_posterior_prob=apply(posterior, 1, which.max)

  #   cutoff_vals = unlist(lapply(seq(k-1), function(i){
  #     # print(range(edge_norm_occurences[which(higher_posterior_prob == i)]))
  #     .5 * (max(edge_norm_occurences[which(higher_posterior_prob == i)]) + min(edge_norm_occurences[which(higher_posterior_prob == i+1)]))
  #     }))

  #   if(!is.null(debug_plot_basename)){
  #     # png(paste0(debug_plot_basename, '/edge_occurence_clustering.png'))
  #     pdf(paste0(debug_plot_basename, '_edge_occurence_clustering.pdf'))
  #     h <- hist(edge_norm_occurences, breaks=1000, plot=FALSE)

  #     # print(cutoff_vals)
  #     # print(nm.res$mu)
  #     # print(nm.res$sigma)
  #     # print(nm.res$loglik)

  #     cuts <- cut(h$breaks, c(-Inf,cutoff_vals,Inf))
  #     plot(h, col=cuts, lty="blank", main=paste0(nm.all[[kid]]$loglik, 'Edge occurence in aggregated rMSTs (keeping top 2 clusters)'))
  #     graphics.off()
  #   }

  # }
  # rMST_adjmatrix[which(rMST_adjmatrix < cutoff_vals[2])] <- 0

  # Hardcoded thershold
  # ###################
  cutoff_val = 0.1

  if(!is.null(debug_plot_basename)){
    # png(paste0(debug_plot_basename, '/edge_occurence_clustering.png'))
    pdf(paste0(debug_plot_basename, '_edge_occurence_clustering.pdf'))
    h <- hist(edge_norm_occurences, breaks=1000, plot=FALSE)
    cuts <- cut(h$breaks, c(-Inf,cutoff_val,Inf))
    plot(h, col=cuts, lty="blank", main=paste0('Edge occurence in aggregated rMSTs (selecting frequency > ', as.integer(100*cutoff_val), '%'))
    graphics.off()
  }

  rMST_adjmatrix[which(rMST_adjmatrix < cutoff_val)] <- 0
  

  # return a binary matrix
  rMST_adjmatrix <- 1*(rMST_adjmatrix>0)

  return(rMST_adjmatrix)
}


aggregateRMSTs_LowMemory <- function(distmat, Nrand=1000, RandomShotRatio=.1, conv=FALSE, plot=TRUE, allowed_adjacency=NA, numcores=12, version='edge_based'){

  # require(classInt)
  require(parallel)

  Nsample = ncol(distmat)

  # sufficient
  adj0 = matrix(0, ncol=Nsample, nrow=Nsample)
  
  if(conv==FALSE){

    if(numcores > 1){
      num_clusters = min(Nrand, numcores)
      print(paste0('Make ', num_clusters, ' clusters'))
      cl <- parallel::makeCluster(num_clusters, outfile='/dev/null') # remove outfile arguments to stop printing workers' outputs
      
      print("Assigns variables to clusters' environments")
      for(i in seq(num_clusters)){
        parallel::clusterExport(cl[i], c('distmat', 'RandomShotRatio', 'allowed_adjacency', 'getRandomMST_edge_based'), envir=environment())
      }
      print('Register clusters (doSNOW)')
      doSNOW::registerDoSNOW(cl)
      # doMPI::registerDoMPI(cl)

      print(paste0('Generate ', Nrand,' perturbed minimal spanning trees'))


      pb <- txtProgressBar(max = Nrand, style = 3)
      progress <- function(n) setTxtProgressBar(pb, n)
      opts <- list(progress = progress)

      adj <- foreach(x = 1:Nrand, .options.snow = opts, .noexport=ls(), .packages=c('igraph'), .verbose=FALSE, .combine='+', .init=adj0, .inorder=FALSE) %dopar% {
        getRandomMST_edge_based(distmat, RandomShotRatio, allowed_adjacency)
      }

      parallel::stopCluster(cl)
      print('Stop clusters')

    } else {

      for(i in 1:Nrand){
        print(paste0(i, ' / ', Nrand))
        adj0 = adj0 + getRandomMST_edge_based(distmat, RandomShotRatio, allowed_adjacency)
      }

    }

  } else {

    stop('TODO: implement convergence case')

    # batch_size=500
    
    # num_edges_2plusclasses_list = c()
    # Nbatch=1
    # while(TRUE){

    #     print(paste("batch",Nbatch))
        
    #     set.seed(Nbatch) 
    #     # Setting different seeds changes mclapply behaviour at each calling. Yet, the result is reproducible. 
    #     # However, the result depends on the number of cores used. Setting "mc.preschedule = FALSE" would allow reproducibility with different numcores but the computation time increases a lot.
    #     adj_list = parallel::mclapply(1:batch_size, function(x) getRandomMST_edge_based(distmat, RandomShotRatio, allowed_adjacency), mc.cores=numcores, mc.set.seed = TRUE, mc.preschedule = TRUE)

    #     adj <- adj + Reduce('+', adj_list)
    #     rm(adj_list)

    #     adj.norm = adj / (Nbatch*batch_size)

    #     if( (sum(adj.norm[adj.norm>0]) == length(adj.norm[adj.norm>0])) 
    #             & (var(adj.norm[adj.norm>0]) == 0) ){
    #       print("A unique MST is generated.")
    #       break
    #     }

    #     # get class intervals for continuous numerical variables using  Jenks Natural Breaks optimization
    #     # https://en.wikipedia.org/wiki/Jenks_natural_breaks_optimization
    #     cl <- classInt::classIntervals(adj.norm[adj.norm>0], style="fisher")

    #     # first interval boundaries are [cl$brks[[1]], cl$brks[[2]]]
    #     print(paste("     # edges with every class:", as.integer(table(adj.norm > 0)[2])))
    #     num_edges_2plusclasses = as.integer(table(adj.norm > cl$brks[[2]])[2])
    #     print(paste("     # edges wo first class  :", num_edges_2plusclasses,"in interval [",cl$brks[[1]],",",cl$brks[[2]],"]",", num. class:",length(cl$brks)))
    #     if(plot){
    #       pal1 <- c("wheat1", "red3")
    #       plot(cl, pal=pal1)
    #     }
    #     num_edges_2plusclasses_list = c(num_edges_2plusclasses_list, num_edges_2plusclasses)
        
    #     # we need the last 3 batch to have the same number of "valid" edges in order to stop the loop
    #     if(Nbatch > 5){
    #       print(paste("     # variance:",var(num_edges_2plusclasses_list[(Nbatch-2):Nbatch])))
    #       if(var(num_edges_2plusclasses_list[(Nbatch-2):Nbatch]) < 1.5){
    #         break
    #       }
    #     }

    #     Nbatch=Nbatch+1
    # }
  }

  return(adj/(Nrand))
}

getRandomMST_edge_based <- function(dist, RandomShotRatio=.1, allowed_adjacency=NA){

  require(igraph)

  Nsample = dim(dist)[1]

  randmat=matrix(runif(Nsample*Nsample), Nsample, Nsample) # uniform dist. between 0 and 1
  
  if(!identical(NA, allowed_adjacency)) {
    curr_dist = dist + 999999999.0 * ( (randmat < RandomShotRatio) | !allowed_adjacency )
  } else{
     curr_dist = dist + 999999999.0 * ( randmat < RandomShotRatio )
  }

  # mode="upper" -> An undirected graph will be created, only the upper right triangle (including the diagonal) is used for the number of edges.
  g <- igraph::graph_from_adjacency_matrix( curr_dist, weighted=TRUE, mode="upper" ) # undirected weighted graph
  
  rm(curr_dist)
  gc()

  g_mst <- igraph::minimum.spanning.tree(g, algorithm='prim') # Prim's uses the weights
  
  g_mst_adj = igraph::as_adjacency_matrix(g_mst, sparse=FALSE) # this adjacency matrix stores only 0 and 1 (use attr="weight" to store 0 and weights, useless here)
  
  return(g_mst_adj)
}


get_kNN_adjacency_matrix <- function(distmat, k=3, allowed_adjacency=NA){

  if(!any(is.na(allowed_adjacency))){
    distmat[which(allowed_adjacency == 0)] <- NA # only allowed neighbors will be used
  }

  diag(distmat) <- 0 # diagonal elements have null distance

  # verify that each cell has enough potential neighbors
  if(min(apply(distmat, 1, function(r){length(which(!is.na(r)))-1})) < k){
    stop('k is too large (some cells does not have enough allowed neighbors)')
  }

  # structure similar to FNN::get.knn
  knn = list()
  knn$nn.index = do.call(rbind, 
                         lapply(1:nrow(distmat), function(i){
                # order(distmat[i,], na.last=T)[2:(1+k)]
                # remove row id to avoid error with multiple identical samples
                setdiff(order(distmat[i,], na.last=T), i)[1:k]

                          }))

  knn.adjacency = matrix(0, nrow=ncol(distmat), ncol=ncol(distmat))

  for(id in 1:ncol(distmat)){
    knn.adjacency[id, knn$nn.index[id,]] <- 1
  }

  # Knn not symetric
  knn.adjacency <- 1 * ((knn.adjacency+t(knn.adjacency)) > 0)

  return(knn.adjacency)
}


get_graph_from_adjacency_matrix <- function(d, adjmatrix){

  require(igraph)

  gf <- igraph::graph_from_adjacency_matrix(
    Matrix::Matrix(as.matrix(
      as.dist(d) * as.dist(adjmatrix))),
    weighted = TRUE,
    mode = "upper")
  
  w = E(gf)$weight
  w_max = max(w)
  E(gf)$weight = 1 + w_max-w # closer neighbor (min dist) must have the higher attraction coefficient
  E(gf)$width = 1 # for display 
  E(gf)$distance <- w

  return(gf)
}
