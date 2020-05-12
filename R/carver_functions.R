newClusterWeightedDistance <- function(Wcc, compact_data){

  res <- fastEuclideanDist(t(sqrt(Wcc) * compact_data))

  # Normalized distance matrix
  res <- res / sum(res)

  # add a minimal distance between cells
  res <- res + 1e-10
  diag(res) <- 0

  return(res)

}


newScoreFunctionWeightedDistance <- function(
  W,
  objs=c("timecorr", "smoothness"),
  compact_data,
  allow_multiple_roots,
  root_cluster_id,
  cluster_ids_from_cell_id,
  ordering,
  gm_smoothnesses_scaling_factors,
  compact_data_bin){

  # Normalize weights
  Wcc <- W / sum(W)

  cluster_cluster_distance <- newClusterWeightedDistance(Wcc, compact_data)

  # Calculate MST from cluster distance matrix
  all_cluster_graph <- igraph::graph_from_adjacency_matrix(
  	as.matrix(cluster_cluster_distance),
  	weighted = TRUE,
  	mode     = "upper")

  mstree <- igraph::minimum.spanning.tree(all_cluster_graph, algorithm='prim')

  # Evaluate individual
  objs_score <- list()

  if("timecorr" %in% objs){

    # Objective 1: Time correlation
    manifold_cluster_mst_pseudotime <- as.numeric(igraph::distances(
    	mstree,
    	v       = root_cluster_id,
    	weights = igraph::E(mstree)$weight))

    cell_mst_pseudotime <- manifold_cluster_mst_pseudotime[cluster_ids_from_cell_id]

    objs_score[["timecorr"]] <- - cor(
	    cell_mst_pseudotime,
	    ordering,
	    method = 'spearman')
  }

  if("timeMaxDiff" %in% objs){

    manifold_cluster_mst_pseudotime = as.numeric(igraph::distances(mstree, v=root_cluster_id, weights=NA)) # , weights=igraph::E(mstree)$weight))

    mst_pseudotime_norm = (manifold_cluster_mst_pseudotime - min(manifold_cluster_mst_pseudotime)) / (max(manifold_cluster_mst_pseudotime) - min(manifold_cluster_mst_pseudotime))

    cluster_ordering = unlist(lapply(seq(ncol(compact_data)), function(i) mean(ordering[which(cluster_ids_from_cell_id == i)])))

    cluster_ordering_norm = (cluster_ordering - min(cluster_ordering)) / (max(cluster_ordering) - min(cluster_ordering))

    # cell_mst_pseudotime_norm_avg = 

    objs_score[["timeMaxDiff"]] = max(abs(cluster_ordering_norm-mst_pseudotime_norm))
  }

  if("timeMaxDiffCells" %in% objs){

    manifold_cluster_mst_pseudotime = as.numeric(igraph::distances(mstree, v=root_cluster_id, weights=NA)) #igraph::E(mstree)$weight))

    cell_mst_pseudotime = manifold_cluster_mst_pseudotime[cluster_ids_from_cell_id]

    cell_mst_pseudotime_norm = (cell_mst_pseudotime - min(cell_mst_pseudotime)) / (max(cell_mst_pseudotime) - min(cell_mst_pseudotime))

    ordering_norm = (ordering - min(ordering)) / (max(ordering) - min(ordering))

    objs_score[["timeMaxDiffCells"]] = max(abs(ordering_norm-cell_mst_pseudotime_norm))
  }

  if("smoothness" %in% objs){

    edge_list = get.edgelist(mstree, names=FALSE)

    if(allow_multiple_roots){
      # not sure whether edge are ordered, test both vertices just in case...
      excluded_edges = c(which(edge_list[, 2] == root_cluster_id), which(edge_list[, 1] == root_cluster_id))
      edge_list <- edge_list[-excluded_edges,]
    }

    # Objective 2: Smoothness
    # We score using the clustering gms
    gm_smoothnesses = unlist(lapply(seq(nrow(compact_data)), function(i){

        gm_level = compact_data[i, ]

        gm_edge_levels = t(apply(edge_list, 1, function(x){gm_level[x]}))

        gm_mstree_smoothness = mean(abs(gm_edge_levels[,1]-gm_edge_levels[,2])) # we do not need pseudotime ordered edges if we take absolute value...
        
        gm_mstree_smoothness / gm_smoothnesses_scaling_factors[[i]]
      }))

    objs_score[['smoothness']] = sum(Wcc * gm_smoothnesses)

  }
  
  if("coverage" %in% objs){
    # Objective 3: Coverage
    # aim at maximizing the weighted coverage
    # cs = colSums(compact_data_bin * Wcc_mean / sum(Wcc_mean))
    cs = colSums(compact_data_bin * Wcc)
    objs_score[["coverage"]] = 1 / mean(cs)
  }
  
  if("regularization" %in% objs){
    objs_score[["regularization"]] = sum(abs(W))
  }

  if("complexity" %in% objs){
    objs_score[["complexity"]] = length(which(igraph::degree(mstree)==1))    
  }

  y = as.numeric(objs_score[objs])
  # print(y)

  return(y)

}


constraints_approxSumTo1 <- function(
                                    W,
                                    objs=c("timecorr", "smoothness"),
                                    # W_learned,
                                    compact_data,
                                    # optim_gms,
                                    # dataset_gms,
                                    root_cluster_id,
                                    cluster_ids_from_cell_id,
                                    # cell_ids_from_cluster_id,
                                    ordering,
                                    # clustering_num_used_gms,
                                    gm_smoothnesses_scaling_factors,
                                    compact_data_bin
                                    # training_cell_ids,
                                    # forbidden_cluster_edges
                              ){

  c(
      sum(W) - .9, # sum(x) >= 0.8
      -sum(W) +  1.1 # sum(x) <= 1.2
    )
}


# adapted from https://stackoverflow.com/questions/32791911/fast-calculation-of-pareto-front-in-python
which_pareto_efficient = function(inputData){
 
  num_points <- nrow(inputData)
 
  is_efficient <- logical(length = num_points)

  inputData.mat <- as.matrix(inputData)

  for(i in seq(num_points)){

    c <- matrix(rep(inputData.mat[i, ,drop=F], num_points), nrow = num_points, byrow = T)

    is_efficient[i] <- all(apply(c[-i, ,drop=F] < inputData.mat[-i,,drop = F], 1, any))

  }
  return(which(is_efficient))
}


plot_compact_data <- function(
  output_folder,
  compact_data) {

  silent_pheatmap <- pheatmap::pheatmap(
    compact_data,
    cluster_rows      = FALSE,
    cluster_cols      = FALSE,
    legend            = FALSE,
    annotation_legend = FALSE,
    annotation_col    = data.frame("Cluster" = seq(ncol(compact_data))),
    annotation_colors = list("Cluster" = brew_more_colors(seq(ncol(compact_data)), "Set2")),
    show_colnames     = FALSE,
    show_rownames     = TRUE,
    labels_row        = paste0("GM ", seq(nrow(compact_data))),
    silent            = TRUE)

  pdf(
    paste0(output_folder, "Compact_dataset.pdf"),
    width  = 7,
    height = 5)
  # plot(silent_pheatmap$gtable)
  grid::grid.draw(silent_pheatmap$gtable)
  graphics.off()
}

make_init_plots <- function(
	output_folder,
	pf,
	nsga2_res){

 
  if(ncol(pf) >= 2){
    
    # Plot Pareto Front(s)
    pdf(paste0(output_folder, "Optimization_Pareto_front.pdf"))
    pairs(pf, labels=colnames(pf))
    graphics.off()

    pdf(paste0(output_folder, "Optimization_Pareto_front2D.pdf"))
    plot(pf[,1], pf[,2], xlab=colnames(pf)[1], ylab=colnames(pf)[2])
    graphics.off()

    progression = do.call(rbind.data.frame, lapply(seq(length(nsga2_res)), function(i){
        parteo_opt = nsga2_res[[i]]$pareto.optimal
        rg <- cbind.data.frame(
          nsga2_res[[i]]$value[parteo_opt, 1, drop=F],
          nsga2_res[[i]]$value[parteo_opt, 2, drop=F],
          i)
        colnames(rg) <- c(colnames(pf), "Generation")
        rg
      }))

    timecorr_id = which(colnames(pf) == "timecorr")
    if (length(timecorr_id) == 1) {
    	progression[, timecorr_id] <- -progression[, timecorr_id]
    	pf[, timecorr_id] <- -pf[, timecorr_id]
    }

    pdf(paste0(output_folder, "Optimization_Pareto_progression.pdf"),
      height = 4)
    p <- ggplot() +
    	geom_point(data = progression, aes_string(x = colnames(pf)[1], y = colnames(pf)[2], color = "Generation")) + 
    	geom_line(data = progression, aes_string(x = colnames(pf)[1], y = colnames(pf)[2], color = "Generation", group = "Generation")) + 
    	geom_point(data = pf, aes_string(x = colnames(pf)[1], y = colnames(pf)[2]), color = "red") + 
    	geom_line(data = pf, aes_string(x = colnames(pf)[1], y = colnames(pf)[2]), color = "red") + 
    	xlab(colnames(pf)[1]) +
    	ylab(colnames(pf)[2]) +
    	theme_classic()
    print(p)
    graphics.off()

    for(i in seq(ncol(pf))) {
	    pdf(paste0(output_folder, "Optimization_", colnames(pf)[i], "_progression.pdf"))
	    plot(progression[,3], progression[,i], pch=16, cex=.5, xlab="Generation", ylab=colnames(pf)[i])
	    graphics.off()
		}

  }
}


make_all_loop_plots <- function(
  basename,
  expSet,
  pf_ordered_norm,
  pf_ordered_pos_corr,
  num_graphclust,
  num_graphclust_MSTs,
  num_graphclust_Ws,
  mst_graph_clusters,
  W_graph_clusters,
  pf_mst_adj.distances,
  W.diss,
  pf_mst_adj.hc,
  W.hc,
  num_solutions,
  num_obj,
  all_clusters,
  used_clusters_name,
  ps,
  MSTcluster_sizes,
  min_quantile_threshold,
  TopScore_MSTclusters_gms,
  cell_ids_from_cluster_id,
  TopScore_MST_ids_from_MSTcluster_id,
  W_learned,
  optim_gms,
  dataset_gms,
  training_cell_ids,
  compact_data,
  MSTclusters_Ws_average,
  selected_cluster_ids,
  MSTclusters_average_scores,
  forbidden_cluster_edges,
  color_maps,
  extra_pd             = NULL,
  allow_multiple_roots = FALSE,
  vs_factor            = 1000) {

  # MST distance matrix heatmap
  obj_cols = c("blue", "green", "red", "orange", "black")

  pf_integer = apply(pf_ordered_norm, 2, function(x){1+as.integer(100*x)})

  rsc = do.call(cbind, lapply(seq(num_obj), function(i){
        colorRampPalette(c("white", obj_cols[i]))(n = 101)[pf_integer[,i]]
    }))

  rsc_MST=cbind(
    rsc,
    colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(num_graphclust_MSTs)[mst_graph_clusters]
  )

  pdf(paste0(basename, "_mst_distances.pdf"), pointsize=.2) #  height=15, width=15
  heatmap.plus(
    pf_mst_adj.distances,
    Colv=as.dendrogram(pf_mst_adj.hc),
    Rowv=as.dendrogram(pf_mst_adj.hc),
    scale='none',
    labRow=paste0('(', mst_graph_clusters, ') ', seq(nrow(pf_mst_adj.distances))),
    margins=c(10, 10),
    RowSideColors=rsc_MST,
    ColSideColors=rsc_MST[,rev(seq(ncol(rsc_MST)))],
    col=rev(colorRampPalette(RColorBrewer::brewer.pal(11, "RdYlBu"))(n=1000))
    # col=colorRampPalette(c("blue", "white", "red"))(n = 1000)
    # col=colorRampPalette(c("#0464DF", "#FFE800"))(n = 1000)
    )
  legend(x="topright", legend=seq(num_graphclust_MSTs), col=colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(num_graphclust_MSTs), pch=15,
    # trace=TRUE,
    # ncol=1+as.integer(length(legend[['text']])/10),
    # ncol=8,
    # border=FALSE, bty="n", 
    y.intersp = 0.7, cex=.7,
    bg='white',
    box.col='black'
    )
  graphics.off()


  # W distance matrix heatmap
  rsc_W=cbind(
    rsc_MST,
    colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(num_graphclust_Ws)[W_graph_clusters]
  )

  pdf(paste0(basename, "_W_distances.pdf"), pointsize=.2) #  height=15, width=15
  heatmap.plus(
    as.matrix(W.diss),
    Colv=as.dendrogram(W.hc),
    Rowv=as.dendrogram(W.hc),
    scale='none',
    labRow=paste0('(', W_graph_clusters, ') ', seq(num_solutions)),
    margins=c(10, 10),
    RowSideColors=rsc_W,
    ColSideColors=rsc_W[,rev(seq(ncol(rsc_W)))],
    col=rev(colorRampPalette(RColorBrewer::brewer.pal(11, "RdYlBu"))(n=1000)))

  legend(x="topright", legend=seq(num_graphclust_Ws), col=colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(num_graphclust_Ws), pch=15,
    y.intersp = 0.7, cex=.7,
    bg='white',
    box.col='black'
    )

  graphics.off()

  # Alluvial plot to show clusters
  pdf(paste0(basename, '_Cross_clusterings.pdf'))
  test = alluvial::alluvial(
            all_clusters,
            freq=rep(1, num_solutions),
            border=NA, gap.width = 0.2, cw=0.1, blocks=T,
            col=colorRampPalette(c("black","red"))(n = num_solutions) # cw: width of cluster blocks, blocks: show blocks, 
            )
  graphics.off()


  # All solutions colors
  rsc_all = cbind(rsc, do.call(cbind,lapply(seq(ncol(all_clusters)), function(i){

        cramp = c(colorRampPalette(RColorBrewer::brewer.pal(9, "Set1")), colorRampPalette(RColorBrewer::brewer.pal(12, "Set3")), colorRampPalette(RColorBrewer::brewer.pal(8, "Set2")))[[i]](length(unique(all_clusters[,i])))

        cramp[all_clusters[,i]]
    })))


  # Plot all solutions Ws
  pdf(paste0(basename, "_allsolutions_Ws.pdf"), pointsize=.2)

  plot_ordering = order(all_clusters[, used_clusters_name])

  heatmap.plus(
  ps[plot_ordering, ],
  Colv=NA,
  Rowv=NA,
  scale='none',
  labRow=paste0('(', rep(seq(num_graphclust), MSTcluster_sizes), ') ', plot_ordering),
  RowSideColors=rsc_all[plot_ordering, ],
  col=colorRampPalette(c("white","black"))(n = 101),

  )
  graphics.off()


  # Plot Prototype Ws
  pdf(paste0(basename, "_", used_clusters_name, "clust_Top_", min_quantile_threshold,'_GMs_', num_graphclust, '_numClust.pdf'))
  image(TopScore_MSTclusters_gms, col=colorRampPalette(c("white","black"))(n = 101), xlab="GMs", ylab="Top MST", axes=FALSE)
  axis(2, at = seq(0, 1, length=ncol(TopScore_MSTclusters_gms)), labels=seq(ncol(TopScore_MSTclusters_gms))) 
  axis(1, at = seq(0, 1, length=nrow(TopScore_MSTclusters_gms)), labels=seq(nrow(TopScore_MSTclusters_gms)))
  grid(nrow(TopScore_MSTclusters_gms), ncol(TopScore_MSTclusters_gms), lwd = 1)
  graphics.off()

  plot_ordering = order(all_clusters[, used_clusters_name])

  pdf(paste0(basename, "_", used_clusters_name, "clust_Top_", min_quantile_threshold,'_GMs_', num_graphclust, '_numClust_Repeated.pdf'), pointsize=.2)
  heatmap.plus(
    t(TopScore_MSTclusters_gms[, rep(seq(num_graphclust), MSTcluster_sizes)]),#[plot_ordering,],
    Colv=NA,
    Rowv=NA, #as.dendrogram(pf_mst_adj.hc),
    scale='none',
    labRow=paste0('(', rep(seq(num_graphclust), MSTcluster_sizes), ') ', plot_ordering),
    RowSideColors=rsc_all[plot_ordering, ],
    col=colorRampPalette(c("white","black"))(n = 101),
    )
  graphics.off()

  # Plot Cluster Prototypes MSTs

  pdatas = pData(expSet) #[,c('timepoint', 'treatment', 'cells_colors')]

  pdatas$plot_order = unlist(lapply(seq(length(pData(expSet)$treatment)),
         function(i){
            t=pData(expSet)$treatment[[i]]
            if(substr(t, nchar(t)-1, nchar(t))=='RA'){
              pData(expSet)$timepoint[[i]]
            }else{
              -pData(expSet)$timepoint[[i]]
            }
          }))

  samples_cols_categories_l = list()
  vpie_vals_l = list()
  for (pd in c('cells_samples', extra_pd)) {
    samples_cols_categories_l[[pd]] = unique(pdatas[order(pdatas$plot_order), pd])
    vpie_vals_l[[pd]] = lapply(cell_ids_from_cluster_id, function(l){
        samples_cols = pdatas[l, pd]
        res = structure(
        	rep(0, length(samples_cols_categories_l[[pd]])),
        	names = samples_cols_categories_l[[pd]])
        samples_counts             <- table(samples_cols)
        res[names(samples_counts)] <- samples_counts
        res
      })
  }

  vsize = unlist(lapply(cell_ids_from_cluster_id, length))

  if(allow_multiple_roots){
    vsize <- c(vsize, 1)
    for(pd in names(vpie_vals_l)){
      vpie_vals_l[[pd]][[ncol(compact_data)]] <- setNames(rep(1, length(vpie_vals_l[[pd]][[1]])), names(vpie_vals_l[[pd]][[1]]))
    }
  }

  for(id in seq(num_graphclust)){

    i=TopScore_MST_ids_from_MSTcluster_id[id]

    if(!is.na(i)){

      W = ps[i, ]
      
      Wcc <- W / sum(W)
      cluster_cluster_distance = newClusterWeightedDistance(Wcc, compact_data)

      all_cluster_graph = igraph::graph_from_adjacency_matrix(as.matrix(cluster_cluster_distance), weighted=TRUE, mode="upper")
      mst_graph = igraph::minimum.spanning.tree(all_cluster_graph, algorithm='prim')

      set.seed(1)

      # Calculate MST layout
      temp_csg = CellStateGraph$new()
      V(mst_graph)$label <- as.character(seq(vcount(mst_graph)))
      temp_csg$graph = mst_graph
      temp_csg$project(method="gephiForceAtlas2", num_iter=10000, seed=0)
      mst_layout = temp_csg$layout
      rm(temp_csg)

      pdf(paste0(basename, "_prototype_cluster_", id, "_MSTid_", i, if(id %in% selected_cluster_ids) "_selected" else "", ".pdf"))
      par(mfrow=(length(extra_pd)+2) %>% {c(as.integer(sqrt(.)), ceiling(./as.integer(sqrt(.))))},
          mar = rep(.5, 4))
          # mar = rep(1.5, 4))
      for (pd in extra_pd) {
      	if (pd %in% names(color_maps)) {
      		vpc <- color_maps[[pd]][samples_cols_categories_l[[pd]]]
      	} else {
      		# use same default colors as pheatmap
      		vpc <- scales::dscale(
      			factor(1:length(samples_cols_categories_l[[pd]])),
      			scales::hue_pal(l = 75))[rank(samples_cols_categories_l[[pd]])]
      	}

        igraph::plot.igraph(
        	mst_graph,
					layout             = mst_layout,
					vertex.size        = vs_factor * sqrt(vsize/max(vsize))/vcount(mst_graph),
					vertex.frame.color = NA,
					vertex.label.color = 'black',
					vertex.label       = NA,
					vertex.shape       = "pie",
					vertex.pie         = vpie_vals_l[[pd]],
					vertex.pie.color   = list(vpc))
        title(pd, cex.main=.6)
      }
      igraph::plot.igraph(
      	mst_graph,
				layout             = mst_layout,
				vertex.size        = vs_factor*sqrt(vsize/max(vsize))/vcount(mst_graph),
				vertex.frame.color = NA,
				vertex.label.color = 'black',
				vertex.label       = NA,
				vertex.shape       = "circle",
        vertex.color       = brew_more_colors(seq(ncol(compact_data)), "Set2"))
      # title("Clusters", cex.main=.6)
      title(sprintf("C %.03f S %.03f", pf_ordered_pos_corr[i, 1], pf_ordered_pos_corr[i, 2]), cex.main=.5)
      graphics.off()

      pdf(paste0(basename, "_prototype_cluster_", id, "_MSTid_", i, "_GMs", ifelse(id %in% selected_cluster_ids, "_selected", ""), ".pdf"))
      nplot = nrow(compact_data) + 1
      ncol = as.integer(sqrt(nplot))
      nrow = as.integer(ceiling(nplot/ncol))
      par(mfrow = c(nrow, ncol),
        mar = c(0, 0, 2, 0), # space for one row of text at ticks and to separate plots
        xpd = NA)            # allow content to protrude into outer margin (and beyond)
      igraph::plot.igraph(
      	mst_graph,
				layout             = mst_layout,
				vertex.size        = vs_factor*sqrt(vsize/max(vsize))/vcount(mst_graph),
				vertex.frame.color = NA,
				vertex.label.color = 'black',
				vertex.label       = NA,
				vertex.shape       = "pie",
				vertex.pie         = vpie_vals_l[['cells_samples']],
				vertex.pie.color   = list(color_maps[["cells_samples"]][samples_cols_categories_l[["cells_samples"]]]))
				# vertex.pie         = vpie_vals_l[['cells_colors']],
				# vertex.pie.color   = list(samples_cols_categories_l[['cells_colors']]))
      title(sprintf("C %.03f S %.03f", pf_ordered_pos_corr[i, 1], pf_ordered_pos_corr[i, 2]), cex.main=.5)
      for(mod in seq(nrow(compact_data))){
        plotGeneGraph(
					mst_graph,
					mst_layout,
					compact_data[mod, , drop = FALSE],
					palette                    = "zscore",
					vsize                      = vs_factor*sqrt(vsize/max(vsize))/vcount(mst_graph),
					center_rescale             = TRUE,
					colormap.truncate_outliers = FALSE,
					gene.name                  = as.character(mod),
					title                      = FALSE)
        title(paste0("GM ", rownames(compact_data)[mod]), cex.main=.5, col.main="black")
      }
      graphics.off()
    }

  }

  # Average MSTs Weights
  pdf(paste0(basename, "_", used_clusters_name, "clust_Average_", min_quantile_threshold,'_GMs_', num_graphclust, '_numClust.pdf'))
  image(
  	MSTclusters_Ws_average,
  	col  = colorRampPalette(c("white","black"))(n = 101),
  	xlab = "GMs",
  	ylab = "MST cluster id",
  	axes = FALSE)
  axis(2, at = seq(0, 1, length=ncol(MSTclusters_Ws_average)), labels=seq(ncol(MSTclusters_Ws_average))) 
  axis(1, at = seq(0, 1, length=nrow(MSTclusters_Ws_average)), labels=seq(nrow(MSTclusters_Ws_average)))
  grid(nrow(MSTclusters_Ws_average), ncol(MSTclusters_Ws_average), lwd = 1)
  graphics.off()



  plot_ordering = order(all_clusters[, used_clusters_name])

  pdf(paste0(basename, "_", used_clusters_name, "clust_Average_", min_quantile_threshold,'_GMs_', num_graphclust, '_numClust_Repeated.pdf'), pointsize=.2)
  heatmap.plus(
    t(MSTclusters_Ws_average[, rep(seq(num_graphclust), MSTcluster_sizes)]),
		Colv          = NA,
		Rowv          = NA, #as.dendrogram(pf_mst_adj.hc),
		scale         = 'none',
		labRow        = paste0('(', rep(seq(num_graphclust), MSTcluster_sizes), ') ', plot_ordering),
		RowSideColors = rsc_all[plot_ordering, ],
		col           = colorRampPalette(c("white","black"))(n = 101),
    )
  graphics.off()


  # Average MSTs plots
  if(length(selected_cluster_ids) > 0){

    MSTclusters_Ws_average_selected = MSTclusters_Ws_average
    # MSTclusters_Ws_average_selected[is.na(MSTclusters_Ws_average_selected)] <- 0
    MSTclusters_Ws_average_selected[, -selected_cluster_ids] <- 0

    # Plot selected cluster average 
    pdf(paste0(basename, "_", used_clusters_name, "clust_Average_", min_quantile_threshold,'_GMs_', num_graphclust, '_numClust_Selected.pdf'))
    image(
    	MSTclusters_Ws_average_selected,
    	col  = colorRampPalette(c("white","black"))(n = 101),
    	xlab = "GMs",
    	ylab = "MST cluster id",
    	axes = FALSE)
    axis(2, at = seq(0, 1, length=ncol(MSTclusters_Ws_average_selected)), labels=seq(ncol(MSTclusters_Ws_average_selected))) 
    axis(1, at = seq(0, 1, length=nrow(MSTclusters_Ws_average_selected)), labels=seq(nrow(MSTclusters_Ws_average_selected)))
    grid(nrow(MSTclusters_Ws_average_selected), ncol(MSTclusters_Ws_average_selected), lwd = 1)
    graphics.off()

    # Plot selected cluster average 
    pdf(paste0(basename, "_", used_clusters_name, "clust_Average_", min_quantile_threshold,'_GMs_', num_graphclust, '_numClust_Repeated_Selected.pdf'), pointsize=.2)
    heatmap.plus(
      t(MSTclusters_Ws_average_selected[, rep(seq(num_graphclust), MSTcluster_sizes)]),
      Colv=NA,
      Rowv=NA, #as.dendrogram(pf_mst_adj.hc),
      scale='none',
      labRow=paste0('(', rep(seq(num_graphclust), MSTcluster_sizes), ') ', plot_ordering),
      RowSideColors=rsc_all[plot_ordering, ],
      col=colorRampPalette(c("white","black"))(n = 101),
      )
    graphics.off()


    # Plot Cluster Average
  
    for(id in selected_cluster_ids){

      W = MSTclusters_Ws_average_selected[, id]

      Wcc <- W / sum(W)
      cluster_cluster_distance = newClusterWeightedDistance(Wcc, compact_data)

      all_cluster_graph = igraph::graph_from_adjacency_matrix(as.matrix(cluster_cluster_distance), weighted=TRUE, mode="upper")
      mst_graph = igraph::minimum.spanning.tree(all_cluster_graph, algorithm='prim')

      temp_csg = CellStateGraph$new()
      V(mst_graph)$label <- as.character(seq(vcount(mst_graph)))
      temp_csg$graph = mst_graph
      temp_csg$project(method="gephiForceAtlas2", num_iter=10000, seed=0)
      mst_layout = temp_csg$layout
      rm(temp_csg)

      # Sample Ratios

      pdf(paste0(basename, "_average_cluster_", id, ".pdf"))
      par(mfrow=(length(extra_pd)+2) %>% {c(as.integer(sqrt(.)), ceiling(./as.integer(sqrt(.))))},
          mar = rep(.5, 4))
      for(pd in extra_pd) {
      	if (pd %in% names(color_maps)) {
      		vpc <- color_maps[[pd]][samples_cols_categories_l[[pd]]]
      	} else {
      		# use same default colors as pheatmap
      		vpc <- scales::dscale(
      			factor(1:length(samples_cols_categories_l[[pd]])),
      			scales::hue_pal(l = 75))[rank(samples_cols_categories_l[[pd]])]
      	}

        igraph::plot.igraph(
					mst_graph,
					layout             = mst_layout,
					vertex.size        = vs_factor*sqrt(vsize/max(vsize))/vcount(mst_graph),
					vertex.frame.color = NA,
					vertex.label.color = 'black',
					vertex.label       = NA,
					vertex.shape       = "pie",
					vertex.pie         = vpie_vals_l[[pd]],
					vertex.pie.color   = list(vpc))
        title(pd, cex.main=.6)
      }
      igraph::plot.igraph(
				mst_graph,
				layout             = mst_layout,
				vertex.size        = vs_factor*sqrt(vsize/max(vsize))/vcount(mst_graph),
				vertex.frame.color = NA,
				vertex.label.color = 'black',
				vertex.label       = NA,
				vertex.shape       = "circle",
        vertex.color       = brew_more_colors(seq(ncol(compact_data)), "Set2"))
      title(sprintf("C %.03f S %.03f", -MSTclusters_average_scores[id, 1], MSTclusters_average_scores[id, 2]), cex.main=.5)
      graphics.off()

      # Plot with gene patterns

      pdf(paste0(basename, "_average_cluster_", id, "_GMs.pdf"))

      nplot = nrow(compact_data) + 1
      ncol = as.integer(sqrt(nplot))
      nrow = as.integer(ceiling(nplot/ncol))
      
      par(mfrow = c(nrow, ncol),
        mar = c(0, 0, 2, 0), # space for one row of text at ticks and to separate plots
        xpd = NA)            # allow content to protrude into outer margin (and beyond)

      igraph::plot.igraph(
				mst_graph,
				layout             = mst_layout,
				vertex.size        = vs_factor*sqrt(vsize/max(vsize))/vcount(mst_graph),
				vertex.frame.color = NA,
				vertex.label.color = 'black',
				vertex.label       = NA,
				vertex.shape       = "pie",
				vertex.pie         = vpie_vals_l[['cells_samples']],
				vertex.pie.color   = list(color_maps[["cells_samples"]][samples_cols_categories_l[["cells_samples"]]]))
      title(sprintf("C %.03f S %.03f", -MSTclusters_average_scores[id, 1], MSTclusters_average_scores[id, 2]), cex.main=.5)

      for(mod in seq(nrow(compact_data))){
        plotGeneGraph(
        	mst_graph,
        	mst_layout,
        	compact_data[mod, , drop=F],
					palette                    = "zscore",
					vsize                      = vs_factor*sqrt(vsize/max(vsize))/vcount(mst_graph),
					center_rescale             = TRUE,
					colormap.truncate_outliers = FALSE,
					gene.name                  = as.character(mod),
					title                      = FALSE)
        title(paste0("GM ", rownames(compact_data)[mod]), cex.main=.5, col.main="black")
      }

      graphics.off()

    }
  }
}


make_final_plot <- function(
	basename,
	expSet,
	cell_ids_from_cluster_id,
	Final_mst_graph,
	compact_data,
	color_maps,
	extra_pd             = NULL,
	allow_multiple_roots = FALSE,
	vs_factor            = 1000){

  pdatas = pData(expSet) #[,c('timepoint', 'treatment', 'cells_colors')]

  # # remove replicate bias on colors
  # pdatas$cells_colors = generate_cells_color(pData(expSet), replicate_variation=F)

  pdatas$plot_order = unlist(lapply(seq(length(pData(expSet)$treatment)),
         function(i){
            t=pData(expSet)$treatment[[i]]
            if(substr(t, nchar(t)-1, nchar(t))=='RA'){
              pData(expSet)$timepoint[[i]]
            }else{
              -pData(expSet)$timepoint[[i]]
            }
          }))

  samples_cols_categories_l = list()
  vpie_vals_l = list()
  for(pd in c('cells_samples', extra_pd)){
    samples_cols_categories_l[[pd]] = unique(pdatas[order(pdatas$plot_order), pd])
    vpie_vals_l[[pd]] = lapply(cell_ids_from_cluster_id, function(l){
        samples_cols = pdatas[l, pd]
        res = structure(
        	rep(0, length(samples_cols_categories_l[[pd]])),
        	names = samples_cols_categories_l[[pd]])
        samples_counts             <- table(samples_cols)
        res[names(samples_counts)] <- samples_counts
        res
      })
  }

  vsize = unlist(lapply(cell_ids_from_cluster_id, length))

  if(allow_multiple_roots){
    vsize <- c(vsize, 1)
    for(pd in names(vpie_vals_l)){
      vpie_vals_l[[pd]][[ncol(compact_data)]] <- setNames(rep(1, length(vpie_vals_l[[pd]][[1]])), names(vpie_vals_l[[pd]][[1]]))
    }
  }

  # Calculate MST layout
  temp_csg = CellStateGraph$new()
  temp_csg$graph = Final_mst_graph
  temp_csg$project(method="gephiForceAtlas2", num_iter=50000, seed=0)
  Final_mst_layout = temp_csg$layout
  rm(temp_csg)

  n_row <- as.integer(sqrt(length(extra_pd)+2))
  n_col <- ceiling((length(extra_pd)+2) / n_row)

  pdf(
    paste0(basename, "_Final_cluster_MST_graphs.pdf"),
    width = n_col * 3,
    height = n_row * 3)
  par(
    mfrow = c(n_row, n_col),
    mar   = c(.5, .5, 1, 3),
    xpd=TRUE)
  for (pd in extra_pd) {
  	if (pd %in% names(color_maps)) {
  		vpc <- color_maps[[pd]][samples_cols_categories_l[[pd]]]
  	} else {
  		vpc <- scales::dscale(
  			factor(1:length(samples_cols_categories_l[[pd]])),
  			scales::hue_pal(l = 75))[rank(samples_cols_categories_l[[pd]])]
  	}
    igraph::plot.igraph(
    	Final_mst_graph,
			layout             = Final_mst_layout,
			vertex.size        = vs_factor*sqrt(vsize/max(vsize))/vcount(Final_mst_graph),
			vertex.frame.color = NA,
			vertex.label.color = 'black',
			vertex.label       = NA,
			vertex.shape       = "pie",
			vertex.pie         = vpie_vals_l[[pd]],
			vertex.pie.color   = list(vpc))
    title(pd)
  }
  igraph::plot.igraph(
		Final_mst_graph,
		layout             = Final_mst_layout,
		vertex.size        = vs_factor*sqrt(vsize/max(vsize))/vcount(Final_mst_graph),
		vertex.frame.color = NA,
		vertex.label.color = 'black',
		vertex.label       = NA,
		vertex.shape       = "circle",
		vertex.color       = brew_more_colors(seq(ncol(compact_data)), "Set2"))
  title("Carving Clusters")
  graphics.off()

  pdf(paste0(basename, "_Final_cluster_MST_graph_GMs.pdf"))

    nplot = nrow(compact_data) + 1
    ncol = as.integer(sqrt(nplot))
    nrow = as.integer(ceiling(nplot/ncol))
    
    par(mfrow = c(nrow, ncol),
      mar = c(0, 0, 2, 0), # space for one row of text at ticks and to separate plots
      xpd = NA)            # allow content to protrude into outer margin (and beyond)

    igraph::plot.igraph(
			Final_mst_graph,
			layout             = Final_mst_layout,
			vertex.size        = vs_factor*sqrt(vsize/max(vsize))/vcount(Final_mst_graph),
			vertex.frame.color = NA,
			vertex.label.color = 'black',
			vertex.label       = NA,
			vertex.shape       = "pie",
			vertex.pie         = vpie_vals_l[['cells_samples']],
			vertex.pie.color   = list(color_maps[["cells_samples"]][samples_cols_categories_l[["cells_samples"]]]))

    for(mod in seq(nrow(compact_data))){

      plotGeneGraph(
      	Final_mst_graph,
				Final_mst_layout,
				compact_data[mod, , drop   = F],
				palette                    = "zscore",
				vsize                      = vs_factor*sqrt(vsize/max(vsize))/vcount(Final_mst_graph),
				center_rescale             = TRUE,
				colormap.truncate_outliers = TRUE,
				gene.name                  = as.character(mod),
				title                      = FALSE)
      title(paste0("GM ", rownames(compact_data)[mod]), cex.main=.5, col.main="black")
    }

    graphics.off()


}


get_weighted_distance_matrix <- function(
	readcounts,
	W_local,
	W_local_pow    = 1,
	data_transform = "log_scale"){
 
  if(W_local_pow != 1){
    W_local <- W_local ** W_local_pow
    W_local <- apply(W_local, 2, function(x){x/sum(x)})
  }

  # Weighted distance between log-scaled gene levels
  if(data_transform == "log_scale") {
    rc <- t(scale(t(log(readcounts+1)), center=T, scale=T))
  } else if(data_transform == "log") {
    rc <- log(readcounts+1)
  } else if(data_transform == "none") {
    rc <- readcounts 
  } else if(data_transform == "scale") {
    rc <- t(scale(t(readcounts), center=T, scale=T))
  } else {
    stop("'data_transform' must be either 'none', 'log' or 'log_scale'")
  }
  
  # Check whether rc contains NaN
  na_genes <- which(is.na(rowSums(rc)))
  num_na   <- length(na_genes)
  if(num_na > 0){
    print(paste0("Warning: ", num_na, " genes contains NaN values after '", data_transform, "'' transformation (", paste0(rownames(rc)[na_genes], collapse=', '), "). It is most likely null genes which should be removed upstream in the pipeline. These genes are ignored from current cell-cell distance calculation."))
    rc[na_genes, ] <- 0
  }

  res <- fastEuclideanDist(t(sqrt(W_local) * rc))
  
  # Normalized distance matrix
  res <- res / sum(res)

  # add a minimal distance between cells
  res <- res + 1e-10
  diag(res) <- 0

  return(res)
}

plot_SA_history <- function(
  trace.mat,
  optim_objectives,
  output_folder) {

  obj_label <- optim_objectives

  if (optim_objectives == "timecorr") {
    trace.mat[, "function.value"] <- -trace.mat[, "function.value"]
    trace.mat[, "current.minimum"] <- -trace.mat[, "current.minimum"]
    obj_label <- "Correlation"
  }

  pdf(output_folder, height = 4)

  p <- ggplot(data = data.frame(trace.mat)) +
    geom_point(
      aes(x = nb.steps, y = function.value),
      size = .3) +
    geom_point(
      aes(x = nb.steps, y = current.minimum),
      color = "red",
      size = .3) +
    ylab(obj_label) +
    ggtitle("Simulated Annealing History") +
    theme_classic()

  print(p)

  graphics.off()
}