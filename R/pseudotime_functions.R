
# PT_graph vertices, gene level stats
getPTstructures <- function(pt_graph, data, csg, cellsubset=NA, loess_span=0.5){

  # For each pseudotime point, calculate gene levels' mean and std.
  V(pt_graph)$gene_stats <- lapply(seq(igraph::vcount(pt_graph)), function(i){

    cellnames <- V(pt_graph)[[i]]$cellnames

    if(!identical(cellsubset, NA)){
      cellnames <- intersect(cellnames, cellsubset)
    }

    if (length(cellnames)==0) {
      return(matrix(
        NA,
        nrow     = nrow(data),
        ncol     = 3,
        dimnames = list(rownames(data), c("mean", "sd", "w"))))
    }

    gene_stats = do.call(rbind, lapply(seq(nrow(data[, cellnames, drop=F])), function(gid){
        x = data[gid, cellnames]
        return(list("mean" = mean(x), "sd" = sd(x), "w" = length(cellnames)))
      }))
    rownames(gene_stats) <- rownames(data)

    return(gene_stats)
  })

  pt_graph_paths <- Reduce(c, lapply(unique(csg$communities_df$subgraph_id), function(i) {

    start_community_id <- which(csg$communities_df$subgraph_id == i & csg$communities_df$type == "Start")
    end_community_ids <- which(csg$communities_df$subgraph_id == i & csg$communities_df$type == "End")

    start_node <- which(V(pt_graph)$landmark_id == start_community_id)
    end_nodes <- which(V(pt_graph)$landmark_id %in% end_community_ids)

    igraph::shortest_paths(
        pt_graph,
        from = start_node,
        to   = end_nodes)$vpath
    }))

  branched_dataset <- do.call(rbind, lapply(seq(length(pt_graph_paths)), function(i) {
    pt_graph_path <- pt_graph_paths[[i]]
    pt_graph_path_id <- V(pt_graph)$landmark_id[tail(pt_graph_path, n=1)]
    do.call(
      rbind.data.frame,
      lapply(seq(length(pt_graph_path)), function(j){
        v  <- pt_graph_path[[j]]
        gs <- as.data.frame(v$gene_stats)
        gs$gene_name <- rownames(gs)
        data.frame( 
          "name"         = rep(v$name, nrow(data)),
          "pt"           = rep(v$pt, nrow(data)),
          "path_id"      = rep(pt_graph_path_id, nrow(data)),
          "branch_id"    = rep(v$branch_id, nrow(data)),
          "branch_color" = rep(v$branch_color, nrow(data)),
          "landmark_id"  = rep(v$landmark_id, nrow(data)),
          gs)
      }))
  }))
  branched_dataset$mean <- as.numeric(branched_dataset$mean)
  branched_dataset$sd <- as.numeric(branched_dataset$sd)
  branched_dataset$w <- as.numeric(branched_dataset$w)

  # Smooth_levels of each path independently

  branched_dataset$mean.smooth = rep(NA, nrow(branched_dataset))
  branched_dataset$sd.smooth = rep(NA, nrow(branched_dataset))
  
  for (i in unique(branched_dataset$path_id)) {
    
    path_rows = which(branched_dataset$path_id == i)

    # if there is no cell on a path (all pt gene pt value are NA), return zero values
    if (all(is.na(branched_dataset[path_rows, c("mean")]))) {
      branched_dataset[path_rows, "mean.smooth"] <- rep(0, length(path_rows))
      branched_dataset[path_rows, "sd.smooth"]   <- rep(0, length(path_rows))
    } else {

      # mean
      y = reshape2::dcast(
        branched_dataset[path_rows, c("gene_name", "pt", "mean")],
        gene_name ~ pt,
        value.var="mean")
      rownames(y) = y[,1]
      y <- y[,-1]

      # pt weight
      w = reshape2::dcast(
        branched_dataset[path_rows, c("gene_name", "pt", "w")],
        gene_name ~ pt,
        value.var="w")
      rownames(w) = w[,1]
      w <- w[,-1]
      w <- unname(as.integer(w[1, ]))

      # spline fitting
      y.smooth = t(apply(y, 1, function(y2){
                        x = seq(1, length(y2))
                        navals= is.na(y2)
                        ss   <- smooth.spline(x[!navals], y2[!navals], w=w[!navals], df=3)
                        predict(ss, x)$y
                      }))

      y.smooth[which(y.smooth < 0)] <- 0 # some predict value are negative
      y.smooth.molten=reshape2::melt(y.smooth[rownames(data), ])
      branched_dataset[path_rows, "mean.smooth"] <- y.smooth.molten$value

      # sd
      y=reshape2::dcast(branched_dataset[path_rows, c("gene_name", "pt", "sd")], gene_name ~ pt, value.var="sd")
      rownames(y) = y[,1]
      y <- y[,-1]
      y.smooth <- t(apply(
        y,
        1,
        function(x){
          # "loess.control(surface="direct")" allows the fitting of NA pts at the beginning or the end of the time series (with subsetting...)
          y.loess = loess(x ~ t, span=loess_span, data.frame(x=x, t=seq(ncol(y))), na.action=na.omit, control=loess.control(surface="direct"))
          predict(y.loess, data.frame(t=seq(ncol(y))))
        }))
      y.smooth[which(y.smooth < 0)] <- 0 # some predict value are negative
      y.smooth.molten=reshape2::melt(y.smooth[rownames(data), ])
      branched_dataset[path_rows, "sd.smooth"] <- y.smooth.molten$value
    }
  }

  # Average mean.smooth and sd.smooth on overlapping branches
  mean.smooth.avg = aggregate(mean.smooth ~ name + gene_name, data=branched_dataset, FUN="mean")
  colnames(mean.smooth.avg)[3] <- "mean.smooth.avg"
  sd.smooth.avg = aggregate(sd.smooth ~ name + gene_name, data=branched_dataset, FUN="mean")
  colnames(sd.smooth.avg)[3] <- "sd.smooth.avg"

  branched_dataset = merge(x=branched_dataset, y=mean.smooth.avg, by=c("name", "gene_name"), all=FALSE)
  branched_dataset = merge(x=branched_dataset, y=sd.smooth.avg, by=c("name", "gene_name"), all=FALSE)

  list(PT_graph=pt_graph, PT=branched_dataset)

}
