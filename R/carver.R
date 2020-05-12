
# Allow igraph S3 class to be used as field of Carver
setOldClass("igraph")

Carver <- setRefClass(

  Class = "Carver",

  fields = list(
    antler_env            = "environment",
    W_optim               = "numeric",
    optim_res_all_batches = "list",
    plot_folder_path      = "character",
    runVariables          = "list",
    allow_multiple_roots  = "logical"
    ),

  methods = list(

    initialize = function(
      ...,
      antler_env            = new.env(),
      W_optim               = numeric(),
      optim_res_all_batches = list(),
      plot_folder_path      = character(0),
      runVariables          = list(),
      allow_multiple_roots  = logical()) {

        callSuper(
          ...,
          antler_env            = antler_env,
          W_optim               = W_optim,
          optim_res_all_batches = optim_res_all_batches,
          plot_folder_path      = plot_folder_path,
          runVariables          = runVariables,
          allow_multiple_roots  = allow_multiple_roots)
      },

    # Local GA optimizatioin
    run = function(

      gene_modules_name                            = antler_env$gene_modules$names(),
      optim_strategy                               = c("GA", "SA"),
      optim_objectives                             = c("timecorr", "smoothness", "complexity"), # currently only 1 objecitve can be selected for SA
      
      GA_solutions_selection_quantile_threshold = .9,            # only valid for GA
      GA_solutions_selection_objectives         = c("timecorr"), # only valid for GA
      GA_solutions_selection_strategy           = "average",     # "top", "average", only valid for GA
      GA_num_generations                        = seq(10),       # only valid for GA
      GA_popsize                                = 8,             # only valid for GA
      GA_constraint                             = NULL,          # NULL, "SumTo1",  # only valid for GA
      
      SA_control                                   = list(
        smooth          = FALSE,
        maxit           = 500,
        nb.stop.improvement = 200,
        trace.mat       = TRUE,
        simple.function = FALSE,
        verbose         = TRUE,
        seed            = 1), # only valid for SA
      
      seed                                         = NULL, # re-initialize random generator
      num_clusters                                 = NULL,
      size_limit                                   = NULL,
      size_limit_ratio                             = .1,
      
      extra_pd                                     = NULL,
      polish                                       = TRUE,
      allow_multiple_roots                         = FALSE,
      clustering_timepoint_first                   = FALSE,

      plot_folder_path                             = NA){

      set.seed(seed)

      stopifnot(length(antler_env$gene_modules$names()) > 0)

      stopifnot(!identical(antler_env$readcounts_norm, matrix()))

      allow_multiple_roots <<- allow_multiple_roots

      plot_folder_path <<- plot_folder_path
      dir.create(plot_folder_path, showWarnings = FALSE, recursive = TRUE)

      gene_modules_name <- match.arg(gene_modules_name)

      optim_strategy    <- match.arg(optim_strategy)

      ordering          <- pData(antler_env$expressionSet)$timepoint
      
      gene_modules      <- antler_env$gene_modules$get(gene_modules_name)
      
      # Calculate cell clusters -------------------------------------------------------------

      # print('Scale dataset')
      data_logscaled = t(scale(
        t(log(antler_env$readcounts_norm+1)),
        center = T,
        scale  = T))

      # average by gene modules
      cat('Average by gene modules\n')
      data_gmavg <- do.call(
        rbind,
        lapply(gene_modules,
          function(i){
            colMeans( data_logscaled[i, , drop = FALSE], na.rm=TRUE) # some subset gene level may be null
          }
        )
      )

      if (clustering_timepoint_first) {
        cell_ids_from_cluster_id <- lapply(
          sort(unique(ordering)),
          function(i){
            which(ordering == i)
          })
      } else {
        cell_ids_from_cluster_id <- list(seq(antler_env$num_cells()))
      }

      if (is.null(size_limit)) {
        current_size_limit <- as.integer(antler_env$num_cells() * size_limit_ratio)
      } else {
        current_size_limit <- size_limit
      }

      large_clusters = function(x, current_size_limit){
        cluster_sizes <- unlist(lapply(x, length))
        which(cluster_sizes > current_size_limit)
      }

      # Iterative dichotomy to obtain similarly-sized clusters
      while (length(large_clusters(cell_ids_from_cluster_id, current_size_limit)) > 0) {
       
        cluster_to_split          <- large_clusters(cell_ids_from_cluster_id, current_size_limit)[1]
        cluster_to_split_cell_ids <- cell_ids_from_cluster_id[[cluster_to_split]]  
        euclidean.dist            <- fastEuclideanDist(t(data_gmavg[, cluster_to_split_cell_ids]))
        
        hc          <- hclust(as.dist(euclidean.dist), method = "ward.D2")
        cluster_ids <- cutree(hc, k=2)
        cell_ids_from_cluster_id[[cluster_to_split]]                   <- cluster_to_split_cell_ids[which(cluster_ids == 1)]
        cell_ids_from_cluster_id[[length(cell_ids_from_cluster_id)+1]] <- cluster_to_split_cell_ids[which(cluster_ids == 2)]
        # reorder clusters (send last one to its correct position)
        new_order <- c(
          seq(cluster_to_split),
          length(cell_ids_from_cluster_id),
          if (cluster_to_split < length(cell_ids_from_cluster_id) - 1) {
            seq(cluster_to_split + 1, length(cell_ids_from_cluster_id) - 1)
          } else {
            NULL
          })
        cell_ids_from_cluster_id <- cell_ids_from_cluster_id[new_order]          
      }
      
      cluster_ids_from_cell_id <- rep(
        seq(length(cell_ids_from_cluster_id)),
        unlist(lapply(cell_ids_from_cluster_id, length)))[order(unlist(cell_ids_from_cluster_id))]

      num_clusters <- length(cell_ids_from_cluster_id)

      cat(paste0(
        'Number of clusters with iterative dichotomy: ',
        num_clusters,
        ' (max cluster size ',
        current_size_limit,
        ')\n'))

      
      cat('Average expression by cell clusters\n')
      compact_data <- do.call(
        cbind,
        lapply(
          seq(num_clusters),
          function(i){
            rowMeans(data_gmavg[, cell_ids_from_cluster_id[[i]], drop = FALSE])
            }
          ))
      colnames(compact_data) <- seq(num_clusters)
      rownames(compact_data) <- names(gene_modules)
      
      plot_compact_data(
        paste0(plot_folder_path, '/'),
        compact_data)

      # print(paste0("Calculate cluster average ordering"))
      cluster_ordering <- unlist(lapply(
        seq(num_clusters),
        function(i){
          # ordering ordered by batch
          orderings <- ordering[cell_ids_from_cluster_id[[i]]]
          outliers_thres <- quantile(x = orderings, c(.2, .8))
          # outliers_thres = quantile(x=orderings, c(0, 1))
          orderings[which(orderings < outliers_thres[[1]])] <- outliers_thres[[1]]
          orderings[which(orderings > outliers_thres[[2]])] <- outliers_thres[[2]]
          mean(orderings)
        }))

      # print('GM smoothness scaling factors')
      # Need to be the same length than the score objective
      cluster_pairs <- combn(seq(num_clusters), 2)
      gm_smoothnesses_scaling_factors <- unlist(lapply(
        seq(length(gene_modules)),
        function(i){
          gm_pair_levels <- t(apply(
            cluster_pairs,
            2,
            function(x){
              compact_data[i, x]
            }))
          mean(abs(gm_pair_levels[,1] - gm_pair_levels[,2]))  # normalization factors
        }))

      # Add extra root cluster if allow_multiple_roots is TRUE
      if(allow_multiple_roots){

        # select clusters with lowest average ordering value
        orderings <- sort(unique(ordering))
        early_clusters <- which(cluster_ordering < orderings[2])

        if(length(early_clusters) < 2){
          print("dirty fix for early_clusters")
          print(sort(cluster_ordering))
          early_clusters <- order(cluster_ordering)[1:2]
        }

        # add root cluster to compact_data
        root_cluster_id <- num_clusters + 1
        compact_data    <- cbind(
          compact_data,
          rowMeans(compact_data[, early_clusters]))
        colnames(compact_data)[root_cluster_id] <- root_cluster_id
      } else {
        root_cluster_id <- which.min(cluster_ordering)
        # print(paste0("Define root clusters as lowest '", root_cluster_id, "'"))
      }

      if (optim_strategy == "GA") {

        # Run GA Optimization --------------------------------------------------------------

        cat("Start GA optimization\n")

        if(identical(GA_constraint, 'SumTo1')){
          constraints <- constraints_approxSumTo1
          cdim        <- 2
        } else {
          constraints <- NULL
          cdim        <- 0
        }

        GA_res <- mco::nsga2(
          fn                              = newScoreFunctionWeightedDistance,
          idim                            = length(gene_modules),
          odim                            = length(optim_objectives), 
          objs                            = optim_objectives,
          compact_data                    = compact_data,
          allow_multiple_roots            = allow_multiple_roots,
          root_cluster_id                 = root_cluster_id,
          cluster_ids_from_cell_id        = cluster_ids_from_cell_id,
          ordering                        = ordering,
          gm_smoothnesses_scaling_factors = gm_smoothnesses_scaling_factors,
          compact_data_bin                = compact_data_bin,           
          constraints                     = constraints,
          cdim                            = cdim,
          generations                     = GA_num_generations,
          popsize                         = GA_popsize,
          lower.bounds                    = rep(0, length(gene_modules)),
          upper.bounds                    = rep(1, length(gene_modules)))
        
        class(GA_res) <- "list"

        optim_res_all_batches <<- GA_res
          
        # Process Optimization Results --------------------------------------------------------

        print("Process solutions")

        # pareto front from all generations
        aggregated_scores        <- do.call(
          rbind,
          lapply(
            seq_along(GA_res),
            function(i){
              cbind(
                GA_res[[i]]$value[GA_res[[i]]$pareto.optimal, ,drop = FALSE],
                "gen" = i)
              }))
        colnames(aggregated_scores) <- c(optim_objectives, "gen")
        unique_ids               <- !duplicated(aggregated_scores[, -ncol(aggregated_scores)])
        aggregated_scores_unique <- aggregated_scores[unique_ids, -ncol(aggregated_scores), drop = FALSE]

        if (nrow(aggregated_scores_unique) == 0) {
          stop("No pareto optimal solution found. Increase population size / num. generation or remove constraint.")
        }

        aggregated_scores_unique.df   <- as.data.frame(aggregated_scores_unique)

        order_ids                     <- do.call(order, aggregated_scores_unique.df,)
        
        aggregated_scores_ordered     <- aggregated_scores_unique.df[order_ids,,drop = FALSE]

        front_ids                     <- which_pareto_efficient(aggregated_scores_ordered)

        pareto_front                  <- aggregated_scores_ordered[front_ids, , drop = FALSE]
        colnames(pareto_front)        <- optim_objectives

        aggregated_parameters         <- do.call(
          rbind,
          lapply(
            seq_along(GA_res),
            function(i) {
              cbind(
                GA_res[[i]]$par[GA_res[[i]]$pareto.optimal, , drop = FALSE],
                "gen" = i)
            }))
        aggregated_parameters_unique  <- aggregated_parameters[unique_ids, -ncol(aggregated_parameters), drop = FALSE]
        aggregated_parameters_ordered <- aggregated_parameters_unique[order_ids, , drop = FALSE]

        ps                            <- aggregated_parameters_ordered[front_ids, , drop = FALSE]

        num_obj                       <- ncol(pareto_front)

        num_solutions                 <- nrow(pareto_front)

        cat(paste0("\nNum. solutions on Pareto: ", num_solutions, "\n"))

        pf_range <- apply(-pareto_front, 2, range)

        pf_norm <- do.call(
          cbind,
          lapply(
            seq(num_obj),
            function(i){        
              if (pf_range[1, i] != pf_range[2, i]) {
                (-pareto_front[, i] - pf_range[1, i]) / (pf_range[2, i] - pf_range[1, i])
              } else {
                rep(1, num_solutions)
              }
          }))

        selection_objs_ids <- which(optim_objectives %in% GA_solutions_selection_objectives)
        GA_solutions_selection_score <- apply(pf_norm[, selection_objs_ids, drop = FALSE], 1, geo_mean)
        min_threshold <- quantile(
          GA_solutions_selection_score,
          prob = GA_solutions_selection_quantile_threshold)

        make_init_plots(
          paste0(plot_folder_path, '/'),
          pareto_front,
          GA_res)

        if (GA_solutions_selection_strategy == "top" | num_obj == 1 | num_solutions == 1) {

          top_id     <- which.max(GA_solutions_selection_score)
          top_Ws     <- ps[top_id, ]
          top_scores <- pareto_front[top_id, ]

        } else if (GA_solutions_selection_strategy == "average") {

          # Solution clustering via MST topology
          # ####################################

          print("Calculate MST distance")

          pf_mst_adj <- lapply(seq(num_solutions), function(i){

            W                        <- as.numeric(ps[i, , drop = FALSE])
            Wcc                      <- W / sum(W)
            cluster_cluster_distance <- newClusterWeightedDistance(Wcc, compact_data)
            
            all_cluster_graph        <- igraph::graph_from_adjacency_matrix(
              as.matrix(cluster_cluster_distance),
              weighted = TRUE,
              mode     = "upper")

            mst_graph                <- igraph::minimum.spanning.tree(
              all_cluster_graph,
              algorithm = 'prim')

            igraph::as_adjacency_matrix(
              mst_graph,
              sparse = FALSE,
              type = "both",
              attr = "weight")
          })

          # normalize so that all total MST edge distances are equal
          pf_mst_adj <- lapply(pf_mst_adj, function(x){

                        # # keep distance
                        # x/sum(x)

                        # edges length to 1
                        1*(x>0)

                  })

          # MST shortest path from all to all
          pf_rep_list <- lapply(pf_mst_adj, function(x){

                  sps <- igraph::distances(igraph::graph_from_adjacency_matrix(x, weighted=TRUE, mode='undirected'))
                  sps
                  })

          # Manifold distances
          ngraph <- length(pf_rep_list)
          pf_mst_adj.distances <- matrix(0, ngraph, ngraph, dimnames=list(seq(ngraph), seq(ngraph)))

          upper.tri.idxs <- which(upper.tri(pf_mst_adj.distances), arr.ind=T)

          upper_pf_mst_adj.distances <- unlist(mclapply(
            1:nrow(upper.tri.idxs),
            function(x){

              i <- upper.tri.idxs[x, 1]
              j <- upper.tri.idxs[x, 2]

              if(i >= j) 
                stop('i >= j')

              # Euclidean distance (suited for variable edge length)
              return(sqrt(sum((pf_rep_list[[i]]-pf_rep_list[[j]])**2)))
              
              # # Chebyshev distance (suited for constant edge length)
              # return(max(abs(pf_rep_list[[i]]-pf_rep_list[[j]])))

            },
            mc.cores = 1))

          pf_mst_adj.distances[upper.tri(pf_mst_adj.distances, diag=F)] <- upper_pf_mst_adj.distances
          pf_mst_adj.distances[lower.tri(pf_mst_adj.distances)]         <- t(pf_mst_adj.distances)[lower.tri(pf_mst_adj.distances)]

          pf_mst_adj.hc <- hclust(as.dist(pf_mst_adj.distances), method = "complete")

          # Solution clustering via Weights
          # ###############################

          print("Calculate W distance")

          W.diss <- dist(ps)
          W.hc <- hclust(W.diss, method="complete") #, method = "ward.D2")

          # Crossed clusters
          # ################

          print("Evaluate starting cluster number")

          # MST clusters number
          if(sum(pf_mst_adj.distances) != 0){
            num_graphclust_MSTs <- get_optimal_cluster_number(
              pf_mst_adj.distances,
              pf_mst_adj.hc,
              100000000,
              1,
              display = FALSE,
              method  = "min_area")
          } else {
            # in simple cases, all MSTs are identical
            num_graphclust_MSTs <- 1
          }

          # W clusters number
          num_graphclust_Ws <- get_optimal_cluster_number(
            as.matrix(W.diss),
            W.hc,
            100000000,
            1,
            display = FALSE,
            method  = "min_area")

          # Loop starts here
          # ################

          goon    <- TRUE
          loop_id <- 0

          while(goon){

            # MST clusters
            mst_graph_clust_ids <- cutree(pf_mst_adj.hc, k=num_graphclust_MSTs)
            clustid_ordered     <- order(unique(mst_graph_clust_ids[pf_mst_adj.hc$order]))

            mst_graph_clusters  <- clustid_ordered[mst_graph_clust_ids]

            # W clusters
            W_clust_ids         <- cutree(W.hc, k=num_graphclust_Ws)
            W_clustid_ordered   <- order(unique(W_clust_ids[W.hc$order]))

            W_graph_clusters    <- W_clustid_ordered[W_clust_ids]

            print("Cross both cluterings")

            both_cluster_ids    <- cbind(
              "MST" = sprintf("%03d", mst_graph_clusters),
              "W"   = sprintf("%03d", W_graph_clusters))
            both_cluster_names  <- apply(both_cluster_ids, 1, paste0, collapse='_')

            all_clusters        <- data.frame(
            'MST'     = mst_graph_clusters,
            'X'       = both_cluster_names,
            "W"       = W_graph_clusters,
            row.names = seq(num_solutions))#, stringsAsFactors=F)

            # Define solutions clusters
            # #########################

            used_clusters_name <- "X" # "MST", "W"

            used_clusters <- all_clusters[, used_clusters_name]

            unique_clusters <- sort(unique(used_clusters))

            MST_ids_from_MSTcluster_id <- lapply(unique_clusters, function(i){which(used_clusters == i)})
            names(MST_ids_from_MSTcluster_id) <- unique_clusters

            MSTcluster_sizes <- unlist(lapply(MST_ids_from_MSTcluster_id, length))

            num_graphclust <- length(unique(used_clusters))

            # Top scores of MSTclusters tops
            # ##############################
            
            print("Identify best prototype in each cluster")

            TopScore_MST_ids_from_MSTcluster_id <- unlist(lapply(MST_ids_from_MSTcluster_id, function(x){
                                x[which.max(GA_solutions_selection_score[x])]
                              }))

            SelectedTopScore_MST_ids_from_MSTcluster_id <- unlist(lapply(TopScore_MST_ids_from_MSTcluster_id, function(x){
                if( GA_solutions_selection_score[x]>=min_threshold ) {
                  x
                } else {
                  NA
                }
              }))

            print(paste0("   Number of best prototypes: ", sum(!is.na(SelectedTopScore_MST_ids_from_MSTcluster_id))))

            TopScore_MSTclusters_gms <-  t(ps[SelectedTopScore_MST_ids_from_MSTcluster_id,, drop = FALSE])
            colnames(TopScore_MSTclusters_gms) <- unique_clusters

            # Average of cluster weights
            # ##########################

            print("Identify best clusters (by average)")

            # We average the solutions with sufficient scores

            MSTclusters_Ws_average <- do.call(
              cbind, 
              lapply(
                MST_ids_from_MSTcluster_id,
                function(l){
                  l_filtered <- l[which(GA_solutions_selection_score[l] >= min_threshold)]
                  colMeans(ps[l_filtered, , drop = FALSE])}))

            MSTclusters_Ws_average <- apply(MSTclusters_Ws_average, 2, function(x){x/sum(x)})

            # Get objective scores
            MSTclusters_average_scores <- do.call(rbind, lapply(seq(num_graphclust), function(id){

              W <- MSTclusters_Ws_average[, id]

              if(any(is.na(W))){
                rep(NA, length(optim_objectives))
              } else {
                newScoreFunctionWeightedDistance(
                  W,
                  objs=optim_objectives,
                  compact_data=compact_data,
                  allow_multiple_roots=allow_multiple_roots,
                  root_cluster_id=root_cluster_id,
                  cluster_ids_from_cell_id=cluster_ids_from_cell_id,
                  ordering=ordering,
                  gm_smoothnesses_scaling_factors=gm_smoothnesses_scaling_factors,
                  compact_data_bin=compact_data_bin)
              }

            }))

            # Filter averaged Ws by score
            MSTclusters_average_scores_norm <- do.call(cbind, lapply(seq(num_obj), function(i){
                      if(pf_range[1, i]!=pf_range[2, i]){
                        (-MSTclusters_average_scores[, i]-pf_range[1, i]) / (pf_range[2, i]-pf_range[1, i])
                      } else {
                        rep(1, nrow(MSTclusters_average_scores))
                      }
              }))

            MSTclusters_average_selection_score <- apply(MSTclusters_average_scores_norm[, selection_objs_ids, drop = FALSE], 1, geo_mean)

            selected_cluster_ids <- which(MSTclusters_average_selection_score >= min_threshold)

            print(paste0("   Number of best averaged-W clusters: ", length(selected_cluster_ids)))

            if(length(selected_cluster_ids) == 0){

              print("   We increase the number of clusters to reveal best averaged-W clusters")

              num_graphclust_MSTs <- num_graphclust_MSTs + 1

              num_graphclust_Ws <- num_graphclust_Ws + 1

              loop_id <- loop_id + 1

            } else {
              goon <- FALSE
            }

            make_all_loop_plots(
              paste0(plot_folder_path, '/Loop_', loop_id),
              antler_env$expressionSet,
              pf_norm,
              pareto_front,
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
              GA_solutions_selection_quantile_threshold,
              TopScore_MSTclusters_gms,
              cell_ids_from_cluster_id,
              TopScore_MST_ids_from_MSTcluster_id,
              W_optim,
              optim_gms,
              dataset_gms,
              training_cell_ids,
              compact_data,
              MSTclusters_Ws_average,
              selected_cluster_ids,
              MSTclusters_average_scores,
              forbidden_cluster_edges,
              extra_pd             = extra_pd,
              allow_multiple_roots = allow_multiple_roots,
              color_maps           = antler_env$color_maps)

          } # end of num cluster incremental loop

          # Final solution is top cluster score

          final_clust_id <- which.max(MSTclusters_average_selection_score)

          print(paste0("Final cluster is ", final_clust_id))

          top_Ws <- MSTclusters_Ws_average[, final_clust_id]

          top_scores <- MSTclusters_average_scores[final_clust_id, ]

          print(paste0("Top scores from NSGA-II: ", paste0(top_scores, collapse=" / ")))

        } # end of "average" selection strategy
        
      } else if (optim_strategy == 'SA') {

        cat("Start SA optimization\n")

        sa_res <- GenSA::GenSA(
          par                             = NULL,
          fn                              = newScoreFunctionWeightedDistance,
          lower                           = rep(0, length(gene_modules)),
          upper                           = rep(1, length(gene_modules)),
          control                         = SA_control,
          objs                            = optim_objectives,
          compact_data                    = compact_data,
          allow_multiple_roots            = allow_multiple_roots,
          root_cluster_id                 = root_cluster_id,
          cluster_ids_from_cell_id        = cluster_ids_from_cell_id,
          ordering                        = ordering,
          gm_smoothnesses_scaling_factors = gm_smoothnesses_scaling_factors,
          compact_data_bin                = compact_data_bin)

        top_Ws = sa_res$par

        top_scores = sa_res$value

        plot_SA_history(
          sa_res$trace.mat,
          optim_objectives,
          paste0(plot_folder_path, "/SA_history.pdf"))
      }

      # Plot Best MST for current batch
      # ###############################

      top_Ws_norm <- top_Ws / sum(top_Ws)
      Final_cluster_cluster_distance <- newClusterWeightedDistance(top_Ws_norm, compact_data)

      top_all_graph <- igraph::graph_from_adjacency_matrix(
          as.matrix(Final_cluster_cluster_distance),
          weighted = TRUE,
          mode     = "upper")

      top_mst <- igraph::minimum.spanning.tree(top_all_graph, algorithm = 'prim')

      V(top_mst)$label <- as.character(seq(vcount(top_mst))) # fix error in rgexf::write.gexf 

      make_final_plot(
        paste0(plot_folder_path, '/Top'),
        antler_env$expressionSet,
        cell_ids_from_cluster_id,
        top_mst,
        compact_data, 
        extra_pd             = extra_pd,
        allow_multiple_roots = allow_multiple_roots,
        color_maps           = antler_env$color_maps)

      # Polish top solution with gradient descent

      if (polish) {

        print("Polish top solution with gradient descent")
        gd_res <- optim(
          par                             = top_Ws_norm,
          fn                              = newScoreFunctionWeightedDistance,
          objs                            = "timecorr",
          compact_data                    = compact_data,
          allow_multiple_roots            = allow_multiple_roots,
          root_cluster_id                 = root_cluster_id,
          cluster_ids_from_cell_id        = cluster_ids_from_cell_id,
          ordering                        = ordering,
          gm_smoothnesses_scaling_factors = gm_smoothnesses_scaling_factors,
          compact_data_bin                = compact_data_bin,
          method                          = "L-BFGS-B",
          lower                           = 0,
          upper                           = 1
        )

        if (gd_res$convergence == 0) {

          top_scores <- gd_res$value
          print(paste0("Top scores from L-BFGS-B: ", paste0(top_scores, collapse=" / ")))

          top_Ws_norm <- gd_res$par / sum(gd_res$par)
          Final_cluster_cluster_distance <- newClusterWeightedDistance(top_Ws_norm, compact_data)

          top_all_graph <- igraph::graph_from_adjacency_matrix(
              as.matrix(Final_cluster_cluster_distance),
              weighted=TRUE,
              mode="upper")

          top_mst <- igraph::minimum.spanning.tree(top_all_graph, algorithm='prim')

          V(top_mst)$label <- as.character(seq(vcount(top_mst))) # fix error in rgexf::write.gexf 

          make_final_plot(
            paste0(plot_folder_path, '/TopGD'),
            antler_env$expressionSet,
            cell_ids_from_cluster_id,
            top_mst,
            compact_data, 
            extra_pd             = extra_pd,
            allow_multiple_roots = allow_multiple_roots,
            color_maps           = antler_env$color_maps)
        } else {
          print(gd_res$message)
        }
      }

      # Store run variables
      runVariables[["cluster_ids_from_cell_id"]] <<- cluster_ids_from_cell_id
      runVariables[["cell_ids_from_cluster_id"]] <<- cell_ids_from_cluster_id
      runVariables[["compact_data"]]             <<- compact_data
      runVariables[["cluster_ordering"]]         <<- cluster_ordering
      runVariables[["top_scores"]]               <<- top_scores
      runVariables[["top_mst"]]                  <<- top_mst

      if (optim_strategy == "GA") {
        runVariables[["GA_pareto_scores"]]      <<- aggregated_scores
        runVariables[["GA_pareto_parameters"]]  <<- aggregated_parameters
      } else if (optim_strategy == "SA") {
        # ...        
      }

      # Update learned Weights
      # ######################

      W_optim <<- top_Ws_norm
      
      cat(paste0("\nFinal score(s): ", paste0(top_scores, collapse=' '), '\n'))

      cat(paste0("\nCarving completed.\n"))

    }

  )
)

