
# Allow igraph S3 class to be used as field of CellStateGraph
setOldClass("igraph")

CellStateGraph <- setRefClass(
  Class = "CellStateGraph",

  fields = list(
    antler_env                 = "environment",
    dist_matrix                = "matrix",
    graph                      = "igraph",
    graph_adjacency_matrix     = "matrix",
    build_method               = "character",
    settings                   = "list",
    layout                     = "matrix",
    layout_method              = "character",
    community_detection_method = "character",
    communities                = "numeric",
    communities_number         = "integer",
    score_MSTness              = "numeric",
    communities_graph          = "igraph",
    communities_df             = "data.frame",
    communities_cleantree      = "igraph",
    communities_mstree         = "igraph",
    root_community_id          = "integer"),

  methods = list(

    initialize = function(
      ...,
      antler_env                 = new.env(),
      dist_matrix                = matrix(0),
      graph                      = igraph::make_ring(0),
      graph_adjacency_matrix     = matrix(0),
      build_method               = '',
      settings                   = list(),
      layout                     = matrix(0),
      layout_method              = character(0),
      community_detection_method = character(0),
      communities                = numeric(0),
      communities_number         = integer(0),
      score_MSTness              = numeric(0),
      communities_graph          = igraph::make_ring(0),
      communities_df             = data.frame(),
      communities_cleantree      = igraph::make_ring(0),
      communities_mstree         = igraph::make_ring(0),
      root_community_id          = integer(0)) {

        callSuper(
          ...,
          antler_env                 = antler_env,
          dist_matrix                = dist_matrix,
          graph                      = graph,
          graph_adjacency_matrix     = graph_adjacency_matrix,
          build_method               = build_method,
          settings                   = settings,
          layout                     = layout,
          layout_method              = layout_method,
          community_detection_method = community_detection_method,
          communities                = communities,
          communities_number         = communities_number,
          score_MSTness              = score_MSTness,
          communities_graph          = communities_graph,
          communities_df             = communities_df,
          communities_cleantree      = communities_cleantree,
          communities_mstree         = communities_mstree,
          root_community_id          = root_community_id)
      },

    build = function(
      gene_modules_name   = NULL,
      use_distance        = FALSE,
      data_status         = c('Raw', 'Normalized', 'Smoothed'),
      method              = c('kNN', 'rMST'),
      numRand             = NA,
      maxDegree           = NA,
      RandomShotRatio     = NA,
      removedClass        = 1,
      allowed_adjacency   = NA,
      numcores            = antler_env$num_cores,
      k                   = NA,
      debug_plot_basename = NULL){

      build_method <<- match.arg(method)

      if (!is.null(gene_modules_name))
        gene_modules_name <- match.arg(
          gene_modules_name,
          antler_env$gene_modules$names())

      if ((is.null(gene_modules_name) & !use_distance) | (!is.null(gene_modules_name) & use_distance))
        stop(
          "The cell state graph is calculated via two mutually exclusive modes:\n",
          "1. from the readcounts by specifying the gene module list used to calculate cell-cell distances ('gene_modules_name' parameter), or\n", 
          "2. directly from the cell-cell distance matrix with 'use_distance' set to TRUE. The distance matrix must be specified beforehand.")

      # calculate distance matrix if gene module list is provided
      if (!use_distance) {

        used_genes.list <- unlist(antler_env$gene_modules$get(gene_modules_name))

        if (length(used_genes.list) == 0)
          stop("The specified gene modules list is empty.")

        cell_cell_distance.mat <- fastEuclideanDist(
          t(log(antler_env$read_counts(data_status = data_status) + 1)[used_genes.list, ]))

        cell_cell_distance.mat <- cell_cell_distance.mat / sum(cell_cell_distance.mat)

        dist_matrix <<- cell_cell_distance.mat

        antler_env$processing_state <<- c(
          antler_env$processing_state,
          paste0(
            'Calculated cell-cell distance matrix with Euclidean distance applied on the log-levels of genes belonging to "', gene_modules_name,
            '"" (length: ', length(used_genes.list), ')'))
      }

      if (ncol(dist_matrix) == 0 | nrow(dist_matrix) == 0)
        stop("missing or invalid input distance matrix")

      settings$max_degree   <<- as.integer(maxDegree)

      if (build_method == "rMST") {
        
        settings$rMST_randomShotRatio <<- if(is.na(RandomShotRatio)){.2}else{RandomShotRatio}
        settings$rMST_num             <<- as.integer(numRand)
        settings$rMST_removedClass    <<- as.integer(removedClass)

        graph_adjacency_matrix <<- get_rMST_adjacency_matrix_edge_based(
          dist_matrix,
          RandomShotRatio     = settings$rMST_randomShotRatio,
          numRand             = settings$rMST_num,
          plot                = FALSE,
          removedClass        = settings$rMST_removedClass,
          allowed_adjacency   = allowed_adjacency,
          numcores            = numcores,
          debug_plot_basename = debug_plot_basename)

        antler_env$processing_state <<- c(
          antler_env$processing_state,
          paste0(
            fdate(), ' Create cell state graph (method: ', method,
            ', numRand: ', numRand, ', maxDegree: ', maxDegree,
            ', RandomShotRatio: ', RandomShotRatio,
            ', removedClass: ', removedClass, ') '))

      } else if (build_method == 'kNN') {

        settings$kNN_k <<- if(is.na(k)){3}else{as.integer(k)}
        
        graph_adjacency_matrix <<- get_kNN_adjacency_matrix(
          dist_matrix,
          k                 = settings$kNN_k,
          allowed_adjacency = allowed_adjacency)

        antler_env$processing_state <<- c(
          antler_env$processing_state,
          paste0(
            fdate(), 'Created cell state graph (method: ', method, ', k: ',
            k, ') '))
      }

      if (!is.na(settings$max_degree)) {
        graph_adjacency_matrix <<- pruneAdjacencyMatrix(
          graph_adjacency_matrix,
          MaxDegree = settings$max_degree)
      }

      graph <<- get_graph_from_adjacency_matrix(dist_matrix, graph_adjacency_matrix)
    },

    project = function(
      method         = c("forceatlas2", "fruchterman_reingold", "gephiForceAtlas2"),
      num_iter       = 2000,
      stop_mean      = 1,
      seed           = 0,
      gephi_jar_path = system.file(
        "java",
        "LayoutGephi.jar",
        package="Antler"),
      verbose        = FALSE) {

      method <- match.arg(method)

      set.seed(seed)

      layout_method <<- method

      if (method == 'forceatlas2') {
        layout <<- layout.forceatlas2(
          graph,
          directed   = TRUE,
          iterations = num_iter,
          linlog     = FALSE,          # default
          pos        = NULL,           # default
          nohubs     = FALSE,          # default
          k          = 1000,
          gravity    = 0,
          ks         = 0.1,            # default
          ksmax      = 10,             # default
          delta      = .1,
          center     = NULL,           # default
          tolerance  = 0.1,            # default
          dim        = 2,              # default
          plotstep   = 0,              # > 0 for live plots
          plotlabels = FALSE)
      } else if (method == 'fruchterman_reingold') {
        layout <<- layout.fruchterman.reingold(
          graph,
          niter   = num_iter,
          weights = rep(1, length(E(graph))))
      } else if (method == 'gephiForceAtlas2') {

        tdir          <- tempdir()
        gexf_filepath <- paste0(tdir, '/current_graph.gexf')
        csv_filepath  <- paste0(tdir, '/current_graph.csv')

        # print('Write gexf file')
        export_graph_as_gexf(graph, filepath = gexf_filepath)

        # print('Run java code')
        args <- paste0(
          "-jar ", gephi_jar_path, " ", gexf_filepath, " ",
          as.integer(num_iter), " ", stop_mean)

        if (verbose) {
          stdout <- ""
        } else {
          stdout <- NULL
        }

        system2(
          command = "java",
          args    = args,
          stdout  = stdout,
          stderr  = stdout)

        # print('Read layout csv')
        layout <<- as.matrix(read.table(
          csv_filepath, header = FALSE, sep = ','))

        # print('Remove temp files')
        file.remove(gexf_filepath)
        file.remove(csv_filepath)
      }
    },

    show = function(){},

    remove_cells = function(ids = NULL){
      
      dist_matrix            <<- dist_matrix[-cell_ids, -cell_ids]
      graph_adjacency_matrix <<- graph_adjacency_matrix[-cell_ids, -cell_ids]
      graph                  <<- igraph::delete_vertices(graph, cell_ids)
      layout                 <<- layout[-cell_ids,] # TODO: check if rownames matters
      
      initial_comm_ids   <- sort(unique(communities))
      communities        <<- communities[-cell_ids]
      remaining_comm_ids <- sort(unique(communities))
  
      # if some communities have disappeared, we re-order all of them

      if (communities_number != length(remaining_comm_ids)) {

        print(paste0(
          "Community(ies) ",
          paste0(setdiff(initial_comm_ids, remaining_comm_ids), collapse=', '),
          " are now empty. All communities are relabelled between 1 and ",
          length(remaining_comm_ids)))

        communities_number <<- length(remaining_comm_ids)

        comm_names         <- names(communities)
        communities        <<- setNames(
          seq(communities_number),
          sort(unique(communities)))[as.character(communities)]
        names(communities) <<- comm_names # keep old community labels in case they are used elsewhere (normally cell names)

        root_community_id  <<- setNames(
          seq(communities_number),
          sort(unique(communities)))[as.character(root_community_id)]

        communities_mstree <<- delete_vertices(
          communities_mstree, setdiff(initial_comm_ids, remaining_comm_ids))
        communities_graph  <<- delete_vertices(
          communities_graph, setdiff(initial_comm_ids, remaining_comm_ids))
        communities_cleantree <<- delete_vertices(
          communities_cleantree, setdiff(initial_comm_ids, remaining_comm_ids))

        communities_df <<- communities_df[remaining_comm_ids, , drop=F]
      }

    },


    plot_from_color_vector = function(
      fullpath           = paste0(
        paste0(sample(letters, 6), collapse=''),
        '_coloredPlot_', build_method, '_', layout_method, '.pdf'),
      sample_colors      = NA,
      labels             = NA,
      vertex.size        = 3,
      vertex.shape       = "circle",
      vertex.frame.color = 'black',
      edge.width         = vertex.size/5,
      edge.color         = 'grey',
      width              = 7,
      height             = 7,
      save_pdf            = TRUE,
      shuffle            = FALSE,
      vertex.plot.order  = NA) {

        if(save_pdf)
          pdf(fullpath, width = width, height = height)
        op <- par(mar = c(0,0,0,0))
        plotColorGraph2(
          graph,
          layout,
          shape              = vertex.shape,
          colors             = sample_colors,
          vertex.size        = vertex.size,
          label              = labels,
          vertex.frame.color = vertex.frame.color,
          edge.width         = edge.width,
          edge.color         = edge.color,
          shuffle            = shuffle,
          vertex.plot.order  = vertex.plot.order)
        par(op)
        if(save_pdf)
          silent <- dev.off()

      },

    detect_communities = function(
      method                       = 'spin_glass_no_weight',
      n_communities                = NULL,
      checkCommunitiesConnectivity = TRUE,
      seed                         = 1234) {

      # Some community detection methods does not work with multi-component graphs. We apply the methods on each subgraph independently.
      graph_comp <- components(graph)

      communities_l <- list()
      n_comm_total <- 0

      for (i in seq(graph_comp$no)) {

        subgraph <- delete.vertices(
          graph,
          setdiff(seq(vcount(graph)), which(graph_comp$membership == i)))

        communities_l[[i]] <- n_comm_total + get_communities_from_graph(
            subgraph,
            method                       = method,
            communities_number           = n_communities,
            checkCommunitiesConnectivity = checkCommunitiesConnectivity,
            seed                         = seed)

        n_comm_total <- max(communities_l[[i]])
      }

      communities                <<- unlist(communities_l)[V(graph)$name]
      communities_number         <<- length(unique(communities))
      community_detection_method <<- method
      communities_df             <<- data.frame(
        ID        = seq(length(unique(communities))),
        row.names = seq(length(unique(communities))))

      # print(paste0("Final number of communities: ", length(unique(communities))))

      # Store community id and tree id in expressionSet
      pData(antler_env$expressionSet)[, "subgraph_id"] <<- graph_comp$membership
      pData(antler_env$expressionSet)[, "community_id"] <<- communities
      
      communities_df$cell_ids <<- lapply(
        seq(communities_number),
        function(i) {which(pData(antler_env$expressionSet)$community_id == i)})
      
      communities_df$subgraph_id <<- lapply(
        communities_df$cell_ids,
        function(l) {unique(pData(antler_env$expressionSet)$subgraph_id[l])})

      antler_env$add_color_map(
        name = "community_id",
        content = brew_more_colors(seq(communities_number), "Set1"))
    },

    # requires dist_matrix, graph_adjacency_matrix and communities
    build_community_graph = function(
      min_occurence = 1) {

      # To define local manifolds, we need a cleaned community graph (no self loop, trimmed edges).
      # From the cell state graph, we calculate the number of edges relating each communities.
      communities_adjmat <- do.call(
        cbind,
        lapply(
          seq(communities_number),
          function(i){
            rowSums(
              graph_adjacency_matrix[, which(communities==i), drop=F])
          }))

      communities_adjmat <- do.call(
        rbind,
        lapply(
          seq(communities_number),
            function(i){
              colSums(
                communities_adjmat[which(communities==i), , drop=F])
            }))

      # ignore diagonal (within community edges)
      diag(communities_adjmat) <- 0

      # Optionally, low occurenc edges are removed
      communities_adjmat_trim_edges <- communities_adjmat
      communities_adjmat_trim_edges[communities_adjmat < min_occurence] <- 0

      # Set positive adjacency to 1
      communities_adjmat_single_edges = 1*(communities_adjmat_trim_edges > 0)

      # The "clean" graph is undirected
      communities_cleantree <<- igraph::graph_from_adjacency_matrix(
        communities_adjmat_single_edges,
        mode = "undirected")

      # The community "graph" stores the number of edges in-between the communities
      communities_graph <<- igraph::graph_from_adjacency_matrix(communities_adjmat,
        weighted = NULL,
        mode     = "upper")

      # The community "mstree" is calculated from the community adjacency matrix only, not from the actual community states. The weight between two communities is the inverse of the number of edges relating them, so that well connected communities are favored
      communities_graph_weighted <- igraph::graph_from_adjacency_matrix(
        communities_adjmat,
        weighted = TRUE,
        mode     = "upper")

      communities_mstree <<- igraph::minimum.spanning.tree(
        communities_graph_weighted,
        algorithm = 'prim',
        weights   = 1/E(communities_graph_weighted)$weight)

      # Also, calculate score_MSTness used in plots (for now ?)
      calculate_mstness()
    },

    # Count edges not in community MST / total edges in community graph. Do not take into account inner community edge counts.
    # requires communities_mstree, graph_adjacency_matrix and communities
    calculate_mstness = function(){

      communities_adjmat <- do.call(
        cbind,
        lapply(
          seq(communities_number),
          function(i){
            rowSums(
              graph_adjacency_matrix[, which(communities==i), drop=F])
          }))

      communities_adjmat <- do.call(
        rbind,
        lapply(
          seq(communities_number),
            function(i){
              colSums(
                communities_adjmat[which(communities==i), , drop=F])
            }))

      # ignore diagonal (within community edges)
      diag(communities_adjmat) <- 0

      communities_mstree_adjmat <- igraph::as_adjacency_matrix(
        communities_mstree,
        sparse = FALSE)

      num_edge_out_mst <- .5 * sum((1-communities_mstree_adjmat) * communities_adjmat)
      num_edge_in_mst <- .5 * sum(communities_mstree_adjmat * communities_adjmat)

      score_MSTness <<- num_edge_out_mst / num_edge_in_mst
    }

  )
)


#' Plot cell state graph 
#' 
#' @name cell_state_graph_plot
#' @param cell_colors character vector indicating the metadata to use to color the vertices. Names can be either one of the column names of the phenotypic metadata structure (see \code{colnames(pData(antler$expressionSet))}), names of a cluster entries, or gene names. The element of the character vector can be named themselves, in that case the element name will serve as a title for phenotypic metadata or cluster entries. If gene names are specified, their vector names must indicate the type of read counts that will be used to render the expression levels: "Raw", "Normalized" or "Smoothed". If no gene name's name is specified raw level will be rendered.
#' @param labels character vector. Labels to draw on top of the cells. Default to NA.
#' @param vertex.size numeric vector indicating the size of the vertices representing the cells.
#' @param vertex.shape character vector indicating the shape of the vertices. Can be either "circle" (default), "square" or "triangle".
#' @param vertex.frame.color character vector indicating the color of the frame of each vertex. Default to NA.
#' @param edge.width numeric value indicating the width of graph edges. Default to a fifth of the vertex size.
#' @param edge.color character string indicating the color of the edges. Default to "grey".
#' @param width numeric value indicating the width (in inches) of each individual subplots. Default to 3.
#' @param height numeric value indicating the height (in inches) of each individual subplots. Default to 3.
#' @param title logical value indicating whether to show each subplots' title. Default to TRUE.
#' @param save_pdf logical value indicating whether to save plot as a pdf file. Default to TRUE.
#' @param shuffle logical value indicating whether to randomize the z-order of each vertex. Default to FALSE.
#' @param vertex.plot.order integer vector. The z-order used to render the vertices. Default to NA.
#' @param gene_palette character string. The color palette to use if genes are specified as cell colors.
#' @param filename character string. The filename of the output file, if any.
#' @param ... Extra arguments passed to \code{\link[igraph]{plot.igraph}}.
CellStateGraph$methods(
  plot = function(
    cell_colors,
    labels             = NA,
    vertex.size        = 1000 / antler_env$num_cells(),
    vertex.shape       = "circle",
    vertex.frame.color = NA,
    edge.width         = .2 * vertex.size,
    edge.color         = 'grey',
    width              = 3,
    height             = 3,
    title              = TRUE,
    save_pdf           = TRUE,
    shuffle            = FALSE,
    vertex.plot.order  = NA,
    gene_palette       = colorRampPalette(
      c("#0464DF", "#FFE800"))(n = 1000),
    filename           = paste0(
      'Cell_state_graph_',
      build_method, '_',
      layout_method, '_',
      paste0(cell_colors, collapse = "-"),
      '.pdf'),
    ...) {

    if (!(all(
      cell_colors %in% colnames(pData(antler_env$expressionSet)) |
      cell_colors %in% antler_env$cell_clusters$names() |
      cell_colors %in% antler_env$gene_names())))
      stop("'cell_colors' must contain elements of the phenotypic metadata structure (see 'colnames(pData(antler$expressionSet))'), names of a cluster entries (see 'antler$cell_clusters$names()'), or gene names (see 'antler$gene_names()')")

    n_row <- as.integer(sqrt(length(cell_colors)))
    n_col <- ceiling(length(cell_colors) / n_row)

    if (save_pdf)
      pdf(
        file.path(antler_env$output_folder, filename), 
        width = n_col * width,
        height = n_row * height)

    if (length(cell_colors) > 1) {
      par(
        mfrow = c(n_row, n_col),
        mar   = c(.5, .5, 1, 3), xpd=TRUE)
    } else {
      par(
        mar   = c(.5, .5, 1, 3), xpd=TRUE)
    }

    for (i in seq(length(cell_colors))) {

      color_maps <- antler_env$color_maps
      
      values <- unlist(get_color_matrix(
        cell_colors[i],
        antler_env$expressionSet,
        antler_env$cell_clusters,
        antler_env$read_counts(
          ifelse(
            cell_colors[i] %in% antler_env$gene_names() &
              names(cell_colors[i]) %in% c("Raw", "Normalized", "Smoothed"),
            names(cell_colors[i]),
            "Raw"
          ))))

      if (!cell_colors[i] %in% names(color_maps)
          & is.double(values))
        color_maps[[cell_colors[i]]] <- gene_palette

      cols <- colorize(
        cell_colors[i],
        values,
        color_maps)

      plot_color_graph(
        graph,
        layout,
        shape              = vertex.shape,
        colors             = cols,
        vertex.size        = vertex.size,
        label              = labels,
        vertex.frame.color = if (identical(vertex.frame.color, NA)) cols else vertex.frame.color,
        edge.width         = edge.width,
        edge.color         = edge.color,
        shuffle            = shuffle,
        vertex.plot.order  = vertex.plot.order,
        ...)

      if (title) {
        if (is.null(names(cell_colors)[i]) |
            identical(names(cell_colors)[i], "")) {
          title(cell_colors[i])
        } else {
          if (cell_colors[i] %in% antler_env$gene_names()) {
            title(paste0(cell_colors[i], " (", names(cell_colors)[i], ")"))
          } else {
            title(names(cell_colors)[i])
          }
        }
      }

      if (length(unique(values)) <= 20 &
          !is.double(values))
        legend(
          'right',
          inset   = c(-0.1, 0),
          legend  = sort(unique(values)),
          col     = cols[as.character(sort(unique(values)))],
          pt.bg   = cols[as.character(sort(unique(values)))],
          pt.cex  = 0.5,
          box.col = "white",
          cex     = 0.5,
          pch     = 21,
          ncol    = 1)
    }

    if (save_pdf)
      silent <- dev.off()

  }
)
