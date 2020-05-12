
CellClusters <- setRefClass( 

  Class = "CellClusters",

  fields = list(
    antler_env = "environment",
    lists      = "list"
  ),

  methods = list(

    initialize = function(
      ...,
      antler_env = new.env(),
      lists      = list()) {

      callSuper(
        ...,
        antler_env = antler_env,
        lists      = lists)
    },

    # add 'base::' to get function to disambiguate with class get
    copy = function(shallow = FALSE) {
      def <- .refClassDef
      value <- new(def)
      vEnv <- as.environment(value)
      selfEnv <- as.environment(.self)
      for (field in base::names(def@fieldClasses)) {
        current <- base::get(field, envir = selfEnv)
        if (is(current, "envRefClass")) 
            current <- current$copy(FALSE)
        assign(field, current, envir = vEnv)
      }
      value
    },

    get = function(name) {

      if (missing(name))
        stop("The argument 'name' is required")

      if (!is.character(name))
        stop("The argument 'name' must be a character string")

      if (!name %in% base::names(lists))
        stop(paste0(
          "There is no '", name, "' cell clusters. Available name(s): '",
          paste0(base::names(lists), collapse= "', '"), "'"))

      return(lists[[name]][['content']])
    },

    set = function(name, content) {

      if (missing(name))
        stop("The argument 'name' is required")

      if (!is.character(name))
        stop("The argument 'name' must be a character string")

      if (!all(unlist(content) %in% antler_env$gene_names())) {
        stop("Some gene names in 'content' are not available in the Antler structure")
      }

      lists[[name]] <<- list("content" = content)

    },

    names = function() {
      return(base::names(lists))
    }
  )
)


#' Cluster cells from a selected gene module list.
#' 
#' A text file containing the gene modules is written to the Antler output folder.
#' @name cell_clusters_identify
#' @param name character string. The name of the cell clusters object. Multiple cell clusters objects under different names can be stored.
#' @param gene_modules_name character string. The name of the list of gene modules from which the cell-to-cell distance will be calculated.
#' @param method character string. The clustering method to use. Currently only \code{\link[stats]{hclust}} is available (using "ward.D2" as agglomeration method).
#' @param cell_cell_dist numeric matrix. Optionally a matrix containing the cell-to-cell distances can be provided. If NA (default), the matrix is calculated internally using Euclidean distance.
#' @param logscaled logical. If TRUE (default), the gene expression level are log-transformed and z-normalized before calculating the cell-to-cell distance. If FALSE, no transformation is applied to the levels.
#' @param num_clusters an integer value indicating the number of cell clusters. If NULL (default) a heuristic method is used to find the optimal number of clusters.
#' @param data_status character string. Specifies whether the gene expression levels to be used are raw ("Raw"), normalized ("Normalized") or have been imputed ("Smoothed"). Default to "Normalized".
#' @param process_plots a logical value indicating whether to render the heuristic trade-off determining the number clusters. Default to FALSE.
CellClusters$methods(

  identify = function(
    name,
    gene_modules_name  = NULL,
    method             = 'hclust',
    cell_cell_dist     = NA,
    logscaled          = TRUE,
    num_clusters       = NULL,
    data_status        = c('Normalized', 'Raw', 'Smoothed'),
    process_plots      = FALSE) {

    if (missing(name))
      stop("The argument 'name' is required")

    if (!is.character(name))
      stop("The argument 'name' must be a character string")

    lists[[name]] <<- list("method" = method, "content" = list())

    method        <- match.arg(method)
    data_status   <- match.arg(data_status)

    if (! gene_modules_name %in% antler_env$gene_modules$names())
      stop("invalid gene modules name (must be one of the following: '", paste0(antler_env$gene_modules$names(), collapse="', '"), "')")

    data.1 <- antler_env$read_counts(data_status)[unlist(antler_env$gene_modules$get(gene_modules_name)), ]
    gene_modules_name.report <- paste0('the reduced dimensions "', gene_modules_name, '"" (length: ', nrow(data.1), ')')


    if (logscaled){
      data.2 <- t(scale(
        t(log(data.1 + 1)),
        center = TRUE,
        scale = TRUE))

      # Check whether data.2 contains NaN
      na_genes = which(is.na(rowSums(data.2)))
      num_na = length(na_genes)
      if (num_na > 0){
        print(paste0("Warning: ", num_na, " genes contains NaN values after logscaled transformation (", paste0(rownames(data.2)[na_genes], collapse=', '), "). It is most likely null genes which should be removed upstream in the pipeline. These genes are ignored from current cell-cell distance calculation."))
        data.2[na_genes, ] <- 0
      }

    } else {
      data.2 = log(data.1+1)
    }

    if (all(is.na(cell_cell_dist))){
      cell_cell_dist = as.dist(fastEuclideanDist(t(data.2)))
    }

    if (method == 'hclust') {
      hc.cells <- hclust(cell_cell_dist, method = "ward.D2")

      if (is.null(num_clusters)){
        num_clusters <- get_optimal_cluster_number(
          t(data.2),
          hc.cells,
          100000,
          1,
          display       = FALSE,
          pdf_plot_path = if (process_plots) {paste0(antler_env$output_folder, '/Cell_clustering_heuristic_', name, '.pdf')} else {NULL},
          method = "piecewise"
          )
      }

      hc.cells.clust_ids = cutree(hc.cells, k = num_clusters)

      # We order the cluster id from left to right
      clust_ids_ordered = unique(hc.cells.clust_ids[hc.cells$order]) %>% setNames(seq(num_clusters), .)

      lists[[name]][['content']]$res          <<- hc.cells
      lists[[name]][['content']]$cell_ids     <<- hc.cells.clust_ids %>% {setNames(clust_ids_ordered[as.character(.)], base::names(.))}
      lists[[name]][['content']]$gene_modules_name   <<- gene_modules_name
      lists[[name]][['content']]$gene_levels  <<- logscaled
      lists[[name]][['content']]$num_clusters <<- num_clusters
      lists[[name]][['content']]$names        <<- NA

      pData(antler_env$expressionSet)[, name]  <<- lists[[name]][['content']]$cell_ids

      antler_env$add_color_map(
        name = name,
        content = setNames(
          colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(num_clusters),
          seq(num_clusters)))
    }

    antler_env$processing_state <<- c(
      antler_env$processing_state,
      paste0(
        fdate(),
        ' Identify ', num_clusters,
        ' cell clusters (method,', method, ' from the Euclidean distance matrix generated with ',
        ifelse(logscaled, 'the z-scored log-levels of ', ''), gene_modules_name.report))

  }
)
