
GeneModules <- setRefClass( 

  Class = "GeneModules",

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

    # saveToCSV = function(){
    #   stop('TODO: csv export')
    # },

    # whichModuleContains = function(gene_name){

    #   for(used_gm in grep('genemodules', names(.self$.refClassDef@fieldClasses), value = TRUE)) {

    #     if( length(.self[[used_gm]]) > 0 ){

    #       mod_id = Filter(function(i){gene_name %in% .self[[used_gm]][[i]]},
    #              seq(length(.self[[used_gm]]))
    #              )

    #       cat(paste0('Module id in "', used_gm, '" -> ', mod_id, '\n'))
    #     }
    #   }
    # },

    get = function(name) {

      if (missing(name))
        stop("The argument 'name' is required")

      if (!is.character(name))
        stop("The argument 'name' must be a character string")

      if (!name %in% base::names(lists))
        stop(paste0(
          "There is no '", name, "' gene module list. Available names: '",
          paste0(base::names(lists), collapse= "', '"), "'"))

      return(lists[[name]][['content']])
    },

    set = function(name, content) {

      if (missing(name))
        stop("The argument 'name' is required")

      if (!is.character(name))
        stop("The argument 'name' must be a character string")

      if (!all(unlist(content) %in% antler_env$gene_names())) {
        missing_genes <- setdiff(
          unlist(content),
          antler_env$gene_names())
        stop(paste0(
          "The following gene names are not available in the current project: ",
          paste0(missing_genes, collapse = ", ")))
      }

      lists[[name]] <<- list("content" = content)

    },

    names = function() {
      return(base::names(lists))
    }
  )
)


#' Save gene module to file.
#' 
#' A text file containing the gene modules is written to the Antler output folder.
#' @name gene_modules_export
#' @param name character string. The name of the list of gene modules that will be exported.
GeneModules$methods(

  export = function(name = names()) {

    name <- match.arg(name)

    gm_filename = paste0(antler_env$output_folder, '/', name, '.txt')
    
    file.create(gm_filename)
    current_gms <- get(name)

    for(i in seq(length(current_gms))) {

      modname = base::names(current_gms)[i]
      if (is.null(modname)){
        modname = paste0('GM', i)
      }

      write(
        paste0(
          modname,
          ': ',
          paste0(
            current_gms[[i]],
            collapse=', ')),
        file   = gm_filename,
        append = TRUE)
    }
  }
)

#' Identify gene modules.
#' 
#' This method groups genes demonstrating concerted patterns of expression. 
#' TODO details from methods
#' @name gene_modules_identify
#' @param name character string. The name of the list of gene modules. Multiple lists of gene modules under different names can be stored.
#' @param corr_matrix numeric matrix. Specifies an optional correlation matrix. If not provided a Spearman correlation matrix will be calculated by the function.
#' @param corr_t numeric value. Correlation cutoff value or correlation quantile cutoff value used to select genes. Either \code{corr_t} or \code{corr_num_genes} must be specified.
#' @param corr_num_genes integer value. Approximate number of genes to select in the correlation filtering step. Either \code{corr_t} or \code{corr_num_genes} must be specified. 
#' @param corr_quantile logical. If FALSE (default), \code{corr_t} specifies the correlation threshold. If TRUE, \code{corr_t} is the percentile correlation threshold over all \code{corr_matrix} correlation values.
#' @param corr_min integer value. The minimum number of high gene- gene correlation values a gene needs to have to be kept. Default to 3.
#' @param mod_min_cell an integer value indicating the minimum number of expressed cells a gene module needs to have to be selected. Default to 0.
#' @param mod_max_cell an integer value indicating the maximum number of expressed cells a gene module needs to have to be selected. Default to Inf.
#' @param mod_consistency_thres numeric value. If in a gene module the ratio of expressed genes in the "positive" cells (see Details) is less than this value, the gene module is excluded. Default to 0.4.
#' @param min_cell_level numeric value. If in a gene module the average expression level of the "positive" cells is less than this value, the gene module is excluded. Default to 0.
#' @param mod_skewness_thres numeric value. In each gene module the skewness of the cells' average expression level is calculated. If the skewness is lower than this value, the gene module is excluded. Default to -Inf.
#' @param ordering_correlation_thres numeric value. if the Spearman correlation between the cells' average expression level of a gene module and the phenotypic metadata specified by the \code{ordering} parameter (usually the cells' timepoint) is lower than this value. The gene module is excluded. If equal to NA (Default), this filtering step is not applied.
#' @param ordering a character string indicating the phenotypic metadata to use for the ordering correlation filtering step.
#' @param num_final_gms integer value. Constraints the resulting number of gene modules after applying the filtering steps. If NULL (Default), a heuristic method is used to find the optimal number of gene modules.
#' @param num_max_final_gms integer value. Caps the maximum number of resulting gene modules. If NULL (Default), no capping is applied.
#' @param num_min_final_gms integer value. Caps the minimum number of resulting gene modules. If NULL (Default), no capping is applied.
#' @param num_initial_gms integer value. Constraints the initial number of clustered gene modules, before applying the filtering steps. If NULL (Default), a heuristic method is used to find the optimal number of gene modules.
#' @param clustering_method character string. The agglomeration method to be used by the hierarchical clustering algorithm (\code{\link[stats]{hclust}}). Default to "ward.D2".
#' @param display logical. Whether to display intermediate processing plots. Default to FALSE.
#' @param num_cores. an integer value indicating the number of cores to use.
#' @param verbose a logical value indicating whether live processing information should be printed. Default to TRUE.
#' @param verbose_final a logical value indicating whether processing information should be printed just before the function returns. Default to FALSE.
#' @param process_plots a logical value indicating whether to render the heuristic trade-off determining the number of gene modules for each iteration. Default to FALSE.
GeneModules$methods(

  identify = function(
    name,
    corr_matrix                = NA,
    corr_t                     = NA,
    corr_num_genes             = NA,
    corr_quantile              = FALSE,
    corr_min                   = 3,
    mod_min_cell               = 0,
    mod_max_cell               = Inf,
    mod_consistency_thres      = .4,
    min_cell_level             = 0,
    mod_skewness_thres         = -Inf,
    ordering_correlation_thres = NA,
    ordering                   = NA,
    num_final_gms              = NULL,
    num_max_final_gms          = NULL,
    num_min_final_gms          = NULL,
    num_initial_gms            = NULL,
    clustering_method          =  c("ward.D2", "ward.D2", "single", "complete", "average", "mcquitty", 
        "median", "centroid"),
    display                    = FALSE,
    num_cores                  = antler_env$num_cores,
    verbose                    = TRUE,
    verbose_final              = FALSE,
    process_plots              = FALSE) {

    if (missing(name))
      stop("The argument 'name' is required")

    if (!is.character(name))
      stop("The argument 'name' must be a character string")

    if (identical(antler_env$readcounts_norm, matrix()))
      stop("Gene levels must be normalized before running this function")

    lists[[name]] <<- list("content" = list(), "history" = list())

    report <- character(0)

    clustering_method <- match.arg(clustering_method)

    # Identify correlated genes
    ###########################

    # cat('Select correlated genes...')
    if (is.na(corr_num_genes) + is.na(corr_t) != 1)
      stop('either "corr_num_genes" or "corr_t" must be specified, not both.')

    if (identical(corr_matrix, NA)) {
      # if (verbose){cat("Calculate gene-gene correlation matrix\n")}
      used.genes <- antler_env$gene_names()
      
      corr_matrix <- fastCor(t(antler_env$read_counts('Normalized')), method = "spearman")
    
    } else {
      used.genes <- rownames(corr_matrix)
    }

    if(is.na(corr_num_genes) & !is.na(corr_t)){

      if(corr_quantile){
        correlation_thres <- quantile(abs(corr_matrix), prob=corr_t)
      } else {
        correlation_thres <- corr_t
      }
    } else if(!is.na(corr_num_genes) & is.na(corr_t)){
      correlation_thres <- getCorThreshold(corr_matrix, init_corr_t=.2, numgenes_target=corr_num_genes, corr_min=corr_min, numcores=num_cores)
    }

    corr_matrix_thres <- abs(corr_matrix) >= correlation_thres
    corr_genes <- base::names(which(colSums(corr_matrix_thres) > corr_min))

    rm(corr_matrix_thres)
    clean_memory()
    
    verbose.1 <- paste0('Select correlated genes: ', length(corr_genes),' genes out of ', length(used.genes),' (threshold: ', sprintf("%.3f", correlation_thres),')')
    
    if(verbose){cat(verbose.1)}
    report <- paste0(report, verbose.1)

    # 1. optimal number of gene modules
    ###################################

    data.logscaled <- t(scale(t(log(antler_env$read_counts('Normalized')[corr_genes,]+1)), center=T, scale=T))

    clusterGenes_res <- clusterGenes(
      antler_env$read_counts('Normalized')[corr_genes,],
      data.logscaled[corr_genes, ],
      corr_matrix[corr_genes, corr_genes],
      num_gms           = num_initial_gms,
      num_gms_min       = NULL,
      num_gms_max       = NULL,
      numcores          = num_cores,
      display           = display,
      basename          = if (process_plots) {paste0(antler_env$output_folder, '/Gene_module_identification_', name, '_heuristic_loop_0.pdf')} else {NULL},
      verbose           = verbose,
      clustering_method = clustering_method)

    lists[[name]][['content']] <<- clusterGenes_res$clusters
    report <- paste0(report, clusterGenes_res$verbose.log)

    if(length(lists[[name]][['content']])==0){
      lists[[name]][['content']] <<- list()
      return()
    }

    lists[[name]][['history']][['corrgenes']] <<- lists[[name]][['content']]

    # ############################################
    # 2. Apply quality filters on the gene modules
    ##############################################

    # applying gene module filtering on Normalized filters (not on smoothed ones)
    data_normalized.logscaled = t(scale(t(log(antler_env$read_counts('Normalized')[corr_genes,]+1)), center=T, scale=T))
    
    if(mod_consistency_thres > 0){
      # cat("\nStart binarizing... ")
      data_normalized.bin = silent.binarize.array(data_normalized.logscaled) # no need to normalize non-corrgenes
      # cat("done\n")
      } else {data_normalized.bin = NULL}

    goon <- TRUE
    loop_id <- 0

    while(goon){

      genemodules.pre_filters <- lists[[name]][['content']]

      mod_avg = do.call(rbind, lapply(genemodules.pre_filters, function(m){colMeans(data_normalized.logscaled[m,, drop=F])}))

      mod_avg.bin = silent.binarize.array(mod_avg)

      if(display){dev.new(); heatmap.plus(mod_avg.bin, scale='none', Rowv=NA, Colv=NA, labCol = NA)}

      current_gm_id_list = seq(length(genemodules.pre_filters))

      # Criteria: Keep modules expressed in at least N cells in Unsmoothed (ie Normalized) readcounts
      # #############################################################################################

      if(!(identical(mod_min_cell, 0) & identical(mod_max_cell, Inf))) {

        select_res = selectModulesByPositiveCellRange(current_gm_id_list, mod_avg.bin, mod_min_cell, mod_max_cell, verbose)

        current_gm_id_list = select_res$gm_ids

        lists[[name]][['content']] <<- genemodules.pre_filters[current_gm_id_list]
        report <- paste0(report, select_res$verbose.log, ' (', length(unlist(genemodules.pre_filters[current_gm_id_list])), ' genes)')

        if(length(lists[[name]][['content']])==0){
          lists[[name]][['content']] <<- list()
          return()
        }

        lists[[name]][['history']][[paste0("Loop_", loop_id, "_mincells")]] <<- lists[[name]][['content']]
      }
      

      # Criteria: Keep modules with consistent expression in positive cells
      # ###################################################################

      if(mod_consistency_thres > 0){
        select_res = selectModulesByConsistency(genemodules.pre_filters, current_gm_id_list, mod_avg.bin, data_normalized.bin, consistency_thres=mod_consistency_thres, verbose=verbose)
        
        current_gm_id_list = select_res$gm_ids

        lists[[name]][['content']] <<- genemodules.pre_filters[current_gm_id_list]
        report <- paste0(report, select_res$verbose.log, ' (', length(unlist(genemodules.pre_filters[current_gm_id_list])), ' genes)')

        if(length(lists[[name]][['content']])==0){
          lists[[name]][['content']] <<- list()
          return()
        }

        lists[[name]][['history']][[paste0("Loop_", loop_id, "_consistency")]] <<- lists[[name]][['content']]
      }          
      

      # Criteria: Keep modules with sufficient expression in positive cells
      # ###################################################################

      if(min_cell_level > 0){

        select_res = selectModulesByLevelInPositiveCells(genemodules.pre_filters, current_gm_id_list, mod_avg.bin, antler_env$read_counts('Normalized'), min_cell_level=min_cell_level, verbose=verbose)
        
        current_gm_id_list = select_res$gm_ids
        
        lists[[name]][['content']] <<- genemodules.pre_filters[current_gm_id_list]
        report <- paste0(report, select_res$verbose.log, ' (', length(unlist(genemodules.pre_filters[current_gm_id_list])), ' genes)')
    
        if(length(lists[[name]][['content']])==0){
          lists[[name]][['content']] <<- list()
          return()
        }

        lists[[name]][['history']][[paste0("Loop_", loop_id, "_levels")]] <<- lists[[name]][['content']]

      }

      # Criteria: keep modules with sufficient skewness
      # ###############################################

      if(!identical(mod_skewness_thres, -Inf)){

        select_res = selectModulesBySkewness(current_gm_id_list, mod_avg, skewness_thres=mod_skewness_thres, verbose=verbose)
        
        current_gm_id_list = select_res$gm_ids
        
        lists[[name]][['content']] <<- genemodules.pre_filters[current_gm_id_list]
        report <- paste0(report, select_res$verbose.log, ' (', length(unlist(genemodules.pre_filters[current_gm_id_list])), ' genes)')

        if(length(lists[[name]][['content']])==0){
          lists[[name]][['content']] <<- list()
          return()
        }

        lists[[name]][['history']][[paste0("Loop_", loop_id, "_skewed")]] <<- lists[[name]][['content']]

      }

      # # Criteria: keep modules with feature correlation
      # # ###############################################

      if(!identical(ordering_correlation_thres, NA) & !identical(ordering, NA)) {

        select_res = selectModulesByOrderingCorrelation(current_gm_id_list, mod_avg, order_corr_thres=ordering_correlation_thres, ordering = ordering, verbose=verbose)
        
        current_gm_id_list = select_res$gm_ids
        
        lists[[name]][['content']] <<- genemodules.pre_filters[current_gm_id_list]
        report <- paste0(report, select_res$verbose.log, ' (', length(unlist(genemodules.pre_filters[current_gm_id_list])), ' genes)')

        if(length(lists[[name]][['content']])==0){
          lists[[name]][['content']] <<- list()
          return()
        }

        lists[[name]][['history']][[paste0("Loop_", loop_id, "_OrderingCorrelation")]] <<- lists[[name]][['content']]

      }

      # 3. Recluster gene modules without removed genes
      ##################################################
      
      # Stop loop if no gene module has been excluded during iteration, and 
      # all gm size criteria are fullfilled (if any)
      if(
        identical(lists[[name]][['content']], genemodules.pre_filters) &
        ifelse(is.null(num_final_gms), TRUE, length(lists[[name]][['content']]) == num_final_gms) &
        ifelse(is.null(num_min_final_gms), TRUE, length(lists[[name]][['content']]) >= min(0, num_min_final_gms)) &
        ifelse(is.null(num_max_final_gms), TRUE, length(lists[[name]][['content']]) <= max(0, num_max_final_gms))
        ) {

        lists[[name]][['history']][["Final"]] <<- lists[[name]][['content']]

        goon <- FALSE
        report <- paste0(report, '\nStop iterations as gene modules are not affected by filters.\n')
      } else {

        # Recluster genes
        data.logscaled <- data.logscaled[unlist(lists[[name]][['content']]),]
      
        clusterGenes_res <- clusterGenes(
          antler_env$read_counts('Normalized')[unlist(lists[[name]][['content']]),],
          data.logscaled[unlist(lists[[name]][['content']]),],
          corr_matrix[unlist(lists[[name]][['content']]), unlist(lists[[name]][['content']])],
          num_gms     = num_final_gms,
          num_gms_min = num_min_final_gms,
          num_gms_max = num_max_final_gms,
          numcores    = num_cores,
          display     = display,
          basename    = if (process_plots) {paste0(antler_env$output_folder, '/Gene_module_identification_', name, '_heuristic_loop_', (loop_id+1),'.pdf')} else {NULL},
          verbose           = verbose,
          clustering_method = clustering_method)
        
        report <- paste0(report, clusterGenes_res$verbose.log)

        lists[[name]][['content']] <<- setNames(clusterGenes_res$clusters, seq(length(clusterGenes_res$clusters)))
        lists[[name]][['history']][[paste0("Loop_", (loop_id+1), "_reclustered")]] <<- lists[[name]][['content']]

        loop_id = loop_id + 1
        report <- paste0(report, '\nIteration ', loop_id)
      }

    }

    if(verbose_final){
      cat(report)
    }

    antler_env$processing_state <<- c(
      antler_env$processing_state,
      paste0(
        fdate(),
        ' TopCorr gene modules identification:'), 
        paste0('    ', strsplit(report, split='\n')[[1]]))

  }
)


#' Select gene modules from indices
#' 
#' @name gene_modules_select_from_ids
#' @param name_to character string. The name of the new list of gene modules.
#' @param name_from character string. The name of the original list from which gene modules are selected.
#' @param ids an integer vector. The indices of the gene modules to select.
GeneModules$methods(

  select_from_ids = function(
    name_to,
    name_from,
    ids) {

    lists[[name_to]] <<- list(
      "content" = lists[[name_from]][["content"]][ids],
      "history" = "manual")
  }
)


#' Select gene modules from gene names
#' 
#' @name gene_modules_select_from_gene_names
#' @param name_to character string. The name of the new list of gene modules.
#' @param name_from character string. The name of the original list from which gene modules are selected.
#' @param gene_names a character vector. Gene modules containing at least one of these gene names are kept.
GeneModules$methods(

  select_from_gene_names = function(
    name_to,
    name_from,
    gene_names) {

    gm_ids <- Filter(
      function(i) {
        length(intersect(
          lists[[name_from]][["content"]][[i]],
          gene_names
          )) > 0
      },
      seq(length(lists[[name_from]][["content"]]))
      )

    select_from_ids(
      name_to,
      name_from,
      gm_ids)
  }
)


#' Report whether genes belong to modules
#' 
#' @name gene_modules_report_genes
#' @param name character string. The name of the list of gene modules.
#' @param gene_names a character vector. Gene names to report.
GeneModules$methods(

  report_genes = function(
    name = names(),
    gene_names) {

    name <- match.arg(name)

    report <- lapply(
      seq(length(get(name))),
      function(i){
        sort(intersect(
          get(name)[[i]],
          gene_names))      
      })
    names(report) <- seq(length(get(name)))

    report <- Filter(
      function(x) length(x) != 0,
      report)

    not_in_dataset <- setdiff(
      gene_names,
      antler_env$gene_names())

    not_kept <- sort(setdiff(
      setdiff(gene_names, not_in_dataset),
      unlist(report)))

    if (length(not_kept) > 0)
      report[['not in modules']] <- not_kept

    if (length(not_in_dataset) > 0)
      report[['not in dataset']] <- not_in_dataset

    print(report)
  }
)