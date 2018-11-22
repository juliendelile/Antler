
DR <- setRefClass(
    Class = "DR",

    fields = list(
       genemodules="list",
       genemodules.selected="list",
       method="character",
       report="character"
      ),

    method = list(

      initialize = function(..., genemodules=list(), genemodules.selected=list(), method=character(0), report=character(0)){

        callSuper(..., genemodules=genemodules, genemodules.selected=genemodules.selected, method=method, report=report)
      },

      saveToCSV = function(){
        stop('TODO: csv export')
      },

      whichModuleContains = function(gene_name){

        for(used_gm in grep('genemodules', names(.self$.refClassDef@fieldClasses), value=T)){

          if( length(.self[[used_gm]]) > 0 ){
            mod_id = Filter(function(i){gene_name %in% .self[[used_gm]][[i]]},
                   seq(length(.self[[used_gm]]))
                   )

            cat(paste0('Module id in "', used_gm, '" -> ', mod_id, '\n'))
          }
        }

      }
      )
  )

#' @include common.R

TopCorr_DR <- setRefClass(
    Class = "TopCorr_DR",

    fields = list(
        correlation_matrix='matrix',
        correlation_thres='numeric',
        genemodules.history='list',
        reduce_memory='logical'
      ),

    contains = "DR",

    methods = list(

      initialize = function(..., correlation_matrix=matrix(0), correlation_thres=numeric(0), genemodules.history=list(), reduce_memory=TRUE){

        callSuper(..., correlation_matrix=correlation_matrix, correlation_thres=correlation_thres, genemodules.history=genemodules.history, reduce_memory=reduce_memory)
      },

      run = function(readcounts, readcounts_normalized, corr_matrix=NA, corr_t=NA, corrgenes_t=NA, corr_quantile=FALSE, corr_min=3, num_mod_max = 10000, mod_min_cell=0, mod_max_cell=Inf, mod_consistency_thres=.4, min_cell_level = 0, mod_skewness_thres=-Inf, ordering_correlation_thres=NA, ordering=NA, num_final_gms=NULL, num_max_final_gms=NULL, num_min_final_gms=NULL, clustering_method="ward.D2", display=TRUE, num_initial_gms=NULL, numcores=3, verbose=FALSE, verbose.final=FALSE, displayed.genes=rownames(readcounts), filter_genes_from_names=TRUE, basename=NULL){
        
        report <<- character(0)

        method <<- "TopCorr_DR"

        # Identify correlated genes
        ###########################

        # cat('Select correlated genes...')
        if(any(is.na(corr_matrix))) {
          if(verbose){cat("Calculate gene-gene correlation matrix\n")}
          used.genes = rownames(readcounts)
          # correlation_matrix <<- fastCor(t(readcounts[used.genes,]), method="spearman")

          if(
              (is.na(corrgenes_t) & is.na(corr_t)) |
              (!is.na(corrgenes_t) & !is.na(corr_t)) 
              ) {
            stop('either "corrgenes_t" or "corr_t" must be specified, not both.')
          }
          
          corr_matrix = fastCor(t(readcounts), method="spearman")
          
          if(!is.na(corrgenes_t)) {
            correlation_thres <<- getCorThreshold(corr_matrix, init_corr_t=.2, numgenes_target=corrgenes_t, corr_min=corr_min, numcores=numcores)
          } else {
            correlation_thres <<- corr_t
          }

          corr_matrix_thres = abs(corr_matrix) >= correlation_thres
          corr_genes = names(which(colSums(corr_matrix_thres) > corr_min))

        } else {

          used.genes = rownames(corr_matrix)
          # correlation_matrix <<- corr_matrix
          # corr_matrix = corr_matrix

          if(is.na(corrgenes_t) & !is.na(corr_t)){

            if(corr_quantile){
              correlation_thres <<- quantile(abs(corr_matrix), prob=corr_t)
            } else {
              correlation_thres <<- corr_t
            }
          } else if(!is.na(corrgenes_t) & is.na(corr_t)){
            correlation_thres <<- getCorThreshold(corr_matrix, init_corr_t=.2, numgenes_target=corrgenes_t, corr_min=corr_min, numcores=numcores)
          } else {
            stop('either "corrgenes_t" or "corr_t" must be specified, not both.')
          }

          # corr_matrix_thres = abs(fastCor(t(readcounts), method="spearman")) >= correlation_thres
          corr_matrix_thres = abs(corr_matrix) >= correlation_thres
          corr_genes = names(which(colSums(corr_matrix_thres) > corr_min))
        }
       
        # if(is.na(corrgenes_t) & !is.na(corr_t)){
        #   if(corr_quantile){
        #     correlation_thres <<- quantile(abs(local_correlation_matrix), prob=corr_t)
        #   } else {
        #     correlation_thres <<- corr_t
        #   }
        # } else if(!is.na(corrgenes_t) & is.na(corr_t)){
        #   correlation_thres <<- getCorThreshold(local_correlation_matrix, init_corr_t=.2, numgenes_target=corrgenes_t, corr_min=corr_min, numcores=numcores)
        # } else {
        #   stop('either "corrgenes_t" or "corr_t" must be specified, not both.')
        # }

        # corr_matrix_thres = abs(local_correlation_matrix) >= correlation_thres
        # corr_genes = names(which(colSums(corr_matrix_thres) > corr_min)) 

        if(reduce_memory && exists("corr_matrix_thres")){
          rm(corr_matrix_thres)
          for (i in 1:10) gc()
        }
        
        verbose.1 = paste0('\nSelect correlated genes: ', length(corr_genes),' genes out of ', length(used.genes),' (threshold: ', sprintf("%.3f", correlation_thres),')')
        # verbose.1 = paste0('Select correlated genes: ', length(corr_genes),' genes out of ', length(used.genes), ' (r > ', correlation_thres,')\n')
        if(verbose){cat(verbose.1)}
        report <<- paste0(report, verbose.1)

        # 1. optimal number of gene modules
        ###################################

        data.logscaled = t(scale(t(log(readcounts[corr_genes,]+1)), center=T, scale=T))

        clusterGenes_res = clusterGenes(
                                       readcounts[corr_genes,],
                                       data.logscaled[corr_genes, ],
                                       corr_matrix[corr_genes, corr_genes],
                                       num_gms=num_initial_gms,
                                       num_gms_min = NULL,
                                       num_gms_max = NULL,
                                       numcores = numcores,
                                       display=display,
                                       basename=if(is.null(basename))NULL else paste0(basename, '_PreFilters.pdf'),
                                       verbose=verbose,
                                       clustering_method=clustering_method
                                        )
        genemodules <<- clusterGenes_res$clusters
        report <<- paste0(report, clusterGenes_res$verbose.log)

        if(length(genemodules)==0){
          genemodules <<- list()
          return()
        }

        genemodules.history[['corrgenes']] <<- genemodules

        # ############################################
        # 2. Apply quality filters on the gene modules
        ##############################################

        # applying gene module filtering on Normalized filters (not on smoothed ones)
        data_normalized.logscaled = t(scale(t(log(readcounts_normalized[corr_genes,]+1)), center=T, scale=T))
        
        if(mod_consistency_thres > 0){
          cat("\nStart binarizing... ")
          data_normalized.bin = silent.binarize.array(data_normalized.logscaled) # no need to normalize non-corrgenes
          cat("done\n")
          } else {data_normalized.bin = NULL}

        goon = TRUE

        loop_id = 1

        while(goon){

          genemodules.pre_filters <- genemodules

          mod_avg = do.call(rbind, lapply(genemodules.pre_filters, function(m){colMeans(data_normalized.logscaled[m,, drop=F])}))

          mod_avg.bin = silent.binarize.array(mod_avg)

          if(display){dev.new(); heatmap.plus(mod_avg.bin, scale='none', Rowv=NA, Colv=NA, labCol = NA)}

          current_gm_id_list = seq(length(genemodules.pre_filters))

          # Criteria: Keep modules expressed in at least N cells in Unsmoothed (ie Normalized) readcounts
          # #############################################################################################

          if(!(identical(mod_min_cell, 0) & identical(mod_max_cell, Inf))) {

            select_res = selectModulesByPositiveCellRange(current_gm_id_list, mod_avg.bin, mod_min_cell, mod_max_cell, verbose)

            current_gm_id_list = select_res$gm_ids

            genemodules <<- genemodules.pre_filters[current_gm_id_list]
            report <<- paste0(report, select_res$verbose.log, ' (', length(unlist(genemodules.pre_filters[current_gm_id_list])), ' genes)')

            if(length(genemodules)==0){
              genemodules <<- list()
              return()
            }

            genemodules.history[[paste0("Loop_", loop_id, "_mincells")]] <<- genemodules
          }
          

          # Criteria: Keep modules with consistent expression in positive cells
          # ###################################################################

          if(mod_consistency_thres > 0){
            select_res = selectModulesByConsistency(genemodules.pre_filters, current_gm_id_list, mod_avg.bin, data_normalized.bin, consistency_thres=mod_consistency_thres, verbose=verbose)
            
            current_gm_id_list = select_res$gm_ids

            genemodules <<- genemodules.pre_filters[current_gm_id_list]
            report <<- paste0(report, select_res$verbose.log, ' (', length(unlist(genemodules.pre_filters[current_gm_id_list])), ' genes)')

            if(length(genemodules)==0){
              genemodules <<- list()
              return()
            }

            genemodules.history[[paste0("Loop_", loop_id, "_consistency")]] <<- genemodules
          }          
          

          # Criteria: Keep modules with sufficient expression in positive cells
          # ###################################################################

          if(min_cell_level > 0){

            select_res = selectModulesByLevelInPositiveCells(genemodules.pre_filters, current_gm_id_list, mod_avg.bin, readcounts, min_cell_level=min_cell_level, verbose=verbose)
            
            current_gm_id_list = select_res$gm_ids
            
            genemodules <<- genemodules.pre_filters[current_gm_id_list]
            report <<- paste0(report, select_res$verbose.log, ' (', length(unlist(genemodules.pre_filters[current_gm_id_list])), ' genes)')
        
            if(length(genemodules)==0){
              genemodules <<- list()
              return()
            }

            genemodules.history[[paste0("Loop_", loop_id, "_levels")]] <<- genemodules

          }

          # Criteria: keep modules with sufficient skewness
          # ###############################################

          if(!identical(mod_skewness_thres, -Inf)){

            select_res = selectModulesBySkewness(current_gm_id_list, mod_avg, skewness_thres=mod_skewness_thres, verbose=verbose)
            
            current_gm_id_list = select_res$gm_ids
            
            genemodules <<- genemodules.pre_filters[current_gm_id_list]
            report <<- paste0(report, select_res$verbose.log, ' (', length(unlist(genemodules.pre_filters[current_gm_id_list])), ' genes)')

            if(length(genemodules)==0){
              genemodules <<- list()
              return()
            }

            genemodules.history[[paste0("Loop_", loop_id, "_skewed")]] <<- genemodules

          }

          # # Criteria: keep modules with time correlation
          # # ############################################

          if(!identical(ordering_correlation_thres, NA) & !identical(ordering, NA)) {

            select_res = selectModulesByOrderingCorrelation(current_gm_id_list, mod_avg, order_corr_thres=ordering_correlation_thres, ordering = ordering, verbose=verbose)
            
            current_gm_id_list = select_res$gm_ids
            
            genemodules <<- genemodules.pre_filters[current_gm_id_list]
            report <<- paste0(report, select_res$verbose.log, ' (', length(unlist(genemodules.pre_filters[current_gm_id_list])), ' genes)')

            if(length(genemodules)==0){
              genemodules <<- list()
              return()
            }

            genemodules.history[[paste0("Loop_", loop_id, "_OrderingCorrelation")]] <<- genemodules

          }

          # 3. Recluster gene modules without removed genes
          ##################################################
          
          # Stop loop if no gene module has been excluded during iteration, and 
          # all gm size criteria are fullfilled (if any)
          if(
              identical(genemodules, genemodules.pre_filters) &
              !(!is.null(num_final_gms) & identical(length(genemodules), num_final_gms)) &
              !(!is.null(num_min_final_gms) & length(genemodules) < min(0, num_min_final_gms)) &
              !(!is.null(num_max_final_gms) & length(genemodules) > max(0, num_max_final_gms))
            ) {

            genemodules.history[["Final"]] <<- genemodules

            goon = FALSE
            report <<- paste0(report, '\nStop iterations as gene modules are not affected by filters.\n')
          } else {

            # Recluster genes
            data.logscaled = data.logscaled[unlist(genemodules),]
          
            clusterGenes_res = clusterGenes(
                                             readcounts[unlist(genemodules),],
                                             data.logscaled[unlist(genemodules),],
                                             corr_matrix[unlist(genemodules), unlist(genemodules)],
                                             num_gms=num_final_gms,
                                             num_gms_min = num_min_final_gms,
                                             num_gms_max = num_max_final_gms,
                                             numcores = numcores,
                                             display=display,
                                             basename=if(is.null(basename))NULL else paste0(basename, '_Loop_', loop_id, '.pdf'),
                                             # basename=if(is.null(basename))NULL else paste0(basename, '_PostFilters.pdf'),
                                             verbose=verbose,
                                             clustering_method=clustering_method
                                              )
            
            # report <<- paste0(report, clusterGenes_res$verbose.log, '\n')
            report <<- paste0(report, clusterGenes_res$verbose.log)

            genemodules <<- clusterGenes_res$clusters
            genemodules.history[[paste0("Loop_", loop_id, "_reclustered")]] <<- genemodules

            loop_id = loop_id + 1
            report <<- paste0(report, '\nIteration ', loop_id)
          }

          # # Test stop
          # if(identical(sort(unlist(genemodules)), sort(unlist(clusterGenes_res$clusters)))){
          #   goon = FALSE
          #   report <<- paste0(report, '\nStop iterations as gene modules are unchanged.\n')
          # } else {
          #   loop_id = loop_id + 1
          #   report <<- paste0(report, '\nIteration ', loop_id)
          # }

          # current_gm_id_list = seq(length(genemodules))

        }

        if(verbose.final){
          cat(report)
        }

        }
      )
  )

# Allow prcomp S3 class to be used as field of PCA_DR
setOldClass("prcomp")

PCA_DR <- setRefClass(
    Class = "PCA_DR",

    fields = list(
        pca="prcomp"
      ),

    contains = "DR",

    methods = list(

      initialize = function(..., pca=prcomp(0)){
        callSuper(..., pca=pca)
      },

      run = function(readcounts){

        method <<- "PCA_DR"

        readcounts_scaled=t(scale(t(readcounts),center=TRUE,scale=TRUE))
        
        pca <<- prcomp(readcounts_scaled)

        },

      getCellCoordinates = function(){
        return(pca$rotation)
      },

      extractGeneModules = function(num_pcs=NULL, prop_variance_explained=.90, num_genes_per_module = 30, basename=NULL){

        if(is.null(num_pcs)) {
          # variance explained per PC
          std_dev <- pca$sdev
          pr_var <- std_dev^2
          #proportion of variance explained
          prop_varex <- pr_var/sum(pr_var)
          cum_prop_varex <- cumsum(prop_varex)

          if(!is.null(basename)) {
            pdf(paste0(basename, '_variance_explained_per_PCs.pdf'))
            graphics::plot(prop_varex, xlab = "Principal Component", ylab = "Proportion of Variance Explained", type = "b")
            dev.off()
            pdf(paste0(basename, '_cumulative_explained_variance.pdf'))
            graphics::plot(cum_prop_varex, xlab = "Principal Component",
                ylab = "Cumulative Proportion of Variance Explained",
                type = "b")
            dev.off()
          }

          num_pcs = length(which(cum_prop_varex<prop_variance_explained))

        }

        genemodules <<- list()
        for(i in seq(num_pcs)){

          pc_scores = pca$x[,i]
          pc_scores_sorted = pc_scores[rev(order(abs(pc_scores)))]

          genemodules[[paste0('PC',i,'_plus')]] <<- names(pc_scores_sorted[pc_scores_sorted > 0])[1:num_genes_per_module]
          genemodules[[paste0('PC',i,'_minus')]] <<- names(pc_scores_sorted[pc_scores_sorted < 0])[1:num_genes_per_module]

        }

        report <<- paste0('Number of PCs: ', num_pcs, ' / Number of genes: ', length(unlist(genemodules)), '\n')
      }
    )
  )


Dispersion_DR <- setRefClass(
    Class = "Dispersion_DR",

    fields = list(
      ),

    contains = "DR",

    methods = list(

      initialize = function(...){
        callSuper(...)
      },

      run = function(readcounts, zscore_threshold=2.0, mean_threshold=10, numcores=1, clustering_mode=TRUE){
        
        method <<- "Dispersion_DR"
        
        dispersed_genes = intersect(getHighGenes(readcounts, mean_threshold=mean_threshold), getDispersedGenes(readcounts, zscore_threshold))

        if(clustering_mode){

          readcounts.logscaled =  t(scale(t(log(readcounts[dispersed_genes,]+1)), center=T, scale=T))

          data.eucldistance = as.dist(fastEuclideanDist(readcounts.logscaled))

          hc.eucldistance = hclust(data.eucldistance, method="ward.D2")

          num_modules = getOptimalModuleNumber_HC(
                            readcounts.logscaled, hc.eucldistance, 100000, numcores,
                            display=FALSE,
                            pdf_plot_path=NULL
                            )

          hc.mod_ids = cutree(hc.eucldistance, k=num_modules)

          genemodules <<- lapply(1:num_modules,function(i){sort(names(which(hc.mod_ids==i)))})
        } else {
          genemodules <<- list()
          genemodules[[1]] <<- dispersed_genes
        }

        report <<- paste0('Number of dispersed genes: ', length(unlist(genemodules)), '\n', 'Number of modules: ', length(genemodules))

      }
    )
  )

