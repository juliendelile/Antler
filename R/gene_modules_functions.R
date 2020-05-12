

numCorGenes <- function(corr, corr_t, corr_min){
    corr_thres = abs(corr) > corr_t
    length(which(colSums(corr_thres) > corr_min))  
}

getCorThreshold <- function(corr, init_corr_t=.2, numgenes_target=1500, corr_min=3, numcores=3){
  corr_t <- init_corr_t
  while( numCorGenes(corr, corr_t, corr_min=corr_min) >= numgenes_target ){corr_t <- corr_t + .1}
  corr_t <- corr_t - .1
  while( numCorGenes(corr, corr_t, corr_min=corr_min) >= numgenes_target ){corr_t <- corr_t + .01}
  corr_t <- corr_t - .01
  return(corr_t)
}


# Dimension Reduction (TopCorr_DR) - helper functions
clusterGenes <- function(
  readcount,
  readcount_logscaled,
  corr_matrix,
  num_gms           = NULL,
  num_gms_min       = NULL,
  num_gms_max       = NULL,
  numcores          = 1,
  display           = TRUE,
  basename          = NULL,
  verbose           = TRUE,
  clustering_method = "ward.D2") {

  # if(verbose){cat('\nRun hierarchical clustering (dissimilarity from euclidean distance between z-scored log-leveled gene expressions)')}
  # data.distance <- as.dist(fastEuclideanDist(readcount_logscaled))
  # hc <- hclust(data.distance, method=clustering_method)

  data.distance <- as.dist(1-corr_matrix)
  hc            <- hclust(data.distance, method=clustering_method)

  if(!is.null(num_gms)){
    verbose.log = paste0('\nGene HC / number of gene modules given a priori: ', num_gms)
    if(verbose){cat(verbose.log)}
  } else {

    # if(verbose){cat('Estimate number of gene modules\n')}
    log <- ''

    num_gms <- get_optimal_cluster_number(
      readcount_logscaled,
      hc,
      1e10,
      numcores,
      display       = display,
      pdf_plot_path = basename,
      method = "min_area")

    if(!is.null(num_gms_max)){
      if(num_gms > num_gms_max){log = ', bound to max'}
      num_gms = min(num_gms_max, num_gms)
    }

    if(!is.null(num_gms_min)){
      if(num_gms < num_gms_min){log = ', bound to min'}
      num_gms = max(num_gms_min, num_gms)
    }

    verbose.log = paste0('\nGene HC / number of gene modules estimated: ', num_gms, log)
    if(verbose){cat(verbose.log)}
  }

  if(num_gms > nrow(corr_matrix)){
    num_gms <- nrow(corr_matrix)
    verbose.log <- paste0(verbose.log, '\nNumber of GM set to number of genes: ', num_gms)
  }

  hc.mod_ids = cutree(hc, k=num_gms)
  mod_ids_dendrogram_order = unique(hc.mod_ids[order.dendrogram(as.dendrogram(hc))])

  return(
    list("clusters"=lapply(mod_ids_dendrogram_order, function(i){sort(names(which(hc.mod_ids==i)))}),
         "verbose.log"=verbose.log
      )
    )
}

# Dimension Reduction (TopCorr_DR) - helper functions
selectModulesByPositiveCellRange <- function(gm_id_list, mod_avg.bin, min_cell=3, max_cell=1e10, verbose=FALSE){

  row_sums = rowSums(mod_avg.bin[gm_id_list, , drop=F])

  gm_ids = gm_id_list[which(row_sums >= min_cell & row_sums <= max_cell)]

  verbose.log = paste0('\nSelect gene modules expressed in at least ', min_cell, ' cells', if(max_cell==1e10){''}else{paste0(' and at most ', max_cell, ' cells')}, ': ', length(gm_ids), ' remaining')
  
  if(verbose){cat(verbose.log)}
  
  return(

    list("gm_ids"=gm_ids, 
      "verbose.log"=verbose.log
      )
  )
}

# Dimension Reduction (TopCorr_DR) - helper functions
selectModulesByConsistency <- function(gms, gm_id_list, mod_avg.bin, data_normalized.bin, consistency_thres=.4, verbose=FALSE){

  mod.consistency = unlist(lapply(gm_id_list,
              function(i){
                pos_cells_from_avg = which(mod_avg.bin[i,] == 1)
                pos_neg = table(data_normalized.bin[gms[[i]], pos_cells_from_avg])
                pos_neg[2] / sum(pos_neg)
              }))

  gm_ids = gm_id_list[which(mod.consistency > consistency_thres)]

  verbose.log = paste0('\nSelect gene modules with consistency higher than ', consistency_thres, ': ', length(gm_ids), ' remaining')

  if(verbose){cat(verbose.log)}
  
  return(

    list("gm_ids"=gm_ids, 
      "verbose.log"=verbose.log
      )
  )
}

# Dimension Reduction (TopCorr_DR) - helper functions
selectModulesByLevelInPositiveCells <- function(gms, gm_id_list, mod_avg.bin, readcounts, min_cell_level=0, verbose=FALSE){

  mod.cell.level = unlist(lapply(gm_id_list, function(i){
                pos_cells_from_avg = which(mod_avg.bin[i,] == 1)
                mean(readcounts[(gms[[i]]), pos_cells_from_avg])
          }))

  gm_ids = gm_id_list[which(mod.cell.level >= min_cell_level)]

  verbose.log = paste0('\nSelect gene modules with expression level of least ', min_cell_level, ' gene unit per cell: ', length(gm_ids), ' remaining')

  if(verbose){cat(verbose.log)}
  
  return(

    list("gm_ids"=gm_ids, 
      "verbose.log"=verbose.log
      )
  )
}

# Dimension Reduction (TopCorr_DR) - helper functions
selectModulesBySkewness <- function(gm_id_list, mod_avg, skewness_thres=0, verbose=FALSE){

  genemodules.skewness = apply(mod_avg[gm_id_list, ,drop=F], 1, moments::skewness)


  gm_ids = gm_id_list[which(genemodules.skewness >= skewness_thres)]

  verbose.log = paste0('\nSelect gene modules with skewness higher than ', skewness_thres, ': ', length(gm_ids), ' remaining')

  if(verbose){cat(verbose.log)}
  
  return(

    list("gm_ids"=gm_ids, 
      "verbose.log"=verbose.log
      )
  )
}

# Dimension Reduction (TopCorr_DR) - helper functions
selectModulesByOrderingCorrelation <- function(gm_id_list, mod_avg, order_corr_thres=0, ordering=NA, verbose=FALSE){

  genemodules.corr = apply(mod_avg[gm_id_list, ,drop=F], 1,function(x){abs(cor(x, ordering, method="spearman"))})


  gm_ids = gm_id_list[which(genemodules.corr >= order_corr_thres)]

  verbose.log = paste0('\nSelect gene modules with ordering correlation higher than ', order_corr_thres, ': ', length(gm_ids), ' remaining')

  if(verbose){cat(verbose.log)}
  
  return(

    list("gm_ids"=gm_ids, 
      "verbose.log"=verbose.log
      )
  )
}
