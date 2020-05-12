clean_memory <- function(n=10) { 
	for (i in 1:n) gc()
}


#' Convert GO terms to gene names.
#' 
#' This function download or load Genes/GO terms association maps from BioMart.
#' @name genes_from_GOterms
#' @param go_ids a character vector containing the GO terms to convert. Each term must be written with the 'GO:' prefix, e.g. 'GO:0007399'.
#' @param gofilename character string. The file path of a csv file contining a mapping between gene names (first column) and GO terms (second column). If NA, the mapping is downloaded from ensembl's bioMart.
#' @param biomart_dataset a character string specifying the BioMart dataset to use for the conversion. The list of available datasets can be listed with \code{\link[biomaRt]{listDatasets}} from the "biomaRt" package. Default to 'mmusculus_gene_ensembl'.
#' @export genes_from_GOterms
genes_from_GOterms <- function(
	go_ids,
	biomart_dataset    = 'mmusculus_gene_ensembl',
	gofilename = system.file(
		"extdata",
		"Annotations/gene_goterms.txt",
		package="Antler")) {

  if(is.na(gofilename)){
    gofilename <- '/tmp/gene_goterms.txt'
    # Run system command from R
    system(paste0("wget -O ",gofilename," 'http://www.ensembl.org/biomart/martservice?query=<?xml version=\"1.0\" encoding=\"UTF-8\"?><!DOCTYPE Query><Query  virtualSchemaName = \"default\" formatter = \"CSV\" header = \"0\" uniqueRows = \"1\" count = \"\" datasetConfigVersion = \"0.6\" ><Dataset name = \"", biomart_dataset, "\" interface = \"default\" ><Filter name = \"go_parent_term\" value = \"",
      paste(c(go_ids), collapse=","),
      "\"/><Attribute name = \"external_gene_name\" /><Attribute name = \"name_1006\" /></Dataset></Query>'"))
  } else if(!file.exists(gofilename)){
    print('Input file does not exist.')
    return(0)
  }

  # Load genes 
  biomart.df <- read.csv(
  	gofilename,
  	header           = FALSE,
  	col.names        = c('gene_name', 'GO_term'),
  	stringsAsFactors = FALSE)

  return(sort(unique(biomart.df$gene_name)))
}


elbow_pos_piecewise_heuristic <- function(
  x,
  y,
  N             = 10,
  plot          = FALSE,
  save.pdf.path = NULL){

  x.interp = seq(min(x), max(x), 1)
  y.interp = approx(x,y,x.interp)$y
  a.start = (y.interp[N] - y.interp[1]) / (x.interp[N] - x.interp[1])
  b.start = (tail(y.interp, n=N)[1] - tail(y.interp,n=1)) / (tail(x.interp, n=N)[1] - tail(x.interp,n=1)) 
    
  # First heuristic

  f1 <- function(x, a,b){
    af1 = a*(x-x.interp[1])+y.interp[1]
    af2 = b*(x-tail(x.interp,n=1))+tail(y.interp,n=1)
    ifelse(af1 > af2, af1, af2)
  }
  
  nls.fit1 <- NULL
  try(nls.fit1 <- nls(
          y~f1(x,a,b),
          algorithm = 'default',
          data.frame(x=x.interp, y=y.interp),
          start=list(a=a.start,b=b.start)
          ), silent=TRUE)

  if (!is.null(nls.fit1)) {

    nls.fit1.coef = coef(nls.fit1)
    x_intersection = (nls.fit1.coef[['b']]*tail(x.interp,n=1) + y.interp[1] - nls.fit1.coef[['a']]*x.interp[1] - tail(y.interp,n=1)) / (nls.fit1.coef[['b']]-nls.fit1.coef[['a']])
    y_intersection = nls.fit1.coef[['b']] * (x_intersection-tail(x.interp,n=1)) + tail(y.interp,n=1)#nls.fit1.coef[['c']]

    num_clusters = as.integer(x_intersection)
    
    custom_plot <- function(){
      ggplot2::ggplot() +
        geom_line(
          data = data.frame("x" = x.interp, "y" = y.interp),
          aes(x, y), color = "gray50") +
        geom_point(
          data = data.frame("x" = x, "y" = y),
          aes(x, y), size = .8) +
        geom_line(
          data = data.frame(
            "x" = x.interp,
            "y" = nls.fit1.coef[['a']] * (x.interp-x.interp[1]) + y.interp[1]),
          aes(x, y), color = "grey", linetype = "dashed") +
        geom_line(
          data = data.frame(
            "x" = x.interp,
            "y" = nls.fit1.coef[['b']] * (x.interp-tail(x.interp,n=1)) + tail(y.interp,n=1)), 
          aes(x, y), color = "grey", linetype = "dashed") +
        geom_line(
          data = data.frame(
            "x" = c(x.interp[1], x_intersection),
            "y" = c(y.interp[1], y_intersection)),
          aes(x, y), color = "orangered") +
        geom_line(
          data = data.frame(
            "x" = c(tail(x.interp, n=1), x_intersection),
            "y" = c(tail(y.interp, n=1), y_intersection)),
          aes(x, y), color = "orangered") +
        geom_point(
          data = data.frame(
            "x" = c(x.interp[1], x_intersection, tail(x.interp, n=1)),
            "y" = c(y.interp[1], y_intersection, tail(y.interp, n=1))),
          aes(x, y), size = 2, color = "orangered", alpha = 1) +
        labs(
          subtitle = "Piecewise heuristic method",
          title = paste0("Optimal number of clusters: ", as.integer(num_clusters))) +
        xlab("Number of candidate clusters") +
        ylab("Distortion") +
        scale_x_continuous(
          breaks = c(pretty(x.interp), x_intersection), #, tail(x.interp2, n=1)),
          labels = c(pretty(x.interp), num_clusters)) + #, tail(x.interp2, n=1))) +
        coord_cartesian(ylim = c(0, 1.1 * y.interp[1])) +
        theme_minimal() +
        theme(
          panel.grid.minor = element_blank())
    }
    if(plot){
      dev.new()
      custom_plot()
    }
    if(!is.null(save.pdf.path)){
      pdf(save.pdf.path, width = 6, height = 4, useDingbats = FALSE)
      print(custom_plot())
      graphics.off()
    }

    return(num_clusters)
  } else {
      stop("The piecewise heuristic method to estimate the optimal number of clusters did not converge. Try the minimal area heuristic or set number of clusters manually.")
  }
}

elbow_pos_area_heuristic <- function(
  x,
  y,
  plot=FALSE,
  save.pdf.path=NULL) {

  x.interp = seq(min(x), max(x), 1)

  # normalize y to match x range
  norm_coeff = (max(x)-1) / max(y)
  y2 <- norm_coeff * y
  y.interp = approx(x, y2, x.interp)$y
  num_points = length(x.interp)

  euc.dist <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2))

  segments_length = c(0, unlist(lapply(seq(num_points-1), function(i){euc.dist(c(x.interp[i], y.interp[i]), c(x.interp[i+1], y.interp[i+1]))})))

  cumulative_segments_length = cumsum(segments_length)

  num_points2 = 100

  x.interp2 = approx(cumulative_segments_length, 
           seq(length(cumulative_segments_length)),
           seq(0, max(cumulative_segments_length), length.out=num_points2)
           )$y

  y.interp2 = approx(x,y2,x.interp2)$y

  # Point M (xm, ym) and line passing through A (xa, ya) and B (xb, yb)
  point_line_distance <- function(xa, ya, xb, yb, xm, ym){

    nx = -(yb-ya)
    ny = (xb-xa)

    n_norm = sqrt(nx^2+ny^2)

    d = abs(nx * (xa-xm) + ny *(ya-ym)) / n_norm
    return(d)
  }

  all_distances = unlist(lapply(seq(2, num_points2-1), function(k){

      dist_1 = 0
      dist_2 = 0

      for(i in seq(k)) {
        dist_1 = dist_1 + point_line_distance(x.interp2[1], y.interp2[1], x.interp2[k], y.interp2[k], x.interp2[i], y.interp2[i])
      }

      for(i in seq(k+1, num_points2)) {
        dist_2 = dist_2 + point_line_distance(x.interp2[k], y.interp2[k], tail(x.interp2, n=1), tail(y.interp2, n=1), x.interp2[i], y.interp2[i])
      }
      return(dist_1 / k + dist_2 / (num_points2 - k))
    }))

  optimal_point_id = 1 + which.min(all_distances)
  
  num_clusters = as.integer(x.interp2[optimal_point_id])

  custom_plot <- function(){
    ggplot2::ggplot() +
      geom_polygon(
        data = data.frame(
          "x" = c(x.interp2, x.interp2[optimal_point_id]),
          "y" = c(y.interp2, y.interp2[optimal_point_id])),
        aes(x, y), fill = "moccasin", color = "gray80", alpha = .5) +
      geom_line(
        data = data.frame("x" = x.interp2, "y" = y.interp2),
        aes(x, y), color = "gray50") +
      geom_point(
        data = data.frame("x" = x, "y" = y2),
        aes(x, y), size = .8) +
      geom_line(
        data = data.frame(
          "x" = c(x.interp2[optimal_point_id], x.interp2[optimal_point_id]),
          "y" = c(y.interp2[optimal_point_id], 0)),
        aes(x, y), linetype = "dashed"
        ) +
      geom_point(
        data = data.frame(
          "x" = x.interp2[optimal_point_id],
          "y" = y.interp2[optimal_point_id]),
        aes(x, y), color = "orangered", size = 2, alpha = .8) +
      labs(
        subtitle = "Minimal area heuristic method",
        title = paste0("Optimal number of clusters: ", as.integer(num_clusters))) +
      xlab("Number of candidate clusters") +
      ylab("Distortion") +
      scale_x_continuous(
        breaks = c(pretty(x.interp2), x.interp2[optimal_point_id]), #, tail(x.interp2, n=1)),
        labels = c(pretty(x.interp2), num_clusters)) + #, tail(x.interp2, n=1))) +
      theme_minimal() +
      theme(
        panel.grid.minor = element_blank())
  }
  if(plot){
    dev.new()
    custom_plot()
  }
  if(!is.null(save.pdf.path)){
    pdf(save.pdf.path, width = 6, height = 4, useDingbats = FALSE)
    print(custom_plot())
    graphics.off()
  }

  return(num_clusters)
}

# use within-cluster sum of squares as distortion measure
getDistortion <- function(k, hc, dataset){

  hc.mod_ids <- cutree(hc, k=k)
  distortion <- 0

  for(ik in 1:k){
    clust      <- dataset[which(hc.mod_ids==ik),,drop=F]
    clust_mean <- colMeans(clust)
    distortion <- distortion + sum(apply(clust, 1, function(r){ rc = r-clust_mean; return(sum(rc^2))}))
  }

  return(distortion)
}

get_optimal_cluster_number <- function(
  dataset,
  hc,
  num_mod_max,
  numcores,
  display       = FALSE,
  pdf_plot_path = NULL,
  method = c("min_area", "piecewise")){

  method <- match.arg(method)

  krange = c(seq(1,99,10), seq(100, 999, 100), seq(1000, 100000, 1000))
  krange <- krange[which(krange < nrow(dataset))]
  krange <- krange[which(krange <= num_mod_max)] # avoid useless calculation
  krange <- c(krange, nrow(dataset))

  if(numcores == 1){
    distortion.vals <- unlist(lapply(krange, function(k){getDistortion(k, hc, dataset)}))
  } else {

    #  outfile="/dev/null" : remove outfile arguments to stop printing workers' outputs
    #  methods=FALSE : do not load packages
    cl <- parallel::makeCluster(
            min(numcores, length(krange)),
            outfile = "/dev/null",
            methods = FALSE) 

    parallel::clusterExport(cl, c('hc', 'dataset', 'krange'), envir = environment())
    # parallel::clusterExport(cl, c('krange'), envir = pryr::where("krange"))
    # parallel::clusterExport(cl, c('hc'), envir = pryr::where('hc'))
    # parallel::clusterExport(cl, c('dataset'), envir = pryr::where('dataset'))

    doSNOW::registerDoSNOW(cl)

    cat('\nRun parallel distortion jobs\n')
    pb <- txtProgressBar(max = length(krange), style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)

    distortion.vals = unlist(foreach(
            k = krange,
            .options.snow = opts,
            .noexport=ls(),
            .packages=c('Matrix')
            ) %dopar% {
            
      if(any(is.na(dataset))){
        print(paste(k, " contains NA values"))
      }

      hc.mod_ids = cutree(hc, k=k)
      distortion=0
      for(ik in 1:k){
        clust = dataset[which(hc.mod_ids==ik),,drop=F]
        clust_mean = colMeans(clust)
        distortion <- distortion + sum(apply(clust, 1, function(r){ rc = r-clust_mean; return(sum(rc^2))}))
      }

      print(paste(k, "before"))
      print(paste(k, ":", distortion))
      return(distortion)
    })

    parallel::stopCluster(cl)
  }
  
  distortion.df = data.frame(
    x = krange,
    y = distortion.vals
    )
  
  num_mod_max_idx <- length(which(distortion.df$x <= num_mod_max))

  if (method == "piecewise"){
    num_modules <- as.integer(elbow_pos_piecewise_heuristic(
      distortion.df$x[1:num_mod_max_idx],
      distortion.df$y[1:num_mod_max_idx],
      plot = display,
      save.pdf.path = pdf_plot_path)[1])

  } else if (method == "min_area") {
    num_modules <- as.integer(elbow_pos_area_heuristic(
      distortion.df$x[1:num_mod_max_idx],
      distortion.df$y[1:num_mod_max_idx],
      plot = display,
      save.pdf.path = pdf_plot_path))
  }

  return(num_modules)
}

# Remove the annoying report and deal with single row matrix binarization
silent.binarize.array <- function(x){

  # require(ArrayBin)

  nr = nrow(x)

  x <- rbind(x, x[1,,drop=FALSE])

  invisible(capture.output(x.bin <- binarize.array(x)))

  return(x.bin[1:nr, , drop=F])
}

seurat_dispersion_zscore <- function(data, nBin = 20){
  # For each gene, calculate log(mean(x)+1) with x distribution of the gene level over all cells
  data_x=apply(data, 1, function(x) log(mean(x)) ) 
  # For each gene, calculate log(var(x)/mean(x)) with x distribution of the gene level over all cells
  data_y=apply(data, 1, function(x) log(var(x)/mean(x) )) 
  # Determine nBin intervals between min(data_x) and max(data_x), and associate each gene to the interval it belongs to
  data_x_bin=cut(data_x, nBin)
  # Within each bin, compute the mean of the dispersion (data_y) of all genes
  mean_y=tapply(data_y,data_x_bin,mean)
  # Within each bin, compute the standard deviation of the dispersion (data_y) over all genes
  # !! If a bin contains a single gene, sd returns NA !!
  sd_y=tapply(data_y,data_x_bin,sd)
  # For each gene, calculate the normalized dispersion with the mean and std corresponding to the gene level bin
  data_norm_y = (data_y - mean_y[as.numeric(data_x_bin)]) / sd_y[as.numeric(data_x_bin)]
  rownames(data_norm_y) = rownames(data)

  return(data_norm_y)
}

geo_mean <- function(x){

  if(any(is.na(x))) {
    return(NA)
  } else if(any(x==0)){
    return(0)
  } else{
    return(exp(mean(log(x))))
  }
}

quantile_mean <- function(x, quantile_threshold=0){
    x_thres = quantile(x=x, c(quantile_threshold, 1-quantile_threshold))
    x[which(x < x_thres[[1]])] <- x_thres[[1]]
    x[which(x > x_thres[[2]])] <- x_thres[[2]]
    mean(x)
}


fastCor <- function(x, method='spearman', subdims){
  if(method=='spearman'){
    # rank features
    x = t(apply(t(x), 1, rank))
  } else if(method=='pearson'){
    x = t(x)
  } else {
    stop("Fast correlation works only with 'spearman' or 'pearson'")
  }
  # Center each variable
  x = x - rowMeans(x);
  # Standardize each variable
  x = x / sqrt(rowSums(x^2));

  # Calculate correlations
  if(missing(subdims)){
    return(tcrossprod(x))
  } else {
    return(tcrossprod(x, x[subdims, , drop=F]))
  }

}

testFastCor <- function(n_samples=100, n_features=40){

  randmat = matrix(runif(n_samples * n_features), nrow=n_features)
  
  test_spearman = table( abs(fastCor(randmat, method='spearman') - cor(randmat, method='spearman')) < 1e-7)
  
  stopifnot(test_spearman == (n_samples * n_samples))

  test_pearson = table( abs(fastCor(randmat, method='pearson') - cor(randmat, method='pearson')) < 1e-7)

  stopifnot(test_pearson == n_samples * n_samples)

  print('fastCor test passed.')

}

# http://nonconditional.com/2014/04/on-the-trick-for-computing-the-squared-euclidian-distances-between-two-sets-of-vectors/
fastEuclideanDist <- function(x, y=x){

  xy <- x %*% t(y)

  xx <- matrix(rep(rowSums(x^2), times = nrow(y)), ncol = nrow(y), byrow = F)

  yy <- matrix(rep(rowSums(y^2), times = nrow(x)), nrow = nrow(x), byrow = T)

  return(sqrt(abs(xx + yy - 2 * xy)))
}

compareEuclideanDistSpeed <- function(x, nc){
  t1=Sys.time()
  d1=as.matrix(dist(t(x), method='euclidean', diag=T, upper=T))
  t2=Sys.time()
  print(paste0('Default distance time: ', t2-t1))

  t3=Sys.time()
  d2=fastEuclideanDist(t(x))
  t4=Sys.time()
  print(paste0('Fast distance time: ', t4-t3))

  t5=Sys.time()
  d3=as.matrix(
              amap::Dist(
                          t(x),
                          method='euclidean',
                          nbproc=if(missing(nc)) as.integer(.5*parallel::detectCores()) else nc
                          )
              )
  t6=Sys.time()
  print(paste0('Parallel distance time: ', t6-t5))
}


fdate <- function(){
  format(Sys.time(), "%a %b %d %X")
  # format(Sys.time(), "%a %b %d %X %Y")
}

