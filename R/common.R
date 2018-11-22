
getTFs <- function(){
  read.table(system.file("extdata", "Annotations/TF_mmusculus_150722.dat", package="Antler"), row.names=1, header=TRUE, stringsAsFactors=F)[,'gene_name']
}

fdate <- function(){
  format(Sys.time(), "%a %b %d %X")
  # format(Sys.time(), "%a %b %d %X %Y")
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

  # print(dim(d1))
  # print(dim(d2))p
  # print(dim(d3))

}

generateCellColors <- function(phenoD, method='hsv', seed=10, hue_shift=.0, replicate_variation=TRUE){

  if(method=='hsv'){

    sample_t = unique(phenoD$timepoint)
    sample_rid = unique(phenoD$replicate_id)
    # sample_treatment = lapply(unique(phenoD$treatment), function(st){paste0(tail(st, n=1)[[1]], collapse='')})
    sample_treatment = lapply(unique(as.character(phenoD$treatment)), function(st){
            tail(strsplit(st, split='_')[[1]], n=1)
            })

    low.t = min(unlist(sample_t))
    hig.t = max(unlist(sample_t))
    low.r = min(unlist(sample_rid))
    hig.r = max(unlist(sample_rid))
    all.ldt = list(unique(unlist(sample_treatment)))[[1]] # first 'list' in case there is a unique treatment in the dataset

    cells_colors <- unlist(lapply(seq(nrow(phenoD)), function(i) {
          # browser()
          # sample day
          d = phenoD$timepoint[i]
          # sample replicate
          r = phenoD$replicate_id[i]
          # sample treatment
          # t = phenoD$treatment[i]
          t = tail(strsplit(as.character(phenoD$treatment[i]), split='_')[[1]], n=1)

          # Time -> Saturation between .1 and 1
          if(low.t == hig.t){
            d.scaled=.5
          } else {
            d.scaled = .1 + .9 * ((d - low.t) / (hig.t-low.t))
          }

          # Replicates -> Value between .35 and .65
          if(low.r == hig.r | !replicate_variation){
            r.scaled=1
          } else {
            r.scaled = .8 + .2 * ((r - low.r) / (hig.r-low.r))
          }

          # Hue -> last day's treatment
          t.scaled = (hue_shift + ( which(all.ldt == t) / length(all.ldt))) %% 1

          # print(paste0(c(t.scaled, d.scaled, r.scaled), collapse=" "))
          
          hsv(t.scaled, d.scaled, r.scaled, 1)
      }))

    return(cells_colors)

  } else if(method=='random'){
    
    set.seed(seed)
    samples_fullname = unique(phenoD$cells_samples_fullname)
    df = data.frame(cols = sample(colours(), length(samples_fullname)), row.names=samples_fullname)
    return(as.character(df[phenoD$cells_samples_fullname, 'cols']))

  }
}

generateCellSampleNames <- function(phenoD){
  apply(phenoD, 1, function(x){paste0('day', x["timepoint"], ' rep', x["replicate_id"], ' ', tail(strsplit(as.character(x["treatment"]), split='_')[[1]], n=1), collapse='')})
}

generateCellSampleFullNames <- function(phenoD){
  apply(phenoD, 1, function(x){paste0('day', x["timepoint"], '_rep', x["replicate_id"], '_', x["treatment"], collapse='')})
}

generateCellFullTreatments <- function(phenoD){
  apply(phenoD, 1, function(x){paste0('day', x["timepoint"], '_', x["treatment"], collapse='')})
}

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

geo_mean <- function(x){

  if(any(is.na(x))) {
    return(NA)
  } else if(any(x==0)){
    return(0)
  } else{
    return(exp(mean(log(x))))
  }
  
}

# Calculate knee-point
# version post 171128
knee_pos_piecewise_heuristic <- function(x, y, N=10, plot=FALSE, save.pdf.path=NULL){

  # save(x,y,N,plot,save.pdf.path, file='tmp.Rdata')

  # load('tmp.Rdata')


  x.interp = seq(min(x), max(x), 1)
  
  # x.interp <- x.interp / max(x)
  # x <- x / max(x)
  # y <- y / max(y)

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

  if(!is.null(nls.fit1)){

    nls.fit1.coef = coef(nls.fit1)
    x_intersection = (nls.fit1.coef[['b']]*tail(x.interp,n=1) + y.interp[1] - nls.fit1.coef[['a']]*x.interp[1] - tail(y.interp,n=1)) / (nls.fit1.coef[['b']]-nls.fit1.coef[['a']])
    y_intersection = nls.fit1.coef[['b']] * (x_intersection-tail(x.interp,n=1)) + tail(y.interp,n=1)#nls.fit1.coef[['c']]

    num_clusters = x_intersection

    # angle_a = atan(nls.fit1.coef[['a']])
    # angle_b = atan(nls.fit1.coef[['b']])
    # new_slope = tan(.5*(angle_a + angle_b))

    # # new_slope = .5*(nls.fit1.coef[['a']] + nls.fit1.coef[['b']])

    # # distance between bisectrix and distortion curve
    # k_dist = unlist(lapply(seq(length(x.interp)), function(i){

    #     x0 = x.interp[i]
    #     y0 = y.interp[i]

    #     A = new_slope
    #     B = -1
    #     C = y_intersection - new_slope * x_intersection

    #     return(abs(A*x0+B*y0+C)/sqrt(A*A+B*B))
    #   }))

    # num_clusters = which.min(k_dist)

    nls.fit1.plot <- function(){
      plot(x,y, main=paste0('1st heuristic (', as.integer(num_clusters), ')'), cex=.8)
      abline(v=0, col='grey')
      abline(h=0, col='grey')
      lines(x.interp, y.interp, col='red')
      lines(x.interp, nls.fit1.coef[['a']] * (x.interp-x.interp[1]) + y.interp[1])
      lines(x.interp, nls.fit1.coef[['b']] * (x.interp-tail(x.interp,n=1)) + tail(y.interp,n=1))
      points(x_intersection, y_intersection, col='green', lwd=2)
      # points(x_af1_intersect, 0, col='green', lwd=2)
      lines(x.interp, abs(y.interp - f1(x.interp, nls.fit1.coef[['a']], nls.fit1.coef[['b']])))
      # lines(x.interp, new_slope * (x.interp - x_intersection) + y_intersection, col='blue', lwd=2) 
    }

    if(plot){
      dev.new()
      nls.fit1.plot()
    }

    if(!is.null(save.pdf.path)){
      pdf(save.pdf.path)
      nls.fit1.plot()
      graphics.off()
    }

  } else {

    # Second heuristic in case first heuristic fails
    f2 <- function(x, ax, ay, bx, by){

      unlist(lapply(x, function(t){
              if(t < ax){
                (ay - y.interp[1]) * (t - x.interp[1]) / (ax - x.interp[1]) + y.interp[1]
              } else if(t < bx){
                (by - ay) * (t - ax) / (bx - ax) + ay
              } else {
                (tail(y.interp, n=1) - by) * (t - bx) / (tail(x.interp, n=1) - bx) + by
              }
        }))
    }


    ax.start = x.interp[as.integer(length(x.interp)/30)]
    ay.start = y.interp[as.integer(length(x.interp)/30)]
    bx.start = x.interp[as.integer(length(x.interp)/2)]
    by.start = y.interp[as.integer(length(x.interp)/2)]

    nls.fit2 <- NULL
    try(
        nls.fit2 <- nls(
                y~f2(x,ax,ay,bx,by),
                algorithm = 'default',
                data.frame(x=x.interp, y=y.interp),
                start=list(ax=ax.start,ay=ay.start,bx=bx.start,by=by.start),
                trace=F,
                control=nls.control(maxiter=2000)
                ),
        silent=TRUE
        )

    if(!is.null(nls.fit2)){

      nls.fit2.coef = coef(nls.fit2)

      num_clusters = 2 * nls.fit2.coef[['ax']]

      nls.fit2.plot <- function(){

        plot(x,y,type='p', cex=.8, main=paste0('2nd heuristic (', num_clusters, ')'))
        abline(v=0, col='grey')
        abline(h=0, col='grey')
        lines(x.interp, y.interp, col='red')

        x.interp.part1 = x.interp[x.interp < nls.fit2.coef[['ax']]]
        lines(
              x.interp.part1, 
              (nls.fit2.coef[['ay']] - y.interp[1]) * (x.interp.part1 - x.interp[1]) / (nls.fit2.coef[['ax']] - x.interp[1]) + y.interp[1]
              )
        x.interp.part2 = x.interp[x.interp >= nls.fit2.coef[['ax']] & x.interp < nls.fit2.coef[['bx']]]
        lines(
              x.interp.part2,
              (nls.fit2.coef[['by']] - nls.fit2.coef[['ay']]) * (x.interp.part2 - nls.fit2.coef[['ax']]) / (nls.fit2.coef[['bx']] - nls.fit2.coef[['ax']]) + nls.fit2.coef[['ay']]
          )
        x.interp.part3 = x.interp[x.interp >= nls.fit2.coef[['bx']]]
        lines(
              x.interp.part3,
              (tail(y.interp, n=1) - nls.fit2.coef[['by']]) * (x.interp.part3 - nls.fit2.coef[['bx']]) / (tail(x.interp, n=1) - nls.fit2.coef[['bx']]) + nls.fit2.coef[['by']]
          )
        abline(v=num_clusters, col='green')
      }

      if(plot){
        dev.new()
        nls.fit2.plot()
      }

      if(!is.null(save.pdf.path)){
        pdf(save.pdf.path)
        nls.fit2.plot()
        graphics.off()
      }
    } else {
      # no heuristic
      num_clusters = length(x)

      no.fit.plot <- function(){
        plot(x,y,type='p', cex=.8, main=paste0('No fit (', num_clusters, ')'))
        abline(v=0, col='grey')
        abline(h=0, col='grey')
        lines(x.interp, y.interp, col='red')
        abline(v=num_clusters, col='green')
      }

      if(plot){
        dev.new()
        no.fit.plot()
      }

      if(!is.null(save.pdf.path)){
        pdf(save.pdf.path)
        no.fit.plot()
        graphics.off()
      }
    }
 
  }

  diff = abs(x - num_clusters)
  dist_score = y[which.min(diff)]
  
  return(c(num_clusters, dist_score))

}

knee_pos_piecewise_heuristic_v2 <- function(x, y, plot=FALSE, save.pdf.path=NULL){

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
    plot(x,y2, main=paste0('Recent heuristic (', as.integer(num_clusters), ' clusters)'), cex=.8, asp=1)
    abline(v=0, col='grey')
    abline(h=0, col='grey')
    lines(x.interp, y.interp, col='red')
    points(x.interp2, y.interp2, col='blue', pch=4, cex=.3)
    segments(x.interp2[1], y.interp2[1], x.interp2[optimal_point_id], y.interp2[optimal_point_id], col='green')
    segments(x.interp2[optimal_point_id], y.interp2[optimal_point_id], tail(x.interp2, n=1), tail(y.interp, n=1), col='blue')
  }

  if(plot){
    dev.new()
    custom_plot()
  }

  if(!is.null(save.pdf.path)){
    pdf(save.pdf.path)
    custom_plot()
    graphics.off()
  }

  return(num_clusters)
}

# use within-cluster sum of squares as distortion measure
getDistortion <- function(k, hc, dataset){
  hc.mod_ids = cutree(hc, k=k)
  distortion=0
  for(ik in 1:k){
    clust = dataset[which(hc.mod_ids==ik),,drop=F]
    clust_mean = colMeans(clust)
    # distortion <- distortion + sum(apply(clust, 1, function(r){ rc = r-clust_mean; return(sqrt(sum(rc^2)))}))
    distortion <- distortion + sum(apply(clust, 1, function(r){ rc = r-clust_mean; return(sum(rc^2))}))
  }
  return(distortion)
}

getOptimalModuleNumber_HC <- function(dataset, hc, num_mod_max, numcores, display=FALSE, pdf_plot_path=NULL, old_heuristic=TRUE){

  # krange = c(seq(1,49,1), seq(50, 499, 10), seq(500, 40000, 100))
  krange = c(seq(1,99,10), seq(100, 999, 100), seq(1000, 100000, 1000))
  # krange = c(1,10,100,1000,10000,100000)
  krange <- krange[which(krange < nrow(dataset))]
  krange <- krange[which(krange <= num_mod_max)] # avoid useless calculation
  krange <- c(krange, nrow(dataset))
  
  if(numcores == 1){
    distortion.vals = unlist(lapply(krange, function(k){getDistortion(k, hc, dataset)}))
  } else {

    #  outfile="/dev/null" : remove outfile arguments to stop printing workers' outputs
    #  methods=FALSE : do not load packages
    cl <- parallel::makeCluster(
            min(numcores, length(krange)),
            outfile="/dev/null",
            methods=FALSE) 

    parallel::clusterExport(cl, c('hc', 'dataset', 'krange'), envir=environment())

    doSNOW::registerDoSNOW(cl)

    cat('Run parallel distortion jobs\n')
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
        # distortion <- distortion + sum(apply(clust, 1, function(r){ rc = r-clust_mean; return(sqrt(sum(rc^2)))}))
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
  
  num_mod_max_idx = length(which(distortion.df$x <= num_mod_max))
    
  if(old_heuristic){

    num_modules = as.integer(knee_pos_piecewise_heuristic(distortion.df$x[1:num_mod_max_idx], distortion.df$y[1:num_mod_max_idx], plot=display, save.pdf.path = pdf_plot_path)[1])
  } else {
    num_modules = as.integer(knee_pos_piecewise_heuristic_v2(distortion.df$x[1:num_mod_max_idx], distortion.df$y[1:num_mod_max_idx], plot=display, save.pdf.path = pdf_plot_path))
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

plotHeatmapCellsGeneModules <- function(mat, genemodules, filename, cells_colors, displayed.geneset=rownames(mat), heatmap.palette=colorRampPalette(c("#191d73", "white", "#ed7901"))(n = 1000), genemodules.palette=colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(length(genemodules)), genemodules.extra_colors=NA, genemodules.text_colors=NA, genes.extra_colors=NA, zscore.thres = 2, legend=NULL, width=15, height=15, dendrogram=NULL, heatmap3=T,rect_overlay, pretty.params=list("size_factor"=1, "ngenes_per_lines" = 8, "side.height.fraction"=.3), main = NULL) {

  require(heatmap.plus)

  mat.capped = mat
  cap_val = zscore.thres
  mat.capped[which(mat > cap_val)] <- cap_val
  mat.capped[which(mat < - cap_val)] <- - cap_val

  if(tools::file_ext(filename) == 'pdf'){
    pdf(filename, width=width, height=height)
  } else if(tools::file_ext(filename) == 'png'){
    png(filename, width=width, height=height, res=width/10)
  } else if(tools::file_ext(filename) == 'jpeg'){
    jpeg(filename, width=width, height=height, pointsize = 40)
  } else {
    stop('wrong filetype.')
  }

  rlab_de_ids = unlist(lapply(1:length(genemodules), function(i){rep(i, length(genemodules[[i]]))}))

  rlab_de = as.matrix(unlist(lapply(genemodules.palette[rlab_de_ids], as.character)))

  if(!identical(genemodules.extra_colors, NA)){
    rlab_de <- cbind(rlab_de, genemodules.extra_colors[rlab_de_ids, ])
  }

  if(!identical(genes.extra_colors, NA)){
    rlab_de <- cbind(rlab_de, genes.extra_colors)
  }

  if(!identical(genemodules.text_colors, NA)){
    genemodules.text_colors = rep(genemodules.text_colors, lapply(genemodules, length))
  }

  # Heatmap.plus requires at least two columns for ColSideColors
  if(class(cells_colors)=='character' | ncol(cells_colors)==1) {
    ccol = cbind(cells_colors, cells_colors)
  } else {
    ccol = cells_colors
  }

  # Heatmap.plus requires at least two columns for RowSideColors
  if(class(rlab_de)=='character' | ncol(rlab_de)==1) {
    rlab_de = cbind(rlab_de, rlab_de)
  }

  # colnames(ccol) <- NULL

  mod_name = if(is.null(names(genemodules))){seq(length(genemodules))}else{names(genemodules)}
  if(identical(displayed.geneset, NA)){
    
    lr=paste(
          rownames(mat[unlist(genemodules),]),
          paste0('(', rep(mod_name, unlist(lapply(genemodules, length))), ')')
          )
    rowsep=c()
  } else if(displayed.geneset == "naked"){
    lr = rownames(mat[unlist(genemodules),])
  } else {
    
    modulenames = lapply(1:length(genemodules), function(i) {
                      l=genemodules[[i]]
                      list = intersect(l, displayed.geneset)
                      list_split = split(list, ceiling(seq_along(list)/pretty.params$ngenes_per_lines))
                      # list_newline = lapply(list_split, function(ll){paste0(c(ll,'\n'), collapse=', ')})
                      # paste(mod_name[i], ': ', paste0(Reduce(c, list_newline), collapse=','))
                      paste(mod_name[i], ':', paste0(lapply(list_split, function(ll){paste0(ll, collapse=' ')}), collapse='\n'))
                      })
    lr = rep("", length(unlist(genemodules)))
    mod.halflength = unlist(lapply(genemodules, function(l){as.integer(length(l)/2)}))
    mod.firstrow = cumsum(c(1, unlist(lapply(genemodules, length))))
    mod.firstrow = mod.firstrow[1:(length(mod.firstrow)-1)]
    lr[mod.firstrow + mod.halflength] <- modulenames
    rowsep=mod.firstrow-1
  }

  if(heatmap3){
    # browser()
    # pdf(filename, width=width, height=height)
    heatmap.3(
      as.matrix(mat.capped[unlist(genemodules),]),
      Colv = dendrogram,
      Rowv = FALSE, 
      dendrogram=if(!is.null(dendrogram)) "column" else "none",
      col=heatmap.palette,
      scale='none',
      margins=c(0.2,0.2),
      ColSideColors=ccol,
      RowSideColors=t(rlab_de),
      labRow=lr,
      colRow=genemodules.text_colors,
      labCol=F,
      cexRow = pretty.params$size_factor*(0 + .5/log10(length(unlist(genemodules)))),
      useRaster=if(tools::file_ext(filename) == 'pdf') TRUE else FALSE,
      # lmat = rbind(c(6, 0,5,0), c(0,0,2,0), c(4, 1, 3,0), c(0,0,0,0)),
      # lhei = c(3, 3, 10, 2),
      # lwid = c(3, 1.5, 10, 5),
      lmat = rbind(c(6, 5, 0,0), c(0,2,0, 0), c(4, 3,1, 0)),
      lhei = if(!identical(NA, dendrogram)) c(3, 3*pretty.params$side.height.fraction, 10) else c(.1, 3*pretty.params$side.height.fraction, 10),
      lwid = c(3, 10, .3, 5),
      key=F,
      # trace="row",
      rowsep=NULL,
      sepcolor = "gray70",
      sepwidth = c(0.001, 0.001),
      rect_overlay=rect_overlay,
      main = main
        )
    # graphics.off()
    hm=0
  } else {

    hm = heatmap.plus(as.matrix(mat.capped[unlist(genemodules),]),
      scale='none',
      Colv = dendrogram,
      Rowv = NA, 
      col=heatmap.palette,
      labRow = lr,
      labCol = FALSE,
      RowSideColors=rlab_de,
      ColSideColors=ccol,
      margins=c(1,15),
      keep.dendro=TRUE
    )

  }

  if(!is.null(legend)){

    par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
    plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
    legend(x="left", 
      legend=legend[['text']],
      col=legend[['colors']],
      pch=15,  
      y.intersp = .7, cex=.6,
      bg='transparent',
      box.col='transparent'
        )



    # legend(x="topright", legend=legend[['text']], col=legend[['colors']], pch=15,
    #   # trace=TRUE,
    #   # ncol=1+as.integer(length(legend[['text']])/10),
    #   # ncol=8,
    #   # border=FALSE, bty="n", 
    #   y.intersp = 0.7, cex=.7,
    #   # bg='white',
    #   bg="transparent",
    #   box.col='black'
    #   )
  }

  graphics.off()

  # if(!is.null(legend)){
  #   slegend <- function(){
  #     plot(0, xaxt = "n", bty = "n", yaxt = "n", type = "n", xlab = "", 
  #       ylab = "")
  #     legend(x="topleft", legend=legend[['text']], col=legend[['colors']], pch=15,
  #     # y.intersp = 0.7, 
  #     cex=.4,
  #     # box.col='white',
  #     bty="n"
  #     # bg='white',
  #     )
  #   }
  # }

  # source()

  # # png(paste0(m$plot_folder, '/test.png'), width=2500, height=5000)
  # pdf(paste0(m$plot_folder, '/test2.pdf'))#, width=5, height=15)
  # heatmap3(
  #     as.matrix(mat.capped[unlist(genemodules),]),
  #     Colv = dendrogram,
  #     Rowv = NA, 
  #     ColSideLabs=colnames(ccol),
  #     showColDendro=T,
  #     col=heatmap.palette,
  #     scale='none',
  #     ColSideWidth=.5,
  #     margins=c(0.2,0.2),
  #     ColSideColors=ccol,
  #     # RowSideColors=rlab_de,
  #     labRow=lr,
  #     labCol=F,
  #     keep.dendro=T,
  #     useRaster=T,
  #     legendfun=slegend,
  #     cexRow=.2
  #       )
  # graphics.off()

  return(hm)
}

# Download and load Genes/GO terms association from BioMart
getGenesFromGOterms <- function(goids_list, dataset='mmusculus_gene_ensembl', gofilename='./inst/extdata/Annotations/gene_goterms.txt'){

  require(dplyr) # summarise_each

  if(is.na(gofilename)){
    gofilename = '/tmp/gene_goterms.txt'
    # Run system command from R
    system(paste0("wget -O ",gofilename," 'http://www.ensembl.org/biomart/martservice?query=<?xml version=\"1.0\" encoding=\"UTF-8\"?><!DOCTYPE Query><Query  virtualSchemaName = \"default\" formatter = \"CSV\" header = \"0\" uniqueRows = \"1\" count = \"\" datasetConfigVersion = \"0.6\" ><Dataset name = \"", dataset, "\" interface = \"default\" ><Filter name = \"go_parent_term\" value = \"",
      paste(c(goids_list), collapse=","),
      "\"/><Attribute name = \"external_gene_name\" /><Attribute name = \"name_1006\" /></Dataset></Query>'"))
  } else if(!file.exists(gofilename)){
    print('Input file does not exist.')
    return(0)
  }

  # Load genes 
  biomart.df = read.csv(gofilename, header=F, col.names=c('gene_name', 'GO_term'),stringsAsFactors = F)

  return(sort(unique(biomart.df$gene_name)))
}

seuratDispersionZScore <- function(data, nBin = 20){
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

getDispersedGenes = function(readcounts, zscore_threshold=-Inf, nBin=20){
  dispZscore = seuratDispersionZScore(readcounts, nBin = nBin)
  return(rownames(readcounts)[which(dispZscore >= zscore_threshold)])
}

getHighGenes = function(readcounts, mean_threshold=0, mean_log_threshold=-Inf){
  return(rownames(readcounts)[which(rowMeans(readcounts) >= mean_threshold & rowMeans(log(readcounts+1)) >= mean_log_threshold)])
}

plotRCPerSample <- function(eSet, filename='/tmp/readcounts.pdf', width=4, height=7){

  # check read counts per samples

  pData(eSet)$colSums = colSums(exprs(eSet))

  pdf(filename, width=width, height=height)
  p <- ggplot(pData(eSet), aes(factor(cells_samples), colSums)) + 
      geom_violin()+geom_boxplot(width=.1) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
    labs(x = NULL, y = "Total read counts per cell",
           fill = NULL) +
    coord_flip(ylim=c(.9*min(pData(eSet)$colSums), 1.1*max(pData(eSet)$colSums)))
  print(p)
  graphics.off()

  pData(eSet)$colSums = log10(1+colSums(exprs(eSet)))

  pdf(paste0(tools::file_path_sans_ext(filename), '_log10.', tools::file_ext(filename)), width=width, height=height)
  p <- ggplot(pData(eSet), aes(factor(cells_samples), colSums)) + 
      geom_violin()+geom_boxplot(width=.1) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
    labs(x = NULL, y = "Total read counts per cell (log10)",
           fill = NULL) +
    coord_flip(ylim=c(.9*min(pData(eSet)$colSums), 1.1*max(pData(eSet)$colSums)))
  print(p)
  graphics.off()

  pData(eSet)$colSums <- NULL

}

# from a comment in SO
# http://stackoverflow.com/questions/19610805/r-memory-usage-of-each-variable
list_obj_sizes<-function(
    list_obj=ls(envir=.GlobalEnv)){
  sizes<-sapply(list_obj,
          function(n){
            object.size(get(n))
          },
          simplify=FALSE)
  print(sapply(sizes[order(-as.integer(sizes))],function(s)format(s,unit='auto')))
}

getCommitShortName <- function(repository_path=find.package("Antler")) {


  repo <- git2r::repository(repository_path)

  if(length(git2r::status(repo)$unstaged) > 0){
    if(any(unlist(git2r::status(repo)$unstaged) %in% list.files('R', full.names=T))){
      stop('some changes have not been commited. Commit first!')
    }
  }

  fullcommit = git2r::commits(repo)[[1]]@sha
  shortcommit = substr(fullcommit, 0, 7)

  return(shortcommit)
}

getClusterColors <- function(v=1){
  if(v==1){
    c(RColorBrewer::brewer.pal(8, "Set2"), RColorBrewer::brewer.pal(12, "Set3"), RColorBrewer::brewer.pal(9, "Set1"), '#EBFF00','#8F00FF','#14FF00','#00FFFF','#24FF00','#FF00E6','#FF000F','#FF3D00','#3300FF','#00B2FF','#1400FF','#0075FF','#FF7A00','#9E00FF','#FF8A00','#CCFF00','#0066FF','#FF0F00','#FF9900','#05FF00','#FF004C','#FFF500','#FF005C','#0047FF','#FFD600','#00A3FF','#2400FF','#FF00D6','#00FF57','#00FF47','#FF001F','#FF00B8','#ADFF00','#00FF75','#70FF00','#FF0099','#00FFA3','#00FF29','#FFE500','#0500FF','#FF6B00','#00E0FF','#00FF1A','#0085FF','#00FFF0','#FF007A','#BD00FF','#0057FF','#8000FF','#00FFE0','#00FF85','#FF00C7','#FFA800','#0094FF','#FFC700','#DBFF00','#EB00FF','#00FFC2','#BDFF00','#33FF00','#0038FF','#FF0000','#FA00FF','#FF003D','#00C2FF','#4200FF','#52FF00','#00FFD1','#9EFF00','#00FF94','#8FFF00','#FF5C00','#FFB800','#FF002E','#000AFF','#FF00A8','#CC00FF','#00FF38','#00FF66','#00FFB3','#00FF0A','#FF00F5','#00D1FF','#00F0FF','#5200FF','#0019FF','#80FF00','#7000FF','#0029FF','#FF1F00','#FF008A','#FF4D00','#AD00FF','#FAFF00','#DB00FF','#42FF00','#6100FF','#FF006B','#FF2E00','#61FF00')
  } else if(v==2){
    c(RColorBrewer::brewer.pal(9, "Set1"), RColorBrewer::brewer.pal(12, "Set3"), '#EBFF00','#8F00FF','#14FF00','#00FFFF','#24FF00','#FF00E6','#FF000F','#FF3D00','#3300FF','#00B2FF','#1400FF','#0075FF','#FF7A00','#9E00FF','#FF8A00','#CCFF00','#0066FF','#FF0F00','#FF9900','#05FF00','#FF004C','#FFF500','#FF005C','#0047FF','#FFD600','#00A3FF','#2400FF','#FF00D6','#00FF57','#00FF47','#FF001F','#FF00B8','#ADFF00','#00FF75','#70FF00','#FF0099','#00FFA3','#00FF29','#FFE500','#0500FF','#FF6B00','#00E0FF','#00FF1A','#0085FF','#00FFF0','#FF007A','#BD00FF','#0057FF','#8000FF','#00FFE0','#00FF85','#FF00C7','#FFA800','#0094FF','#FFC700','#DBFF00','#EB00FF','#00FFC2','#BDFF00','#33FF00','#0038FF','#FF0000','#FA00FF','#FF003D','#00C2FF','#4200FF','#52FF00','#00FFD1','#9EFF00','#00FF94','#8FFF00','#FF5C00','#FFB800','#FF002E','#000AFF','#FF00A8','#CC00FF','#00FF38','#00FF66','#00FFB3','#00FF0A','#FF00F5','#00D1FF','#00F0FF','#5200FF','#0019FF','#80FF00','#7000FF','#0029FF','#FF1F00','#FF008A','#FF4D00','#AD00FF','#FAFF00','#DB00FF','#42FF00','#6100FF','#FF006B','#FF2E00','#61FF00')
  } else if(v==3){
    c(RColorBrewer::brewer.pal(9, "Set1"), RColorBrewer::brewer.pal(12, "Set3"), RColorBrewer::brewer.pal(8, "Set2"), '#EBFF00','#8F00FF','#14FF00','#00FFFF','#24FF00','#FF00E6','#FF000F','#FF3D00','#3300FF','#00B2FF','#1400FF','#0075FF','#FF7A00','#9E00FF','#FF8A00','#CCFF00','#0066FF','#FF0F00','#FF9900','#05FF00','#FF004C','#FFF500','#FF005C','#0047FF','#FFD600','#00A3FF','#2400FF','#FF00D6','#00FF57','#00FF47','#FF001F','#FF00B8','#ADFF00','#00FF75','#70FF00','#FF0099','#00FFA3','#00FF29','#FFE500','#0500FF','#FF6B00','#00E0FF','#00FF1A','#0085FF','#00FFF0','#FF007A','#BD00FF','#0057FF','#8000FF','#00FFE0','#00FF85','#FF00C7','#FFA800','#0094FF','#FFC700','#DBFF00','#EB00FF','#00FFC2','#BDFF00','#33FF00','#0038FF','#FF0000','#FA00FF','#FF003D','#00C2FF','#4200FF','#52FF00','#00FFD1','#9EFF00','#00FF94','#8FFF00','#FF5C00','#FFB800','#FF002E','#000AFF','#FF00A8','#CC00FF','#00FF38','#00FF66','#00FFB3','#00FF0A','#FF00F5','#00D1FF','#00F0FF','#5200FF','#0019FF','#80FF00','#7000FF','#0029FF','#FF1F00','#FF008A','#FF4D00','#AD00FF','#FAFF00','#DB00FF','#42FF00','#6100FF','#FF006B','#FF2E00','#61FF00')
  }
}

# Dimension Reduction (TopCorr_DR) - helper functions
clusterGenes <- function(
                         readcount,
                         readcount_logscaled,
                         corr_matrix,
                         num_gms=NULL,
                         num_gms_min = NULL,
                         num_gms_max = NULL,
                         numcores = 1,
                         display=TRUE,
                         basename=NULL,
                         verbose=TRUE,
                         clustering_method="ward.D2"
                         ){

    if(verbose){cat('\nRun hierarchical clustering (dissimilarity from euclidean distance between z-scored log-leveled gene expressions)')}

    # data.eucldistance = as.dist(fastEuclideanDist(readcount_logscaled))
    # hc = hclust(data.eucldistance, method=clustering_method)

    data.distance = as.dist(1-corr_matrix)
    hc = hclust(data.distance, method=clustering_method)

    if(!is.null(num_gms)){
      verbose.log = paste0('\nGene HC / number of gene modules given a priori: ', num_gms)
      if(verbose){cat(verbose.log)}
    } else {

      # if(verbose){cat('Estimate number of gene modules\n')}
      log = ''

      num_gms = getOptimalModuleNumber_HC(
                        readcount_logscaled, hc, 1e10, numcores,
                        display=display,
                        pdf_plot_path=basename,
                        old_heuristic=FALSE
                        )
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


