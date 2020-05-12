
normalize_feature_values <- function(feature_values) {

  unique_values <- sort(unique(feature_values))

  if (length(unique_values) == 1) {
    return(rep(1, length(feature_values)))
  }

  if (is.numeric(unique_values)) {

    low  <- min(unique_values)
    high <- max(unique_values)

    return(unlist(lapply(
      feature_values,
      function(x) {
        (x - low) / (high - low)
      })))

  } else if (is.character(unique_values)) {
    
    return(unlist(lapply(
      feature_values,
      function(x) {
        (which(unique_values == x) - 1) / (length(unique_values) - 1)
      })))
  }
}

generate_sample_hsv_map <- function(
  pheno_features,
  seed                = 10,
  hue_shift           = .0) {

  pf <- pheno_features[!duplicated(pheno_features), , drop = FALSE]

  normalized_pheno_features <- apply(pf[, -4], 2, normalize_feature_values)

  if (!is.matrix(normalized_pheno_features)) {
    normalized_pheno_features <- t(as.matrix(normalized_pheno_features))
  }

  sample_colors <- setNames(
    apply(normalized_pheno_features, 1, function(x){
      hsv(x[1], x[2], x[3], 1)
    }),
    pf[, 4])

  return(sample_colors[sort(names(sample_colors))])
}

get_color_matrix <- function(
  cell_colors,
  expressionSet,
  cell_clusters,
  read_counts) {

  if (!identical(cell_colors, "none")) {
    cell_colors_matrix <- do.call(cbind.data.frame, lapply(
      cell_colors,
      function(x) {
        if (x %in% colnames(pData(expressionSet)))
          return(as.factor(pData(expressionSet)[, x]))
        if (x %in% cell_clusters$names()) 
          return(as.factor(cell_clusters$get(x)$cell_ids))
        if (x %in% rownames(read_counts))
          return(as.double(read_counts[x, ]))
        return(rep(NA, nrow(expressionSet)))
      }))
    colnames(cell_colors_matrix) <- cell_colors
    rownames(cell_colors_matrix) <- colnames(expressionSet)
    cell_colors_matrix <- cell_colors_matrix[, colnames(cell_colors_matrix), drop = FALSE]
  } else {
    cell_colors_matrix <- NA
  }
  return(cell_colors_matrix) 
}

colorize <- function(
  label,
  values,
  color_maps,
  outlier_quantile_t = NA) {

  if (label %in% names(color_maps)) {
  
    if (label %in% names(color_maps) & is.double(values)) {

      if (identical(outlier_quantile_t, NA)) {
        values <- values / max(values)
      } else {

        q <- quantile(x = value, outlier_quantile_t)
        values <- values / q
        values[which(values > 1)] <- 1
      }
      values <- as.integer(1+(length(color_maps[[label]])-1)*values)
    }
    cols <- color_maps[[label]][values]
 
  } else {
    # use same default colors as pheatmap
    cols <- scales::dscale(
      factor(1:length(unique(values))),
      scales::hue_pal(l = 75))[values]
  }
  setNames(
    cols,
    values)
}

plot_color_graph <- function(
  gf,
  lay,
  colors             = NA,
  label              = NA,
  shape              = "circle",
  vertex.size        = 8,
  rescale            = TRUE,
  vertex.frame.color = 'black',
  edge.width         = 1,
  edge.color         = 'grey',
  shuffle            = FALSE,
  vertex.plot.order  = NA,
  ...) {

  require(igraph)

  V(gf)$color=colors
 
  labels=NA
  if(!is.na(label[1])){
    labels=label
  }

  if(length(shape)==1){
    shape <- rep(shape, vcount(gf))
  }

  if(length(vertex.frame.color)==1){
    vertex.frame.color <- rep(vertex.frame.color, vcount(gf))
  }

  if(length(vertex.size)==1){
    vertex.size <- rep(vertex.size, vcount(gf))
  }

  # no rescaling for the graph, so the coordinates must belong to [-1, 1]**2
  lay <- t(apply(lay, 1, function(r){r - colMeans(lay)}))
  lay <- lay / max(abs(lay))

  # Change vertex z-order if needed
  if(!identical(vertex.plot.order, NA) | shuffle){

    if(shuffle)
      new_order = sample(vcount(gf))
    if(!identical(vertex.plot.order, NA))
      new_order = vertex.plot.order

    adj_mat = igraph::as_adjacency_matrix(gf, type="both", attr="weight", names=TRUE, sparse=TRUE)
    adj_mat <- adj_mat[new_order, new_order]
    gf <- igraph::graph_from_adjacency_matrix(adj_mat, mode="undirected", weighted="weight")

    lay <- lay[new_order, , drop=F]
    colors <- colors[new_order]
    vertex.size <- vertex.size[new_order]
    label <- label[new_order]
    shape <- shape[new_order]
    vertex.frame.color <- vertex.frame.color[new_order]    
  }

  igraph::plot.igraph(
    gf,
    layout             = lay,
    vertex.color       = colors,
    vertex.size        = vertex.size,
    vertex.label       = label,
    # vertex.label.dist  = .0,
    # vertex.label.color = 'grey',
    # vertex.label.cex   = .2,
    vertex.shape       = shape,
    vertex.frame.color = vertex.frame.color,
    edge.width         = edge.width,
    edge.color         = edge.color,
    rescale            = rescale,
    ...)
}

#' Build and name color palette vectors
#' 
#' @name brew_more_colors
#' @param ids character or integer vector. The names that will be associated with the colors.
#' @param palette_name character string. The name of the RColorBrewer palette that will be extrapolated if needed. See \code{RColorBrewer::brewer.pal.info} for list of available colors.
#' @export brew_more_colors
brew_more_colors <- function(
  ids,
  palette_name = rownames(
    RColorBrewer::brewer.pal.info)) {

  palette_name <- match.arg(palette_name)

  palette <- RColorBrewer::brewer.pal(
    RColorBrewer::brewer.pal.info[palette_name, "maxcolors"],
    palette_name)

  n <- length(ids)

  if (n <= length(palette) &
      RColorBrewer::brewer.pal.info[palette_name, "category"] == "qual") {
    return(setNames(
      palette[seq(n)],
      ids))
  } else {
    return(setNames(
      colorRampPalette(palette)(n),
      ids))
  }
}