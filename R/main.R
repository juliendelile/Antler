#' Interface of the Antler package (reference class).
#'
#' @export Antler
Antler <- setRefClass(

  Class = "Antler",

  fields = list(
      expressionSet       = "ExpressionSet",
      readcounts_raw      = "matrix",
      readcounts_norm     = "matrix",
      readcounts_smoothed = "matrix",
      gene_modules        = "GeneModules",
      cell_clusters       = "CellClusters",
      cellClusters        = "list",
      output_folder       = 'character',
      log_filepath        = 'character',
      num_cores           = 'numeric',
      processing_state    = 'character',
      write_report        = 'logical',
      favorite_genes      = "character",
      cell_state_graph    = "CellStateGraph",
      PT_graph            = "igraph",
      carver              = "Carver",
      color_maps          = "list"
      ),  

  methods = list(

    initialize = function(
      ...,
      expressionSet       = ExpressionSet(),
      readcounts_raw      = matrix(),
      readcounts_norm     = matrix(),
      readcounts_smoothed = matrix(),
      gene_modules        = GeneModules$new(
        antler_env = as.environment(.self)),
      cell_clusters       = CellClusters$new(
        antler_env = as.environment(.self)),
      cellClusters        = list('hclust' = NULL),
      output_folder       = './',
      log_filepath        = paste0(output_folder, '/analysis.log'),
      num_cores           = 1,
      processing_state    = character(0),
      write_report        = TRUE,
      favorite_genes      = character(0),
      cell_state_graph    = CellStateGraph$new(antler_env = as.environment(.self)),
      PT_graph            = igraph::make_ring(0),
      carver              = Carver$new(
        antler_env = as.environment(.self)),
      color_maps          = list()) {

        # Enable methods' autocompletion in R session
        methods <- getRefClass(class(.self))$methods()
        eval(parse(text=paste0(".self$", methods)))

        callSuper(
          ...,
          expressionSet       = expressionSet,
          readcounts_raw      = readcounts_raw,
          readcounts_norm     = readcounts_norm,
          readcounts_smoothed = readcounts_smoothed,
          gene_modules        = gene_modules,
          cell_clusters       = cell_clusters,
          cellClusters        = cellClusters,
          output_folder       = output_folder,
          log_filepath        = log_filepath,
          num_cores           = num_cores,
          processing_state    = processing_state,
          write_report        = write_report,
          favorite_genes      = favorite_genes,
          cell_state_graph    = cell_state_graph,
          PT_graph            = PT_graph,
          carver              = carver,
          color_maps          = color_maps)

        # Random number generator ("L'Ecuyer-CMRG") 
        # compatible with parallel::mclapply (see ?parallel::mclapply)
        RNGkind("L'Ecuyer-CMRG")

        # Create plot directory
        if (length(output_folder) == 1){
          dir.create(
            output_folder,
            showWarnings = FALSE,
            recursive    = TRUE)
        }
      },

    copy = function() {

      def <- .refClassDef
      value <- new(def)
      vEnv <- as.environment(value)
      selfEnv <- as.environment(.self)
      for (field in names(def@fieldClasses)) {
        current <- get(field, envir = selfEnv)
        if (is(current, "envRefClass")) 
            current <- current$copy(FALSE)
        assign(field, current, envir = vEnv)
      }
      value$gene_modules$antler_env     <- vEnv
      value$cell_clusters$antler_env    <- vEnv
      value$carver$antler_env           <- vEnv
      value$cell_state_graph$antler_env <- vEnv
      value
    },

    write_report_file = function() {
    write(
      paste0(
        "# Antler report file (last modification on ", 
        fdate(),
        ')\n\n',
        paste0(unlist(lapply(processing_state, function(r){paste0(r, '\n')})), collapse='')),
      file = log_filepath)
    }
  )
)

#' Print the current analysis log in the console.
#'
#' @name report
Antler$methods(
  report = function() {
    cat(paste0(unlist(lapply(processing_state, function(r) {paste0(r, '\n')})), collapse=''))
  }
)


#' Save dataset to file.
#'
#' Export the gene expression matrix and the table of phenotypic sample information as csv files.
#' @name save_dataset
#' @param folder_path character string. The path of the directory in which the csv files are written.
#' @param data_status character string. Specifies whether the exported expression levels are raw ("Raw"), normalized ("Normalized") or have been imputed ("Smoothed"). Default to "Raw".
Antler$methods(

  save_dataset = function(
    folder_path,
    data_status = c('Raw', 'Normalized', 'Smoothed')) {

    data_status <- match.arg(data_status)

    dir.create(
      folder_path,
      showWarnings = FALSE)

    write.table(
      read_counts(data_status = data_status),
      file      = file.path(folder_path, "assayData.csv"),
      quote     = FALSE,
      sep       = '\t',
      row.names = TRUE)
    write.table(
      pData(expressionSet),
      file      = file.path(folder_path, "phenoData.csv"),
      quote     = FALSE,
      sep       = '\t',
      row.names = TRUE)
    })

#' Import a RNA-seq dataset into Antler.
#' 
#' @name load_dataset
#' @param folder_path character string. The path of the directory in which the dataset files are stored. if NULL (default), the "assayData" and "phenoData" arguments must be specified.
#' @param assayData_filename a character string specifying the name of the file containing the expression values.
#' @param phenoData_filename character string. The name of the file containing the cell metadata.
#' @param featData_filename  a character string specifying the name of the file containing the gene metadata.
#' @param assayData Matrix. A matrix storing the expression values.
#' @param phenoData Matrix. A matrix storing the sample metadata.
#' @param data_status character string. Specifies whether the loaded expression are raw ("Raw"), normalized ("Normalized") or have been imputed ("Smoothed"). Default to "Raw".
#' @param default_gene_names character string. Label describing the source of the gene names. Default to "ensembl_gene_id".
Antler$methods(
  load_dataset = function(
    folder_path         = NULL,
    phenoData_filename = "phenoData.csv",
    assayData_filename = "assayData.csv",
    featData_filename  = "featureData.csv",
    assayData          = NULL,
    phenoData          = NULL,
    data_status        = c('Raw', 'Normalized', 'Smoothed'),
    default_gene_names = 'ensembl_gene_id') {

    data_status <- match.arg(data_status)

    if (is.null(folder_path)) {
      if (is.null(assayData) & is.null(phenoData)) {
        stop('input dataset is missing')
      }
      expressionSet <<- ExpressionSet(
        assayData = assayData,
        phenoData = phenoData
        )
      processing_state <<- c(processing_state, paste0(fdate(), ' Dataset loaded from provided data.frames (', data_status, ')'))
    } else {
      expressionSet <<- ExpressionSet(
        assayData = as.matrix(read.table(
          paste0(folder_path, '/', assayData_filename),
          header           = TRUE,
          sep              = "\t",
          row.names        = 1,
          as.is            = TRUE,
          stringsAsFactors = FALSE,
          check.names      = FALSE)),
        phenoData = new(
          "AnnotatedDataFrame",
          data=read.table(
            paste0(folder_path, '/', phenoData_filename),
            header           = TRUE,
            sep              = "\t",
            row.names        = 1,
            as.is            = TRUE,
            stringsAsFactors = FALSE,
            check.names      = FALSE))
      )

      if (file.exists(file.path(folder_path, featData_filename))) {
        fd <- read.table(
          file        = file.path(folder_path, featData_filename),
          sep         = '\t',
          row.names   = 1,
          header      = TRUE,
          check.names = FALSE)
        for(n in colnames(fd)){
          fData(expressionSet)[,n] <<- fd[,n]
        }
      }

      processing_state <<- c(
        processing_state,
        paste0(fdate(), " '", data_status, "' dataset loaded from '", folder_path, "'"))
      if (write_report) write_report_file()
    }

    pData(expressionSet)$cells_samples <<- rep('sample', ncol(exprs(expressionSet)))

    fData(expressionSet)[, default_gene_names]  <<- rownames(exprs(expressionSet))
    fData(expressionSet)$current_gene_names     <<- fData(expressionSet)[, default_gene_names]

    if (data_status == "Raw") {
      readcounts_raw      <<- exprs(expressionSet)
      readcounts_norm     <<- matrix()
      readcounts_smoothed <<- matrix()
    } else if (data_status == "Normalized") {
      readcounts_norm     <<- exprs(expressionSet)
      readcounts_raw      <<- matrix()
      readcounts_smoothed <<- matrix()
    } else if (data_status == "Smoothed") {
      readcounts_smoothed <<- exprs(expressionSet)
      readcounts_raw      <<- matrix()
      readcounts_norm     <<- matrix()
    }

    if (write_report) write_report_file()
  }
)

#' Define the cells sample names and map the cell colors to some of their metadata.
#' 
#' @name define_cell_samples
#' @param sample_feature a character vector of cell features. Must be included in the column names of the cell metadata matrix.
#' @param color_mapping a character vector having the same length as 'sample_features' and whose values are all different and either 'Hue', 'Saturation' or 'Value'. If NULL (default) the colors are chosen randomly.
#' @param seed an integer setting the seed for random color generation.
Antler$methods(
  define_cell_samples = function(
    sample_features,
    color_mapping   = NULL,
    seed            = 0) {

    if (!all(sample_features %in% colnames(pData(expressionSet)))) {
      
      stop(paste0(
        "Some sample feature(s) are not included in the input dataset: '",
        paste0(setdiff(sample_features, colnames(pData(expressionSet))), collapse = "', '"),
        "'. Check feature names and/or add missing feature(s) to the sample meta-data structure pData(...) <- ...\n"))
    } else {

      pData(expressionSet)$cells_samples <<- apply(
        pData(expressionSet),
        1,
        function(x) {
          paste0(
            lapply(sample_features, function(s) {
              paste(substr(s, 0, 1), x[s]) 
            }),
            collapse = " / "
          )
        })

      if (is.null(color_mapping)) {

        samples <- sort(unique(pData(expressionSet)$cells_samples))
        set.seed(seed)
        color_maps[["cells_samples"]] <<- setNames(
          sample(colours(), length(samples), replace = FALSE),
          samples)

        processing_state <<- c(
          processing_state,
          paste0(fdate(), " Set 'cells_samples' colors using random colors"))
        if (write_report) write_report_file()

      } else {

        if (!is.character(color_mapping) |
          length(unique(color_mapping)) != length(sample_features) |
          any(!color_mapping %in% c('Hue', 'Saturation', 'Value'))) {

          stop("'color_mapping' must be a character vector having the same length as 'sample_features' and whose values are all different and either 'Hue', 'Saturation' or 'Value'. For example, \"sample_features = c('feature_1', 'feature_2')\" and \"color_mapping = c('Hue', 'Value')\"")
        }

        pheno_features <- data.frame(
          "Hue"        = rep(1, num_cells()),
          "Saturation" = rep(1, num_cells()),
          "Value"      = rep(1, num_cells()))

        pheno_features[, color_mapping] <- pData(expressionSet)[, sample_features]

        pheno_features <- cbind(
          pheno_features,
          "cells_samples" = pData(expressionSet)$cells_samples)

        color_maps[["cells_samples"]] <<- generate_sample_hsv_map(pheno_features)

        processing_state <<- c(
          processing_state,
          paste0(fdate(), " Set 'cells_samples' colors using the '", sample_features, "' feature"))
        if (write_report) write_report_file()
      }
    }
  }
)

#' Add a new color map.
#' 
#' @name add_color_map
#' @param name a character string specifying the name of the new color map.
#' @param content a named vector of colors stored as character strings.
Antler$methods(
  add_color_map = function(name, content) {

    if (missing(name))
      stop("The argument 'name' is required")

    if (!is.character(name))
      stop("The argument 'name' must be a character string")

    color_maps[[name]] <<- content
  }
)

#' Store a list of favorite gene names.
#' 
#' @name set_favorite_genes
#' @param gene_list a character containing gene names.
Antler$methods(
  set_favorite_genes = function(gene_list=NULL){

    unrecognized_genes <- gene_list[which(!gene_list %in% gene_names())]

    if (length(unrecognized_genes) > 0){
      cat(paste0("Some genes are not included in the dataset: ", paste0(unrecognized_genes, collapse=', '), "\n"))
    }

    favorite_genes <<- unique(gene_list[which(gene_list %in% gene_names())])
  }
)

#' Returns a gene expression matrix.
#'
#' @name read_counts
#' @param data_status a character string specifying whether the returned expression levels are raw ("Raw"), normalized ("Normalized") or have been imputed ("Smoothed"). Default to "Raw".
#' @return a Matrix containing the gene expression levels with genes as rows and cells as columns.
Antler$methods(
  read_counts = function(data_status = c('Raw', 'Normalized', 'Smoothed')) {

    data_status <- match.arg(data_status)

    if (data_status == "Raw") {
      readcounts = readcounts_raw
    } else if (data_status == "Normalized") {
      readcounts = readcounts_norm
    } else if (data_status == "Smoothed") {
      readcounts = readcounts_smoothed
    }

    if (identical(readcounts, matrix())){
      stop(paste0(data_status, " dataset has not been loaded or calculated."))
    }

    rownames(readcounts) <- fData(expressionSet)$current_gene_names

    return(readcounts)
  }
)

#' Returns the current gene names.
#' 
#' @name gene_names
#' @return a character vector.
Antler$methods(
  gene_names = function(){
    return(fData(expressionSet)$current_gene_names)
    }
)

#' Returns the current number of genes.
#'
#' @name num_genes
#' @return an integer number.
Antler$methods(
  num_genes = function(){
    return(nrow(exprs(expressionSet)))
  }
)

#' Returns the current number of cells.
#' 
#' @name num_cells
#' @return an integer number.
Antler$methods(
num_cells = function(){
  return(ncol(exprs(expressionSet)))
  }
)

#' Returns the current cell names.
#' 
#' @name cell_names
#' @return a character vector.
Antler$methods(
  cell_names = function(){
    return(colnames(exprs(expressionSet)))
  }
)

#' Returns a contingency table with the counts per type of cell samples.
#' 
#' @name sample_size
#' @return a table.
Antler$methods(
  sample_size = function(){
    return(table(pData(expressionSet)$cells_samples))
  }
)

#' Render various plots indicating the tallied gene and cell count distributions.
#' 
#' @name plot_QC
#' @param data_status a character string specifying whether the considered expression levels should be raw ("Raw"), normalized ("Normalized") or have been imputed ("Smoothed"). Default to "Raw".
#' @param basename a character string specifying the prefix of the rendered plots. Defaults to NULL.
#' @param reads_type a character string indicating whether the internal dataset stores generic "read" (default) or "UMI" elements.
#' @param feature_1 an optional character string indicating a metadata feature used to aggregate the cells.
#' @param feature_2 an optional character string indicating a metadata feature used to aggregate the cells.
#' @return write plots as pdf files in \strong{Antler}'s \code{output_folder} output directory.
Antler$methods(
  plot_QC = function(
    data_status    = c('Raw', 'Normalized', 'Smoothed'),
    basename       = NULL,
    reads_type     = c("read", "UMI"),
    feature_1      = NULL,
    feature_2      = NULL) {

    data_status <- match.arg(data_status)
    reads_type  <- match.arg(reads_type)

    if (!is.null(feature_1) &
      !feature_1 %in% colnames(pData(expressionSet))) {
      stop(paste0("'", feature_1, "' is not a valid feature."))
    }

    if (!is.null(feature_2)) {
      if (!feature_2 %in% colnames(pData(expressionSet))) {
        stop(paste0("'", feature_2, "' is not a valid feature."))
      }
    }

    plotnames <- c()

    stat_df <- data.frame(
      totalread_per_cell           = colSums(read_counts(data_status = data_status)),
      totalexpressedgenes_per_cell = colSums(read_counts(data_status = data_status) > 0))

    plotname  <- paste0(output_folder, "/", basename, ifelse(is.null(basename), "", "_"), "library_sizes.pdf")
    plotnames <- c(plotnames, plotname)
    pdf(plotname, width = 4, height = 4, useDingbats = FALSE)
    p <- ggplot(stat_df, aes(x = totalread_per_cell)) +
      geom_histogram(bins = 50) +
      # scale_x_log10() +
      ylab("Cell count") +
      xlab(paste0("Number of ", reads_type, "s per cell")) +
      ggtitle(paste0("Library size distribution")) + 
      theme_classic()
    print(p)
    graphics.off()

    plotname  <- paste0(output_folder, "/", basename, ifelse(is.null(basename), "", "_"), "expressed_genes_per_cell.pdf")
    plotnames <- c(plotnames, plotname)
    pdf(plotname, width = 4, height = 4, useDingbats = FALSE)
    p <- ggplot(stat_df, aes(x = totalexpressedgenes_per_cell)) +
      geom_histogram(bins = 50) +
      # scale_x_log10() +
      ylab("Cell count") +
      xlab(paste0("Number of expressed genes per cell")) +
      ggtitle(paste0("Expressed genes distribution")) + 
      theme_classic()
    print(p)
    graphics.off()

    if (!is.null(feature_1)) {
      stat_df[, feature_1] <- factor(
        pData(expressionSet)[, feature_1],
        levels = sort(unique(pData(expressionSet)[, feature_1])))

      if (!is.null(feature_2)) {
        stat_df[, feature_2] <- factor(
          pData(expressionSet)[, feature_2],
          levels = sort(unique(pData(expressionSet)[, feature_2]))) 
      } else {
        stat_df[, feature_2] <- factor(rep("", num_cells()))
      }

      rplots <- lapply(
        sort(unique(stat_df[, feature_1])),
        function(t){
          ggplot(stat_df[which(stat_df[, feature_1] == t), ]) +
            geom_freqpoly(
              aes_string(
                x        = "totalread_per_cell",
                linetype = feature_2),
              position = "identity",
              bins     = 30) +
            # scale_x_log10() +
            theme(legend.position="none") +
            xlab(paste0("Total ", reads_type, " count")) +
            ylab("Cell count") + 
            ggtitle(paste(feature_1, t)) +
            scale_colour_gradient(limits=c(.5,1)) + 
            labs(linetype = feature_2) +
            theme_classic()
        })

      if (nlevels(stat_df[, feature_2]) > 1) {

        g_legend <- function(a.gplot){
          pdf(file=NULL)
          tmp <- ggplot_gtable(ggplot_build(a.gplot))
          graphics.off()
          leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
          legend <- tmp$grobs[[leg]]
          return(legend)}

        lgd <- g_legend(rplots[[1]])

        rplots <- lapply(rplots, function(x) {
          x + theme(legend.position="none")
          })

        rplots[[length(rplots) + 1]] <- lgd
      } else {
        rplots <- lapply(rplots, function(x) {
          x + theme(legend.position="none")
          })
      }

      plotname  <- paste0(
        output_folder, "/", basename,
        ifelse(is.null(basename), "", "_"),
        "library_sizes_I_by_",
        paste0(c(feature_1, feature_2), collapse = "-"),
        ".pdf")
      plotnames <- c(plotnames, plotname)
      
      pdf(plotname, width = 7, height = 4, useDingbats = FALSE)
      numplots = length(rplots)
      numcol = min(3, numplots)
      numrow = ceiling(numplots / 3)
      gridExtra::grid.arrange(
        grobs         = rplots,
        layout_matrix = matrix(
          seq(numcol * numrow),
          ncol = numcol,
          byrow = T),
        top           = grid::textGrob(
          paste0("Library size distribution"),
          x = 0,
          hjust = 0))
      graphics.off()


      fill_feature <- ifelse(is.null(feature_2), feature_1, feature_2)

      if (fill_feature %in% c(character(0), names(color_maps))) {
        sfm <- scale_fill_manual(
          values = color_maps[[fill_feature]],
          name = fill_feature)
      } else {
        sfm <- scale_fill_manual(
          values = getOption("antler.feature.colors"),
          name = fill_feature)
      }

      p <- stat_df %>%
        count(!!!syms(c(feature_1, feature_2))) %>%
        ggplot() +
        geom_bar(
          aes_string(
            feature_1,
            weight = "n",
            fill   = fill_feature),
          color = "black",
          size  = .3) +
        ylab("Number of cells") +
        xlab(tools::toTitleCase(feature_1)) +
        ggtitle("Cell counts") +
        theme_classic() +
        sfm

      if (fill_feature == feature_1 | nlevels(stat_df[, feature_2]) == 1)
        p <- p + theme(legend.position="none")

      plotname  <- paste0(
        output_folder, "/", basename,
        ifelse(is.null(basename), "", "_"),
        "cell_count_by_",
        paste0(c(feature_1, feature_2), collapse = "-"),
        ".pdf")
      plotnames <- c(plotnames, plotname)
      pdf(plotname, width=5, height=3, useDingbats=FALSE)
      print(p)
      graphics.off()

      p <- stat_df %>% 
        ggplot() +
        geom_violin(
          aes_string(
            x    = feature_1,
            y    = "totalread_per_cell",
            fill = fill_feature)) +
        ylab(paste0(reads_type, " per cell")) +
        xlab(tools::toTitleCase(paste0(c(feature_1, feature_2), collapse = " x "))) +
        ggtitle("Library size distribution") +
        theme_classic() +
        sfm

      plotname  <- paste0(
        output_folder, "/", basename,
        ifelse(is.null(basename), "", "_"),
        "library_sizes_II_by_",
        paste0(c(feature_1, feature_2), collapse = "-"),
        ".pdf")
      plotnames <- c(plotnames, plotname)
      pdf(plotname, width=5, height=5, useDingbats=FALSE)
      print(p)
      graphics.off()

      
      p <- stat_df %>% 
        ggplot() +
        geom_violin(
          aes_string(
            x    = feature_1,
            y    = "totalexpressedgenes_per_cell",
            fill = fill_feature)) +
        ylab("Gene per Cell") +
        xlab(tools::toTitleCase(paste0(c(feature_1, feature_2), collapse = " x "))) +
        ggtitle("Expressed genes distribution") +
        theme_classic() +
        sfm

      if (fill_feature == feature_1 | nlevels(stat_df[, feature_2]) == 1)
        p <- p + theme(legend.position="none")

      plotname  <- paste0(
        output_folder, "/", basename,
        ifelse(is.null(basename), "", "_"),
        "expressed_genes_count_by_",
        paste0(c(feature_1, feature_2), collapse = "-"),
        ".pdf")
      plotnames <- c(plotnames, plotname)
      pdf(plotname, width=5, height=5, useDingbats=FALSE)
      print(p)
      graphics.off()

    }

    processing_state <<- c(processing_state, paste0(fdate(), ' Save ', reads_type, ' count statistics plots (', paste0(plotnames, collapse=', '), ')'))
    if (write_report) write_report_file()
  }
)

#' Replace ensembl gene ids with gene names
#' 
#' When the dataset is loaded with genes specified by their ensembl ids, this function maps these ids to their associated gene names. 
#' @name convert_ensembl_id_to_name
#' @param biomart_dataset a character string specifying the BioMart dataset to use for the conversion. The list of available datasets can be listed with \link[biomaRt]{listDatasets} from the "biomaRt" package. Default to 'mmusculus_gene_ensembl'.
#' @param gene_ids_mapping_file a character string. If NULL, the conversion is performed from a gene ids/names mapping obtained using the biomaRt package (Default). This requires an internet connection. If not NULL, this argument must specify the location of a csv file with header containing the ensembl gene ids in the first column and the gene names in the second column.
Antler$methods(

  convert_ensembl_id_to_name = function(
    biomart_dataset     = 'mmusculus_gene_ensembl',
    gene_ids_mapping_file = NULL){

    if (is.null(gene_ids_mapping_file)){

      # It happens quite often that the Biomart server is not reachable.
      ensembl = tryCatch({
          useMart("ensembl", dataset=biomart_dataset)
      }, warning = function(w) {
          NA
      }, error = function(e) {
          cat('Biomart server can not be reached, please use the local file version.')
          NA
      })

      if (identical(ensembl, NA)) {
        cat('Biomart server can not be reached, please use the local file version with the "gene_ids_mapping_file" argument.')
        return()
      }
      # listMarts()    # to see which database options are present
      # ensembl=useMart("ensembl")
      # ensembl=useMart("ensembl", dataset="mmusculus_gene_ensembl")
      # listDatasets(ensembl)     # function to see which datasets are present in ensembl
      # ensembl=useDataset("mmusculus_gene_ensembl",mart=ensembl) 
      # listFilters(ensembl)  # check which filters are available
      # listAttributes(ensembl) # check attributes are available to select.More information on ensembl data base
      gene_names_ids.df <- biomaRt::getBM(
        attributes = c('ensembl_gene_id', 'external_gene_name'),
        mart= ensembl)

      # some genomes does not have annotated "external_gene_name" (e.g. 'ggallus')
      missing_attr <- which(gene_names_ids.df[['external_gene_name']] == '')
      gene_names_ids.df[missing_attr, 'external_gene_name'] <- gene_names_ids.df[missing_attr, 'ensembl_gene_id']

    } else {

      gene_names_ids.df <- read.csv(
        gene_ids_mapping_file,
        stringsAsFactors = FALSE,
        header           = TRUE)

      colnames(gene_names_ids.df) <- c('ensembl_gene_id', 'external_gene_name')
    }
    
    # ! first column must have unique names (default: unique "ensembl_id")

    renamed_genes_ids = which(fData(expressionSet)$ensembl_gene_id %in% gene_names_ids.df$ensembl_gene_id)

    renamed_genes = fData(expressionSet)$ensembl_gene_id[renamed_genes_ids]

    rownames(gene_names_ids.df) <- gene_names_ids.df$ensembl_gene_id

    fData(expressionSet)$current_gene_names[renamed_genes_ids] <<- make.unique(gene_names_ids.df[renamed_genes ,'external_gene_name'])

    fData(expressionSet)[, 'external_gene_name'] <<- fData(expressionSet)$current_gene_names

    processing_state <<- c(
          processing_state,
          paste0(fdate(), ' Load new gene names ("external_gene_name" from ', biomart_dataset, ')'))
    if (write_report) write_report_file()
  }
)

#' Exclude cells from the current dataset
#' 
#' @name remove_cells
#' @param ids an integer vector storing the indices of the  cells to remove.
Antler$methods(

  remove_cells = function(ids = NULL){

    if (length(ids) > 0){

      expressionSet <<- expressionSet[, -ids]

      if (!identical(readcounts_raw, matrix()))
        readcounts_raw <<- readcounts_raw[, -ids]
      if (!identical(readcounts_norm, matrix()))
        readcounts_norm <<- readcounts_norm[, -ids]
      if (!identical(readcounts_smoothed, matrix()))
        readcounts_smoothed <<- readcounts_smoothed[, -ids]

      if (!all(dim(cell_state_graph$dist_matrix) == c(1,1))) {
        cell_state_graph$remove_cells(ids)
      }

      processing_state <<- c(processing_state, paste0(fdate(), ' Removed ', length(ids), ' cells (', num_cells(), ' left)'))
      if (write_report) write_report_file()
    }
  }
)


#' Exclude cells belonging to one or more clusters
#' 
#' If some genes were expressed only in the excluded cells, they will also be removed from the dataset.
#' @name remove_clusters
#' @param cell_clusters_name character string. The name of the cell clusters.
#' @param ids an integer vector storing the indices of the clusters to remove.
Antler$methods(

  remove_clusters = function(
    cell_clusters_name = cell_clusters$names(),
    ids
    ) {

    cell_clusters_name <- match.arg(cell_clusters_name)

    if (!all(ids %in% unique(cell_clusters$get(cell_clusters_name)$cell_ids)))
      stop("Provided indices are not compatible with the specified clusters.")

    excluded_ids <- which(cell_clusters$get(cell_clusters_name)$cell_ids %in% ids)

    remove_cells(excluded_ids)

    processing_state <<- c(
      processing_state,
      paste0(fdate(), ' Exclude cells belonging to "', length(ids), '" cluster(s) (', length(excluded_ids), ' cells)'))
    if (write_report) write_report_file()

    exclude_unexpressed_genes(
      min_cells   = 1,
      data_status = "Raw")
  }
)


#' Exclude genes from the current dataset.
#' 
#' @name remove_genes 
#' @param ids an integer vector specifying the indices of the genes to remove.
#' @param names a character vector specifying the names of the genes to remove.
Antler$methods(

  remove_genes = function(ids = NULL, names = NULL) {

    if ((identical(ids, NULL) + identical(names, NULL)) %in% c(0, 2))
      stop("Either 'ids' or 'names' must be specified.")

    if (identical(ids, NULL))
      ids <- which(gene_names() %in% names)

    if (length(ids) > 0){
    
      expressionSet <<- expressionSet[-ids, , drop=F]

      if (!identical(readcounts_raw, matrix()))
        readcounts_raw <<- readcounts_raw[-ids, , drop=F]
      if (!identical(readcounts_norm, matrix()))
        readcounts_norm <<- readcounts_norm[-ids, , drop=F]
      if (!identical(readcounts_smoothed, matrix()))
        readcounts_smoothed <<- readcounts_smoothed[-ids, , drop=F]
      if (length(favorite_genes) > 0)
        favorite_genes <<- favorite_genes[favorite_genes %in% gene_names()]

      processing_state <<- c(processing_state, paste0(fdate(), ' Removed ', length(ids), ' genes (', num_genes(), ' left)'))
      if (write_report) write_report_file()
    }
  }
)

#' Print highly-dispersed genes.
#' 
#' The genes are allocated to a bin based on their average log-level of expression, then with each bin the dispersion (variance over mean of the log-levels) is z-scored. Genes having a sufficiently high dispersion z-scores are returned.
#' This method replicates \code{\link[Seurat]{FindVariableFeatures}} from the Seurat package.
#' @name dispersed_genes
#' @param zscore_threshold a numeric value indicating the zcored dispersion threshold above which the gene names are returned. Default to 0.
#' @param num_bins a integer value indicating the number of bins used to calculate z-score into. Default to 20.
#' @param data_status character string. Specifies whether the gene expression levels used for calculation are raw ("Raw"), normalized ("Normalized") or have been imputed ("Smoothed"). Default to "Raw".
#' @param invert logical. If FALSE (default) genes with z-scored dispersion higher than the threshold are returned. If FALSE, the complementary set is returned.
#' @return a character vector of genes showing a sufficient level of dispersion.
Antler$methods(

  dispersed_genes = function(
    zscore_threshold = 0,
    num_bins         = 20,
    data_status      = c('Raw', 'Normalized', 'Smoothed'),
    invert           = FALSE) {

    data_status <- match.arg(data_status)

    dispZscore <- seurat_dispersion_zscore(
      read_counts(data_status),
      nBin = num_bins)

    dispersed_genes_ids <- which(dispZscore >= zscore_threshold)

    if (invert) {
      dispersed_genes_ids <- setdiff(seq(num_genes()), dispersed_genes_ids)
    }

    return(gene_names()[dispersed_genes_ids])
  }
)

#' Select highly-dispersed genes.
#' 
#' The genes are allocated to a bin based on their average log-level of expression, then with each bin the dispersion (variance over mean of the log-levels) is z-scored. Genes not having a sufficiently high dispersion z-scores are excluded from the dataset.
#' This method replicates \code{\link[Seurat]{FindVariableFeatures}} from the Seurat package.
#' @name select_dispersed_genes
#' @param zscore_threshold a numeric value indicating the zcored dispersion threshold above which the gene names are returned. Default to 0.
#' @param num_bins a integer value indicating the number of bins used to calculate z-score into. Default to 20.
#' @param data_status character string. Specifies whether the gene expression levels used for calculation are raw ("Raw"), normalized ("Normalized") or have been imputed ("Smoothed"). Default to "Raw".
Antler$methods(
  
  select_dispersed_genes = function(
    zscore_threshold = 0,
    num_bins         = 20,
    data_status      = c('Raw', 'Normalized', 'Smoothed')) {

    remove_genes(
      names = dispersed_genes(
        zscore_threshold = zscore_threshold,
        num_bins         = num_bins,
        data_status      = data_status,
        invert           = TRUE))
  }
)

#' Exclude genes with low expression.
#' 
#' Genes must be expressed in a given number of cell with a certain level to be kept.
#' @name exclude_unexpressed_genes
#' @param min_cells a integer indicating the minimal number of cells a gene must be expressed in. Default to 3.
#' @param min_level a numeric value indicating the minimal level of expression required. Default to 0.
#' @param verbose a logical indicating whether to print the list of excluded genes. Default to FALSE.
#' @param data_status character string. Specifies whether the gene expression levels used for calculation are raw ("Raw"), normalized ("Normalized") or have been imputed ("Smoothed"). Default to "Raw".
Antler$methods(

  exclude_unexpressed_genes = function(
    min_cells   = 3,
    min_level   = 0,
    verbose     = FALSE,
    data_status = c('Raw', 'Normalized', 'Smoothed')){

    readcounts = read_counts(data_status=data_status)

    unexpressed_genes_ids = which(rowSums(readcounts > min_level) < min_cells)

    log = paste0("Genes expressed in less than ", min_cells, " cells: ", paste(fData(expressionSet)$current_gene_names[unexpressed_genes_ids], collapse=" "))
    if (verbose){cat(paste0(log, '\n'))}
    processing_state <<- c(processing_state, log)
    if (write_report) write_report_file()

    remove_genes(unexpressed_genes_ids)

  }
)

#' Remove outlier genes and cells.
#' 
#' This method excludes cells with a total read counts outside of a specified interval, cells having less than a certain number of expressed genes, and genes expressed in less than a given number of cells.
#' @name remove_outliers
#' @param lowread_thres numeric value. Cells will be excluded if their total expression level is below this threshold. Default to -Inf.
#' @param lowread_high numeric value. Cells will be excluded if their total expression level is higher than this threshold. Default to Inf.
#' @param min_cells integer value. Genes expressed in less than this value will be excluded.
#' @param min_genes integer value. If a cell has a total number of expressed genes lower than this value, it is excluded.
#' @param verbose a logical indicating whether to print the dataset dimension at each filtering step. Default to TRUE.
#' @param data_status character string. Specifies whether the gene expression levels used for calculation are raw ("Raw"), normalized ("Normalized") or have been imputed ("Smoothed"). Default to "Raw".
Antler$methods(

  remove_outliers = function(
    lowread_thres  = -Inf,
    highread_thres = Inf,
    min_cells      = 10,
    min_genes      = 1000,
    verbose        = TRUE,
    data_status    = c('Raw', 'Normalized', 'Smoothed')) {

    log.0 = paste0('Pre-filtering (dataset dimension:', nrow(exprs(expressionSet)), ' ', ncol(exprs(expressionSet)), ')')
    if (verbose){cat(paste0(log.0, '\n'))}
    processing_state <<- c(processing_state, log.0)
    
    # Cell filtering: select cells counting at least X reads
    if (!is.infinite(lowread_thres)){
      remove_cells(ids =which(apply(read_counts(data_status=data_status), 2, sum) <= lowread_thres))
    }
    log.1 = paste0('Filter I - select cells with at least ', lowread_thres, ' reads (dataset dimension:', nrow(exprs(expressionSet)), ' ', ncol(exprs(expressionSet)), ')')
    if (verbose){cat(paste0(log.1, '\n'))}
    processing_state <<- c(processing_state, log.1)

    # Cell filtering: select cells counting less than X reads
    if (!is.infinite(highread_thres)){
      remove_cells(ids =which(apply(read_counts(data_status=data_status), 2, sum) >= highread_thres))
    }
    log.1b = paste0('Filter II - select cells with less than ', highread_thres, ' reads (dataset dimension:', nrow(exprs(expressionSet)), ' ', ncol(exprs(expressionSet)), ')')
    if (verbose){cat(paste0(log.1b, '\n'))}
    processing_state <<- c(processing_state, log.1b)

    # Gene filtering: select genes expressed in at least N cells
    if (min_cells > 0){
      remove_genes(ids=which(rowSums(read_counts(data_status=data_status)>0) < min_cells))
    }
    log.2 = paste0('Filter III - select genes expressed in at least ', min_cells,' cells (dataset dimension:', nrow(exprs(expressionSet)), ' ', ncol(exprs(expressionSet)), ')')
    if (verbose){cat(paste0(log.2, '\n'))}
    processing_state <<- c(processing_state, log.2)

    # Cell filtering: select cells expressing at least N genes
    if (min_genes > 0){
      remove_cells(ids = which(colSums(read_counts(data_status=data_status)>0) < min_genes))
    }
    log.3 = paste0('Filter IV - select cells expressing at least ', min_genes, ' genes (dataset dimension:', nrow(exprs(expressionSet)), ' ', ncol(exprs(expressionSet)), ')')
    if (verbose){cat(paste0(log.3, '\n'))}
    processing_state <<- c(processing_state, log.3)
    if (write_report) write_report_file()
  }
)

#' Normalize gene expression levels.
#' 
#' The raw gene expression level are normalized according to one the following methods:
#' \itemize{
#'   \item 'CPM' or 'Count-Per-Million': The raw gene expression levels are divided by
#' the total read count of the cell and multiplied by 1e6. This method corrects for 
#' sequencing depth.
#'   \item 'FPKM': The CPM levels are divided by the gene length in kilobase. This
#' method corrects for sequencing depth and gene length.
#'   \item 'MR' or 'Median-of-ratios': The raw gene expression levels are divided by cell size
#' factors equals to the median ratio of gene counts relative to geometric mean per gene.
#' This method corrects for library size and RNA composition bias and is a
#' reimplementation of DESeq2' `estimateSizeFactors()` with the "poscounts" estimator
#' which deals with genes with some zeros. 
#' }
#' @name normalize
#' @param method a character string. The method to use for normalization: either 'CPM', 'Count-Per-Million', 'FPKM', 'MR', 'Median-of-ratios' (Default to 'CPM').
#' @param gene_length if the normalization method is 'FPKM', a numeric vector containing the length of each genes (Default to NULL).
Antler$methods(

  normalize = function(
    method      = c('Count-Per-Million', 'CPM', 'FPKM', 'Median-of-ratios', 'MR'),
    gene_length = NULL) {

    method <- match.arg(method)

    if (identical(readcounts_raw, matrix())){
      stop("'normalize' operates on 'Raw' datasets. Please provide one.")
    }

    if (method == 'Count-Per-Million' | method == "CPM") {
      readcounts_norm <<- 1e6 * sweep(readcounts_raw, 2, colSums(readcounts_raw), FUN="/")
    } else if (method == "FPKM") {
      readcounts_norm <<- 1e6 * sweep(readcounts_raw, 2, colSums(readcounts_raw), FUN="/")
      readcounts_norm <<- 1000 * sweep(readcounts_norm, 1, gene_length, "/")
    } else if (method %in% c("Median-of-ratios", "MR")) {
      
      if (any(rowSums(readcounts_raw) == 0))
        stop("'Median-of-ratios' does not work with null genes. Use remove_outliers() to exclude null genes first.")
      # pseudo-referance sample (geometric mean per gene)
      geom_means <- apply(
        readcounts_raw,
        1,
        function(x) {
          if (all(x == 0)) { 0 } else { exp( sum(log(x[x > 0])) / length(x) ) }
        })

      # medians of ratios
      log_geom_means <- log(geom_means)

      medians <- apply(
        readcounts_raw,
        2, 
        function(x) {
          exp(
            median(
              (log(x) - log_geom_means)[is.finite(log_geom_means) & x > 0]
            )
          )
        })
       medians <- medians/exp(mean(log(medians)))

      # divide row by medians
      readcounts_norm <<- t(t(readcounts_raw) / medians)
    }
    
    processing_state <<- c(processing_state, paste0(fdate(), ' Gene levels normalized (', method, ')'))
    if (write_report) write_report_file()
  }
)

#' Exclude genes from their name pattern.
#' 
#' @name exclude_exotic_genes
#' @param pattern a character vector listing the patterns of gene names to exclude. Default to 'Gm[0-9].*', '.*-ps.*', 'ENSMUSG.*', '.*Rik' and '.*-rs*'.
#' @param verbose a logical indicating whether to print the list of excluded genes. Default to TRUE.
Antler$methods(

  exclude_exotic_genes = function(
    exotic_patterns = c('Gm[0-9].*', '.*-ps.*', 'ENSMUSG.*', '.*Rik', '.*-rs*'),
    verbose         = TRUE) {

    exotic_gene_ids = unique(unlist(lapply(exotic_patterns, function(ep){grep(ep, fData(expressionSet)$current_gene_names)})))

    remove_genes(ids=exotic_gene_ids)

    log = paste0('Excluded ', length(exotic_gene_ids), ' "exotic" gene(s) whose name matches the following patterns: ', paste0(exotic_patterns, collapse=', '), ' (dataset dimension: ', nrow(exprs(expressionSet)), ' ', ncol(exprs(expressionSet)), ')')

    if (verbose){cat(paste0(log, '\n'))}
    processing_state <<- c(processing_state, log)
    if (write_report) write_report_file()

  }
)

#' Heatmap of the transcriptomic profile.
#' 
#' Once the gene modules and optionally the cell clusters have been calculated, this function produces a heatmap showing the gene expression levels.
#' @name plot_transcriptome_summary
#' @param gene_module_name. character string. The name of the gene module list shown along the rows. 
#' @param cell_clusters_name character string. The name of the cell clusters shown along the columns. If 'none' (default), the columns will be ordered as in the current dataset.
#' @param cell_side_colors character vector indicating the metadata to show in the column annotation bar. Names can be either one of the column names of the phenotypic metadata structure (see \code{pData(antler$expressionSet)}), names of a cluster entries, or gene names. Default to 'cells_samples'.
#' @param cell_side_colors_hide_legend character vector. The name of the phenotypic metadata to hide in the legend. Default to "none".
#' @param gene_side_colors character vector. The name of the feature metadata to show in the row annotation bar. Default to "none".
#' @param data_status character string. Specifies whether the gene expression levels displayed are raw ("Raw"), normalized ("Normalized") or have been imputed ("Smoothed"). Default to "Normalized".
#' @param gene_level_func character string. Specifies whether the gene expression levels are rendered as is ('none'), log-transformed ('log') or z-scored log-transformed ('logscaled', default).
#' @param heatmap_palette_centered a character vector storing a color map to use if the \code{gene_level_func} is set to 'none' or 'log'.
#' @param heatmap_palette a character vector storing a color map to use if the \code{gene_level_func} is set to logscaled'.
#' @param color_limit_value numeric value. The upper limit of the color map if \code{gene_level_func} is set to 'none' or 'log' (default NULL for no set limit), the absolute value of the upper and lower limits of the color map if \code{gene_level_func} is set to 'logscaled' (default to 2).
#' @param gene_labels_filter character vector. Either a list of genes whose intersection with the gene modules will be displayed as paragraph, or 'all' for rendering each individual gene name on each row, or 'all_gms' for rendering each individual gene name plus the gene module id on each row, or 'naked' for no row labels at all.
#' @param num_gene_per_line if \code{gene_labels_filter} is a list of genes, the maximum number of gene names to render on a single paragraph row. Default to 8.
#' @param fontsize integer value indicating the plot font size. Default to 10.
#' @param subsampling_ratio numeric value. The ratio of cells to render, see parameter \code{subsampling_min_cells} for an exception. Default to 1.
#' @param subsampling_min_cells integer value. If a cell cluster object is specified with \code{cell_clusters_name}, the minimum of cells to render per cell cluster. Default to 20.
#' @param file_settings a list of lists indicating the type and dimensions of the output plots. For each plot type, a list containing three arguments must be specified: 'type' (either 'pdf' or 'png'), 'width' and 'height', and a fourth one, 'res', can be set to define the resolution of png files. The unit of the pdf device dimension is inches and pixels for png device. The png resolution is 200 ppi (default) so the appearance of a w10xh15 pdf file is the same as a w2000xh3000 png file Default to \code{list(list(type='pdf', width=15, height=15))}  
#' @param suffix a character string that will be appended to the file name. Default to NULL. 
#' @param ... Extra arguments passed to \code{\link[pheatmap]{pheatmap}}.
Antler$methods(

  plot_transcriptome_summary = function(
    gene_modules_name            = gene_modules$names(),
    cell_clusters_name           = c('none', cell_clusters$names()),
    cell_side_colors             = "cells_samples",
    cell_side_colors_hide_legend = "none",
    gene_side_colors             = "none",
    data_status                  = c('Normalized', 'Raw', 'Smoothed'),
    gene_level_func              = c('logscaled', 'log', 'none'),
    heatmap_palette_centered     = colorRampPalette(
      c("#191d73", "white", "#ed7901"))(n = 1000),
    heatmap_palette              = colorRampPalette(
      c("#0464DF", "#FFE800"))(n = 1000),
    color_limit_value            = if (gene_level_func == 'logscaled') 2.0 else NULL,
    gene_labels_filter           = gene_names(),
    num_gene_per_line            = 8,
    fontsize                     = 10,
    subsampling_ratio            = 1,
    subsampling_min_cells        = 20,
    file_settings                = list(list(type='pdf', width=15, height=15)),
    suffix                       = NULL,
    ...) {

    gene_modules_name  <- match.arg(gene_modules_name, several.ok = TRUE)
    cell_clusters_name <- match.arg(cell_clusters_name, several.ok = TRUE)

    # Define cell colors
    cell_side_colors_choices <- c(
        "none",
        colnames(pData(expressionSet)),
        cell_clusters$names())
    cell_side_colors   <- match.arg(
      cell_side_colors,
      choices    = cell_side_colors_choices,
      several.ok = TRUE)
    cell_side_colors_hide_legend   <- match.arg(
      cell_side_colors_hide_legend,
      choices    = cell_side_colors_choices,
      several.ok = TRUE)
    # reorder cell annotations if hiding legend
    if (!identical(cell_side_colors_hide_legend, "none"))
      cell_side_colors <- c(
        setdiff(cell_side_colors,cell_side_colors_hide_legend),
        cell_side_colors_hide_legend)

    if (!identical(cell_side_colors, "none")) {
      cell_side_colors_matrix <- get_color_matrix(
        cell_side_colors,
        expressionSet,
        cell_clusters,
        read_counts(data_status))
      cell_side_colors_matrix <- cell_side_colors_matrix[, rev(seq(ncol(cell_side_colors_matrix))), drop = FALSE] # first row on top
    } else {
      cell_side_colors_matrix <- NA
    }

    # Define gene colors
    gene_side_colors   <- match.arg(
      arg        = gene_side_colors,
      choices    = c(
        "none",
        colnames(fData(expressionSet)),
        gene_modules$names()),
      several.ok = TRUE)
    data_status        <- match.arg(data_status, several.ok = TRUE)

    if (!identical(gene_side_colors, "none")) {
      gene_side_colors_matrix <- do.call(cbind.data.frame, lapply(
        gene_side_colors,
        function(x) {
          if (x %in% colnames(fData(expressionSet))) {
            if (typeof(fData(expressionSet)[, x]) == "double") {
              return(fData(expressionSet)[, x])
            } else {
              return(as.factor(fData(expressionSet)[, x]))
            }
          }
          if (x %in% gene_modules$names()) 
            return(factor(
              setNames(
                rep(
                  names(gene_modules$get(gmn)),
                  unlist(lapply(gene_modules$get("unbiasedGMs"), length))),
                unlist(gene_modules$get(gmn)))[gene_names()],
              levels = names(gene_modules$get(gmn))))
          }))
      colnames(gene_side_colors_matrix) <- gene_side_colors
      rownames(gene_side_colors_matrix) <- gene_names()
      gene_side_colors_matrix <- gene_side_colors_matrix[, rev(colnames(gene_side_colors_matrix)), drop = FALSE]
    } else {
      gene_side_colors_matrix <- NA
    }

    for(ds in data_status){

      readcounts = read_counts(data_status=ds)
      
      filenames = c()
      for(norm in gene_level_func){

        if (norm == 'logscaled') {
          plot_data  <- t(scale(t(log(readcounts+1)), center = TRUE, scale = TRUE))
          hm_palette <- heatmap_palette_centered
          breaks     <- if (is.null(color_limit_value)) NA else seq(-color_limit_value, color_limit_value, length.out = length(hm_palette) + 1)
        } else {              
          hm_palette <- heatmap_palette
          breaks     <- if (is.null(color_limit_value)) NA else seq(0, color_limit_value, length.out = length(hm_palette) + 1)

          if (norm == "log") {
            plot_data <- log(readcounts+1)
          } else if (norm == "none") {
            plot_data <- readcounts
          }
        }

        for(gmn in gene_modules_name){

          # Define row names  
          mod_names = if (is.null(names(gene_modules$get(gmn)))) seq(length(gene_modules$get(gmn))) else names(gene_modules$get(gmn))
          
          show_rownames <- TRUE
          labels_row <- NULL
          if (identical(gene_labels_filter, "all_gms")) {
            labels_row <- paste(
              unlist(gene_modules$get(gmn)),
              paste0('(', rep(mod_names, unlist(lapply(gene_modules$get(gmn), length))), ')')
              )
          } else if (identical(gene_labels_filter, "all")) {
            labels_row <- unlist(gene_modules$get(gmn))
          } else if (identical(gene_labels_filter, "naked")) {
            show_rownames <- FALSE
          } else {

            labels_row  <- rep("", length(unlist(gene_modules$get(gmn))))
            gm_halfsize <- unlist(lapply(
              gene_modules$get(gmn),
              function(l) as.integer(length(l)/2)))
            gm_index    <- head(
              cumsum(c(1, unlist(lapply(gene_modules$get(gmn), length)))),
              -1)

            labels_row[gm_index + gm_halfsize] <- lapply(
              seq(length(gene_modules$get(gmn))),
              function(i) {
                list       <- intersect(gene_modules$get(gmn)[[i]], gene_labels_filter)
                list_split <- split(list, ceiling(seq_along(list)/num_gene_per_line))
                paste(mod_names[i], ':', paste0(lapply(list_split, function(ll){paste0(ll, collapse=' ')}), collapse='\n'))
              })
          }

          for(ccn in cell_clusters_name){

            for(fs in file_settings){

              if (subsampling_ratio == 1) {
                
                silent_pheatmap <- pheatmap::pheatmap(
                  mat               = plot_data[unlist(gene_modules$get(gmn)), , drop = FALSE],
                  color             = hm_palette,
                  breaks            = breaks,
                  cluster_cols      = if(identical(ccn, 'none')) FALSE else cell_clusters$get(ccn)$res,
                  cutree_cols       = if(identical(ccn, 'none')) NA else cell_clusters$get(ccn)$num_clusters,
                  annotation_row    = gene_side_colors_matrix,
                  annotation_col    = cell_side_colors_matrix,
                  show_colnames     = FALSE,
                  cluster_rows      = FALSE,
                  gaps_row          = head(cumsum(unlist(lapply(gene_modules$get(gmn), length))), -1),
                  annotation_colors = color_maps,
                  fontsize          = fontsize,
                  width             = fs$width,
                  height            = fs$height,
                  labels_row        = labels_row,
                  show_rownames     = show_rownames,
                  silent            = TRUE,
                  ...)
              } else {

                if (identical(ccn, 'none')) {
                  cell_ids <- sort(sample(
                    seq(num_cells()),
                    subsampling_ratio * num_cells(),
                    replace = FALSE))
                } else {
                  # The hclust clusters are resampled and the matrix is
                  # rendered by reordering the columns, without dendrogram
                  # but with calculated gaps.
                  cell_ids <- lapply(
                    seq(cell_clusters$get(ccn)$num_clusters),
                    function(i) {
                      # cell_ids <- which(cell_clusters$get(ccn)$cell_ids == i)
                      cell_ids <- cell_clusters$get(ccn)$res$order[which(
                        cell_clusters$get(ccn)$cell_ids[cell_clusters$get(ccn)$res$order] == i)]
                      new_size = subsampling_ratio * length(cell_ids)
                      if (new_size < subsampling_min_cells)
                        new_size <- min(subsampling_min_cells, length(cell_ids))
                      cell_ids[sort(sample(seq(length(cell_ids)), new_size, replace = FALSE))]
                      })  
                }

                silent_pheatmap <- pheatmap::pheatmap(
                  mat               = plot_data[unlist(gene_modules$get(gmn)), unlist(cell_ids), drop = FALSE],
                  color             = hm_palette,
                  breaks            = breaks,
                  cluster_cols      = FALSE,
                  annotation_row    = gene_side_colors_matrix,
                  annotation_col    = cell_side_colors_matrix[unlist(cell_ids), , drop = FALSE],
                  show_colnames     = FALSE,
                  cluster_rows      = FALSE,
                  gaps_row          = head(
                    cumsum(unlist(lapply(gene_modules$get(gmn), length))),
                    -1),
                  gaps_col          = head(cumsum(unlist(lapply(cell_ids, length))), -1),
                  annotation_colors = color_maps,
                  fontsize          = fontsize,
                  width             = fs$width,
                  height            = fs$height,
                  labels_row        = labels_row,
                  show_rownames     = show_rownames,
                  silent            = TRUE,
                  ...)
              }

              filename  <- paste0(
                output_folder,
                '/Transcriptome_summary_', gmn, '_', ccn, '_', ds, '_', norm,
                if (is.null(suffix)) NULL else "_", suffix, 
                '.',
                ifelse(fs$type == "cairo_pdf", "cairo.pdf", fs$type))
              filenames <- c(filenames, filename)

              if (fs$type == 'pdf') {
                pdf(
                  filename,
                  width       = fs$width,
                  height      = fs$height,
                  useDingbats = FALSE)
              } else if (fs$type == 'png') {
                png(
                  filename,
                  width  = fs$width,
                  height = fs$height,
                  res    = if(is.null(fs$res)) 200 else fs$res,
                  type   = "cairo-png")
              } else if (fs$type == 'cairo_pdf') {
                cairo_pdf(
                  filename,
                  width  = fs$width,
                  height = fs$height,
                  fallback_resolution = if(is.null(fs$res)) 200 else fs$res)
              } else {
                stop("'Wrong file type in 'file_settings' argument.")
              }

              gtable = silent_pheatmap$gtable

              if (!identical(cell_side_colors_hide_legend, "none")) {
                annotation_legend_grob <- which(gtable$layout$name == "annotation_legend")
                grob_names <- names(gtable$grobs[[annotation_legend_grob]]$childrenOrder)
                for (n in cell_side_colors_hide_legend) {
                  for (gn in grep(n, grob_names, value = TRUE)) {
                    gtable$grobs[[annotation_legend_grob]] <- grid::removeGrob(
                      gtable$grobs[[annotation_legend_grob]],
                      gtable$grobs[[annotation_legend_grob]]$childrenOrder[gn])
                  }
                }
              }

              grid::grid.draw(gtable)
              # grid::grid.draw(silent_pheatmap$gtable)
              # plot(silent_pheatmap$gtable)
              graphics.off()
            }
          }
        }
      }
    }

    processing_state <<- c(processing_state, paste0(fdate(), ' Plot transcriptome summary heatmap(s) (', paste0(filenames, collapse=', '), ')'))

    if (write_report) write_report_file()
  }
)

#' Export internal expression set.
#' 
#' The exported dataset is written in the Antler object's output directory.
#' @name export_expressionset
#' @param exported_dataset. character string. The name of the directory containing the exported dataset. Default to 'exported_dataset'.
Antler$methods(

  export_expressionset = function(
    name = 'exported_dataset') {

    dir.create(file.path(output_folder, name), showWarnings = FALSE)

    readcounts <- read_counts(data_status='Raw')
    dimnames(readcounts) <- dimnames(exprs(expressionSet))
    exprs(expressionSet) <<- readcounts

    write.table(
      x         = readcounts,
      file      = file.path(output_folder, name, 'assayData.csv'), 
      sep       = '\t',
      row.names = TRUE,
      quote     = FALSE,
      col.names = NA)

    write.table(
      x         = pData(expressionSet),
      file      = file.path(output_folder, name, 'phenoData.csv'), 
      sep       = '\t',
      row.names = TRUE,
      quote     = FALSE,
      col.names = NA)     
  }
)

#' Reconstruct cell state lineage trajectories 
#' 
#' See carving cell trajectories vignette.
#' @name carve
#' @param gene_modules_name character string. The name of the list of gene modules whose weights will be optimized by this method.
#' @param optim_strategy character string. Specifies the algorithm used to optimize the gene module weights. Either "SA" for simulated annealing (default) or "GA" for genetic algorithm. Different processing plots will be generated depending of the selected strategy.
#' @param optim_objectives a character vector containing a subset of of either:
#' \describe{
#'    \item{'timecorr'}{to favor a population of candidate cell lineages maximizing their Spearman correlation with the 'timepoint' phenotypic metadata.}
#'    \item{'complexity'}{to favor a population of candidate cell lineages minimizing the number of leaves of the reconstructed tree.}
#'    \item{'smoothness'}{to favor a population of candidate cell lineages minimizing the irregularity of the gene expression pattern along the differentiation trajectories.}
#' }
#' Multiple elements can be specified if the optimization strategy is set to 'GA' but a single one 'SA'. Default to \code{c('timecorr')}.
#' @param GA_selection_objectives character vector. If 'GA' is selected as optimization strategy, the final objectives used to select the optimal cell lineage. Any of the \code{optim_objectives} options can be selected. If more than one objective is specified, the geometric average of these objectives will be used to select the optimal solution. Default to "timecorr".
#' @param GA_num_generations integer value. If 'GA' is selected as optimization strategy, the number of generation of the genetic algorithm. Default to 200.
#' @param GA_pop_size integer value. If 'GA' is selected as optimization strategy, the size of population of the genetic algorithm. Default to 200.
#' @param SA_control List. If 'SA' is selected as optimization strategy, the list controlling the simulated annealing algorithm behavior. See documentation of the underlying \code{\link[GenSA]{GeneSA}} function.
#' @param SA_maxit interger value, the maximum number of iteration of the simulated annealing algorithm. Default to 500.
#' @param cluster_size_limit integer value. The maximum cell cluster size. If NULL (default) the \code{cluster_size_limit_ratio} argument is used instead.
#' @param cluster_size_limit_ratio numeric value. Specifies the maximum cell cluster size as a ratio of the total number of cells. Default to 0.02.
#' @param clustering_timepoint_first logical. Whether to pre-cluster the cells by 'timepoint' or not (default).
#' @param canalize logical. Whether to constraint the cell state graph using a gaussian kernel obtained from the optimal cluster level graph. Default to TRUE.
#' @param kernel_sigma numeric value. Defines the standard deviation of the gaussian kernel used to canalize the trajectories. Default to 1.
#' @param allow_multiple_roots (Experimental) logical. Whether to allow reconstructed cell lineages to start from multiple root states. Default to FALSE.
#' @param output_name character string. The name of the directory storing the processing plots. This directory will be created in the Antler's object output directory.
#' @param cell_colors character string. The name of the phenotypic metadata to show on the cell state graph plots. Default to 'cells_samples'.
#' @param seed integer value. Set the random number generator seed for reproducibility.

Antler$methods(

  carve = function(
    gene_modules_name,
    optim_strategy                       = c('SA', 'GA'),
    optim_objectives                     = c("timecorr"),
    GA_selection_objectives              = c("timecorr"),
    GA_num_generations                   = 200,
    GA_pop_size                          = 200,
    SA_maxit                             = 500,
    SA_control                           = list(
        smooth          = FALSE,
        maxit           = SA_maxit,
        nb.stop.improvement = 200,
        trace.mat       = TRUE,
        simple.function = FALSE,
        verbose         = TRUE,
        seed            = 1),
    cluster_size_limit                   = NULL,
    cluster_size_limit_ratio             = .02,
    clustering_timepoint_first           = FALSE,
    canalize                             = TRUE,
    kernel_sigma                         = 1,
    allow_multiple_roots                 = FALSE,
    output_name                          = paste0("carvings_", optim_strategy),
    cell_colors                          = "cells_samples",
    seed                                 = 1234) {

    cell_colors   <- match.arg(
      cell_colors,
      choices    = c(
        "none",
        colnames(pData(expressionSet))),
      several.ok = TRUE)

    optim_strategy <- match.arg(optim_strategy)

    carver$run(
      gene_modules_name                            = gene_modules_name,
      
      optim_strategy                               = optim_strategy,
      optim_objectives                             = optim_objectives, # currently only 1 objecitve can be selected for SA
      
      GA_solutions_selection_quantile_threshold = .8, # only valid for GA
      GA_solutions_selection_objectives         = GA_selection_objectives, # only valid for GA
      GA_solutions_selection_strategy           = "average", # "top", "average", only valid for GA
      
      GA_num_generations                        = seq(GA_num_generations), # only valid for GA
      GA_popsize                                = GA_pop_size, # only valid for GA
      GA_constraint                             = NULL, # NULL, "SumTo1",  # only valid for GA
      
      SA_control                                   = SA_control,
      seed                                         = seed,
      
      size_limit                                   = cluster_size_limit,
      size_limit_ratio                             = cluster_size_limit_ratio,
      
      extra_pd                                     = cell_colors,
      allow_multiple_roots                         = allow_multiple_roots,
      clustering_timepoint_first                   = clustering_timepoint_first,

      plot_folder_path                             = file.path(output_folder, output_name))

    transfer_carver_weights(
      gene_modules_name = gene_modules_name,
      from = "best")

    carving_dist_to_cell_dist(
      allow_multiple_roots = allow_multiple_roots,
      canalize             = canalize,
      kernel_sigma         = kernel_sigma,
      plot_folder_path     = file.path(output_folder, output_name))

    cell_state_graph$build(
      use_distance = TRUE,
      method       = "kNN",
      k            = 4)

    cell_state_graph$detect_communities(
      method             = "spin_glass_no_weight",
      n_communities      = NULL,
      seed               = seed)

    cell_state_graph$build_community_graph(
      min_occurence = 1)

    cell_state_graph$project(
      method    = 'gephiForceAtlas2',
      num_iter  = 50000,
      stop_mean = 5)

    if (!identical(cell_colors, "none")) {

      cell_state_graph$plot(
        filename           = paste0(
          output_name, '_cell_state_graph.pdf'),
        cell_colors        = cell_colors,
        labels             = NA,
        vertex.size        = 1000 / num_cells(),
        save_pdf           = TRUE,
        shuffle            = TRUE)
    }


    if (write_report) write_report_file()
  }
)

#' Smooth gene expression using the cell state graph
#' 
#' This function proposes two methods to smooth the gene expression levels: a simple average over the cells' local neighborhood or a more advanced diffusion-based imputation inspired by the "Rmagic" package.
#' @name smooth
#' @param method a character string among:
#' \describe{
#'    \item{'average'}{to smooth by averaging gene expression over the cell neighbors of the cell state graph (default).}
#'    \item{'diffusion'}{to smooth by diffusing the gene expression levels on the cell state graph.}
#' }
#' @param from character string. Specifies whether the original gene expression levels used are raw ("Raw"), normalized ("Normalized") or have been imputed ("Smoothed"). Default to "Normalized".
Antler$methods(

  smooth = function(
    method  = c("average", "diffusion"),
    from    = c("Normalized", "Raw", "Smoothed"),
    ...) {
    
    method <- match.arg(method)
    from   <- match.arg(from)

    if (igraph::vcount(cell_state_graph$graph) !=
        num_cells())
      stop("The cell state graph must be calculated first.")

    if (method == "average") {
      readcounts_smoothed <<- getGraphAverage(
        graph     = cell_state_graph$graph,
        data_orig = read_counts(data_status = from),
        ...)
    } else if (method == "diffusion") {
      readcounts_smoothed <<- getMagicDiffusionData(
        weighted_adjacency_matrix = cell_state_graph$dist_matrix *
          igraph::as_adjacency_matrix(cell_state_graph$graph, sparse = F),
        data_orig                 = read_counts(data_status = from),
        ...)
    }

    processing_state <<- c(processing_state, paste0(fdate(), "'", from, "' gene levels smoothed (", method, ")"))
    if (write_report) write_report_file()
  }
)

#' Plot cell state graph over different gaussian kernel widths
#'
#' This function must be called after the trajectories have been carved (see \code{carve}).
#' @name test_kernel_widths
#' @param kernel_sigmas a numeric vector specifying the candidate standard deviations of the gaussian kernel used to canalize the trajectories. Default to \code{c(1, 5, 10, 25, 50, 1000)}.
#' @param cell_colors character string. The name of the phenotypic metadata to show on the cell state graph plots. Default to 'cells_samples'.
#' @param plot_name character string. The filename of the output plot file. Default to 'candidate_kernel_params.pdf'.
Antler$methods(
  test_kernel_widths = function(
    kernel_sigmas        = c(.02, .05, .05, .08, .1, .2, .5, .8, 1),
    cell_colors          = "cells_samples",
    plot_name            = "candidate_kernel_params.pdf") {

    # works on object copy to prevent modifying current
    # state
    antler_copy = .self$copy()

    pdf(
      file.path(output_folder, plot_name),
      width = 10,
      height = 10)

    n_row <- as.integer(sqrt(length(kernel_sigmas)))
    n_col <- ceiling(length(kernel_sigmas) / n_row)
    par(mfrow = c(n_row, n_col))

    for (kernel_sigma in kernel_sigmas) {

      antler_copy$carving_dist_to_cell_dist(
        allow_multiple_roots = carver$allow_multiple_roots,
        canalize             = TRUE,
        kernel_sigma         = kernel_sigma)

      antler_copy$cell_state_graph$build(
        use_distance = TRUE,
        method       = "kNN",
        k            = 4)

      # antler_copy$cell_state_graph$detect_communities(
      #   method             = "spin_glass_no_weight",
      #   n_communities      = NULL)

      # antler_copy$cell_state_graph$build_community_graph(
      #   min_occurence = 1)

      antler_copy$cell_state_graph$project(
        method    = 'gephiForceAtlas2',
        num_iter  = 50000,
        stop_mean = 5)

      antler_copy$cell_state_graph$plot(
        cell_colors        = setNames(cell_colors, paste("Sigma: ", kernel_sigma)),
        labels             = NA,
        vertex.size        = 1000 / antler_copy$num_cells(),
        save_pdf           = FALSE,
        shuffle            = TRUE)
    }
    graphics.off()
  }
)

#' Reconstruct pseudotime tree
#' 
#' This function create a pseudotime tree structure and position each cell on one of its branches.
#' @name build_pseudotime_tree
#' @param ordering character string. The name of the phenotypic metadata use to define the pseudotime ordering. Default to 'timepoint'.
#' @param process_plots logical value. Whether to show processing plots. Default to TRUE.
#' @param basename a character string indicating the path and prefix of the processing plots. Default to 'Pt_tree'.
Antler$methods(

  build_pseudotime_tree = function(
    ordering              = "timepoint",
    process_plots         = TRUE,
    basename              = "PT_tree") {

    ordering_vals <- pData(expressionSet)[, ordering]

    communities_ordering.df <- data.frame(
      communities_id = pData(expressionSet)$community_id,
      ordering_vals  = ordering_vals)

    # "average" community ordering value
    cell_state_graph$communities_df$mean_ordering <<- plyr::ddply(
        communities_ordering.df,
        ~communities_id,
        summarise,
        mean=quantile_mean(ordering_vals, .2))$mean

    # Step 1: Identify landmark cells #####################################################

    cat("Identify landmark cells\n")

    # Process each subgraph independently

    num_subgraphs <- igraph::components(cell_state_graph$graph)$no

    excluded_community_ids <- c()

    cell_state_graph$communities_df$type             <<- rep(NA, nrow(cell_state_graph$communities_df))
    cell_state_graph$communities_df$landmark_cell_id <<- rep(NA, nrow(cell_state_graph$communities_df))

    for (subgraph_id in seq(num_subgraphs)){
      subgraph_cell_ids <- which(pData(expressionSet)$subgraph_id == subgraph_id)

      subgraph_community_ids <- which(cell_state_graph$communities_df$subgraph_id == subgraph_id)

      subgraph <- delete_vertices(cell_state_graph$graph, cell_names()[-subgraph_cell_ids])

      # Cells and communities are ordered according to the "ordering" vector
      # (usually cell timepoint is used)

      starting_community_id = subgraph_community_ids[which(
        cell_state_graph$communities_df$mean_ordering[subgraph_community_ids] == min(cell_state_graph$communities_df$mean_ordering[subgraph_community_ids])
        )]

      if (length(starting_community_id) > 1){
        # select community with smallest degree
        starting_community_id = starting_community_id[which.min(degree(cell_state_graph$communities_mstree)[starting_community_id])]
      }

      community_mst_pseudotime <- as.numeric(distances(
        cell_state_graph$communities_mstree,
        v = starting_community_id,
        to = subgraph_community_ids,
        weights=E(cell_state_graph$communities_mstree)$weight))

      # identify start cells from earliest smoothed ordering values
      smooth_ordering = getGraphAverage(
        graph     = subgraph,
        data_orig = matrix(
          ordering_vals[subgraph_cell_ids],
          nrow     = 1,
          dimnames = list("ordering", V(subgraph)$name)),
        depth     = 1,
        rescale   = .99
        )[1,]

      # Start cells -------------------------------------------------------------------------

      # Select start cell as the one having the lowest "smoothed" ordering value.
      subgraph_start_cellname <- names(which.min(smooth_ordering))
      start_cell_id <- which(cell_names() == subgraph_start_cellname)

      # End cells ---------------------------------------------------------------------------

      # Identify end cells as the farthest cells from the start cell on the subgraph (higher excentricity).
      # In the community minimum spanning tree, each leaf community contains a single end cell (except from the root community)
      # We define a community path as the path relating the root community to one of the leaf communities.
      # The end cells are defined as the farthest cells on the graph associated with each community path.
      # If in a path the end cells do not belong to the leaf community, this path is removed from the analysis (malformed trajectories)
      # When a leaf community is removed, the parent community may become a leaf community and define a new community path. We 
      # redo the end cell identification after each community removal, until convergence (no excluded community).

      changes <- TRUE

      while (changes) {

        changes <- FALSE

        # If only one community remains, we discard it
        if (length(setdiff(subgraph_community_ids, excluded_community_ids)) == 1) {
          excluded_community_ids <- unique(c(excluded_community_ids, subgraph_community_ids))
          start_cell_id <- NULL
          break
        }

        # The leaf communities are the communities associated of degree 1 which have not been excluded
        remaining_community_mstree <- igraph::delete.vertices(
          cell_state_graph$communities_mstree,
          v = excluded_community_ids)
        V(remaining_community_mstree)$name <- setdiff(seq(cell_state_graph$communities_number), excluded_community_ids)

        leaf_community_ids <- setdiff(
          intersect(
            as.integer(names(which(degree(remaining_community_mstree)==1))),
            subgraph_community_ids),
          c(starting_community_id, excluded_community_ids))

        # The community paths relate the root community to the leaf communities
        community_paths <- igraph::shortest_paths(
          remaining_community_mstree,
          from = as.character(starting_community_id),
          to   = as.character(leaf_community_ids))$vpath
        
        end_cell_ids <- lapply(seq(length(community_paths)), function(i){

          community_path = as.integer(community_paths[[i]]$name)

          # Get graph associated with the current community path
          subgraph_community_path <- igraph::delete.vertices(cell_state_graph$graph,
                                            v=setdiff(
                                                seq(num_cells()),
                                                unlist(cell_state_graph$communities_df$cell_ids[community_path])
                                              )
                                            )

          # id of start cell in previous graph
          subgraph_start_cells_id <- which(V(subgraph_community_path)$name == subgraph_start_cellname)

          # Calculate shortest paths between start cell and all cells in the current graph
          sps <- shortest_paths(
            subgraph_community_path,
            from    = subgraph_start_cells_id,
            mode    = c("all"),
            to      = V(subgraph_community_path),
            weights = E(subgraph_community_path)$distance)$vpath
          
          sps_length <- unlist(lapply(sps, length))

          # use vertex names instead of id to relate back to original graph
          subgraph_end_cellnames <- V(subgraph_community_path)$name[which(sps_length == max(sps_length))]

          # Get end cells id in global dataset
          curr_end_cell_ids <- which(cell_names() %in% subgraph_end_cellnames)
          
          # If one of the end cell does not belong to the leaf community associated with the current community path, 
          # we exclude that path from the analysis
          if (any(pData(expressionSet)$community_id[curr_end_cell_ids] != tail(community_path, n=1))) {
            
            print(paste0(
              "Excluding path ", paste0(community_path, collapse=' '),
              " because at least one of the end cells do not belong to the leaf community (",
              tail(community_path, n=1), ")"))

            end_cell_id <- NULL
          } else {

            end_cell_id <- curr_end_cell_ids

            # Pick a single end cell randomly if more than one
            if (length(curr_end_cell_ids) > 1){
              end_cell_id <- sample(curr_end_cell_ids, 1)
            }
          }

          end_cell_id
        })
        names(end_cell_ids) <- as.character(leaf_community_ids)

        # if any leaf has been removed, we rerun the end cells / communities identification loop
        
        if (any(unlist(lapply(end_cell_ids, is.null)))){
          changes <- TRUE
        }

        # Update list of excluded communities for next iteration
        excluded_community_ids <- c(
          excluded_community_ids,
          unlist(lapply(which(unlist(lapply(end_cell_ids, is.null))), function(i){
            tail(as.integer(community_paths[[i]]$name), n=1)
          })))
      }

      # Central cells -----------------------------------------------------------------------

      # In each transition community (communities that are neithe root or leaf), we define a "central" cell
      # that will act as a way point in the differentiation trajectories.
      # Central cells are defined as the one having the highest closeness score (inverse of the average length
      # of the shortest paths to/from the other vertices in the graph)

      central_cell_ids = list()

      for (i in setdiff(
        subgraph_community_ids,
        c(leaf_community_ids, excluded_community_ids, starting_community_id))) {

        subgraph_community <- igraph::delete.vertices(
          cell_state_graph$graph,
          v = setdiff(
            seq(num_cells()),
            cell_state_graph$communities_df$cell_ids[[i]]))

        closeness_scores <- closeness(
          subgraph_community,
          vids = V(subgraph_community),
          mode = c("all"),
          weights = NULL,
          normalized = TRUE)

        central_cell_name <- names(which.max(closeness_scores))

        central_cell_id <- which(cell_names() %in% central_cell_name)

        central_cell_ids[[as.character(i)]] <- central_cell_id
      }

      # Aggregate landmark ------------------------------------------------------------------

      if (!is.null(start_cell_id)) {
        cell_state_graph$communities_df$type[starting_community_id]             <<- "Start"
        cell_state_graph$communities_df$landmark_cell_id[starting_community_id] <<- start_cell_id

        for (l in names(end_cell_ids)) {
          cell_state_graph$communities_df$type[as.integer(l)]             <<- "End"
          cell_state_graph$communities_df$landmark_cell_id[as.integer(l)] <<- end_cell_ids[[l]]
        }
        
        if (length(central_cell_ids) > 0) {
          for (c in names(central_cell_ids)) {
            cell_state_graph$communities_df$type[as.integer(c)]             <<- "Central"
            cell_state_graph$communities_df$landmark_cell_id[as.integer(c)] <<- central_cell_ids[[c]]
          }
        }
      }
    }

    cell_state_graph$communities_df$excluded <<- rep(FALSE, nrow(cell_state_graph$communities_df))
    cell_state_graph$communities_df$excluded[excluded_community_ids] <<- TRUE

    # Plot landmarks ----------------------------------------------------------------------

    if (process_plots) {

      vert_cols <- structure(rep("grey", num_cells()), names=cell_names())
      vert_cols[unlist(cell_state_graph$communities_df$cell_ids[excluded_community_ids])] <- "white"

      vert_shapes <- structure(rep("circle", num_cells()), names=cell_names())

      # vertices frame color
      vfc <- color_maps$community_id[pData(expressionSet)$community_id]
      # vfc <- brew_more_colors(seq(cell_state_graph$communities_number), "Set1")[pData(expressionSet)$community_id]
      # vfc = getClusterColors(v=3)[pData(expressionSet)$community_id]
      vfc[na.omit(cell_state_graph$communities_df$landmark_cell_id)] <- "black"

      for (i in which(!cell_state_graph$communities_df$excluded)) {
        vert_cols[cell_state_graph$communities_df$landmark_cell_id[i]] <- color_maps$community_id[i]
      }

      for (i in which(cell_state_graph$communities_df$type == "Start")) {
        vert_cols[cell_state_graph$communities_df$landmark_cell_id[i]] <- "black"
        vert_shapes[cell_state_graph$communities_df$landmark_cell_id[i]] <- "square"
      }

      for (i in which(cell_state_graph$communities_df$type == "End")) {
        vert_shapes[cell_state_graph$communities_df$landmark_cell_id[i]] <- "square"
      }

      cell_state_graph$plot_from_color_vector(
        fullpath           = paste0(output_folder, '/', basename, '_cell_state_graph_landmarks.pdf'),
        sample_colors      = vert_cols,
        vertex.size        = 1,
        vertex.frame.color = vfc,
        edge.width         = 3/5,
        vertex.shape       = vert_shapes,
        vertex.plot.order  = c(setdiff(seq(num_cells()), na.omit(cell_state_graph$communities_df$landmark_cell_id)), na.omit(cell_state_graph$communities_df$landmark_cell_id)))
    }

    # Step 2: Get trajectories for all community edges ####################################

    cat("Get community trajectories\n")

    # The community MST is not storing directed edges. We modify it to capture the 
    # orderinig of the (non-excluded) communities.

    community_mstree_edges <- igraph::get.edgelist(cell_state_graph$communities_mstree)

    # Order edges by community mst ordering, starting from the start communities

    cell_state_graph$communities_df$mst_pt <<- rep(NA, cell_state_graph$communities_number)

    for (subgraph_id in seq(num_subgraphs)) {
    
      subgraph_cell_ids <- which(pData(expressionSet)$subgraph_id == subgraph_id)
      subgraph_community_ids <- sort(unique(pData(expressionSet)[subgraph_cell_ids, "community_id"]))

      starting_community_id <- intersect(
        subgraph_community_ids,
        which(cell_state_graph$communities_df$type == "Start"))

      if (length(starting_community_id) == 1) {

        cell_state_graph$communities_df$mst_pt[subgraph_community_ids] <<- as.numeric(distances(
          cell_state_graph$communities_mstree,
          v = starting_community_id,
          to = subgraph_community_ids,
          weights=E(cell_state_graph$communities_mstree)$weight))
      }
    }

    # Order edges tail to head, by mst ordering
    community_mstree_edges_ordered <- as.data.frame(t(apply(
      community_mstree_edges,
      1,
      function(x){
        d <- cell_state_graph$communities_df$mst_pt[x[1]] - cell_state_graph$communities_df$mst_pt[x[2]]
        if (is.na(d)) {
          c(NA, NA)
        } else {
          if (d < 0) {
            x
          } else {
            rev(x)
          }
        }
      })))
    rownames(community_mstree_edges_ordered) <- apply(community_mstree_edges_ordered[,1:2], 1, paste0, collapse='_')

    # Exclude community MST edges relating excluded communities 
    excluded_edges <- which(is.na(community_mstree_edges_ordered[,2]))
    if (length(excluded_edges) > 0) {
      community_mstree_edges_ordered <- community_mstree_edges_ordered[-excluded_edges, ]
    }
    excluded_edges.2 <- which(
      community_mstree_edges_ordered[, 1] %in% excluded_community_ids |
      community_mstree_edges_ordered[, 2] %in% excluded_community_ids)
    if (length(excluded_edges.2) > 0) {
      community_mstree_edges_ordered <- community_mstree_edges_ordered[-excluded_edges.2, ]
    }

    # Generate ordered community MST 
    community_mstree_ordered = igraph::graph_from_edgelist(
      apply(community_mstree_edges_ordered, 2, as.character),
      directed=TRUE)


    # Step 3: Associate cells to branches #################################################

    cat("Map cells to branches\n")

    # Because some communities stand at the bifurcation point of the differentiation
    # trajectories, some cells belonging to the same community appears to be standing
    # before or after the bifuraction (in one of the downstream paths). To capture this, we
    # define the "branch" structure which store the set of cells lying roughly in-between
    # the graph landmark ponits. A cell branch is community MST edge relating its community
    # to the closest neighbor community with distances defined as the average distance
    # between the cell and all the neighbor community cells. We prefer using this average
    # distance instead of the neighbor community landmark to improve the robustness of the
    # algorithm.

    # Store cell-cell distances from cell state graph
    cell_cell_graph_dist <- igraph::distances(cell_state_graph$graph, weights=NA)

    # Store the neighbor community ids for each community
    community_neighbor_ids <- setNames(
      lapply(
        ego(community_mstree_ordered, order=1, mode="all", mindist=1),
        function(x){
          as.integer(V(community_mstree_ordered)$name[x])
        }),
      V(community_mstree_ordered)$name)

    # Average distance between a cell and its the neighboring community
    cell_comm_avg_dist <- matrix(NA, nrow=num_cells(), ncol=cell_state_graph$communities_number)

    for(cn in seq(num_cells())) {
      comm_id = pData(expressionSet)$community_id[[cn]]
      if (!comm_id %in% excluded_community_ids){
        neighb_comm_ids <- community_neighbor_ids[[as.character(comm_id)]]
        for(nid in neighb_comm_ids){
          cell_comm_avg_dist[cn, nid] <- mean(cell_cell_graph_dist[cn, cell_state_graph$communities_df$cell_ids[[nid]]])
        }
      }
    }

    # For each cell, the branch relates the cell's community to the closest neighbor community
    closest_neigh_comm <- unlist(apply(cell_comm_avg_dist, 1, function(x){if(all(is.na(x))){NA}else{which.min(x)}}))

    cell_branch <- rbind(pData(expressionSet)$community_id, closest_neigh_comm)

    cell_branch_names <- apply(
      cell_branch,
      2,
      function(x){
        if (is.na(x[2])){
          NA
        } else {
          if (cell_state_graph$communities_df$mst_pt[as.integer(x[1])]-cell_state_graph$communities_df$mst_pt[as.integer(x[2])] < 0){
            paste0(x, collapse="_")
          } else {
            paste0(rev(x), collapse="_")
          }
        }
    })

    # Step 4: PT structure ################################################################

    cat("Build PT structure\n")

    # We define the differentation trajectories' pseudotime length as the length of each branch of the community MST. Each branch length is the diameter of the graph composed by the cell of the branch. The diameter is the length of the longest geodesic, ie the length of the longest shortest path between any pair of vertices in the graph.

    community_mstree_edges_ordered$diameter <- unlist(lapply(
      rownames(community_mstree_edges_ordered),
      function(e) {
        igraph::diameter(
          igraph::delete_vertices(cell_state_graph$graph, v = which(cell_branch_names != e)),
          directed    = FALSE,
          unconnected = TRUE,
          weights     = NA)
      }))

    # We construct a "pseudotime graph" where each node is a pseudotime point along a differentation trajectory.

    pt_graph_edges = do.call(
      rbind,
      lapply(
        seq(nrow(community_mstree_edges_ordered)),
        function(i){

          tail      <- community_mstree_edges_ordered[i, 1]

          pt.length <- community_mstree_edges_ordered[i, "diameter"]
         
          df = data.frame(
            "from" = paste0('e', i, '_v', seq(pt.length-1)),
            "to"   = paste0('e', i, '_v', seq(2, pt.length)))
         
          # Connect branch to its parent node
          ancestor_edge = which(community_mstree_edges_ordered[,2] == tail)

          # Only the root edge(s) has no ancestor 
          if (length(ancestor_edge) == 1) {
            df <- rbind(df, 
              data.frame(
                "from" = paste0('e', ancestor_edge, '_v', community_mstree_edges_ordered[ancestor_edge, "diameter"]),
                "to"   = paste0('e', i, '_v1'))
              )
          }

          if (length(ancestor_edge) == 0) {
            df <- rbind(
              df,
              data.frame(
                "from" = paste0("root_", cell_state_graph$communities_df$subgraph_id[[tail]]),
                "to"   = paste0('e', i, '_v1')))
          }
          return(df)
        }))

    pt_graph_vertices <- do.call(
      rbind,
      lapply(
        seq(nrow(community_mstree_edges_ordered)),
        function(i){

          tail      <- community_mstree_edges_ordered[i, 1]
          head      <- community_mstree_edges_ordered[i, 2]
          pt.length <- community_mstree_edges_ordered[i, "diameter"]
          lids      <- c(rep(NA, pt.length-1), head)
     
          data.frame(
              "name"         = paste0('e', i, '_v', seq(pt.length)),
              "branch_id"    = i,
              "branch_color" = colorRampPalette(c(color_maps$community_id[tail], color_maps$community_id[head]))(n = (pt.length+1))[2:(pt.length+1)],
              "landmark_id"  = lids,
              "subgraph_id"  = cell_state_graph$communities_df$subgraph_id[[tail]]
              )
        }))

    pt_graph_vertices <- rbind(
      pt_graph_vertices,
      do.call(
        rbind,
        lapply(
          which(cell_state_graph$communities_df$type == "Start"),
          function(i) {
            data.frame(
              "name"         = paste0("root_", cell_state_graph$communities_df$subgraph_id[i]),
              "branch_id"    = 0,
              "branch_color" = color_maps$community_id[i],
              "landmark_id"  = i,
              "subgraph_id"  = cell_state_graph$communities_df$subgraph_id[[i]]
              )
          })))

    PT_graph <<- igraph::graph_from_data_frame(
      pt_graph_edges,
      directed = TRUE,
      vertices = pt_graph_vertices)

    # Store pt_graph vertex PT

    V(PT_graph)$pt <<- do.call(
      cbind,
      lapply(
        unique(pt_graph_vertices$subgraph_id),
        function(i) {
          vert_ids <- which(pt_graph_vertices$subgraph_id == i)
          igraph::distances(
            PT_graph,
            v       = intersect(vert_ids, which(pt_graph_vertices$branch_id == 0)),
            to      = vert_ids,
            weights = NA)}))[,V(PT_graph)$name]

    # We map cells to the nodes of the PT graph

    # Store communities' upstream communities
    comm_upstream_list <- setNames(
      lapply(
        ego(
          community_mstree_ordered,
          order   = vcount(community_mstree_ordered),
          mode    = "in",
          mindist = 1),
        function(x) {
          as.integer(V(community_mstree_ordered)$name[x])
        }),
      V(community_mstree_ordered)$name)

    # Store communities' downstream communities
    comm_downstream_list <- setNames(
      lapply(
        ego(
          community_mstree_ordered,
          order   = vcount(community_mstree_ordered),
          mode    = "out",
          mindist = 1),
        function(x) {
          as.integer(V(community_mstree_ordered)$name[x])}),
      V(community_mstree_ordered)$name) 

    # For each cell, store mean distance with all upstream-community cells and all downstream-community cells on the cell-cell graph
    cell_flow_avg_dist = matrix(
      NA,
      nrow = num_cells(),
      ncol = 2,
      dimnames = list(cell_names(), c('upstream', 'downstream')))

    for (cn in seq(num_cells())) {
      
      comm_id           <- pData(expressionSet)$community_id[[cn]]
      upstream_comm_ids <- comm_upstream_list[[as.character(comm_id)]]
      
      if (length(upstream_comm_ids) > 0) {
        cell_flow_avg_dist[cn, 1] <-  mean(
          cell_cell_graph_dist[cn, unlist(cell_state_graph$communities_df$cell_ids[upstream_comm_ids])])
      }

      downstream_comm_ids <- comm_downstream_list[[as.character(comm_id)]]

      if (length(downstream_comm_ids) > 0){
        cell_flow_avg_dist[cn, 2] <-  mean(
          cell_cell_graph_dist[cn, unlist(cell_state_graph$communities_df$cell_ids[downstream_comm_ids])])
      }
    }

    cell_flow_avg_dist[is.na(cell_flow_avg_dist)] <- 0

    # Map using the rank of ... in each branch
    for(i in seq(nrow(community_mstree_edges_ordered))) {

      pt.length <- community_mstree_edges_ordered[i, "diameter"]
      
      edge_cell_ids <- which(cell_branch_names == rownames(community_mstree_edges_ordered)[i])
      
      # cells from excluded communities do not have PT coordinate 
      if (length(edge_cell_ids) > 0) {
        cellrank <- rank(cell_flow_avg_dist[edge_cell_ids, "upstream"]) -
          rank(cell_flow_avg_dist[edge_cell_ids, "downstream"])
        # cellrank = - rank(cell_flow_avg_dist[edge_cell_ids, "downstream"])
        # cellrank = rank(cell_flow_avg_dist[edge_cell_ids, "upstream"])
        cellrank_cut <- cut(cellrank, pt.length)

        cellnames_list <- lapply(
          levels(cellrank_cut),
          function(l){
            cell_names()[edge_cell_ids][which(cellrank_cut==l)]
          })

        V(PT_graph)[which(V(PT_graph)$branch_id==i)]$cellnames <<- cellnames_list
      }
    }

    excluded_cell_names = cell_names()[which(pData(expressionSet)$community_id %in% excluded_community_ids)]

    # Transfer pseudotime and branch ids to the cells' expressionSet
    pData(expressionSet)$Pseudotime <<- c(
      setNames(
        rep(V(PT_graph)$pt, lapply(V(PT_graph)$cellnames, length)),
        unlist(V(PT_graph)$cellnames)),
      setNames(
        rep(NA, length(excluded_cell_names)),
        excluded_cell_names))[cell_names()]

    pData(expressionSet)$Branch_id <<- c(
      setNames(
        rep(V(PT_graph)$branch_id, lapply(V(PT_graph)$cellnames, length)),
        unlist(V(PT_graph)$cellnames)),
      setNames(
        rep(NA, length(excluded_cell_names)),
        excluded_cell_names))[cell_names()]

    # Plot pt graph

    if (process_plots) {

      leaves_edges <- unlist(lapply(
        which(!cell_state_graph$communities_df$excluded &
              cell_state_graph$communities_df$type == "End"),
        function(i) {
          which(community_mstree_edges_ordered[,2]==i)
        }))
      
      leaves_names = paste0("e", leaves_edges, "_v", community_mstree_edges_ordered[leaves_edges, "diameter"])

      temp_csg = CellStateGraph$new()
      temp_csg$graph = PT_graph
      temp_csg$graph <- igraph::delete_vertex_attr(temp_csg$graph, "cellnames")
      temp_csg$project(method="gephiForceAtlas2", num_iter=100000, seed=0, stop_mean=.5)
      
      pdf(paste0(output_folder, '/', basename, '_forceatlas_communities.pdf'))
      plot(
        PT_graph,
        layout             = temp_csg$layout,
        vertex.label       = V(PT_graph)$landmark_id,
        vertex.label.color = "black",
        vertex.size        = 8,
        edge.arrow.size    = 0,
        vertex.frame.color = "white",
        vertex.color       = V(PT_graph)$branch_color)
      graphics.off()

      leaf_communities <- which(cell_state_graph$communities_df$type == "End")
      pdf(
        paste0(output_folder, '/', basename, '_communities.pdf'),
        width  = 3 + length(leaf_communities),
        height = 0 + .5 * max(V(PT_graph)$pt))
      plot(
        PT_graph,
        layout             = layout_as_tree(
          PT_graph,
          root = which(pt_graph_vertices$branch_id == 0)),
        vertex.label       = V(PT_graph)$landmark_id,
        vertex.label.color = "black",
        vertex.label.dist  = -0.08,
        vertex.size        = 10,
        vertex.size2       = 10,
        edge.arrow.mode    = 0,
        edge.arrow.size    = 0,
        vertex.frame.color = "white",
        vertex.color       = V(PT_graph)$branch_color,
        asp                = 0)
      graphics.off()

      pdf(paste0(output_folder, '/', basename, '_forceatlas_pt.pdf'))
      plot(
        PT_graph,
        layout             = temp_csg$layout,
        vertex.label       = V(PT_graph)$landmark_id,
        vertex.label.color = "white",
        vertex.size        = 8,
        edge.arrow.size    = 0,
        vertex.frame.color = "white",
        vertex.color       = colorRampPalette(c('gray90','black'))(100)[as.numeric(cut(V(PT_graph)$pt, breaks=100))])
      graphics.off()

      leaf_communities <- which(cell_state_graph$communities_df$type == "End")
      pdf(
        paste0(output_folder, '/', basename, '_pt.pdf'),
        width  = 3 + length(leaf_communities),
        height = 0 + .5 * max(V(PT_graph)$pt))
      plot(
        PT_graph,
        layout             = layout_as_tree(
          PT_graph,
          root = which(pt_graph_vertices$branch_id == 0)),
        vertex.label       = V(PT_graph)$landmark_id,
        vertex.label.color = "white",
        vertex.label.dist  = -0.08,
        vertex.size        = 10,
        vertex.size2       = 10,
        edge.arrow.mode    = 0,
        edge.arrow.size    = 0,
        vertex.frame.color = "gray10",
        vertex.color       = colorRampPalette(c('gray90','black'))(100)[as.numeric(cut(V(PT_graph)$pt, breaks=100))],
        asp                = 0)
      graphics.off()

      rm(temp_csg)

      # Plot branch <-> cell mapping        
      cell_state_graph$plot_from_color_vector(
        fullpath           = paste0(output_folder, '/', basename, '_cell_state_graph_pt.pdf'),
        sample_colors      = colorRampPalette(c('gray90','black'))(max(V(PT_graph)$pt))[pData(expressionSet)$Pseudotime],
        vertex.size        = 3,
        vertex.frame.color = brew_more_colors(seq(length(unique(pData(expressionSet)$Branch_id))), "Dark2")[factor(pData(expressionSet)$Branch_id)],
        edge.width         = 3/5,
        vertex.shape       = vert_shapes,
        vertex.plot.order  = c(
          setdiff(
            seq(num_cells()),
            na.omit(cell_state_graph$communities_df$landmark_cell_id)),
          na.omit(cell_state_graph$communities_df$landmark_cell_id)))

      # Plot branch and landmarks
      cell_state_graph$plot_from_color_vector(
        fullpath           = paste0(output_folder, '/', basename, '_cell_state_graph_branches.pdf'),
        sample_colors      = vert_cols,
        vertex.size        = 3,
        vertex.frame.color = brew_more_colors(seq(length(unique(pData(expressionSet)$Branch_id))), "Dark2")[factor(pData(expressionSet)$Branch_id)],
        # vertex.frame.color = brew_more_colors(seq(length(unique(cell_branch_names))), "Dark2")[factor(cell_branch_names)],
        edge.width         = 3/5,
        vertex.shape       = vert_shapes,
        vertex.plot.order  = c(
          setdiff(
            seq(num_cells()),
            na.omit(cell_state_graph$communities_df$landmark_cell_id)),
          na.omit(cell_state_graph$communities_df$landmark_cell_id)))
    }
  }
)

#' Visualize the smoothed gene expression dynamics along specific branches of the pseudotime tree.
#'
#' @name plot_pseudotime_dynamics
#' @param gene_list character vector. The names of the genes to render.
#' @param leaf_ids integer vector. Specifies the indices of the differentiation trajectories to visualize. If NULL (default), all the trajectories are rendered.
#' @param title character string. The title of the plot. Default to the rendered gene names. 
#' @param subtitle character string. The subtitle of the plot. Default to NULL.
#' @param data_status character string. Specifies whether the gene expression levels to be rendered are raw ("Raw"), normalized ("Normalized") or have been imputed ("Smoothed"). Default to "Smoothed".
#' @param show_cells logical. Whether to show the gene expression levels for each individual cells. Default to FALSE.
#' @param basename character string. The name of the output file, without extension.
Antler$methods(

  plot_pseudotime_dynamics = function(
    gene_list,
    leaf_ids    = NULL,
    title       = paste0(gene_list, collapse="-"),
    subtitle    = NULL,
    data_status = "Smoothed",
    show_cells  = FALSE,
    basename    = paste0(c("PT_dynamics", gene_list, data_status), collapse="_")) {

    # Calculate gene expression along each trajectories
    cat("Generate smoothed gene expression levels along each trajectories\n")
    ptstruct_res <- getPTstructures(
      PT_graph,
      data = read_counts(data_status=data_status)[gene_list, , drop = F],
      csg  = cell_state_graph)

    if (identical(NULL, leaf_ids)) {
      endpoints <- unique(ptstruct_res$PT$path_id)
    } else {
      endpoints <- leaf_ids
      # endpoints <- sort(unique(ptstruct_res$PT$path_id))[leaf_ids]
    }

    pt_max <- max(ptstruct_res$PT$pt)

    data_subset <- ptstruct_res$PT[which(
      ptstruct_res$PT$gene_name %in% gene_list &
      ptstruct_res$PT$path_id %in% endpoints
      # !grepl("root", ptstruct_res$PT$name)
      ), ]

    # use the averaged smoothed trajectories
    mean_type <- "mean.smooth.avg" # "mean.smooth" for non-averaged
    sd_type <- "sd.smooth.avg"     # "sd.smooth" for non-averaged
    
    # rename mean / sd for ggplot call
    colnames(data_subset)[which(colnames(data_subset) == mean_type)] <- "cur_mean"
    colnames(data_subset)[which(colnames(data_subset) == sd_type)]   <- "cur_sd"

    # If the trajectories are averaged we plot only one instance of each vertex
    data_subset_ribbon <- data_subset[, c('name', 'pt', 'cur_mean', 'cur_sd', 'gene_name', 'branch_id')]
    
    if (mean_type == "mean.smooth.avg" & sd_type == "sd.smooth.avg") {

      # each branch is missing the connection between its first vertex and its 
      # ancestor (landmark).
      for (branch_id in setdiff(unique(data_subset$branch_id), 0)) {

        data_branch <- data_subset[which(data_subset$branch_id == branch_id), , drop = FALSE]
        first_vert <- data_branch$name[which.min(data_branch$pt)]
        parent_vert <- igraph::ego(
          PT_graph,
          order = 1,
          nodes = which(V(PT_graph)$name == first_vert),
          mode  = "in")[[1]][2]
        new_rows <- data_subset_ribbon[which(data_subset_ribbon$name == parent_vert$name), ]
        new_rows <- new_rows[!duplicated(new_rows$gene_name), ]
        new_rows$branch_id <- rep(branch_id, nrow(new_rows))

        data_subset_ribbon <- rbind(data_subset_ribbon, new_rows)
      } 

      nonredondant_ids <- which(
        !duplicated(data_subset_ribbon[, c('gene_name', 'branch_id', 'pt')]) &
        data_subset_ribbon$branch_id != 0
        )
    
      data_subset_ribbon <- data_subset_ribbon[nonredondant_ids, ]
    }

    # Build cell level dataset if needed
    if (show_cells) {
      vname_color <- data_subset[, c("name", "branch_color")]
      vname_color <- vname_color[!duplicated(vname_color), ]
      cell_df <- do.call(
        rbind.data.frame,
        lapply(
          seq(nrow(vname_color)),
          function(i) {
            vname <- as.character(vname_color$name[i])
            cellnames <- V(ptstruct_res$PT_graph)[[vname]]$cellnames
            pt <- V(ptstruct_res$PT_graph)[[vname]]$pt
            branch_color <- V(ptstruct_res$PT_graph)[[vname]]$branch_color
            if (length(cellnames) > 0) {
              do.call(
                rbind.data.frame,
                lapply(
                  gene_list,
                  function(gn) {
                    data.frame(
                        "level"        = read_counts(data_status = data_status)[
                          gn,
                          cellnames],
                        "gene_name"    = gn,
                        "PT"           = pt,
                        "branch_color" = branch_color)
              }))
             } else {
              data.frame(
                "level"        = numeric(),
                "gene_name"    = character(),
                "PT"           = integer(),
                "branch_color" = character())
            }
        }))

    }

    # Rescale gene level to max
    cat("Normalize gene expression levels\n")
    for(gn in gene_list){

      gene_row_ids <- which(data_subset$gene_name == gn)

      level_max <- max(data_subset[gene_row_ids, "cur_mean"], na.rm=TRUE)

      data_subset[gene_row_ids, "cur_mean"] <- data_subset[gene_row_ids, "cur_mean"] / level_max
      data_subset[gene_row_ids, "cur_sd"]   <- data_subset[gene_row_ids, "cur_sd"] / level_max

      gene_row_ids <- which(data_subset_ribbon$gene_name == gn)          

      data_subset_ribbon[gene_row_ids, "cur_mean"] <- data_subset_ribbon[gene_row_ids, "cur_mean"] / level_max
      data_subset_ribbon[gene_row_ids, "cur_sd"]   <- data_subset_ribbon[gene_row_ids, "cur_sd"] / level_max

      if (show_cells) {
        cell_df[cell_df$gene_name == gn, "level"] <- cell_df[cell_df$gene_name == gn, "level"] / level_max
      }
    }

    data_subset_bifurcation = data_subset[which(!is.na(data_subset$landmark_id)),]

    bifcols <- unique(as.character(data_subset_bifurcation$branch_color))
    names(bifcols) <- bifcols

    genecols <- RColorBrewer::brewer.pal(9, "Set1")[seq(length(gene_list))]
    names(genecols) <- gene_list

    allcols <- c(bifcols, genecols)
    max_level <- 1

    if (show_cells) {
      allcols <- sort(c(
        allcols,
        setNames(
          unique(as.character(cell_df$branch_color)),
          unique(cell_df$branch_color))))
      max_level <- max(max_level, max(cell_df$level))
    }

    pdf(paste0(output_folder, '/', basename, '.pdf'), width=7, height=4)
    p <- ggplot() +
      geom_ribbon(
        data = data_subset_ribbon,
        aes(
          x = as.numeric(pt),
          ymax = as.numeric(cur_mean) + as.numeric(cur_sd),
          ymin = as.numeric(cur_mean) - as.numeric(cur_sd),
          fill = gene_name,
          group = interaction(gene_name, branch_id)
          ),
        alpha = .1,
        color = NA,
        show.legend=TRUE) +
      {
      if (show_cells)
        geom_point(
            data = cell_df,
            aes(
              PT, level,
              col = branch_color),
              # fill = branch_color),
            # pch = 21,
            size = .8)
      } +
      geom_line(
        data = data_subset,
        aes(
          x     = as.numeric(pt),
          y     = as.numeric(cur_mean),
          color = gene_name,
          group = interaction(gene_name, path_id)),
        size        = .5,
        show.legend = FALSE) +
      scale_x_continuous(
        breaks = seq(0, pt_max, length.out = 5),
        labels = seq(0, 1, length.out=5)) +
      xlab("Pseudotime") +
      ylab(paste0(
        "Expression (", 
        data_status,
        ")")) +
      ggtitle(title, subtitle) + 
      geom_point(
        data = data_subset_bifurcation,
        aes(
          x     = as.numeric(pt),
          y     = as.numeric(cur_mean),
          color = branch_color),
        size        = 3,
        show.legend = FALSE) +
      scale_color_manual(values = allcols, guide = FALSE) +
      scale_fill_manual(values = allcols, breaks = names(genecols), name = "") +
      theme(
        panel.border      = element_blank(), 
        axis.ticks.length = unit(.1, "inch"),
        panel.grid.major  = element_blank(),
        axis.line.x       = element_line(
          colour = "black", size = .5, linetype = "solid"),
        axis.line.y       = element_line(
          colour = "black", size = .5, linetype = "solid")) +
      coord_cartesian(ylim = c(0, max_level), xlim = c(0, pt_max)) +
      theme_classic()
    print(p)
    graphics.off()
  }
)

# KEEPING BUT NOT DOCUMENTING (HIDING ?) +++++++++++++++++++++++++++++++++++++

Antler$methods(
  transfer_carver_weights = function(
    gene_modules_name,
    from = "best",
    seed = 0) {

    fData(expressionSet)$carver_weights <<- rep(NA, nrow(fData(expressionSet)))

    if (identical(from, "best")) {
      fData(
        expressionSet)[unlist(gene_modules$get(gene_modules_name)), "carver_weights"] <<- carver$W_optim[rep(
          seq(length(gene_modules$get(gene_modules_name))),
          unlist(lapply(gene_modules$get(gene_modules_name), length)))]
    } else {
      gen_sols <- which(carver$runVariables$GA_pareto_parameters[, "gen"] == from$generation)

      if (is.integer(from$which)) {
        gen_sol <- gen_sols[[from$which]]
      } else if (identical(from$which, "random")) {
        set.seed(0)
        gen_sol <- sample(gen_sols, 1)
      } else if (identical(from$which, "top_timecorr")) {
        top_id <- which.min(carver$runVariables$GA_pareto_scores[gen_sols, "timecorr"])
        gen_sol <- gen_sols[top_id]
      } else if (identical(from$which, "top_smoothness")) {
        top_id <- which.min(carver$runVariables$GA_pareto_scores[gen_sols, "smoothness"])
        gen_sol <- gen_sols[top_id]
      } else if (is.numeric(from$which)) {
        gen_sol <- gen_sols[as.integer(from$which)]
      }

      weights <- carver$runVariables$GA_pareto_parameters[
        gen_sol,
        -ncol(carver$runVariables$GA_pareto_parameters)]

      fData(
        expressionSet)[unlist(gene_modules$get(gene_modules_name)), "carver_weights"] <<- weights[rep(
          seq(length(gene_modules$get(gene_modules_name))),
          unlist(lapply(gene_modules$get(gene_modules_name), length)))]
    }
  },

  carving_dist_to_cell_dist = function(
    allow_multiple_roots,
    canalize         = TRUE,
    kernel_sigma     = .1,
    seed             = 0,
    plot_folder_path = NA,
    suffix           = ""
    ) {

    pos_carve_w_gene_ids <- which(!is.na(fData(expressionSet)[, "carver_weights"]))

    cell_state_graph$dist_matrix <<- get_weighted_distance_matrix(
      read_counts(data_status = "Normalized")[
        pos_carve_w_gene_ids, ],
      fData(expressionSet)[
        pos_carve_w_gene_ids, "carver_weights"],
      1,
      data_transform='log_scale')

    if (canalize) {

      ## retrieve top cluster mst
      gm_carver_weights <- carver$W_optim

      gm_carver_weights <- gm_carver_weights / sum(gm_carver_weights)

      cluster_cluster_distance <- newClusterWeightedDistance(
        gm_carver_weights,
        carver$runVariables$compact_data)

      graph <- igraph::graph_from_adjacency_matrix(
        as.matrix(cluster_cluster_distance),
        weighted = TRUE,
        mode     = "upper")

      mst <- igraph::minimum.spanning.tree(
        graph,
        algorithm = 'prim')

      # cluster distance on mst
      mst_distances <- igraph::distances(mst)

      mst_proba <- exp(-mst_distances / kernel_sigma ** 2)

      # expcell: clust to cell
      mst_cell_prob <- matrix(
        0,
        ncol     = num_cells(),
        nrow     = num_cells(),
        dimnames = list(cell_names(), cell_names()))

      for (i in seq(ncol(carver$runVariables$compact_data))) {
        for (j in seq(ncol(carver$runVariables$compact_data))) {
          mst_cell_prob[
            carver$runVariables$cell_ids_from_cluster_id[[i]],
            carver$runVariables$cell_ids_from_cluster_id[[j]]
            ] <- mst_proba[i, j]
        }
      }

      cell_state_graph$dist_matrix <<- cell_state_graph$dist_matrix / mst_cell_prob

      cell_state_graph$dist_matrix <<- cell_state_graph$dist_matrix / sum(cell_state_graph$dist_matrix)

      if (!identical(plot_folder_path, NA)) {

        pdf(file.path(
          plot_folder_path, 
          paste0("canal_cluster_cluster_distance", suffix, ".pdf")))
        par(mar = c(.5, .5, 2, .5))
        image(
          cluster_cluster_distance,
          axes = FALSE,
          asp = 1,
          col = colorRampPalette(c("#0464DF", "#FFE800"))(n = 101),
          main = "Carver cluster-cluster distance")
        graphics.off()

        pdf(file.path(
          plot_folder_path,
          paste0("canal_top_MST_distance", suffix, ".pdf")))
        par(mar = c(.5, .5, 2, .5))
        image(
          mst_distances,
          axes = FALSE,
          asp = 1,
          col = colorRampPalette(c("#0464DF", "#FFE800"))(n = 101),
          main = "Carver top MST distance")
        graphics.off()

        pdf(file.path(
          plot_folder_path,
          paste0("canal_top_MST_transition", suffix, ".pdf")))
        par(mar = c(.5, .5, 2, .5))
        image(
          mst_proba,
          axes = FALSE,
          asp = 1,
          col = colorRampPalette(c("#0464DF", "#FFE800"))(n = 101),
          main = '"Carver top MST transition "probability"')
        graphics.off()

      }
    }
  }
)

# MAYBES +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Antler$methods(

  excludeCellPopulations = function(cluster_markers, used_gm, basename=NULL, data_status='Check'){

    cell_to_exclude_ids = identifyCellPopulations(cluster_markers=cluster_markers, used_gm=used_gm, basename=basename, data_status=data_status)

    remove_cells(ids =unlist(cell_to_exclude_ids))

    exclude_unexpressed_genes(min_cells=3, data_status=data_status)

    processing_state <<- c(processing_state, paste0('Exclude genes expressed in less than 3 cells (dataset dimension:', nrow(exprs(expressionSet)), ' ', ncol(exprs(expressionSet)), ')'))

    if (write_report) write_report_file()

    },

  removeGenesFromRatio = function(
    threshold       = 0.1,
    candidate_genes = NA,
    verbose         = TRUE,
    plot_ratios     = FALSE,
    file_path       = NULL,
    data_status     = 'Raw') {

    if (identical(candidate_genes, NA)){
      stop("candidate genes are missing.")
    }

    gene_ratio <- colSums(read_counts(data_status=data_status)[candidate_genes,]) / colSums(read_counts(data_status=data_status))
    
    if (plot_ratios | !is.null(file_path)){
      if (!is.null(file_path)){
        pdf(file_path)
      }
      plot(colSums(read_counts(data_status=data_status)[candidate_genes, ]) / colSums(read_counts(data_status=data_status)), ylab="Candidate Gene Counts Ratio")
      if (!is.null(file_path)){
        graphics.off()
      }
    }

    high_cells <- which(gene_ratio > threshold)
    
    remove_cells(ids = high_cells)

    remove_genes(ids = which(fData(expressionSet)$current_gene_names %in% candidate_genes))

    log <- paste0('Removed ', length(high_cells),' cells containing more than ', round(threshold * 100), '% of candidate genes reads and then excluded all mitochondrial read counts (dataset dimension:', nrow(exprs(expressionSet)), ' ', ncol(exprs(expressionSet)), ', candidate genes: ', paste0(candidate_genes, collapse=', '), ')')

    if (verbose){cat(paste0(log, '\n'))}
    processing_state <<- c(processing_state, log)
    if (write_report) write_report_file()

    },

  # Exclude gene with low level in expressing cells
  removeLowlyExpressedGenes = function(expression_threshold = 1,
    selection_theshold   = 10,
    verbose              = TRUE,
    data_status          = 'Check'){

    # multiply gene level by binary mask to get only "positive" cells
    pos_cells_counts = read_counts(data_status=data_status) * (read_counts(data_status=data_status) > expression_threshold)
    pos_cells_counts[which(pos_cells_counts == 0)] <- NA
    pos_cells_mean_level = apply(pos_cells_counts, 1, mean, na.rm=T)
    
    excluded_genes_ids = which(pos_cells_mean_level < selection_theshold)

    remove_genes(ids = excluded_genes_ids)

    log = paste0('Removed ', length(excluded_genes_ids),' genes with an average expression level lower than ', selection_theshold, ' in "positive" cells, i.e. expression threshold higher than ', expression_threshold, ' (performed on ', data_status, ' dataset, dimension:', nrow(exprs(expressionSet)), ' ', ncol(exprs(expressionSet)), ')')

    if (verbose){cat(paste0(log, '\n'))}
    processing_state <<- c(processing_state, log)
    if (write_report) write_report_file()

    },

  identifyCellPopulations = function(cluster_markers, used_gm, basename=NULL, data_status='Check'){

    if (! used_gm %in% getGeneModuleList())
      stop("invalid gene module(s) (must be one of the following: '", paste0(getGeneModuleList(), collapse="', '"), "')")

    gm_split = parseGMname(used_gm)
    gms = .self[[gm_split$dr_name]][[gm_split$gm_name]]

    if (length(gms) == 0)
      stop("the gene module list is empty. Run gene modules identification first.")
    
    # Remove pop / Trim markers if marker absent from gene modules
    cluster_markers <- unlist(lapply(seq(length(cluster_markers)), function(i){    
        l = cluster_markers[[i]]
        l_out = l[!l %in% unlist(gms)]
        if (length(l_out)){
          print(paste0(names(cluster_markers)[[i]], " is missing: ", paste0(l_out, collapse="', ")))
        }
        l_in=list(l[l %in% unlist(gms)])
        names(l_in) = names(cluster_markers)[[i]]
        if (length(l_in[[1]])==0) NULL else l_in
        }), recursive=F)

    # Number of clusters is the number of excluded populations plus the remaining cells cluster
    num_clusters = length(cluster_markers) + 1

    # For each excluded population, gather the modules containing specified marker genes
    selected.genemodules = lapply(cluster_markers, function(pm){
              unlist(Filter(function(i){length(intersect(i, pm)) > 0}, gms))
      })

    # For safety purpose, assure that proper gene names are used
    selected.genemodules <- lapply(selected.genemodules, function(l){intersect(l, fData(expressionSet)$current_gene_names)})

    data.logscaled = t(scale(t(log(read_counts(data_status=data_status)[unlist(selected.genemodules),]+1)), center=T, scale=T))

    mod_avg = do.call(rbind ,
                      lapply(selected.genemodules, function(m){
                         colMeans(data.logscaled[m, , drop=F])
                       }))

    # distance over all genes
    # cell_cell_distance = as.dist(fastEuclideanDist(t(data.logscaled)))

    # distance over gm average <- better if tiny gm describing cell population
    cell_cell_distance = as.dist(fastEuclideanDist(t(mod_avg)))

    hc.cells = hclust(cell_cell_distance, method = "ward.D2")

    hc.cells.clust_ids = cutree(hc.cells, k=num_clusters)

    excluded_clust_id = list()
    excluded_clust_count = list()
    for(i in 1:num_clusters){
        clust.mean = rowMeans(mod_avg[, which(hc.cells.clust_ids == i), drop=F])
        exclude_clust=Reduce('|', clust.mean  > 0)
        if (exclude_clust){
          clust_name=names(which.max(clust.mean))
          excluded_clust_id[clust_name] = i
          excluded_clust_count[clust_name] = length(which(hc.cells.clust_ids %in% i))
        }
    }

    if (!is.null(basename)){ 

        clust.colors = c(RColorBrewer::brewer.pal(9, "Set1"), RColorBrewer::brewer.pal(12, "Set3"), RColorBrewer::brewer.pal(8, "Set2"))

        hm = plot_transcriptome_summary_heatmap(mat=data.logscaled, genemodules=selected.genemodules, filename=paste0(output_folder, '/', basename, '_heatmap.pdf'), cells_colors=cbind(clust.colors[hc.cells.clust_ids], pData(expressionSet)$cells_colors), displayed.geneset=rownames(data.logscaled), zscore.thres=2, legend=list('text'=names(excluded_clust_id), 'colors'=clust.colors[unlist(excluded_clust_id)]), dendrogram=as.dendrogram(hc.cells))
    }

    processing_state <<- c(processing_state, paste0('Identify cells belong to the following clusters: ', paste0(paste0(names(cluster_markers), ' (', unlist(excluded_clust_count[names(cluster_markers)]), ' cells)'), collapse=', '), ' (dataset dimension:', nrow(exprs(expressionSet)), ' ', ncol(exprs(expressionSet)), ')'))
    
    if (write_report) write_report_file()
    
    pop_ids = lapply(excluded_clust_id, function(i){
            which(hc.cells.clust_ids %in% i)
      })

    # return(pop_ids[names(cluster_markers)])
    return(pop_ids)

    },

  identifyCellPopulationsFromExistingCellClusters = function(cluster_markers, used_gm, basename=NULL, data_status='Check', existing_clusters="hclust", ratio_thres=.5){

    if (! used_gm %in% getGeneModuleList())
      stop("invalid gene module(s) (must be one of the following: '", paste0(getGeneModuleList(), collapse="', '"), "')")

    gm_split = parseGMname(used_gm)

    if (length(.self[[gm_split$dr_name]][[gm_split$gm_name]]) == 0)
      stop("the gene module list is empty. Run gene modules identification first.")

    #

    # Number of clusters is the number of excluded populations plus the remaining cells cluster
    # num_clusters = length(cluster_markers) + 1

    # For each excluded population, gather the modules containing specified marker genes
    selected.genemodules = lapply(cluster_markers, function(pm){
              unlist(Filter(function(i){length(intersect(i, pm)) > 0}, .self[[gm_split$dr_name]][[gm_split$gm_name]]))
      })

    # For safety purpose, assure that proper gene names are used
    selected.genemodules <- lapply(selected.genemodules, function(l){intersect(l, fData(expressionSet)$current_gene_names)})

    data.logscaled = t(scale(t(log(read_counts(data_status=data_status)[unlist(selected.genemodules),]+1)), center=T, scale=T))

    mod_avg = do.call(rbind ,
                      lapply(selected.genemodules, function(m){
                         colMeans(data.logscaled[m, , drop=F])
                       }))

    mod_avg.bin = silent.binarize.array(mod_avg)

    # for each population GM, select cell cluster with ratio of positive cells higher than threshold
    excluded_clust_id = lapply(seq(nrow(mod_avg)), function(i){
      Filter(
        function(x){
          clust_cell_ids = which(cellClusters[[existing_clusters]]$cell_ids == x)
          sum(mod_avg.bin[i, clust_cell_ids]) / length(clust_cell_ids) > ratio_thres
        },
        seq(cellClusters[[existing_clusters]]$num_clusters))
      })
    names(excluded_clust_id) <- names(cluster_markers)
    
    excluded_clust_count = lapply(excluded_clust_id, function(l){
          length(which(cellClusters[[existing_clusters]]$cell_ids %in% l))
      })

    processing_state <<- c(processing_state, paste0('Identify cells belong to the following clusters: ', paste0(paste0(names(cluster_markers), ' (', unlist(excluded_clust_count[names(cluster_markers)]), ' cells)'), collapse=', '), ' (dataset dimension:', nrow(exprs(expressionSet)), ' ', ncol(exprs(expressionSet)), ')'))
    
    if (write_report) write_report_file()
    
    return(excluded_clust_id)

    },

  # Monocle 2 (Install first with http://cole-trapnell-lab.github.io/monocle-release/getting-started/)
  runMonocle2DETestOnClusters = function(
    clustering_method    = 'hclust',
    ref_clusters         = c("1"),
    alt_clusters         = unique(cellClusters[[clustering_method]]$cell_ids),
    dispersion_selection = TRUE) {

    closeAllConnections()

    df_pheno <- data.frame(
      row.names       = colnames(read_counts(data_status='Normalized')),
      "cluster_id"    = cellClusters[[clustering_method]]$cell_ids,
      "test_clusters" = 1 * (cellClusters[[clustering_method]]$cell_ids %in% ref_clusters)
    )

    used_genes <- if (dispersion_selection) {getDispersedGenes(read_counts(data_status = 'Normalized'), zscore_threshold = 0.0)} else {nrow(num_genes())}

    df_feature <- data.frame(row.names = used_genes, "gene_short_name" = used_genes)

    HSMM <- monocle::newCellDataSet(
      as.matrix(read_counts(data_status = 'Normalized')[used_genes, ]),
      phenoData        = new("AnnotatedDataFrame", data = df_pheno),
      featureData      = new('AnnotatedDataFrame', data = df_feature),
      expressionFamily = VGAM::tobit()) # "tobit() is appropriate for FPKM, TPM gene expression data."

    diff_test_res <- monocle::differentialGeneTest(
      HSMM[, which(df_pheno$cluster_id %in% unique(c(ref_clusters, alt_clusters)))],
      fullModelFormulaStr    = "~test_clusters",
      reducedModelFormulaStr = "~1",
      cores                  = num_cores)

    return(diff_test_res)
  },

  # Monocle 2 (Install first with http://cole-trapnell-lab.github.io/monocle-release/getting-started/)
  runMonocle2DETest = function(
    fullModelFormulaStr,
    reducedModelFormulaStr = "~1",
    expressionFamily=VGAM::negbinomial.size(),
    tested_genes=gene_names()) {

    closeAllConnections()

    df_pheno = pData(expressionSet)
    df_feature = cbind(fData(expressionSet), 'gene_short_name'=fData(expressionSet)$current_gene_names)
    rownames(df_feature) <- df_feature$gene_short_name

    HSMM <- monocle::newCellDataSet(
          as.matrix(read_counts(data_status='Raw')),
          phenoData = new("AnnotatedDataFrame", data = df_pheno),
          featureData = new('AnnotatedDataFrame', data = df_feature),
          expressionFamily=expressionFamily
          )

    HSMM <- estimateSizeFactors(HSMM)

    if (expressionFamily@vfamily=="negbinomial.size"){
      HSMM <- estimateDispersions(HSMM) # only with VGAM::negbinomial.size()
    }

    test = monocle::differentialGeneTest(
      HSMM[tested_genes,],
      fullModelFormulaStr=fullModelFormulaStr,
      reducedModelFormulaStr = reducedModelFormulaStr,
      cores=num_cores,
      relative_expr=F, 
      verbose=T)
    return(test)
  },

  equalizeSampleSize = function(pheno_category='timepoint', verbose=FALSE, seed=0, sample_size=NA){

    set.seed(seed)

    # if (pheno_category %in% c('timepoint', 'treatment', "cells_fulltreatment")) {
    if (pheno_category %in% colnames(pData(expressionSet))) {

      if (is.na(sample_size)){
        sample_size = min(table(pData(expressionSet)[, pheno_category]))
      }

      if (verbose){
        print("Cell counts per condition:")
        print(table(pData(expressionSet)[, pheno_category]))
        print(paste0("Selecting at most ", sample_size, " cells per condition"))
      }

      cells_ids_same_category = unlist(lapply(
              unique(pData(expressionSet)[, pheno_category]),
              function(d){
                cat_ids = which(pData(expressionSet)[, pheno_category]==d)
                if (length(cat_ids) <= sample_size) {
                  return(cat_ids)
                } else {
                  sample(
                      cat_ids,
                      sample_size,
                      replace=F)
                }
              }))

      remove_cells(ids =setdiff(1:ncol(exprs(expressionSet)), cells_ids_same_category))

      processing_state <<- c(processing_state, paste0('Exclude cells to equalize sample size by ', pheno_category, ' (dataset dimension:', nrow(exprs(expressionSet)), ' ', ncol(exprs(expressionSet)), ')'))

      # subsampling may leave some genes unexpressed in the retained dataset
      exclude_unexpressed_genes(min_cells=1, data_status="Raw", verbose=TRUE)

      processing_state <<- c(processing_state, paste0('Exclude genes unexpressed in remaining cells (dataset dimension:', nrow(exprs(expressionSet)), ' ', ncol(exprs(expressionSet)), ')'))

    } else {
      stop("Only the 'timepoint', 'treatment' and 'cells_fulltreatment' categories have been implemented")
    }
  }

)
