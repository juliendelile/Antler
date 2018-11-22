
#' Antler Reference Class
#'
#' @field expressionSet test
#' @field dR test
#' @include dimension_reduction.R

Antler <- setRefClass(
    Class = "Antler",

    fields = list(
        expressionSet="ExpressionSet",
        readcounts_raw="matrix",
        readcounts_norm="matrix",
        readcounts_smoothed="matrix",
        dR="DR",
        topCorr_DR="TopCorr_DR",
        pca_DR="PCA_DR",
        dispersion_DR="Dispersion_DR",
        cellClusters="list",
        plot_folder='character',
        num_cores='numeric',
        processing_state='character',
        write_report='logical',
        favorite_genes="character"
        ),  

    methods = list(

      initialize = function(..., expressionSet=ExpressionSet(), readcounts_raw=matrix(), readcounts_norm=matrix(), readcounts_smoothed=matrix(), dR=DR$new(), topCorr_DR=TopCorr_DR$new(), pca_DR=PCA_DR$new(), dispersion_DR=Dispersion_DR$new(), cellClusters=list('hclust'=NULL), plot_folder='./', num_cores=1, processing_state=character(0), write_report=TRUE, favorite_genes=character(0)) {

          # hack to obtain methods autocompletion in R session
          methods <- getRefClass(class(.self))$methods()
          eval(parse(text=paste0(".self$", methods)))

          callSuper(..., expressionSet=expressionSet, readcounts_raw=readcounts_raw, readcounts_norm=readcounts_norm, readcounts_smoothed=readcounts_smoothed, dR=dR, topCorr_DR=topCorr_DR, pca_DR=pca_DR, dispersion_DR=dispersion_DR, cellClusters = cellClusters, plot_folder=plot_folder, num_cores=num_cores, processing_state=processing_state, write_report=write_report, favorite_genes=favorite_genes)

          # Select random number generator ("L'Ecuyer-CMRG") 
          # compatible with parallel::mclapply (see ?parallel::mclapply)
          RNGkind("L'Ecuyer-CMRG")

          # Create plot directory
          if(length(plot_folder) == 1){
            dir.create(plot_folder, showWarnings=FALSE, recursive=TRUE)
          }
        },
    
      report = function(){cat(paste0(unlist(lapply(processing_state, function(r){paste0(r, '\n')})), collapse=''))},

      write_report_file = function(filepath=paste0(plot_folder, '/report.txt')) {
        write(paste0( "# Antler report file (last modification on ", 
                      fdate(),
                      ')\n\n',
                      paste0(unlist(lapply(processing_state, function(r){paste0(r, '\n')})), collapse='')),
              file=filepath)
      },

      loadDataset = function(folderpath=NULL, assayData=NULL, phenoData=NULL, data_status='Raw', default_gene_names='ensembl_gene_id'){

        if(is.null(folderpath)) {
          if(is.null(assayData) & is.null(phenoData)){
            stop('input dataset is missing')
          }
          expressionSet <<- ExpressionSet(
            assayData=assayData,
            phenoData=phenoData
            )
          processing_state <<- c(processing_state, paste0(fdate(), ' Dataset loaded from provided data.frames (', data_status, ')'))
        } else {
          expressionSet <<- ExpressionSet(
            assayData=as.matrix(read.table(paste0(folderpath, '/assayData.csv'), header=TRUE, sep="\t", row.names=1, as.is=TRUE, stringsAsFactors=F, check.names=FALSE)),
            phenoData=new("AnnotatedDataFrame", data=read.table(paste0(folderpath, '/phenoData.csv'), header=TRUE, sep="\t", row.names=1, as.is=TRUE, stringsAsFactors=F, check.names=FALSE))
          )

          if(file.exists(paste0(folderpath, '/featureData.csv'))){
            fd = read.table(file=paste0(folderpath, '/featureData.csv'), sep='\t', row.names=1, header=T, check.names=FALSE)
            for(n in colnames(fd)){
              fData(expressionSet)[,n] <<- fd[,n]
            }
          }

          processing_state <<- c(processing_state, paste0(fdate(), ' Dataset loaded from ', folderpath, ' (', data_status, ')'))
        }

        pData(expressionSet)$cells_samples <<- generateCellSampleNames(pData(expressionSet))

        pData(expressionSet)$cells_samples_fullname <<- generateCellSampleFullNames(pData(expressionSet))

        pData(expressionSet)$cells_fulltreatment <<- generateCellFullTreatments(pData(expressionSet))

        pData(expressionSet)$cells_colors <<- generateCellColors(pData(expressionSet))

        fData(expressionSet)[, default_gene_names] <<- rownames(exprs(expressionSet))
        
        fData(expressionSet)$current_gene_names <<- fData(expressionSet)[, default_gene_names]

        if(data_status=="Raw"){
          readcounts_raw <<- exprs(expressionSet)
          readcounts_norm <<- matrix()
          readcounts_smoothed <<- matrix()
        } else if(data_status=="Normalized"){
          readcounts_norm <<- exprs(expressionSet)
          readcounts_raw <<- matrix()
          readcounts_smoothed <<- matrix()
        } else if (data_status=="Smoothed"){
          readcounts_smoothed <<- exprs(expressionSet)
          readcounts_raw <<- matrix()
          readcounts_norm <<- matrix()
        } else {
          stop(paste0("The status of the dataset ('data_status') must be either 'Raw', 'Normalized' or 'Smoothed' (not '", data_status, "')."))
        }

        if(write_report)
          write_report_file()

      },

      setFavoriteGenes = function(genelist=NULL){

        unrecognized_genes = genelist[which(!genelist %in% getGeneNames())]

        if(length(unrecognized_genes) > 0){
          print(paste0("Some genes are not included in the dataset: ", paste0(unrecognized_genes, collapse=', ')))
        }

        favorite_genes <<- unique(genelist[which(genelist %in% getGeneNames())])
      },

      generateCellSampleColors = function(method='hsv', seed=10, hue_shift=.0){
        pData(expressionSet)$cells_colors <<- generateCellColors(pData(expressionSet), method=method, seed=seed, hue_shift=hue_shift)
        },

      getReadcounts = function(data_status=NULL){

        if(is.null(data_status)){
          stop("During code refactoring phase, 'data_status' must always be specified to avoid troubles later on...")
        }

        if(data_status=="Raw"){
          readcounts = readcounts_raw
        } else if(data_status=="Normalized"){
          readcounts = readcounts_norm
        } else if (data_status=="Smoothed"){
          readcounts = readcounts_smoothed
        } else {
          stop(paste0("The status of the dataset ('data_status') must be either 'Raw', 'Normalized' or 'Smoothed' (not '", data_status,"')."))
        }
        if(identical(readcounts, matrix())){
          stop(paste0(data_status, " dataset has not been loaded or calculated."))
        }
        rownames(readcounts) <- fData(expressionSet)$current_gene_names
        return(readcounts)
        },

      getGeneNames = function(){
        return(fData(expressionSet)$current_gene_names)
        },

      getNumberOfGenes = function(){
        return(nrow(exprs(expressionSet)))
        },

      getNumberOfCells = function(){
        return(ncol(exprs(expressionSet)))
        },

      getCellsNames = function(){
        return(colnames(exprs(expressionSet)))
        },

      getSampleSizes = function(){
        return(table(pData(expressionSet)$cells_samples))
        },

      plotReadCountsPerSample = function(basename=NULL, width=4, height=7, data_status='Raw'){

        readcounts = getReadcounts(data_status=data_status)
        dimnames(readcounts) = dimnames(exprs(expressionSet))
        exprs(expressionSet) <<- readcounts

        plotRCPerSample(expressionSet, filename=paste0(c(plot_folder, '/', paste0(c(basename, 'readcounts.pdf'), collapse='_')), collapse=''), width=width, height=height)

        processing_state <<- c(processing_state, paste0(fdate(), ' Plot read count per sample (', basename, 'readcounts.pdf',')'))

        if(write_report)
          write_report_file()
      },

      plotReadcountStats = function(data_status='Raw', basename=NULL, reads_name="UMI", by="timepoint", category=NA){

        plotnames = c()
        
        plotname = paste0(plot_folder, "/", basename, ifelse(is.null(basename), "", "_"), reads_name, "_statistics_all.pdf")
        plotnames <- c(plotnames, plotname)
        pdf(plotname, width=10, height=4, useDingbats=FALSE)
        par(mfrow = c(1,3))
        hist(log10(colSums(getReadcounts(data_status=data_status))), breaks=200, main=paste0("Total ", reads_name, " count per cell (log10)"), xlab="")
        hist(log10(colSums(getReadcounts(data_status=data_status)>0)), breaks=200, main=paste0("Expressed ", reads_name, " per cell (log10)"), xlab="")
        hist(log10(colMeans(getReadcounts(data_status=data_status))), breaks=200, main=paste0("Average ", reads_name, " count per cell (log10)"), xlab="")
        graphics.off()

        if(!identical(NA, category) & category %in% colnames(pData(expressionSet))) {

          stat.df = cbind.data.frame(
                      "TotalReadcounts"=colSums(getReadcounts(data_status=data_status)),
                      "Gene_count"=colSums(1*(getReadcounts(data_status=data_status)>0)),
                      "by"=pData(expressionSet)[, by],
                      "category"=pData(expressionSet)[, category]
                      ) %>%
                dplyr::group_by(by) %>%
                dplyr::mutate(category = factor(as.numeric(factor(category)))) %>%
                dplyr::ungroup()

          rplots <- lapply(unique(stat.df$by), function(t){
                ggplot(stat.df[which(stat.df$by==t),]) +
                  geom_freqpoly(aes(x=TotalReadcounts, linetype=factor(category)), position="identity", bins=30) +
                  scale_x_log10() +
                  theme(legend.position="none") + xlab(paste0("Total ", reads_name, " count")) +
                  ylab("Count") + 
                  ggtitle(paste0("Stage ",t)) + ylim(0, 800) +
                  scale_colour_gradient(limits=c(.5,1))
            })

          plotname = paste0(plot_folder, "/", basename, ifelse(is.null(basename), "", "_"), reads_name, "_statistics_", category, "_by_", by, ".pdf")
          plotnames <- c(plotnames, plotname)
          pdf(plotname, width=7, height=4, useDingbats=FALSE)
          gridExtra::grid.arrange(grobs=rplots, layout_matrix=matrix(seq(6), ncol=3, byrow=T))
          graphics.off()

          cat_colors = RColorBrewer::brewer.pal(9, "Set1")[seq(length(unique(pData(expressionSet)[, category])))]

          p <- stat.df[, c('by', 'category')] %>%
                dplyr::group_by(by, category) %>%
                dplyr::summarize(count=n()) %>%
                dplyr::ungroup() %>%
                ggplot() + geom_bar(aes(factor(by), weight=count, fill=factor(category)), color="black", size=.3) +
                ylab("Cells") + xlab(by) + scale_fill_manual(values=cat_colors, name = category) +
                theme_classic()

          plotname = paste0(plot_folder, "/", basename, ifelse(is.null(basename), "", "_"), "statistics_cellNumber_", category, "_by_", by, ".pdf")
          plotnames <- c(plotnames, plotname)
          pdf(plotname, width=5, height=3, useDingbats=FALSE)
          print(p)
          graphics.off()

          p <- stat.df %>% 
              ggplot() + geom_violin(aes(x=interaction(category, by), y=TotalReadcounts, fill=category)) +
                scale_y_log10() +
                ylab(paste0(reads_name, " per Cell")) + xlab(paste0(by, " X ", category)) + scale_fill_manual(values=cat_colors, name = category) +
                theme_minimal()

          plotname = paste0(plot_folder, "/", basename, ifelse(is.null(basename), "", "_"), "statistics_counts_", category, "_by_", by, ".pdf")
          plotnames <- c(plotnames, plotname)
          pdf(plotname, width=5, height=5, useDingbats=FALSE)
          print(p)
          graphics.off()

          p <- stat.df %>% 
              ggplot() + geom_violin(aes(x=interaction(category, by), y=Gene_count, fill=category)) +
                scale_y_log10() +
                ylab("Gene per Cell") + xlab(paste0(by, " X ", category)) + scale_fill_manual(values=cat_colors, name = category) +
                theme_minimal()

          plotname = paste0(plot_folder, "/", basename, ifelse(is.null(basename), "", "_"), "statistics_geneCounts_", category, "_by_", by, ".pdf")
          plotnames <- c(plotnames, plotname)
          pdf(plotname, width=5, height=5, useDingbats=FALSE)
          print(p)
          graphics.off()

        }

        processing_state <<- c(processing_state, paste0(fdate(), ' Plot ', reads_name, ' count statistics (', paste0(plotnames, collapse=', '), ')'))

        if(write_report)
          write_report_file()
      },

      # bug in biomart with mirror repository
      # see https://support.bioconductor.org/p/94725/
      # 2.31.6 has a fix but version not available for bioconductor 3.4
      setCurrentGeneNames = function(biomart_attribute='external_gene_name', biomart_dataset='mmusculus_gene_ensembl',  geneID_mapping_file=NULL){

        accepted.names <- c('external_gene_name', '...')
        if (!is.null(biomart_attribute))
          if(! biomart_attribute %in% accepted.names) 
            stop(paste0("invalid gene names category, it must be one of the biomaRt attributes ('", paste0(accepted.names, collapse="' or '"), "')"))

        if(is.null(geneID_mapping_file)){

          # It happens quite often that the Biomart server is not reachable.
          ensembl = tryCatch({
              useMart("ensembl", dataset=biomart_dataset)
          }, warning = function(w) {
              NA
          }, error = function(e) {
              cat('Biomart server can not be reached, please use the local file version.')
              NA
          })

          if(identical(ensembl, NA))
            return()
          # listMarts()    # to see which database options are present
          # ensembl=useMart("ensembl")
          # ensembl=useMart("ensembl", dataset="mmusculus_gene_ensembl")
          # listDatasets(ensembl)     # function to see which datasets are present in ensembl
          # ensembl=useDataset("mmusculus_gene_ensembl",mart=ensembl) 
          # listFilters(ensembl)  # check which filters are available
          # listAttributes(ensembl) # check attributes are available to select.More information on ensembl data base
          gene_names_ids.df=getBM(attributes=c(colnames(fData(expressionSet))[1], biomart_attribute), mart= ensembl)

          # some genomes does not have annotated "external_gene_name" (e.g. 'ggallus')
          missing_attr = which(gene_names_ids.df[[biomart_attribute]] == '')
          gene_names_ids.df[missing_attr, biomart_attribute] <- gene_names_ids.df[missing_attr, colnames(fData(expressionSet))[1]]

        } else {

          gene_names_ids.df=read.csv(geneID_mapping_file, stringsAsFactors=FALSE, header=TRUE)
          colnames(gene_names_ids.df) <- c(colnames(fData(expressionSet))[1], biomart_attribute)

        }
        
        # ! first column must have unique names (default: unique "ensembl_id")

        renamed_genes_ids = which(fData(expressionSet)$ensembl_gene_id %in% gene_names_ids.df$ensembl_gene_id)

        renamed_genes = fData(expressionSet)$ensembl_gene_id[renamed_genes_ids]

        rownames(gene_names_ids.df) <- gene_names_ids.df$ensembl_gene_id

        fData(expressionSet)$current_gene_names[renamed_genes_ids] <<- make.unique(gene_names_ids.df[renamed_genes ,'external_gene_name'])

        fData(expressionSet)[, biomart_attribute] <<- fData(expressionSet)$current_gene_names

        processing_state <<- c(processing_state, paste0('New gene names loaded (', biomart_attribute, ' from ', biomart_dataset, ')'))

        if(write_report){write_report_file()}
        },

      excludeCellsFromIds = function(cell_ids=NULL){

        if(length(cell_ids) > 0){

          expressionSet <<- expressionSet[, -cell_ids]

          if(!identical(readcounts_raw, matrix()))
            readcounts_raw <<- readcounts_raw[, -cell_ids]
          if(!identical(readcounts_norm, matrix()))
            readcounts_norm <<- readcounts_norm[, -cell_ids]
          if(!identical(readcounts_smoothed, matrix()))
            readcounts_smoothed <<- readcounts_smoothed[, -cell_ids]

          processing_state <<- c(processing_state, paste0(fdate(), ' Remove ', length(cell_ids), ' cells specified by user '))

          if(write_report){write_report_file()}

        }

        },

      excludeGenesFromIds = function(gene_ids=NULL){
        if(length(gene_ids) > 0){
          expressionSet <<- expressionSet[-gene_ids, , drop=F]

          if(!identical(readcounts_raw, matrix()))
            readcounts_raw <<- readcounts_raw[-gene_ids, , drop=F]
          if(!identical(readcounts_norm, matrix()))
            readcounts_norm <<- readcounts_norm[-gene_ids, , drop=F]
          if(!identical(readcounts_smoothed, matrix()))
            readcounts_smoothed <<- readcounts_smoothed[-gene_ids, , drop=F]
          if(length(favorite_genes) > 0)
            favorite_genes <<- favorite_genes[favorite_genes %in% getGeneNames()]
        }
      },

      excludeUnexpressedGenes = function(min.cells=3, min.level = 0, verbose=FALSE, data_status='Check'){
  
        readcounts = getReadcounts(data_status=data_status)

        unexpressed_genes_ids = which(rowSums(readcounts > min.level) < min.cells)

        log = paste0("Genes expressed in less than ", min.cells, " cells: ", paste(fData(expressionSet)$current_gene_names[unexpressed_genes_ids], collapse=" "))
        if(verbose){cat(paste0(log, '\n'))}
        processing_state <<- c(processing_state, log)
        if(write_report){write_report_file()}

        excludeGenesFromIds(unexpressed_genes_ids)

      },

      removeOutliers = function(lowread_thres = -Inf, highread_thres = Inf, cellmin = 10, genesmin = 1000, verbose=T, data_status='Check'){

        log.0 = paste0('Pre-filtering (dataset dimension:', nrow(exprs(expressionSet)), ' ', ncol(exprs(expressionSet)), ')')
        if(verbose){cat(paste0(log.0, '\n'))}
        processing_state <<- c(processing_state, log.0)
        
        # Cell filtering: select cells counting at least X reads
        if(!is.infinite(lowread_thres)){
          excludeCellsFromIds(cell_ids=which(apply(getReadcounts(data_status=data_status), 2, sum) <= lowread_thres))
        }
        log.1 = paste0('Filter I - select cells with at least ', lowread_thres, ' reads (dataset dimension:', nrow(exprs(expressionSet)), ' ', ncol(exprs(expressionSet)), ')')
        if(verbose){cat(paste0(log.1, '\n'))}
        processing_state <<- c(processing_state, log.1)

        # Cell filtering: select cells counting less than X reads
        if(!is.infinite(highread_thres)){
          excludeCellsFromIds(cell_ids=which(apply(getReadcounts(data_status=data_status), 2, sum) >= highread_thres))
        }
        log.1b = paste0('Filter II - select cells with less than ', highread_thres, ' reads (dataset dimension:', nrow(exprs(expressionSet)), ' ', ncol(exprs(expressionSet)), ')')
        if(verbose){cat(paste0(log.1b, '\n'))}
        processing_state <<- c(processing_state, log.1b)

        # Gene filtering: select genes expressed in at least N cells
        if(cellmin > 0){
          excludeGenesFromIds(gene_ids=which(rowSums(getReadcounts(data_status=data_status)>0) < cellmin))
        }
        log.2 = paste0('Filter III - select genes expressed in at least ', cellmin,' cells (dataset dimension:', nrow(exprs(expressionSet)), ' ', ncol(exprs(expressionSet)), ')')
        if(verbose){cat(paste0(log.2, '\n'))}
        processing_state <<- c(processing_state, log.2)

        # Cell filtering: select cells expressing at least N genes
        if(genesmin > 0){
          excludeCellsFromIds(cell_ids= which(colSums(getReadcounts(data_status=data_status)>0) < genesmin))
        }
        log.3 = paste0('Filter IV - select cells expressing at least ', genesmin, ' genes (dataset dimension:', nrow(exprs(expressionSet)), ' ', ncol(exprs(expressionSet)), ')')
        if(verbose){cat(paste0(log.3, '\n'))}
        processing_state <<- c(processing_state, log.3)
        if(write_report){write_report_file()}
        },
  
      removeGenesFromRatio = function(threshold=.1, candidate_genes=NA, verbose=F, plot_ratios=FALSE, file_path=NULL, data_status='Raw') {

        if(identical(candidate_genes, NA)){
          stop("candidate genes are missing.")
        }

        gene_ratio = colSums(getReadcounts(data_status=data_status)[candidate_genes,]) / colSums(getReadcounts(data_status=data_status))
        
        if(plot_ratios | !is.null(file_path)){
          if(!is.null(file_path)){
            pdf(file_path)
          }
          plot(colSums(getReadcounts(data_status=data_status)[candidate_genes, ]) / colSums(getReadcounts(data_status=data_status)), ylab="Candidate Gene Counts Ratio")
          if(!is.null(file_path)){
            graphics.off()
          }
        }

        high_cells = which(gene_ratio > threshold)
        
        excludeCellsFromIds(cell_ids = high_cells)

        excludeGenesFromIds(gene_ids = which(fData(expressionSet)$current_gene_names %in% candidate_genes))

        log = paste0('Removed ', length(high_cells),' cells containing more than ', round(threshold * 100), '% of candidate genes reads and then excluded all mitochondrial read counts (dataset dimension:', nrow(exprs(expressionSet)), ' ', ncol(exprs(expressionSet)), ', candidate genes: ', paste0(candidate_genes, collapse=', '), ')')

        if(verbose){cat(paste0(log, '\n'))}
        processing_state <<- c(processing_state, log)
        if(write_report){write_report_file()}

        },

      # Exclude gene with low level in expressing cells
      removeLowlyExpressedGenes = function(expression_threshold=1, selection_theshold=10, verbose=TRUE, data_status='Check'){

        # multiply gene level by binary mask to get only "positive" cells
        pos_cells_counts = getReadcounts(data_status=data_status) * (getReadcounts(data_status=data_status) > expression_threshold)
        pos_cells_counts[which(pos_cells_counts == 0)] <- NA
        pos_cells_mean_level = apply(pos_cells_counts, 1, mean, na.rm=T)
        
        excluded_genes_ids = which(pos_cells_mean_level < selection_theshold)

        excludeGenesFromIds(gene_ids = excluded_genes_ids)

        log = paste0('Removed ', length(excluded_genes_ids),' genes with an average expression level lower than ', selection_theshold, ' in "positive" cells, i.e. expression threshold higher than ', expression_threshold, ' (performed on ', data_status, ' dataset, dimension:', nrow(exprs(expressionSet)), ' ', ncol(exprs(expressionSet)), ')')

        if(verbose){cat(paste0(log, '\n'))}
        processing_state <<- c(processing_state, log)
        if(write_report){write_report_file()}

        },

      excludeLowSummaryReadcounts = function(used_gm='topCorr_DR.genemodules', nsigma=2, verbose=TRUE, data_status='Check'){

        if(! used_gm %in% getGeneModuleList())
          stop("invalid gene module(s) (must be one of the following: '", paste0(getGeneModuleList(), collapse="', '"), "')")

        gm_split = parseGMname(used_gm)

        if(length(.self[[gm_split$dr_name]][[gm_split$gm_name]]) == 0)
          stop("the gene module list is empty. Run gene modules identification first.")

        summary_readcounts = colSums(t(scale(t(log(getReadcounts(data_status=data_status)[unlist(.self[[gm_split$dr_name]][[gm_split$gm_name]]),]+1)), center=T, scale=T)))

        excluded_cell_ids = which(summary_readcounts < -nsigma*sd(summary_readcounts))

        excludeCellsFromIds(excluded_cell_ids)

        excludeUnexpressedGenes(min.cells=1, verbose=T, data_status=data_status)

        log = paste0('Removed ', length(excluded_cell_ids),' cells with a total of log-scaled readcounts over ', used_gm, ' genes farther than ', nsigma, ' sigma(s) from the mean (dataset dimension:', nrow(exprs(expressionSet)), ' ', ncol(exprs(expressionSet)), ')')

        if(verbose){cat(paste0(log, '\n'))}
        processing_state <<- c(processing_state, log)
        if(write_report){write_report_file()}

        },

      excludeLowSummaryFeatureCounts = function(used_gm='topCorr_DR.genemodules', nsigma=2, verbose=TRUE, data_status='Check'){

        if(! used_gm %in% getGeneModuleList())
          stop("invalid gene module(s) (must be one of the following: '", paste0(getGeneModuleList(), collapse="', '"), "')")

        gm_split = parseGMname(used_gm)

        if(length(.self[[gm_split$dr_name]][[gm_split$gm_name]]) == 0)
          stop("the gene module list is empty. Run gene modules identification first.")

        genemodules = .self[[gm_split$dr_name]][[gm_split$gm_name]]

        data.logscaled = t(scale(t(log(getReadcounts(data_status=data_status)[unlist(genemodules),]+1)), center=T, scale=T))

        mod_avg = do.call(rbind ,
                          lapply(genemodules, function(m){
                             colMeans(data.logscaled[m, , drop=F])
                           }))

        mod_avg.bin = silent.binarize.array(mod_avg)

        summary_features = colSums(mod_avg.bin)

        summary_features <- summary_features - mean(summary_features)

        excluded_cell_ids = which(summary_features < -nsigma*sd(summary_features))

        excludeCellsFromIds(excluded_cell_ids)

        excludeUnexpressedGenes(min.cells=1, verbose=T, data_status=data_status)

        log = paste0('Removed ', length(excluded_cell_ids),' cells with a total number of ', used_gm, ' features lower than ', nsigma, ' sigma(s) from the mean (dataset dimension:', nrow(exprs(expressionSet)), ' ', ncol(exprs(expressionSet)), ')')

        if(verbose){cat(paste0(log, '\n'))}
        processing_state <<- c(processing_state, log)
        if(write_report){write_report_file()}

        # return(summary_features)
        },

      normalize = function(method='Count-Per-Million', gene_length=NULL){

        if(identical(readcounts_raw, matrix())){
          stop("'normalize' operates on 'Raw' datasets. Please provide one.")
        }

        if(method=='Count-Per-Million' | method=="CPM"){
          # Normalize reads (for each cell, reads are divided by the sum of reads in the cell)
          readcounts_norm <<- 1e6 * sweep(readcounts_raw, 2, colSums(readcounts_raw), FUN="/")
        } else if(method=="FPKM") {
          readcounts_norm <<- 1e6 * sweep(readcounts_raw, 2, colSums(readcounts_raw), FUN="/")

          readcounts_norm <<- 1000 * sweep(readcounts_norm, 1, gene_length, "/")
        } else if(method=="geometric_mean_sizeFactors"){
          cell_count <- apply(readcounts_raw, 2, sum)
          geom_mean = exp(mean(log(cell_count)))
          sizefactors <- cell_count / geom_mean
          sizefactors[is.na(sizefactors)] <- 1
          readcounts_norm <<- t(t(readcounts_raw) / sizefactors)
        }
        else {
          stop('Only Count-Per-Million and FPKM normalization is implemented.')
        }
        
        processing_state <<- c(processing_state, paste0(fdate(), ' Gene levels normalized (', method, ')'))
        if(write_report){write_report_file()}
        },

      excludeExoticGenes = function(exotic_patterns = c('Gm[0-9].*', '.*-ps.*', 'ENSMUSG.*', '.*Rik', '.*-rs*'), verbose=TRUE){

        exotic_gene_ids = unique(unlist(lapply(exotic_patterns, function(ep){grep(ep, fData(expressionSet)$current_gene_names)})))

        excludeGenesFromIds(gene_ids=exotic_gene_ids)

        log = paste0('Excluded ', length(exotic_gene_ids), ' "exotic" gene(s) whose name matches the following patterns: ', paste0(exotic_patterns, collapse=', '), ' (dataset dimension: ', nrow(exprs(expressionSet)), ' ', ncol(exprs(expressionSet)), ')')

        if(verbose){cat(paste0(log, '\n'))}
        processing_state <<- c(processing_state, log)
        if(write_report){write_report_file()}

      },

      getGeneModuleList = function(){
        return(c(
            'topCorr_DR.genemodules',
            'topCorr_DR.genemodules.selected',
            'pca_DR.genemodules',
            'pca_DR.genemodules.selected',
            'dispersion_DR.genemodules',
            'dispersion_DR.genemodules.selected',
            'dR.genemodules',
            'dR.genemodules.selected'
            ))
      },

      parseGMname = function(gm_fullname){
        gm.split = strsplit(gm_fullname, split='\\.')[[1]]
        dr_name = gm.split[1]
        gm_name = paste0(gm.split[-1], collapse='.')
        return(list(dr_name=dr_name, gm_name=gm_name))
      },

      identifyGeneModules = function(method, corrgenes_t=NA, corr_t=NA, corr_quantile=FALSE, corr=NA, topcorr_corr_min=3, topcorr_mod_min_cell=0,  topcorr_mod_max_cell=Inf, topcorr_min_cell_level=0, topcorr_mod_consistency_thres=.4, topcorr_ordering_correlation_thres=NA, topcorr_ordering=NA, topcorr_num_final_gms=NULL, topcorr_num_max_final_gms=NULL, topcorr_num_min_final_gms=NULL, topcorr_num_initial_gms=NULL, topcorr_mod_skewness_thres=-Inf, topcorr_clustering_method="ward.D2", display=FALSE, verbose=TRUE, verbose.final=FALSE, debug_plot_basename=NULL, num_pcs=NULL, prop_variance_explained=.90, width=30, height=30, zscore_threshold=2.0, mean_threshold=10, clustering_mode=TRUE, data_status='Check') {

        if(method=='TopCorr_DR'){
         
          topCorr_DR <<- TopCorr_DR$new()

          topCorr_DR$run(readcounts=getReadcounts(data_status=data_status), readcounts_normalized=getReadcounts(data_status='Normalized'), corr_matrix=corr, corrgenes_t = corrgenes_t, corr_t=corr_t, corr_quantile=corr_quantile, corr_min=topcorr_corr_min, mod_min_cell=topcorr_mod_min_cell, mod_max_cell=topcorr_mod_max_cell, min_cell_level= topcorr_min_cell_level, mod_consistency_thres=topcorr_mod_consistency_thres, mod_skewness_thres=topcorr_mod_skewness_thres, ordering_correlation_thres=topcorr_ordering_correlation_thres, ordering=topcorr_ordering, num_final_gms=topcorr_num_final_gms, num_max_final_gms=topcorr_num_max_final_gms,  num_min_final_gms=topcorr_num_min_final_gms, num_initial_gms=topcorr_num_initial_gms, clustering_method=topcorr_clustering_method, display=display, numcores=num_cores, verbose=verbose, verbose.final=verbose.final, basename=if(is.null(debug_plot_basename))NULL else paste0(plot_folder, '/', debug_plot_basename))

          processing_state <<- c(processing_state, paste0(fdate(), ' TopCorr gene modules identification:'), paste0('    ', strsplit(topCorr_DR$report, split='\n')[[1]]))

          # dR <<- topCorr_DR

        } else if(method == 'PCA_DR') {

          pca_DR <<- PCA_DR$new()

          pca_DR$run(getReadcounts(data_status=data_status))

          pca_DR$extractGeneModules(basename=if(is.null(debug_plot_basename))NULL else paste0(plot_folder, '/', debug_plot_basename), num_pcs=num_pcs, prop_variance_explained=prop_variance_explained)

          processing_state <<- c(processing_state, paste0(fdate(), ' PCA gene modules identification:'), paste0('    ', strsplit(pca_DR$report, split='\n')[[1]]))

          if(!is.null(debug_plot_basename)){

            projected_coordinates = pca_DR$getCellCoordinates()
            
            for(i in seq(min(10, .5*length(pca_DR$genemodules)))) {
              x1=paste0('PC', i)
              x2=paste0('PC', i+1)
              pdf(paste0(plot_folder, '/', debug_plot_basename, '_PCAplot_PC', i, '_PC', i+1, '.pdf'), width=width, height=height)
                plot2D(m, projected_coordinates[, x1], projected_coordinates[, x2], xlabel=x1, ylabel=x2, cellname=FALSE)
              graphics.off()
            }

            processing_state <<- c(processing_state, paste0(fdate(), ' Save PCA plots (', debug_plot_basename, '_PCAplot_PCx_PCy.pdf', ')'))
          }

          # dR <<- pca_DR

        } else if(method=='Dispersion_DR'){

          dispersion_DR <<- Dispersion_DR$new()

          dispersion_DR$run(getReadcounts(data_status=data_status), zscore_threshold=zscore_threshold, mean_threshold=mean_threshold, numcores=num_cores, clustering_mode=clustering_mode)

          processing_state <<- c(processing_state, paste0(fdate(), ' Dispersion gene modules identification (', ifelse(clustering_mode, 'with', 'without'), ' clustering mode): '), paste0('    ', strsplit(dispersion_DR$report, split='\n')[[1]]))

          # dR <<- dispersion_DR

        } else {
          cat(paste0('"', method, '" is not a valid dimension reduction technique (Please use "TopCorr_DR", "PCA_DR" or "Dispersion_DR")\n'))
        }

        if(write_report){write_report_file()}
      },

      identifyCellClusters = function(method='hclust', clust_name=method, used_genes=NULL, cell_cell.euclidean.dist=NA, logscaled=TRUE, numclusters=NULL, data_status='Check'){

        accepted.clustering_method <- c('hclust')
        if(! method %in% accepted.clustering_method) 
          stop(paste0("invalid clustering method (must be ", paste0(accepted.clustering_method, collapse=' or '), ')'))

        if(length(used_genes) > 1){
          data.1 = getReadcounts(data_status=data_status)[unlist(used_genes), ]
          used_genes.report=paste0('a provided gene list (length: ', length(unlist(used_genes)),')')
        } else if(is.character(used_genes)){

          if(! used_genes %in% getGeneModuleList())
            stop("invalid gene module(s) (must be one of the following: '", paste0(getGeneModuleList(), collapse="', '"), "')")

          gm_split = parseGMname(used_genes)

          if(length(.self[[gm_split$dr_name]][[gm_split$gm_name]]) == 0)
            stop("the gene module list is empty. Run gene modules identification first.")

          data.1 = getReadcounts(data_status=data_status)[unlist(.self[[gm_split$dr_name]][[gm_split$gm_name]]), ]
          used_genes.report=paste0('the reduced dimensions "', used_genes, '"" (length: ', nrow(data.1), ')')

        } else {
          data.1 = getReadcounts(data_status=data_status)
          used_genes.report=paste0("all available genes (length ", nrow(data.1), ')') 
        }

        if(logscaled){
          data.2 = t(scale(t(log(data.1+1)), center=TRUE, scale=TRUE))

          # Check whether data.2 contains NaN
          na_genes = which(is.na(rowSums(data.2)))
          num_na = length(na_genes)
          if(num_na > 0){
            print(paste0("Warning: ", num_na, " genes contains NaN values after logscaled transformation (", paste0(rownames(data.2)[na_genes], collapse=', '), "). It is most likely null genes which should be removed upstream in the pipeline. These genes are ignored from current cell-cell distance calculation."))
            data.2[na_genes, ] <- 0
          }

        } else {
          data.2 = log(data.1+1)
        }

        if(all(is.na(cell_cell.euclidean.dist))){
          cell_cell.euclidean.dist = as.dist(fastEuclideanDist(t(data.2)))
        }

        if(method=='hclust'){
          hc.cells = hclust(cell_cell.euclidean.dist, method = "ward.D2")

          if(is.null(numclusters)){
            numclusters = getOptimalModuleNumber_HC(
              t(data.2), hc.cells, 100000, 1,
              display=FALSE,
              pdf_plot_path=NULL #paste0(plot_folder, '/clusters_NumModules_DistortionPlot.pdf')
              )
          }

          hc.cells.clust_ids = cutree(hc.cells, k=numclusters)

          # We order the cluster id from left to right
          clust_ids_ordered = unique(hc.cells.clust_ids[hc.cells$order]) %>% setNames(seq(numclusters), .)

          cellClusters[[clust_name]]$res <<- hc.cells
          cellClusters[[clust_name]]$cell_ids <<- hc.cells.clust_ids %>% {setNames(clust_ids_ordered[as.character(.)], names(.))}
          cellClusters[[clust_name]]$used_genes <<- used_genes
          cellClusters[[clust_name]]$gene_levels <<- logscaled
          cellClusters[[clust_name]]$num_clusters <<- numclusters
          cellClusters[[clust_name]]$names <<- NA
        }

        processing_state <<- c(processing_state, paste0(fdate(), ' Identify ', numclusters, ' cell clusters (method,', method, ' from the Euclidean distance matrix generated with ', ifelse(logscaled, 'the z-scored log-levels of ', ''), used_genes.report))
        if(write_report){write_report_file()}
        },

      setCellClustersNames = function(names, method="hclust"){
        cellClusters[[method]]$names <<- names 
      },

      # dendrogram NA (nothing), 'auto' (default or NULL), cluster names
      plotGeneModules = function(basename='gene_modules', displayed.geneset=fData(expressionSet)$current_gene_names, displayed.gms = 'all', file_settings=list(list(type='pdf', width=15, height=15)), logscaled=TRUE, use.dendrogram='auto', display.clusters=NULL, extra_colors=NA, extra_legend=NA, cell.ordering=NA, genemodules.palette, genemodules.extra_colors=NA, genemodules.text_colors=NA, genes.extra_colors=NA, display.legend=TRUE, zscore.thres=NULL, gene_transformations=c('logscaled', 'log', 'none'), data_status='Check', rect_overlay, curr_plot_folder=plot_folder, pretty.params=list("size_factor"=1, "ngenes_per_lines" = 8, "side.height.fraction"=.3), main = NULL){

        # accepted.genemodules <- setdiff(getGeneModuleList(), c('dR.genemodules', 'dR.genemodules.selected'))
        accepted.genemodules <- getGeneModuleList()

        if (! Reduce('&', displayed.gms %in% c('all', accepted.genemodules))) 
          stop("invalid gene module(s) (either '", paste0(accepted.genemodules, collapse="', '"), "' or 'all' must be used)")

        if(displayed.gms[1] == 'all'){
          displayed.gms = setdiff(accepted.genemodules, c('dR.genemodules', 'dR.genemodules.selected'))
        }

        cells_colors = as.matrix(pData(expressionSet)$cells_colors)

        if(display.legend){

          pd = pData(expressionSet)

          if(any(!is.na(cell.ordering))){
            pd <- pd[cell.ordering, ]
          }

          samples = pd[which(!duplicated(pd[,c('timepoint', 'replicate_id', 'treatment')])), ]

          legend = list(
                'text'=samples$cells_samples,
                'colors'=samples$cells_colors
                )
        } else {
          legend = NULL
        }

        if(!is.null(display.clusters)){

          # if(display.clusters == 'hclust'){
            
          if(!display.clusters %in% names(cellClusters)) {
            stop("'display.clusters' must a valid cell cluster name: ", names(cellClusters))
          }

          if(length(cellClusters[[display.clusters]]$cell_ids) != getNumberOfCells()){
            stop("Each cell needs a cluster id")
          }

          clust.colors = getClusterColors(v=2)

          cells_colors <- cbind(clust.colors[unlist(cellClusters[[display.clusters]]$cell_ids)], cells_colors)

          if(display.legend & FALSE){

            cluster_ids = unique(cellClusters[[display.clusters]]$cell_ids)

            cluster_names = paste0(cellClusters[[display.clusters]]$names[cluster_ids], " (", cluster_ids, ")")

            legend$text <- c(legend$text, cluster_names)
            

            legend$colors <- c(legend$colors, clust.colors[cluster_ids])
          }
        # } else if(!is.null(display.clusters)){
        #   stop("'display.clusters' must be NULL or 'hclust'")
        # }
        }

        if(!identical(extra_colors, NA)){
          cells_colors <- cbind(extra_colors[, rev(seq(ncol(extra_colors)))], cells_colors)
        }

        if(!identical(extra_legend, NA)){
          legend$text <- c(legend$text, extra_legend$text)
          legend$colors <- c(legend$colors, extra_legend$colors)
        }

        if(is.na(use.dendrogram)){
          dendrogram=NA
        } else if(use.dendrogram=='auto'){
          dendrogram=NULL
        } else if(use.dendrogram %in% c('hclust', names(cellClusters))) {
          dendrogram=as.dendrogram(cellClusters[[use.dendrogram]]$res)
        } else {
          stop("'dendrogram' must be either NA, 'auto' or 'hclust'.")
        }

        for(ds in data_status){

          readcounts = getReadcounts(data_status=ds)
          cells_colors_ordered = cells_colors
          
          if(any(!is.na(cell.ordering))){
            readcounts <- readcounts[, cell.ordering]
            cells_colors_ordered <- cells_colors[cell.ordering, , drop=F]
          }

          filenames = c()
          for(norm in gene_transformations){

            if(norm=='logscaled'){
              data.logscaled = t(scale(t(log(readcounts+1)), center=TRUE, scale=TRUE))
              logscaled=TRUE
              hm.palette=colorRampPalette(c("#191d73", "white", "#ed7901"))(n = 1000)
            } else if(norm=='log'){
              data.logscaled = log(readcounts+1)
              logscaled=FALSE
              hm.palette=colorRampPalette(c("#0464DF", "#FFE800"))(n = 1000) # blue -> yellow
            } else if(norm=='none'){
              data.logscaled = readcounts
              logscaled=FALSE
              hm.palette=colorRampPalette(c("#0464DF", "#FFE800"))(n = 1000) # blue -> yellow
            }

            for(gm in displayed.gms){

              gm_split = parseGMname(gm)

              for(fs in file_settings){

                fn = paste0(basename, '_', gm, '_', ds, '_', norm, '.', fs$type)
                filenames = c(filenames, fn)

                if(length(.self[[gm_split$dr_name]][[gm_split$gm_name]]) > 0){
                  hm = plotHeatmapCellsGeneModules(
                      mat=data.logscaled,
                      genemodules=.self[[gm_split$dr_name]][[gm_split$gm_name]],
                      filename=paste0(curr_plot_folder, '/', fn),
                      cells_colors=cells_colors_ordered,
                      genemodules.palette=if(missing(genemodules.palette)) colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(length(.self[[gm_split$dr_name]][[gm_split$gm_name]])) else genemodules.palette,
                      genemodules.extra_colors=genemodules.extra_colors,
                      genemodules.text_colors=genemodules.text_colors,
                      genes.extra_colors=genes.extra_colors,
                      displayed.geneset=displayed.geneset,
                      legend=legend,
                      width=fs$width, height=fs$height,
                      heatmap.palette=hm.palette,
                      zscore.thres = ifelse(is.null(zscore.thres), ifelse(logscaled, 2, 1e10), zscore.thres),
                      dendrogram=dendrogram,
                      rect_overlay=rect_overlay,
                      pretty.params=pretty.params,
                      main=main
                      )
                }
              }
            }
          }

          processing_state <<- c(processing_state, paste0(fdate(), ' Plot ', basename, ' heatmap(s) (', paste0(filenames, collapse=', '), ')'))

        }

        if(write_report){write_report_file()}

      },

      writeGeneModules = function(basename='Genelist', gms='all', folder_path=plot_folder){

        accepted.genemodules <- getGeneModuleList()

        if (! Reduce('&', gms %in% c('all', accepted.genemodules))) 
          stop("invalid gene module(s) (either '", paste0(accepted.genemodules, collapse="', '"), "' or 'all' must be used)")

        if(gms[1] == 'all'){
          gms = setdiff(accepted.genemodules, c('dR.genemodules', 'dR.genemodules.selected'))
        }

        for(gm in gms){

          gm_split = parseGMname(gm)

          if(length(.self[[gm_split$dr_name]][[gm_split$gm_name]]) > 0){
            
            genemodules=.self[[gm_split$dr_name]][[gm_split$gm_name]]

            gm.filename = paste0(folder_path, '/', basename, '_', gm, '.txt')
            file.create(gm.filename)

            for(i in seq(length(genemodules))){

              modname = names(genemodules)[i]
              if(is.null(modname)){
                modname = paste0('G.M.', i)
              }

              write(paste0(modname, ': ', paste0(genemodules[[i]], collapse=', ')), file=gm.filename, append=TRUE)
            }
            
          }
        }

      },

      exportExpressionSet = function(folderpath='./data/new_dataset', data_status='Check'){

        dir.create(file.path(folderpath), showWarnings = FALSE)

        readcounts = getReadcounts(data_status=data_status)
        dimnames(readcounts) = dimnames(exprs(expressionSet))
        exprs(expressionSet) <<- readcounts

        write.table(x=readcounts, file=paste0(folderpath, '/assayData.csv'), sep='\t', row.names=TRUE, quote=FALSE, col.names=NA)

        write.table(x=pData(expressionSet), file=paste0(folderpath, '/phenoData.csv'), sep='\t', row.names=TRUE, quote=FALSE, col.names=NA)
          
        },

      excludeCellFromClusterIds = function(cluster_ids=NULL, used_clusters='hclust', data_status='Check'){
          
        # accepted.clustering_method <- c('hclust')
        accepted.clustering_method <- names(cellClusters)
        if(! used_clusters %in% accepted.clustering_method) 
          stop(paste0("invalid clustering method (must be ", paste0(accepted.clustering_method, collapse=' or '), ')'))

        if(is.null(cellClusters[[used_clusters]])){
          stop(paste0('the "', used_clusters, '" cell clusters must be calculated first with "identifyCellClusters"'))
        }

        excluded_cell_ids = which(cellClusters[[used_clusters]]$cell_ids %in% cluster_ids)

        excludeCellsFromIds(cell_ids=excluded_cell_ids)

        processing_state <<- c(processing_state, paste0('Exclude cells belong to "', used_clusters, '" cluster (', length(excluded_cell_ids), ' cells)'))

        excludeUnexpressedGenes(min.cells=3, data_status=data_status)

        processing_state <<- c(processing_state, paste0('Exclude genes expressed in less than 3 cells (dataset dimension:', nrow(exprs(expressionSet)), ' ', ncol(exprs(expressionSet)), ')'))

        },

      identifyCellPopulations = function(cluster_markers, used_gm, basename=NULL, data_status='Check'){

        if(! used_gm %in% getGeneModuleList())
          stop("invalid gene module(s) (must be one of the following: '", paste0(getGeneModuleList(), collapse="', '"), "')")

        gm_split = parseGMname(used_gm)
        gms = .self[[gm_split$dr_name]][[gm_split$gm_name]]

        if(length(gms) == 0)
          stop("the gene module list is empty. Run gene modules identification first.")
        
        # Remove pop / Trim markers if marker absent from gene modules
        cluster_markers <- unlist(lapply(seq(length(cluster_markers)), function(i){    
            l = cluster_markers[[i]]
            l_out = l[!l %in% unlist(gms)]
            if(length(l_out)){
              print(paste0(names(cluster_markers)[[i]], " is missing: ", paste0(l_out, collapse="', ")))
            }
            l_in=list(l[l %in% unlist(gms)])
            names(l_in) = names(cluster_markers)[[i]]
            if(length(l_in[[1]])==0) NULL else l_in
            }), recursive=F)

        # Number of clusters is the number of excluded populations plus the remaining cells cluster
        num_clusters = length(cluster_markers) + 1

        # For each excluded population, gather the modules containing specified marker genes
        selected.genemodules = lapply(cluster_markers, function(pm){
                  unlist(Filter(function(i){length(intersect(i, pm)) > 0}, gms))
          })

        # For safety purpose, assure that proper gene names are used
        selected.genemodules <- lapply(selected.genemodules, function(l){intersect(l, fData(expressionSet)$current_gene_names)})

        data.logscaled = t(scale(t(log(getReadcounts(data_status=data_status)[unlist(selected.genemodules),]+1)), center=T, scale=T))

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
            if(exclude_clust){
              clust_name=names(which.max(clust.mean))
              excluded_clust_id[clust_name] = i
              excluded_clust_count[clust_name] = length(which(hc.cells.clust_ids %in% i))
            }
        }

        if(!is.null(basename)){ 

            clust.colors = c(RColorBrewer::brewer.pal(9, "Set1"), RColorBrewer::brewer.pal(12, "Set3"), RColorBrewer::brewer.pal(8, "Set2"))

            hm = plotHeatmapCellsGeneModules(mat=data.logscaled, genemodules=selected.genemodules, filename=paste0(plot_folder, '/', basename, '_heatmap.pdf'), cells_colors=cbind(clust.colors[hc.cells.clust_ids], pData(expressionSet)$cells_colors), displayed.geneset=rownames(data.logscaled), zscore.thres=2, legend=list('text'=names(excluded_clust_id), 'colors'=clust.colors[unlist(excluded_clust_id)]), dendrogram=as.dendrogram(hc.cells))
        }

        processing_state <<- c(processing_state, paste0('Identify cells belong to the following clusters: ', paste0(paste0(names(cluster_markers), ' (', unlist(excluded_clust_count[names(cluster_markers)]), ' cells)'), collapse=', '), ' (dataset dimension:', nrow(exprs(expressionSet)), ' ', ncol(exprs(expressionSet)), ')'))
        
        if(write_report){write_report_file()}
        
        pop_ids = lapply(excluded_clust_id, function(i){
                which(hc.cells.clust_ids %in% i)
          })

        # return(pop_ids[names(cluster_markers)])
        return(pop_ids)

        },

      identifyCellPopulationsFromExistingCellClusters = function(cluster_markers, used_gm, basename=NULL, data_status='Check', existing_clusters="hclust", ratio_thres=.5){

        if(! used_gm %in% getGeneModuleList())
          stop("invalid gene module(s) (must be one of the following: '", paste0(getGeneModuleList(), collapse="', '"), "')")

        gm_split = parseGMname(used_gm)

        if(length(.self[[gm_split$dr_name]][[gm_split$gm_name]]) == 0)
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

        data.logscaled = t(scale(t(log(getReadcounts(data_status=data_status)[unlist(selected.genemodules),]+1)), center=T, scale=T))

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
        
        if(write_report){write_report_file()}
        
        return(excluded_clust_id)

        },

      excludeCellPopulations = function(cluster_markers, used_gm, basename=NULL, data_status='Check'){

        cell_to_exclude_ids = identifyCellPopulations(cluster_markers=cluster_markers, used_gm=used_gm, basename=basename, data_status=data_status)

        excludeCellsFromIds(cell_ids=unlist(cell_to_exclude_ids))

        excludeUnexpressedGenes(min.cells=3, data_status=data_status)

        processing_state <<- c(processing_state, paste0('Exclude genes expressed in less than 3 cells (dataset dimension:', nrow(exprs(expressionSet)), ' ', ncol(exprs(expressionSet)), ')'))

        if(write_report){write_report_file()}

        },

      # Monocle 2 (Install first with http://cole-trapnell-lab.github.io/monocle-release/getting-started/)
      runMonocle2DETestOnClusters = function(clustering_method='hclust', ref_clusters=c("1"), alt_clusters=unique(cellClusters[[clustering_method]]$cell_ids), dispersion_selection=TRUE) {

        closeAllConnections()

        df_pheno = data.frame(
          row.names=colnames(getReadcounts(data_status='Normalized')),
          "sample"=pData(expressionSet)$cells_samples_fullname,
          "cluster_id"=cellClusters[[clustering_method]]$cell_ids,
          "test_clusters"=1*(cellClusters[[clustering_method]]$cell_ids %in% ref_clusters)
        )

        used_genes = if(dispersion_selection){getDispersedGenes(getReadcounts(data_status='Normalized'), zscore_threshold=0.0)}else{nrow(getNumberOfGenes())}

        df_feature = data.frame(row.names=used_genes, "gene_short_name"=used_genes)

        HSMM <- monocle::newCellDataSet(as.matrix(getReadcounts(data_status='Normalized')[used_genes,]), phenoData = new("AnnotatedDataFrame", data = df_pheno), featureData = new('AnnotatedDataFrame', data = df_feature), expressionFamily=VGAM::tobit()) # "tobit() is appropriate for FPKM, TPM gene expression data."

        diff_test_res <- monocle::differentialGeneTest(HSMM[, which(df_pheno$cluster_id %in% unique(c(ref_clusters, alt_clusters)))], fullModelFormulaStr="~test_clusters", reducedModelFormulaStr = "~1", cores=num_cores)

        return(diff_test_res)

      },

      # Monocle 2 (Install first with http://cole-trapnell-lab.github.io/monocle-release/getting-started/)
      runMonocle2DETest = function(fullModelFormulaStr, reducedModelFormulaStr = "~1", expressionFamily=VGAM::negbinomial.size(), tested_genes=getGeneNames()) {

        closeAllConnections()

        df_pheno = pData(expressionSet)
        df_feature = cbind(fData(expressionSet), 'gene_short_name'=fData(expressionSet)$current_gene_names)
        rownames(df_feature) <- df_feature$gene_short_name

        HSMM <- monocle::newCellDataSet(
              as.matrix(getReadcounts(data_status='Raw')),
              phenoData = new("AnnotatedDataFrame", data = df_pheno),
              featureData = new('AnnotatedDataFrame', data = df_feature),
              expressionFamily=expressionFamily
              )

        HSMM <- estimateSizeFactors(HSMM)

        if(expressionFamily@vfamily=="negbinomial.size"){
          HSMM <- estimateDispersions(HSMM) # only with VGAM::negbinomial.size()
        }

        test = monocle::differentialGeneTest(HSMM[tested_genes,],
                                    fullModelFormulaStr=fullModelFormulaStr,
                                    reducedModelFormulaStr = reducedModelFormulaStr,
                                    cores=num_cores,
                                    relative_expr=F, 
                                    verbose=T)

        return(test)
      },

      equalizeSampleSize = function(pheno_category='timepoint', data_status='Check', verbose=FALSE, seed=0, sample_size=NA){

        set.seed(seed)

        if(pheno_category %in% c('timepoint', 'treatment', "cells_fulltreatment")) {

          if(is.na(sample_size)){
            sample_size = min(table(pData(expressionSet)[, pheno_category]))
          }

          if(verbose){
            print("Cell counts per condition:")
            print(table(pData(expressionSet)[, pheno_category]))
            print(paste0("Selecting at most ", sample_size, " cells per condition"))
          }

          cells_ids_same_category = unlist(lapply(
                  unique(pData(expressionSet)[, pheno_category]),
                  function(d){
                    cat_ids = which(pData(expressionSet)[, pheno_category]==d)
                    if(length(cat_ids) <= sample_size) {
                      return(cat_ids)
                    } else {
                      sample(
                          cat_ids,
                          sample_size,
                          replace=F)
                    }
                  }))

          excludeCellsFromIds(cell_ids=setdiff(1:ncol(exprs(expressionSet)), cells_ids_same_category))

          processing_state <<- c(processing_state, paste0('Exclude cells to equalize sample size by ', pheno_category, ' (dataset dimension:', nrow(exprs(expressionSet)), ' ', ncol(exprs(expressionSet)), ')'))

          # subsampling may leave some genes unexpressed in the retained dataset
          excludeUnexpressedGenes(min.cells=1, data_status=data_status, verbose=TRUE)

          processing_state <<- c(processing_state, paste0('Exclude genes unexpressed in remaining cells (dataset dimension:', nrow(exprs(expressionSet)), ' ', ncol(exprs(expressionSet)), ')'))

        } else {
          stop("Only the 'timepoint', 'treatment' and 'cells_fulltreatment' categories have been implemented")
        }
      },

      selectGeneModulesFromGenes = function(priorknowledge.gene.selection, original_gms){

        if(! original_gms %in% getGeneModuleList())
          stop("invalid gene module(s) (must be one of the following: '", paste0(getGeneModuleList(), collapse="', '"), "')")

        gm_split = parseGMname(original_gms)

        if(length(.self[[gm_split$dr_name]][[gm_split$gm_name]]) > 0){

          kept_modules_ids = Filter(function(i){
                length(intersect(.self[[gm_split$dr_name]][[gm_split$gm_name]][[i]], priorknowledge.gene.selection)) > 0
           }, 1:length(.self[[gm_split$dr_name]][[gm_split$gm_name]]))

          .self[[gm_split$dr_name]]$genemodules.selected <- .self[[gm_split$dr_name]][[gm_split$gm_name]][kept_modules_ids]
        } else {
          .self[[gm_split$dr_name]]$genemodules.selected <- list()
        }

      },
        
      displayGeneModuleIdsFromGenes = function(priorknowledge.gene.selection=favorite_genes, original_gms="topCorr_DR.genemodules"){

        if(! original_gms %in% getGeneModuleList())
          stop("invalid gene module(s) (must be one of the following: '", paste0(getGeneModuleList(), collapse="', '"), "')")

        gm_split = parseGMname(original_gms)

        if(length(.self[[gm_split$dr_name]][[gm_split$gm_name]]) > 0){

          kept_modules_ids = lapply(seq(length(.self[[gm_split$dr_name]][[gm_split$gm_name]])), function(i){
                   sort(intersect(.self[[gm_split$dr_name]][[gm_split$gm_name]][[i]], priorknowledge.gene.selection))      
            })
          names(kept_modules_ids) <- seq(length(.self[[gm_split$dr_name]][[gm_split$gm_name]]))

          kept_modules_ids <- Filter(function(x){length(x) != 0}, kept_modules_ids)

          not_kept = sort(setdiff(priorknowledge.gene.selection, unlist(kept_modules_ids)))
          if(length(not_kept) > 0)
            kept_modules_ids[['Not in store']] <- not_kept

          # kept_modules_ids = Filter(function(i){
          #       length(intersect(.self[[gm_split$dr_name]][[gm_split$gm_name]][[i]], priorknowledge.gene.selection)) > 0
          #  }, 1:length(.self[[gm_split$dr_name]][[gm_split$gm_name]]))

          print(kept_modules_ids)
        } else {
          print(paste0(original_gms, " is empty..."))
        }

      },

      selectGeneModulesFromIds = function(kept_modules_ids, original_gms){

        if(! original_gms %in% getGeneModuleList())
          stop("invalid gene module(s) (must be one of the following: '", paste0(getGeneModuleList(), collapse="', '"), "')")

        gm_split = parseGMname(original_gms)

        if(length(.self[[gm_split$dr_name]][[gm_split$gm_name]]) > 0){

          .self[[gm_split$dr_name]]$genemodules.selected <- .self[[gm_split$dr_name]][[gm_split$gm_name]][kept_modules_ids]
        } else {
          .self[[gm_split$dr_name]]$genemodules.selected <- list()
        }

        },

      reportFavoriteGenesSelection = function(favgenes = m$favorite_genes){

        remaininggenes = favgenes

        print("Favorite genes excluded from:")

        for(gm_names in names(topCorr_DR$genemodules.history)){

          # gm_split = parseGMname(used_gm)

          if(length(topCorr_DR$genemodules.history[[gm_names]]) == 0)
            stop("the gene module list is empty. Run gene modules identification first.")

          excl_genes = setdiff(remaininggenes, unlist(topCorr_DR$genemodules.history[[gm_names]]))
          print(paste0("      ", gm_names, ": ", paste0(sort(excl_genes), collapse=", ")))

          remaininggenes <- intersect(remaininggenes, unlist(topCorr_DR$genemodules.history[[gm_names]]))
        }

        print(paste0("Favorite genes remaining: ", paste0(sort(remaininggenes), collapse=', ')))

      }

    )
  )
