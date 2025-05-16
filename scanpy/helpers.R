#' Print a Markdown link for a file
#'
#' @param fn Filename
#' @return String ready to be printed
linkify <- function(fn){paste('[', fn, '](', fn, ')')}

#' Print Markdown (frequently used to programatically print headers).
#' Handles header levels and additional newlines. Make sure you're in a chunk
#' with results='asis' for this to work.
#'
#' @param level Markdown header level. 1="#", 2="##", etc
#' @param ... Arguments passed directly to cat
#' @return NULL; as a side effect the text is printed.

mdcat <- function(..., level=NA){
  dots <- paste(...)
  if (!is.na(level)){
    header <- paste0(rep('#', level), collapse="")
    dots <- paste(header, dots)
  }
  cat(paste0('\n\n', dots, '\n\n'), fill=1500)
}

#' Get sample paths from sample table
#'
#' This function takes the path to a sampletable and
#' returns a named vector of paths.
#'
#' Sample table must have the following columns
#' - 'samplename': This is used as 'orig.ident'
#' - 'data_path': absolute path to cellranger output directory
#' - 'group': (optional) This can be used to specify groups. Samples
#'            sharing the same value of 'group' will be pooled.
#'            This is done by adding suffixes, e.g. '_1', '_2', to
#'            the 'group' value and using this as 'orig.ident'.
#'
#' @param f path to sampletable
#' @param suffix suffix added to data paths
#'
get_sample_paths <- function(f, suffix='outs/filtered_feature_bc_matrix'){
  df <- read.table(f, sep='\t', header=TRUE, stringsAsFactors=FALSE)

  reqd_cols <- c('samplename', 'data_path')
  if(!all(reqd_cols %in% colnames(df))){
    msg <- paste0('Sample table must contain the following columns: ',
                  paste(reqd_cols, collapse=', '))
    stop(msg)
  }

  data_paths <- df[,'data_path']
  data_paths <- file.path(data_paths, suffix)

  if(!'group' %in% colnames(df)){
    names(data_paths) <- df[, 'samplename']
  } else {
    # get group sizes
    grp_sizes <- table(df[, 'group'])
    name_vec <- unlist(
                  lapply(unique(df[, 'group']),
                    function(x) paste0(x, '_', 1:grp_sizes[x])
                  )
                )
    names(data_paths) <- name_vec
  }

  # check to make sure the data paths exist
  exists <- dir.exists(data_paths)
  if(sum(!exists) == 0){
    return(data_paths)
  } else {
    msg <- paste0('Some data paths do not exist!\n',
        paste(data_paths[which(!exists)], collapse='\n'), '\n')
    stop(msg)
  }
}

#' Function to generate cliff-knee plot for a sample
#'
#' @param sample_path path to cell-ranger sample output directory
#'        containing 'outs/filtered_feature_bc_matrix'
#'
get_cliff_knee <- function(sample_path){
  # strip trailing 'filtered_feature_bc_matrix' if it exists
  if(grepl('filtered_feature_bc_matrix', sample_path)){
    d <- dirname(sample_path)
  } else {
    d <- sample_path
  }

  # list of paths to matrix mtx file, raw & filtered barcodes.tsv
  f <- list(
         mtx=file.path(d, 'raw_feature_bc_matrix/matrix.mtx.gz'),
         raw_bc=file.path(d, 'raw_feature_bc_matrix/barcodes.tsv.gz'),
         filt_bc=file.path(d, 'filtered_feature_bc_matrix/barcodes.tsv.gz')
       )

  # read the sparse matrix
  if(!file.exists(f[['mtx']])){
    stop(
      paste0('Matrix file: "', f[['mtx']], '" does not exist!')
    )
  }

  if(grepl('\\.gz$', f[['mtx']])){
    sparse.mat <- Matrix::readMM(gzfile(f[['mtx']]))
  } else {
    sparse.mat <- Matrix::readMM(f)
  }

  # read raw & filtered barcodes
  bc_list <- list()
  for(bc in c('raw_bc', 'filt_bc')){
    if(!file.exists(f[[ bc ]])){
      stop(
        paste0('Barcodes file: "', f[[ bc ]], '" does not exist!')
      )
    }

    # NOTE: here we read the barcodes into a data.frame
    #       with a single column, so '$V1' selects that
    #       column
    if(grepl('\\.gz$', f[[ bc ]])){
      bc_list[[ bc ]] <- read.table(gzfile(f[[ bc ]]),
                             header=FALSE, sep='\t')$V1
    } else {
      bc_list[[ bc ]] <- read.table(f[[ bc ]],
                             header=FALSE, sep='\t')$V1
    }
  }

  cnts <- Matrix::colSums(sparse.mat)

  idx <- bc_list[['raw_bc']] %in% bc_list[['filt_bc']]
  filtered <- ifelse(idx, 'cells', 'background')
  cnts_df <- data.frame(umi_counts=cnts,
                        filtered=filtered)

  cnts_df <- cnts_df[order(cnts_df$umi_counts, decreasing=TRUE),]
  cnts_df$barcodes <- 1:nrow(cnts_df)

  p <- ggplot(cnts_df, aes(x=barcodes,
                           y=umi_counts,
                           color=filtered)) +
       geom_line(size=1.5) +
       scale_x_continuous(trans='log10') +
       scale_y_continuous(trans='log10') +
       scale_color_manual(name='',
                          limits=c('cells', 'background'),
                          values=c('cells'='black',
                                   'background'='grey'))
  p + theme_bw()
}

#' Use scDblFinder package to identify doublets.
#'
#' This is done on the raw counts,
#' and for each sample separately. We need to temporarily convert to SCE for the
#' purposes of doublet detection. The resulting SeuratObj has additional columns:
#'
#' - scDblFinder.class
#' - scDblFinder.score
#' - scDblFinder.weighted
#' - scDblFinder.cxds_score
#'
#' The class and score columns are probably most useful.
#'
#' @param x SeuratObj
#' @return SeuratObj, with additional meta.data columns
find_doublets <- function(x){
  sce <- as.SingleCellExperiment(x)
  sce <- scDblFinder(sce)
  so <- as.Seurat(sce)

  # Otherwise the default ident is a fixed "SingleCellExperiment"
  Idents(so) <- 'orig.ident'
  return(so)
}


#' Add percent mitochondrial and percent ribosomal genes for QC.
#'
#' @param obj SeuratObj
#' @param mito_pattern regex pattern used to detect mitochondrial genes
#' @param ribo_pattern regex pattern used to detect ribosomal genes
#' @param mito_label metadata column name for mitochondrial %
#' @param ribo_label metadata column name for ribosomal %
#'
#' @return updated SeuratObj
add_percentage_features <- function(obj,
                                    mito_pattern='^mt-|Mt_|MT-',
                                    ribo_pattern='^Rp[sl]|RP[SL]',
                                    mito_label='percent.mito',
                                    ribo_label='percent.ribo'){

  obj <- PercentageFeatureSet(
    obj,
    pattern=mito_pattern,
    col.name=mito_label
  )

  obj <- PercentageFeatureSet(
    obj,
    pattern=ribo_pattern,
    col.name=ribo_label
  )
  return(obj)
}

#' Find thresholds for QC metrics
#'
#' Use scuttle::isOutlier() to detect cells with MAD > 3 and return in df
#'
#' @param obj SeuratObj
#' @param name identity class for which thresholds are being calculated
#' @param metrics character vector of QC metrics
#'
find_thresholds <- function(obj, name, metrics){
  f <- function(metric){
    res <- scuttle::isOutlier(obj[[metric]][[metric]])
    res %>%
      attr("thresholds") %>%
      t() %>%
      as.data.frame() %>%
      mutate(metric=metric, orig.id=name) %>%
      dplyr::select(orig.id, metric, everything()) %>%
      tibble::remove_rownames()
  }
  m <- purrr::map_df(metrics, f)
}

#' Print a Markdown link for a file
#'
#' @param fn Filename
#' @return String ready to be printed
linkify <- function(fn){paste('[', fn, '](', fn, ')')}

#' Saves ggplot object as PNG and PDF, and then prints out links to files in
#' Markdown.
save_and_show <- function(p, prefix, ...){
  png.fn <- paste0(prefix, '.png')
  pdf.fn <- paste0(prefix, '.pdf')
  ggsave(p, file=png.fn, ...)
  ggsave(p, file=pdf.fn, ...)
  mdcat('- ', linkify(png.fn))
  mdcat('- ', linkify(pdf.fn))
  mdcat('\n\n')
}

#' Get conserved markers for given cluster
#'
#'
#' Given a cluster, return the genes consistently found across samples that are
#' enriched in this cluster vs all others. Modify the output so it has cluster
#' and gene as additional columns, and remove the rownames so we can combine
#' multiple data.frames together
#'
#' @param obj SeuratObj
#' @param cluster.i cluster to find conserved markers for
#' @param findmarkers_test test to use for FindMarkers
#' @param min.pct minimum fraction of cells in either of the two populations expressing
#'        a gene
#' @param logfc.threshold average minimum difference between groups for which
#'        genes will be tested
#' @param min.cells.feature Minimum number of cells expressing the feature in at
#'        least one of the two groups
#' @param min.cells.group minimum number of cells in one of the groups
#' @param cluster_col metadata column used to define clusters
#'
get_conserved <- function(obj, cluster.i, findmarkers_test,
                          min.pct=0.1, logfc.threshold=0.25,
                          min.cells.feature=3, min.cells.group=3,
                          cluster_col='seurat_clusters'){
  # set Idents to 'cluster_col'
  Idents(obj) <- obj@meta.data[, cluster_col]

  # first check to make sure all identity classes have cells
  # in 'cluster.i', if not, skip cluster
  idx <- obj@meta.data[, cluster_col] == cluster.i
  ncells <- table(obj@meta.data[idx, 'orig.ident'])
  idents <- unique(obj@meta.data[, 'orig.ident'])
  if(all(idents %in% names(ncells)) & all(ncells > 3)){
    df <- FindConservedMarkers(
      obj,
      ident.1=cluster.i,
      assay='SCT',
      slot='data',
      test.use=findmarkers_test,
      grouping.var="orig.ident",
      only.pos=TRUE
    ) %>%
      rownames_to_column(var="gene")

    # return df if # rows > 0
    if(nrow(df) > 0)
      return(cbind(cluster=cluster.i, df))
    else NULL
  } else {
    warning(
      paste0('Skipping ', cluster.i, ' since not enough cells in each identity class (', cluster_col, ')\n')
    )
    NULL
  }
}


# Define a function generating subchunks with DT::datatable(dataframe) each
subchunkify <- function(name, input='df', width=7, height=5, ggplotly=TRUE) {
    if (input == 'df') {
        t_deparsed <- paste0("DT::datatable(df)")
        more <- ""
    } else if (input == 'plot') {
        t_deparsed <- paste0("plotly::ggplotly(p)")
        if (!ggplotly) {
            t_deparsed <- paste0("print(p)")
        }
        more <- paste0(", fig.width=", width, ", fig.height=", height)
    } else {
        stop('Incorrect argument: input')
    }
    sub_chunk <- paste0("```{r sub_chunk_", name, ", results='asis', echo=FALSE", more, "}",
        "\n",
        t_deparsed,
        "\n\n```\n\n\n")

    cat(knitr::knit(text = sub_chunk, quiet=TRUE))
}


#' Print a Markdown link for a plot
#'
#' @param fn filepath
#' @return String ready to be printed
link.plot <- function(fn){
    cat(paste('\n\n- Download Plot: [', fn, '](', fn, ')\n\n'))
}

#' Print a Markdown link for a table
#'
#' @param fn filepath
#' @return String ready to be printed
link.table <- function(fn){
    cat(paste('\n\n- Download Table: [', fn, '](', fn, ')\n\n'))
}

#' Run log-normalization and sketch on a merged seurat object
#'
#' @param obj.merged merged seurat object
#' @param variable.features.n number of variable features
#' @return seurat object
norm_sketch <- function(obj.merged, variable.features.n) {
    # Preprocess the obj
    DefaultAssay(obj.merged) <- 'RNA'
    obj.merged <- NormalizeData(
        obj.merged,
        verbose=FALSE) %>%
        FindVariableFeatures(
            verbose=FALSE,
            nfeatures=variable.features.n
            )

    # Run sketch
    obj.merged <- SketchData(
        obj.merged,
        verbose=FALSE,
        sketched.assay='sketch'
    )
    return(obj.merged) 
}

#' Integrate sketched seurat object
#'
#' @param obj.merged merged seurat object, sketched
#' @param variable.features.n number of variable features
#' @return integrated seurat object
sketched_integration <- function(obj.merged, variable.features.n) {

    # Preprocess the sketched obj
    DefaultAssay(obj.merged) <- 'sketch'

    # Retrieve N variable features computed by `SCTransform()`
    features <- VariableFeatures(obj.merged)
    obj.merged <- FindVariableFeatures(
        obj.merged,
        verbose=FALSE,
        nfeuatures=variable.features.n)

    # Scale and reduce dimensions
    obj.merged <- ScaleData(
        obj.merged,
        features=features,
        verbose=FALSE) %>%
        RunPCA(
        assay='sketch',
        features=features,
        reduction.name='sketch.pca',
        reduction.key='sketchpca_',
        verbose=FALSE) 

    # Integrate layers
    obj.integrated <- IntegrateLayers(
        obj.merged,
        method=CCAIntegration,
        orig.reduction='sketch.pca',
        assay="sketch",
        features=features,
        new.reduction='integrated.cca',
        verbose=FALSE
        )

    # Integrate embeddings from the sketched assay
    obj.integrated <- ProjectIntegration(
        obj.integrated,
        verbose=FALSE,
        sketched.assay='sketch',
        method='sketch',
        assay="RNA",
        reduction='integrated.cca') %>%
        # Project high dimensional scRNA expression data from a full dataset
        # onto the lower dimensional embedding of the sketch of the dataset
        ProjectData(
            sketched.assay='sketch',
            assay='RNA',
            normalization.method='SCT',
            sketched.reduction='integrated.cca.full',
            full.reduction='integrated.cca.full',
            verbose=FALSE,
            dims=1:30
        )
    return(obj.integrated)
}

#' Get the OrgDb for the specified organism, using the cached AnnotationHub.
#'
#' @param species Case-sensitive genus and species
#' @param cache Directory in which the AnnotationHub cache is stored
#' @param annotation_key_override If not NA, forces the hub to use this
#'        accession. Use this when you know exactly which OrgDb to use.
#'
#' @return OrgDb object
get.orgdb <- function(species, cache, annotation_key_override=NA){

    # Workaround to allow AnnotationHub to use proxy. See
    # https://github.com/Bioconductor/AnnotationHub/issues/4, and thanks
    # Wolfgang!
    proxy <- Sys.getenv('http_proxy')
    if (proxy == ""){
        proxy <- NULL
    }

    if (!dir.exists(cache)){
        dir.create(cache, recursive=TRUE)
    }

    ah <- AnnotationHub(hub=getAnnotationHubOption('URL'),
             cache=cache,
             proxy=proxy,
             localHub=FALSE)

    find.annotationhub.name <- function(species.name, override.code) { #autodetect ah names based on loaded database
        if (is.na(override.code)) {
        ah.query <- query(ah, "OrgDb")
        ah.query.speciesmatch <- grepl(paste("^", species.name, "$", sep=""), ah.query$species)
        ah.query.which <- which(ah.query.speciesmatch)
        stopifnot(length(ah.query.which) > 0) #require at least one match
        if (length(ah.query.which) > 1) { #warn of duplicate matches
           print("WARNING: found multiple candidate species in AnnotationHub: ");
           print(ah.query.speciesmatch)
        }
        names(ah.query)[ah.query.which[1]]
        } else {
        override.code
        }
    }
    annotation_key <- find.annotationhub.name(annotation_genus_species, annotation_key_override)
    orgdb <- ah[[annotation_key]]
    return(orgdb)
}

#' Writes out original and split clusterprofiler results
#'
#' @param res clusterProfiler results
#' @param cprof.folder Directory in which to save files. Will be created if needed.
#' @param Label to use for the results. Will generate a filename
#' cprof.folder/label.txt and cprof.folder/label_split.txt
#'
#' @return List of the two files that are created on disk.
write.clusterprofiler.results <- function(res, cprof.folder, label){
    dir.create(cprof.folder, showWarnings=FALSE, recursive=TRUE)
    filename.orig <- file.path(cprof.folder, paste0(label, '.txt'))
    write.table(res, file=filename.orig, sep='\t', quote=FALSE, row.names=FALSE)
    filename.split <- file.path(cprof.folder, paste0(label, '_split.txt'))
    res.split <- split.clusterProfiler.results(res)
    write.table(res.split, file=filename.split, sep='\t', quote=FALSE, row.names=FALSE)
    return(list(orig=filename.orig, split=filename.split))
}

#' Split clusterProfiler output into one line per gene
#'
#' @param x Results from clusterProfiler. It is expected that the
#' clusterProfiler enrichment function was called with "readable=TRUE"
#'
#' @return data.frame with genes one per line instead of "/" separated in one
#' line. The rest of the original line is repeated for each gene.
split.clusterProfiler.results <- function(x){
    df <- x@result
    # loop over all rows
    df.parse <- NULL
    for(k in 1:dim(df)[1]){
        g <- strsplit(as.character(df$geneID[k]), "/")[[1]]
        gg <- df[rep(k, length(g)),]
        gg$geneParse <- g
        if(is.null(df.parse)){
            df.parse <- gg
        } else {
            df.parse <- rbind(df.parse, gg)
        }
    }
    return(df.parse)
}


