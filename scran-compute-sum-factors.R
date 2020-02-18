#!/usr/bin/env Rscript 

#computeSumFactors performs a scaling normalization of single-cell RNA-seq data by deconvolving size factors from cell pools. 
#Note: cells should have non-zero library sizes.  
#The assumption is that most genes are not differentially expressed (DE) between cells, such that any differences in expression across the majority of genes represents some technical bias that should be removed.

suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(workflowscriptscommon))

# argument parsing 
option_list = list(
  make_option(
    c("-i", "--input-sce-object"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Path to the input SCE object in rds format.'
  ),
    make_option(
    c("-a", "--assay-type"),
    action = "store",
    default = "logcounts",
    type = 'character',
    help = 'Specify which assay values to use. Default: "logcounts".'
  ),
   make_option(
    c("-s", "--sizes"),
    action = "store",
    default = NULL,
    type = 'numeric',
    help = 'A numeric vector of pool sizes, i.e., number of cells per pool.'
  ),
  make_option(
    c("-c", "--clusters"),
    action = "store",
    default = NULL,
    type = 'factor',
    help = 'An optional factor specifying which cells belong to which cluster, for deconvolution within clusters. For large data sets, clustering should be performed with the quickCluster function before normalization.'
  ),
  make_option(
    c("-r", "--subset-row"),
    action = "store",
    default = NULL,
    type = 'character',
    help = 'Logical, integer or character vector indicating the rows of SCE to use. If a character vector, it must contain the names of the rows in SCE.'
  ),
  make_option(
    c("-g", "--get-spikes"),
    action = "store",
    default = TRUE,
    type = 'logical',
    help = 'If get-spikes = FALSE, spike-in transcripts are automatically removed. If get.spikes=TRUE, no filtering on the spike-in transcripts will be performed.'
  ),
  make_option(
    c("-l", "--scaling"),
    action = "store",
    default = NULL,
    type = 'numeric',
    help = 'A numeric scalar containing scaling factors to adjust the counts prior to computing size factors.'
  ),
  make_option(
    c("-m", "--min_mean"),
    action = "store",
    default = 0,
    type = 'numeric',
    help = 'A numeric scalar specifying the minimum (library size-adjusted) average count of genes to be used for normalization.'
  ),
  make_option(
    c("-o", "--output-sce-object"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Path to the output SCE object containing the vector of size factors in sizeFactors(x).'
  )
)

opt = wsc_parse_args(option_list, mandatory = c("input_sce_object", "output_sce_object"))

#read SCE object
if(!file.exists(opt$input_sce_object)) stop("Input file does not exist.")
suppressPackageStartupMessages(require(SingleCellExperiment))
sce <- readRDS(opt$input_sce_object)

#compute size Factors
suppressPackageStartupMessages(require(scran))
sce <- computeSumFactors(sce, assay.type=opt$assay_type, sizes=opt$sizes, clusters=opt$custers, scaling=opt$scaling, min.mean = opt$min_mean, subset.row = opt$subset_row, get.spikes = opt$get_spikes)

#save SCE object with size Factors
saveRDS(sce, opt$output_sce_object)
