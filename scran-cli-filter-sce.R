#!/usr/bin/env Rscript 

#This is a non SCRAN function script to: Perform basic QC to SCE object: remove low expressed genes across cells, remove cells with low genes expressed,  add randow genes as spike-ins if these not present (for testing purposes), and log transform CPM if specified. 
suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(workflowscriptscommon))

# parse options
option_list = list(
  make_option(
    c("-i", "--input-sce-object"),
    action = "store",
    default = NA,
    type = 'character',
    help = "Path to the input SCE object in rds format."
  ),
  make_option(
    c("-e", "--exprs-values"),
    action = "store",
    default = 'counts',
    type = 'character',
    help= "String indicating which assay contains the count data that should be used to compute log-transformed expression values."
  ),
 make_option( 
    c("g", "--min-genes"),
    action = "store",
    default = 200,
    type = 'integer',
    help = 'Minimum number of genes to be expressed per cell to pass filtering.'
  ), 
 make_option( 
    c("c", "--min-cells"),
    action = "store",
    default = 50,
    type = 'integer',
    help = 'Minimum number of cells for a gene to be expressed to pass filtering.'
  ), 
 make_option( 
    c("-k", "--spikes"),
    action = "store",
    default = "ERCC",
    type = 'character',
    help = 'String or integer scalar specifying the alternative experiment containing the spike-in transcripts. Default; "ERCC"'
  ), 
 make_option( 
    c("-n", "--n-spikes"),
    action = "store",
    default = 25,
    type = 'integer',
    help = 'Integer specifying the number of genes to add as spike-ins in case there is are no spike-ins in the "ERCC" slot. Default: 25.'
  ), 
 make_option( 
    c("-l", "--log"),
    action = "store",
    default = FALSE,
    type = 'logical',
    help = 'Logical specifying wether log-transformation of CPM counts should be performed.'
  ), 
 make_option(
    c("-o", "--output-sce-object"),
    action = "store",
    default = NA,
    type = 'character',
    help = "File name in which to store serialized R object of type 'SingleCellExperiment'."
  )
)

opt <- wsc_parse_args(option_list, mandatory = c('input_sce_object',  'output_sce_object'))

#read SCE object
if(!file.exists(opt$input_sce_object)) stop("Input file does not exist.")
sce <- readRDS(opt$input_sce_object)

#if object has no colnames, append
if(is.null(colnames(sce)))colnames(sce) <- sce$Barcode

#filter out genes and cells based on parameters
sce <- sce[apply(assay(sce, opt$exprs_values), 1, FUN=function(x) sum(x > 0) >= opt$min_cells), apply(assay(sce, opt$exprs_values), 2, FUN=function(x) sum(x > 0) >=opt$min_genes)]

if(nrow(sce) ==0 | ncol(sce) ==0) stop("One of 2 dimensions of SCE object is 0 due to filtering! Change min-genes or min-cells parameters.")

#add some genes as spike-ins, so that scran_computeSpikeFactors.R can run without error
suppressPackageStartupMessages(require(SingleCellExperiment))
if(!opt$spikes %in% spikeNames(sce)){
  print("Spike-ins not present in input object -- Adding random genes as spike-ins for running scran-compute-spike-factors.R")
	if(opt$n_spikes > nrow(sce)) stop("Number of spike-ins to select is higher than the number of genes in the SCE object.")
  	isSpike(sce, "ERCC") <-  c(sample(rownames(sce), 25))
}

if(opt$log ==T){
	#normalize by CPM
	assay(sce, "normcounts") <- apply(assay(sce, "counts"),2, function(x) (x/sum(x))*1000000) 
	#log transform normcounts
	assay(sce, "logcounts") <- log2(assay(sce, "normcounts") + 1)
}
#Output to a serialized R object
saveRDS(sce, file = opt$output_sce_object)
