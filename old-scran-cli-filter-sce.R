#!/usr/bin/env Rscript 

#Perform random downsampling of SCE object.  For testing purposes only. 
suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(workflowscriptscommon))

# parse options
option_list = list(
  make_option(
    c("-i", "--input-sce-object"),
    action = "store",
    default = NA,
    type = 'character',
    help = "File name in which a serialized R SingleCellExperiment object where object matrix found"
  ),
  make_option(
    c("-e", "--exprs-values"),
    action = "store",
    default = 'counts',
    type = 'character',
    help= "String indicating which assay contains the count data that should be used to compute log-transformed expression values."
  ),
  make_option(
    c("-g", "--n-genes"),
    action = "store",
    default = 1200,
    type = 'integer',
    help= "Integer indicating how many genes should we sample from the input object."
  ),
   make_option(
    c("-c", "--n-cells"),
    action = "store",
    default = 100,
    type = 'integer',
    help= "Integer indicating how many cells should we sample from the input object."
  ),
 make_option( 
    c("-k", "--spikes"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'String or integer scalar specifying the alternative experiment containing the spike-in transcripts.'
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

#check number of GENES to sample is not larger that available
if(nrow(sce) > opt$n_genes){
    sce <- sce[sample(rownames(sce), opt$n_genes),]
}    

#check number of CELLS to sample is not larger that available
if(ncol(sce) > opt$n_cells){
    sce <- sce[, sample(colnames(sce), opt$n_cells)]
}

#add some genes as spike-ins, so that scran_computeSpikeFactors.R can run without error
suppressPackageStartupMessages(require(SingleCellExperiment))
if(!opt$spikes %in% spikeNames(sce)){
  print("Spike-ins not present in input object -- Adding random genes as spike-ins for running scran-compute-spike-factors.R")
  isSpike(sce, "ERCC") <-  c(sample(rownames(sce), 25))
}
#Output to a serialized R object
saveRDS(sce, file = opt$output_sce_object)

