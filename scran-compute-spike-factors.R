#!/usr/bin/env Rscript 

#Normalization based on the spike-in counts (Lun et al. 2017). 
#REQUIREMNT: The input SCE object MUST have spike-in assay in altExp(SCE).
#Size factors are computed to scale the counts such that the total coverage of the spike-in transcripts is equal across cells. 
#The main practical difference is that spike-in normalization preserves differences in total RNA content between cells, whereas computeSumFactors and other non-DE methods do not.

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
    c("-t", "--type"),
    action = "store",
    default = "ERCC",
    type = 'character',
    help = 'A character vector specifying which spike-in sets to use. Default: "ERCC".'
  ),
  #make_option( #[argument in SCRAN 1.14, equivalent to current 'type']
  #  c("-s", "--spikes"),
  #  action = "store",
  #  default = "ERCC",
  #  type = 'character',
  #  help = 'String or integer scalar specifying the alternative experiment containing the spike-in transcripts.'
  #),
    make_option(
    c("-a", "--assay-type"),
    action = "store",
    default = "logcounts",
    type = 'character',
    help = 'Specify which assay values to use. Default: "logcounts".'
  ),
  make_option(
    c("-g", "--general_use"),
    action = "store",
    default = TRUE,
    type = 'logical',
    help = 'A logical scalar indicating whether the size factors should be stored for general use by all genes.'
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
sce <- readRDS(opt$input_sce_object)

#check the input SCE has spike-in assays
suppressPackageStartupMessages(require(SingleCellExperiment))
if(length(spikeNames(sce)) == 0){
  stop("Input SCE does not have spike-ins as alternative Experiment assay")
}
#check provided spike-in names are present in the SCE object
if(!opt$type %in% spikeNames(sce)){
  stop("Provided spike-in name is not present in the input SCE")
}
#compute size Factors
suppressPackageStartupMessages(require(scran))
sce <- computeSpikeFactors(sce, type=opt$type, 
			#spikes=opt$spikes,#argument available in SCRAN 1.14 
			assay.type=opt$assay_type, general.use=opt$general_use)

#save SCE object with size Factors
saveRDS(sce, opt$output_sce_object)
