#!/usr/bin/env Rscript 

#Convert a SCESet object into other classes for entry into other analysis pipelines (edgeR, DESeq2, monocle).  

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
    default = NULL,
    type = 'character',
    help = 'A string specifying the analysis for which the object should be prepared. Any of: "edgeR", "DESeq2", "monocle".'
  ),
#  make_option( #argument available in SCRAN 1.14
#    c("-f", "--fData_col"),
#    action = "store",
#    default = NULL,
#    type = 'numeric',
#    help = 'Any set of indices specifying which columns of fData(x) should be retained in the returned object.'
#  ),
#  make_option( #argument available in SCRAN 1.14
#    c("-p", "--pData_col"),
#    action = "store",
#    default = NULL,
#    type = 'numeric',
#    help = 'Any set of indices specifying which columns of pData(x) should be retained.'
#  ),
  make_option( #this argument is "assay" in SCRAN 1.14
    c("-a", "--assay-type"),
    action = "store",
    default = "counts",
    type = 'character',
    help = 'A string specifying which assay of x should be put in the returned object. For edgeR and DESeq2, assay is set to "counts" such that count data is stored in the output object'
  ),
#  make_option(#argument available in SCRAN 1.14
#    c("-n", "--normalize"),
#    action = "store",
#    default = TRUE,
#    type = 'logical',
#    help = 'A logical scalar specifying whether the assay values should be normalized for type="monocle".'
#  ),
  make_option(
    c("-g", "--get-spikes"),
    action = "store",
    default = FALSE,
    type = 'logical',
    help = 'Logical specifying wether to perform spike-ins filtering(FALSE) or not (TRUE).'
  ),
  make_option(
    c("-o", "--output-converted"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Path to the converted object stored as RDS.'
  )
)

opt = wsc_parse_args(option_list, mandatory = c("input_sce_object", "type", "output_converted"))

#read SCE object
if(!file.exists(opt$input_sce_object)) stop("Input file does not exist.")
sce <- readRDS(opt$input_sce_object)

#Compute PCA and denoise it
suppressPackageStartupMessages(require(scran))

#Compute gene pair-correlation
converted_object <- convertTo(sce, type=opt$type, assay.type=opt$assay_type, 
			#fData.col=opt$fData_col, 
			#dData.col=opt$dData_col, 
			# normalize=opt$normalize, 
			get.spikes=opt$get_spikes)

#save converted object to RDS
saveRDS(converted_object, file = opt$output_converted)
