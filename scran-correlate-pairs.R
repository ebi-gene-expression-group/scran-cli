#!/usr/bin/env Rscript 

#Identify pairs of genes that are significantly correlated in their expression profiles, based on Spearman's rank correlation. 
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
    c("-b", "--block"),
    action = "store",
    default = NULL,
    type = 'factor',
    help = 'A factor specifying the blocking levels for each cell in sce, for instance a donor covariate. If specified, correlations are computed separately within each block and statistics are combined across blocks.'
  ),
  make_option(
    c("-d", "--design"),
    action = "store",
    default = NULL,
    type = 'matrix',
    help = 'A numeric design matrix containing uninteresting factors to be ignored.'
  ),
  make_option(
    c("-a", "--assay-type"),
    action = "store",
    default = "counts",
    type = 'character',
    help = 'A character specifying which assay values to use.'
  ),
#  make_option( #[argument available on SCRAN 1.14]
#    c("-e", "--equiweight"),
#    action = "store",
#    default = TRUE,
#    type = 'logical',
#    help = 'A logical scalar indicating whether statistics from each block should be given equal weight. Otherwise, each block is weighted according to its number of cells. Only used if block is specified.'
#  ),
  make_option(
    c("-k", "--iters"),
    action = "store",
    default = 1e+06,
    type = 'logical',
    help = 'Integer scalar specifying the number of iterations to use in correlateNull to build the null distribution.'
  ),
  make_option(
    c("-u", "--use-names"),
    action = "store",
    default = TRUE,
    type = 'logical',
    help = 'A logical scalar specifying whether the row names of x should be used in the output. Alternatively, a character vector containing the names to use.'
  ),
  make_option(
    c("-s", "--subset_row"),
    action = "store",
    default = NULL,
    type = 'logical',
    help = 'Logical, integer or character vector specifying the rows for which to model the variance. Defaults to all genes in x.'
  ),
  make_option(
    c("-g", "--get-spikes"),
    action = "store",
    default = FALSE,
    type = 'logical',
    help = 'Logical specifying wether to perform spike-ins filtering(FALSE) or not (TRUE).'
  ),
  make_option(
    c("-m", "--use-dimred"),
    action = "store",
    default =  FALSE,
    type = 'character',
    help = 'A string specifying whether existing values in reducedDims(x) should be used.'
  ),
  make_option(
    c("-o", "--output-pairs-df"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Path to the output dataframe with one row per gene pair (rows order by increasing p-values) in txt format.'
  )
)

opt = wsc_parse_args(option_list, mandatory = c("input_sce_object", "output_pairs_df"))

#read SCE object
if(!file.exists(opt$input_sce_object)) stop("Input SCE file does not exist.")
sce <- readRDS(opt$input_sce_object)

#Compute gene pair-correlation
suppressPackageStartupMessages(require(scran))
corr_pairs <- correlatePairs(sce, block=opt$block, subset.row=opt$subset_row, assay.type=opt$assay_type, 
                            design=opt$design,  use.names=opt$use_names, get.spikes=opt$get_spikes)

#save gene-pairs table
write.table(corr_pairs, file = opt$output_pairs_df, sep = "\t")
