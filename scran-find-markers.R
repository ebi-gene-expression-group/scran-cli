#!/usr/bin/env Rscript 

#Find candidate marker genes for groups of cells (e.g., clusters) by testing for differential expression between pairs of groups. 
#If x is scale-normalized but not log-transformed, it can be used with test.type="wilcox" and test.type="binom". If x contains raw counts, it can only be used with test.type="binom".
#Output is a named list of DataFrames, each of which contains a sorted marker gene list for the corresponding group.
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
    c("-c", "--clusters"), #This argument is named "groups" on SCRAN 1.14
    action = "store",
    default = NULL,
    type = 'character',
    help = 'A vector of group assignments for all cells.'
  ),
#  make_option( [argument available on SCRAN 1.14]
#    c("-t", "--test-type"),
#    action = "store",
#    default = "t",
#    type = 'character',
#    help = 'String specifying the type of pairwise test to perform - a t-test with "t", a Wilcoxon rank sum test with "wilcox", or a binomial test with "binom".'
#  ),  
  make_option(
    c("-p", "--pvalue-type"),
    action = "store",
    default = "any",
    type = 'character',
    help = 'A character specifying how p-values are to be combined across pairwise comparisons for a given group/cluster. Setting pval.type="all" requires a gene to be DE between each cluster and every other cluster (rather than any other cluster, as is the default with pval.type="any").'
  ),
  make_option(
    c("-s", "--subset_row"),
    action = "store",
    default = NULL,
    type = 'logical',
    help = 'Logical, integer or character vector specifying the rows for which to model the variance. Defaults to all genes in x.'
  ),
  make_option(
    c("-a", "--assay-type"),
    action = "store",
    default = "counts",
    type = 'character',
    help = 'A character specifying which assay values to use.'
  ),
  make_option(
    c("-k", "--get-spikes"),
    action = "store",
    default = FALSE,
    type = 'logical',
    help = 'Logical specifying wether to perform spike-ins filtering(FALSE) or not (TRUE).'
  ),
  make_option(
    c("-o", "--output-markers"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Path to the rds list of DataFrames with a sorted marker gene list per cluster/group.'
  )
)

opt = wsc_parse_args(option_list, mandatory = c("input_sce_object", "clusters", "output_markers"))

#read SCE object
if(!file.exists(opt$input_sce_object)) stop("Input file does not exist.")
sce <- readRDS(opt$input_sce_object)

#check cluster annotation exists
if(!opt$clusters %in% names(sce@colData))stop("Cluster annotation specified not present in SCE")

#Find Marker Genes
suppressPackageStartupMessages(require(scran))
markers <- findMarkers(sce, clusters=sce[[opt$clusters]], subset.row = opt$subset_row, assay.type = opt$assay_type, pval.type=opt$pvalue_type, get.spikes=opt$get_spikes)

#save markers: list of DataFrames, each of which contains a sorted marker gene list for the corresponding group.
saveRDS(sce, opt$output_markers)
