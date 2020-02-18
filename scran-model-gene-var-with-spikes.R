#!/usr/bin/env Rscript 

#Function to perform modelling the per-gene variance. 
#We decompose the total variance of each gene into its biological and technical components by fitting a trend to the endogenous variances(A. T. Lun, McCarthy, and Marioni 2016).

suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(workflowscriptscommon))

# argument parsing 
option_list = list(
  make_option(
    c("-i", "--input-sce-object"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Path to the input SCE object in rds format'
  ),
  make_option(
    c("-b", "--block"),
    action = "store",
    default = NULL,
    type = 'factor',
    help = 'A factor specifying the blocking levels for each cell in sce, for instance a donor covariate. If specified, variance modelling is performed separately within each block and statistics are combined across blocks.'
  ),
  make_option(
    c("-f", "--size-factors"),
    action = "store",
    default = NULL,
    type = 'numeric',
    help = 'A numeric vector of cell-specific size factors. Alternatively NULL, in which case the size factors are extracted or computed from x.'
  ),
  make_option(
    c("-s", "--spikes"),
    action = "store",
    default = "ERCC",
    type = 'character',
    help = 'String or integer scalar specifying the alternative experiment containing the spike-in transcripts.'
  ),
  make_option(
    c("-k", "--spike-size-factors"),
    action = "store",
    default = NULL,
    type = 'numeric',
    help = 'A numeric vector of cell-specific size factors. Alternatively NULL, in which case the size factors are extracted or computed from x.'
  ),
  make_option(
    c("-d", "--design"),
    action = "store",
    default = NULL,
    type = 'matrix',
    help = 'A numeric matrix containing blocking terms for uninteresting factors of variation.'
  ),
  make_option(
    c("-s", "--subset_row"),
    action = "store",
    default = NULL,
    type = 'logical',
    help = 'Logical, integer or character vector specifying the rows for which to model the variance. Defaults to all genes in x.'
  ),
  make_option(
    c("-b", "--block"),
    action = "store",
    default = NULL,
    type = 'factor',
    help = 'A factor specifying the blocking levels for each cell in sce, for instance a donor covariate. If specified, variance modelling is performed separately within each block and statistics are combined across blocks.'
  ),
    make_option(
    c("-d", "--design"),
    action = "store",
    default = NULL,
    type = 'matrix',
    help = 'A numeric matrix containing blocking terms for uninteresting factors of variation.'
  ),
  make_option(
    c("-s", "--subset_row"),
    action = "store",
    default = NULL,
    type = 'logical',
    help = 'Logical, integer or character vector specifying the rows for which to model the variance. Defaults to all genes in x.'
  ),
  make_option(
    c("-f", "--subset-fit"),
    action = "store",
    default = NULL,
    type = 'logical',
    help = 'Logical, integer or character vector specifying the rows to be used for trend fitting. Defaults to subset.row.'
  ),
  make_option(
    c("-a", "--assay-type"),
    action = "store",
    default = "logcounts",
    type = 'character',
    help = 'String or integer scalar specifying the assay containing the log-expression values.'
  ),
  make_option(
    c("-m", "--min-mean"),
    action = "store",
    default = 0.1,
    type = 'numeric',
    help = 'A numeric scalar specifying the minimum mean to use for trend fitting.'
  ),
  make_option(
    c("-p", "--parametric"),
    action = "store",
    default = FALSE,
    type = 'character',
    help = 'A logical scalar indicating whether a parametric fit should be attempted. f parametric=TRUE, a non-linear curve of the form: y = ax/(x^n + b) s fitted to the variances against the means.'
  ),
  make_option(
    c("-e", "--equiweight"),
    action = "store",
    default = TRUE,
    type = 'logical',
    help = 'A logical scalar indicating whether statistics from each block should be given equal weight. Otherwise, each block is weighted according to its number of cells. Only used if block is specified.'
  ),
make_option(
    c("-m", "--method"),
    action = "store",
    default = "fisher",
    type = 'logical',
    help = 'String specifying how p-values should be combined when block is specified, see ?combinePValues.'
  ),
  make_option(
    c("-o", "--output-geneVarSpikes-table"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Path to the table where each row corresponds to a gene in sce, and contains: mean, total var, bio var, tech var, p.value and FDR'
  )
)

opt = wsc_parse_args(option_list, mandatory = c("input_sce_object", "output_geneVarSpikes_table"))

#read SCE object
if(!file.exists(opt$input_sce_object)) stop("Input file does not exist.")
sce <- readRDS(opt$input_sce_object)

#Model Variance
suppressPackageStartupMessages(require(scran))
var_table <- modelGeneVarWithSpikes(sce, block=opt$block, spikes=opt$spikes, size.factors=opt$size_factors, spike.size.factors=opt$spike_size_factors, 
                    assay.type=opt$assay_type, design=opt$design, 
                    subset.row=opt$subset_row, subset.fit=opt$subset_fill, equiweight=opt$equiweight, method=opt$method,
                    parametric=opt$parametric, min.mean=opt$min_mean)

#save table containing modelled variance and its biological and technical estimates
write.table(var_table, opt$output_geneVarSpikes_table, sep="\t")