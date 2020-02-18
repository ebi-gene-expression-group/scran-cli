#!/usr/bin/env Rscript 

#Function available in SCRAN 1.12, substituted in SCRAN 1.14 by fitTrendVar.
#Fit a mean-dependent trend to the variances of the log-normalized expression values derived from count data.

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
    c("-m", "--min-mean"),
    action = "store",
    default = 0.1,
    type = 'numeric',
    help = 'A numeric scalar specifying the minimum mean to use for trend fitting.'
  ),
   make_option(
    c("-t", "--method"),
    action = "store",
    default = "loess",
    type = 'character',
    help = 'A string specifying the algorithm to use for smooth trend fitting.'
  ),
  make_option(
    c("-p", "--parametric"),
    action = "store",
    default = TRUE,
    type = 'logical',
    help = 'A logical scalar indicating whether a parametric fit should be attempted.'
  ),
 make_option(
    c("-a", "--assay-type"),
    action = "store",
    default = "logcounts",
    type = 'character',
    help = 'String or integer scalar specifying the assay containing the log-expression values. Default: "logcounts"'
  ),
  make_option(
    c("-k", "--use-spikes"),
    action = "store",
    default = FALSE,
    type = 'logical',
    help = 'If FALSE only the rows not labelled as spike-in transcripts will be used. If TRUE, nly rows labelled as spike-ins with isSpike(x) will be used.'
),
  make_option(
    c("-o", "--output-trend-var"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Path to the RDS object with named list containing: mean, var, resid.df, block, design, trend, df2.'
  )
)

opt = wsc_parse_args(option_list, mandatory = c("input_sce_object", "output_trend_var"))

#read SCE object
if(!file.exists(opt$input_sce_object)) stop("Input file does not exist.")
sce <- readRDS(opt$input_sce_object)

#Compute trendVar
suppressPackageStartupMessages(require(scran))
fit_trend_var <-trendVar(sce, method=opt$method, min.mean=opt$min_mean, assay.type=opt$assay_type, parametric=opt$parametric, use.spikes=opt$use_spikes) 

#save RDS object
saveRDS(fit_trend_var, opt$output_trend_var)
