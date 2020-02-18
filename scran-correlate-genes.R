#!/usr/bin/env Rscript 

#Compute per-gene correlation statistics by combining results from gene pair correlations. 
#This provides compute a single set of statistics for each gene, rather than for each pair (correlatePairs).

suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(workflowscriptscommon))

# argument parsing 
option_list = list(
  make_option(
    c("-i", "--input-corr-pairs"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Path to the  DataFrame of pairwise correlation statistics, returned by correlatePairs.'
  ),
  make_option(
    c("-o", "--output-corr-genes"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'A DataFrame with one row per unique gene in stats and containing the fields: gene, rho, p.value, FDR, limited.'
  )
)

opt = wsc_parse_args(option_list, mandatory = c("input_corr_pairs", "output_corr_genes"))

#read SCE object
if(!file.exists(opt$input_corr_pairs)) stop("Input file does not exist.")
corr_pairs <- read.table(opt$input_corr_pairs, header = T, sep="\t")

#Compute PCA and denoise it
suppressPackageStartupMessages(require(scran))

#Compute gene pair-correlation
corr_genes <- correlateGenes(stats=corr_pairs)

#save gene-pairs table
write.table(corr_genes, file = opt$output_corr_genes, sep = "\t")
