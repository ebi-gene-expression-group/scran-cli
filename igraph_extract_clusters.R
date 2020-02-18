#!/usr/bin/env Rscript 

#Extract clustering annotation from igraph class object.

suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(workflowscriptscommon))

# argument parsing 
option_list = list(
  make_option(
    c("-i", "--input-igraph-object"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Path to the input igraph object in rds format.'
  ),
  make_option(
    c("-s", "--input-sce-object"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Path to the input SCE object where to add the cluster annotation extracted from the igraph objecti.'
  ),
  make_option(
    c("-o", "--output-sce-object"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Path to the output SCE object in rds format with cluster annotation in $cluster.'
 )
) 

opt = wsc_parse_args(option_list, mandatory = c("input_igraph_object", "input_sce_object", "output_sce_object"))

#read igraph object
if(!file.exists(opt$input_igraph_object)) stop("Input igraph file does not exist.")
graph <- readRDS(opt$input_igraph_object)
#read SCE object
if(!file.exists(opt$input_sce_object)) stop("Input SCE file does not exist.")
sce <- readRDS(opt$input_sce_object)

suppressPackageStartupMessages(require(igraph))
clust_annot <- cluster_walktrap(graph)$membership
print(clust_annot)
sce$cluster <- factor(clust_annot)

#save sce object with assigned clusters as RDS.
saveRDS(sce, opt$output_sce_object)
