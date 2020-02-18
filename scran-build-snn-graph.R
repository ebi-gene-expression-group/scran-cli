#!/usr/bin/env Rscript 

#Build a shared or k-nearest-neighbors graph for cells based on their expression profiles. 
#By default, buildSNNGraph() uses the mode of shared neighbor weighting described by Xu and Su (2015), but other weighting methods (e.g., the Jaccard index) are also available by setting type=X. 
#Outputs an igraph class object. 

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
    c("-s", "--shared"),
    action = "store",
    default = TRUE,
    type = 'logical',
    help = 'Logical specifying wether to compute a Shared NN Graph (if shared=TRUE) or a kNN Graph(if shared=FALSE).'
  ),  
  make_option(
    c("-k", "--k-value"),
    action = "store",
    default = 10,
    type = 'integer',
    help = 'An integer scalar specifying the number of nearest neighbors to consider during graph construction.'
  ),
    make_option(
    c("-d", "--d-value"),
    action = "store",
    default = 50,
    type = 'integer',
    help = 'An integer scalar specifying the number of dimensions to use for the search.'
  ),
  make_option(
    c("-y", "--type"),
    action = "store",
    default = NULL,
    type = 'character',
    help = 'A string specifying the type of weighting scheme to use for shared neighbors.'
  ),
  make_option(
    c("-t", "--transposed"),
    action = "store",
    default = FALSE,
    type = 'logical',
    help = 'A logical scalar indicating whether x is transposed (i.e., rows are cells).'
  ),
  make_option(
    c("-r", "--subset_row"),
    action = "store",
    default = NULL,
    type = 'logical',
    help = 'Logical, integer or character vector specifying the rows for which to model the variance. Defaults to all genes in x.'
  ),
  make_option(
    c("-a", "--assay-type"),
    action = "store",
    default = "logcounts",
    type = 'character',
    help = 'A string specifying which assay values to use. Default: logcounts.'
  ),
  make_option(
    c("-g", "--get-spikes"),
    action = "store",
    default = FALSE,
    type = 'logical',
    help = 'Logical specifying wether to perform spike-ins filtering(FALSE) or not (TRUE).'
  ),
  make_option(
    c("-u", "--use-dimred"),
    action = "store",
    default = NULL,
    type = 'character',
    help = 'A string specifying whether existing values in reducedDims(x) should be used.'
  ),
  make_option(
    c("-o", "--output-igraph-object"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Path to the output igraph object in RDS format.'
  )
)

opt = wsc_parse_args(option_list, mandatory = c("input_sce_object", "shared", "output_igraph_object"))

#read SCE object
if(!file.exists(opt$input_sce_object)) stop("Input file does not exist.")
sce <- readRDS(opt$input_sce_object)

suppressPackageStartupMessages(require(scran))

#Shared NN Graph
if(opt$shared == T){
    print("--- Computing Shared Nearest Neighbour Graph ---")
    graph <- buildSNNGraph(sce, subset.row=opt$subset_row, k=opt$k, d=opt$d, type=opt$type, assay.type=opt$assay_type, get.spikes=opt$get_spikes, use.dimred=opt$use_dimred)
}
#KNN Graph
if(opt$shared == F){
    print("--- Computing K Nearest Neighbour Graph ---")
    graph <- BuildKNNGraph(sce, subset.row=opt$subset_row, k=opt$k, d=opt$d, assay.type=opt$assay_type, get.spikes=opt$get_spikes, use.dimred=opt$use_dimred)
}

#save sce object 
saveRDS(graph, opt$output_igraph_object)
