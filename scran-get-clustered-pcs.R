#!/usr/bin/env Rscript 
#This function is available on SCRAN 1.14
#This function inputs a low-dimensional embedding and adjusts its n_PCs so that they are not less than the number of subpopulations (which are unknown, of course, so we use the number of clusters as a proxy). 

#TODO: Include more detailed documentation/indications. 

suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(workflowscriptscommon))

# argument parsing 
option_list = list(
  make_option(
    c("-i", "--input-sce-object"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Path to the input SCE object including a dimensionality reduction, in rds format.'
  ),
  make_option(
    c("-c", "--function"),
    action = "store",
    default = NULL,
    type = 'function',
    help = 'A clustering function that takes a numeric matrix with rows as cells and returns a vector containing a cluster label for each cell.
    The default is a graph-based clustering method using buildSNNGraph and cluster_walktrap, where arguments in ... are passed to the former.'
  ),
  make_option(
    c("-m", "--min-rank"),
    action = "store",
    default = 5,
    type = 'integer',
    help = 'Integer scalars specifying the minimum number of PCs to retain.'
  ),
    make_option(
    c("-m", "--max-rank"),
    action = "store",
    default = 50,
    type = 'integer',
    help = 'Integer scalars specifying the maximum number of PCs to retain.'
  ),
  make_option(
    c("-b", "--by"),
    action = "store",
    default =  1
    type = 'integer',
    help = 'Integer scalar specifying what intervals should be tested between min.rank and max.rank.'
  ),
  make_option(
    c("-o", "--output-sce-object"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Path to the output SCE object with denoised PC embedding'
  )
)

opt = wsc_parse_args(option_list, mandatory = c("input_sce_object", "output_sce_object"))

#read SCE object
if(!file.exists(opt$input_sce_object)) stop("Input file does not exist.")
sce <- readRDS(opt$input_sce_object)

#check if SCE has PCA embedding
if(length(reducedDimNames(sce)) == 0){ stop("Input SCE provided does not have a dimensionally reduced embedding")}

#denoise PCA embedding based on the number of subpopulations from clustering
suppressPackageStartupMessages(require(scran))
output <- getClusteredPCs(pcs=opt$pcs, FUN=opt$function, min.rank=opt$min_rank, max.rank=opt$max_rank, by=opt$by)
#number of PCs selected based on clustering 
npcs <- metadata(output)$chosen
#subset original embedding
reducedDim(sce, "PCAsub") <- reducedDim(sce, "PCA")[,1:npcs,drop=FALSE]

#save SCE object with denoised PCA components 
saveRDS(sce, opt$output_sce_object)
