# Wrapper scripts for components of the scran package

Scran implements functions for low-level analyses of single-cell RNA-seq data. Methods are provided for normalization of cell-specific biases, assignment of cell cycle phase, detection of highly variable and significantly correlated genes, identification of marker genes, and other common tasks in routine single-cell analysis workflows.

A vignette with the usage of scran can be found [here](https://bioconductor.org/packages/release/bioc/vignettes/scran/inst/doc/scran.html) 

**Note**: the following are only present in Scran version 1.12: 
- scran-trend-var.R (*trendVar()* is deprecated in latest scran for *FitTrendVar()*)

The following scripts are only present in latest version of scran (1.14):
- scran-model-gene-var.R 
- scran-model-gene-var-with-spikes.R

## Install
The simplest way to install package is to use conda. It is recommended to install into clean environment to avoid dependency conflicts.
```
conda install scran-cli
```

## Test installation
To run post-install tests, you will need to install bats testing system into the environment you created:
```
conda install bats
```
Then, run scpred_post_install_tests.sh

## Commands

Currently available scripts are detailed below, each of which has usage instructions available via --help.

### scran-compute-sum-factors.R -- call computeSumFactors()

Scaling normalization of single-cell RNA-seq data by deconvolving size factors from cell pools. 

The assumption is that most genes are not differentially expressed (DE) between cells, 
such that any differences in expression across the majority of genes represents some technical bias that should be removed.

```
scran-compute-sum-factors.R -i <Path to the input SCE object in rds format> -a <Specify which assay to use> -o <Path to the output SCE object containing the vector of size factors in sizeFactors(x)>
```

### scran-compute-spike-factors.R -- call computeSpikeFactors()

Compute size factors based on the coverage of spike-in transcripts. 

Size factors are computed to scale the counts such that the total coverage of the spike-in transcripts is equal across cells. The main practical difference with computeSumFactors and other non-DE methods is that spike-in normalization preserves differences in total RNA content between cells.

```
scran-compute-sum-factors.R -i <Path to the input SCE object in rds format> -s<String or integer scalar specifying the alternative experiment containing the spike-in transcripts. Default:"ERCC"> -o <Path to the output SCE object containing the vector of size factors in sizeFactors(x)>
```
### scran-trend-var.R -- call trendVar()

Fit a mean-dependent trend to the gene-specific variances in single-cell RNA-seq data.

A function that returns the fitted value of the trend at any mean can be found in  output_list$trend, and serves as argument for denoisePCA.
```
scran-trend-var.R -i <Path to the input SCE object in rds format> -o <Path to the RDS object with named list containing: mean, var, resid.df, block, design, trend, df2.>
```
**Note**: trendVar() is deprecated in latest scran version 1.14, now substitueted by *modelGeneVar*,  *modelGeneVarWithSpikes*, or *FitTrendVar*. 

### scran-denoise-pca.R call denoisePCA()

Remove principal components corresponding to technical noise. 

Choice of the number of PCs to discard is based on the estimates of technical variance in technical, which may be computed by *modelGeneVar* or *modelGeneVarWithSpikes*, or deprecated *FitTrendVar*

```
scran-denoise-pca.R -i <Path to the input SCE object in rds format> -t <function that computes the technical component of the variance for a gene, or a vector with this values> -o <Path to the output SCE object with adjusted PC number>
```

### scran-build-snn-graph.R -- call buildSNNGraph()

Build a shared or k-nearest-neighbors graph for cells based on their expression profiles. 

Outputs an igraph class object. 

```
scran-build-snn-graph.R -i <Path to the input SCE object in rds format> -s <'Logical specifying wether to compute a Shared NN Graph (if shared=TRUE) or a kNN Graph(if shared=FALSE)> -o <Path to the output igraph object>
```

### igraph_extract_clusters.R 

Extract clustering annotation from igraph class object.

```
igraph_extract_clusters.R -i <Path to the input igraph object in rds format> -s <Path to the input SCE object where we want to add the cluster annotation from the igraph object> -o <Path to the output SCE object in rds format with cluster annotation as $cluster>
```


### scran-find-markers.R -- call findMarkers()

Find candidate marker genes for groups of cells (e.g., clusters) by testing for differential expression between pairs of groups. 

*Note*: If x is scale-normalized but not log-transformed, it can be used with test.type="wilcox" and test.type="binom". If x contains raw counts, it can only be used with test.type="binom" [available on scran 1.14].

```
scran-find-markers.R -i <Path to the input SCE object in rds format>  -c <A vector of group assignments for all cells> -o <Path to the rds  list of DataFrames with a sorted marker gene list per cluster/group>
```

### scran-correlate-pairs.R -- call correlatePairs()

Identify pairs of genes that are significantly correlated in their expression profiles, based on Spearman's rank correlation. 
```
 scran-correlate-pairs.R -i <Path to the input SCE object in rds format>  -o <Path to the output dataframe with one row per gene pair (rows order by increasing p-values) in txt format.>
```

### scran-correlate-genes.R - call correlateGenes()

Compute per-gene correlation statistics by combining results from gene pair correlations. 
This provides compute a single set of statistics for each gene, rather than for each pair (correlatePairs).

```
scran-correlate-genes.R -i <Path to the  DataFrame of pairwise correlation statistics, returned by correlatePairs> -o <A DataFrame with one row per unique gene in stats and containing the fields: gene, rho, p.value, FDR, limited>
```

### scran-convert-to.R -- call convertTo()

Convert a SCESet object into other classes for entry into other analysis pipelines (edgeR, DESeq2, monocle).  

```
scran-convert-to.R -i <Path to the input SCE object in rds format> -o <A string specifying the analysis for which the object should be prepared. Any of: "edgeR", "DESeq2", "monocle".>
```


