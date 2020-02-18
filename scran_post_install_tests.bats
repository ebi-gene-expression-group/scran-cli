#!/usr/bin/env bats 

# download test sce object from the link provided in package docs
@test "get experiment data" {
    if [ "$use_existing_outputs" = 'true' ] ; then
        skip "exists and use_existing_outputs is set to 'true'"
    fi

    run get_experiment_data.R --accesssion-code $accession_code --expr-data-type $expr_data_type --normalisation-method $normalisation_method --get-sdrf $get_sdrf --get-marker-genes $get_marker_genes --output-dir-name $data_download_dir --exp-data-dir '10x_data'

    echo "status = ${status}" #exit status
    echo "output = ${output}"

    [ "$status" -eq 0 ] 
}

#read downloaded data
@test "read data into serialized SCE object" {
    if [ "$use_existing_outputs" = 'true' ] && [ -f "$sce_object" ]; then
        skip "$sce_object exists and use_existing_outputs is set to 'true'"
    fi

    run rm -f $sce_object && dropletutils-read-10x-counts.R --samples $input_data_dir --output-object-file $sce_object
    echo "status = ${status}" #exit status
    echo "output = ${output}"

    [ "$status" -eq 0 ] 
    [ -f  "$sce_object" ] 

}
#subsample SCE object
@test "subsample SCE object" {
    if [ "$use_existing_outputs" = 'true' ] && [ -f "$sub_sce" ]; then
        skip "$sub_sce exists and use_existing_outputs is set to 'true'"
    fi

    run rm -f $sub_sce &&\
			scran-cli-filter-sce.R\
				--input-sce-object $sce_object\
				--min-genes $filter_min_genes\
				--min-cells $filter_min_cells\
				--n-spikes $n_spikes\
				--log $log_transform\
				--output-sce-object $sub_sce
    echo "status = ${status}" #exit status
    echo "output = ${output}"

    [ "$status" -eq 0 ] 
    [ -f  "$sub_sce" ] 

}

#scran compute sum factors
@test "compute counts factors" {
    if [ "$use_existing_outputs" = 'true' ] && [ -f "$sce_factors" ]; then
        skip "$sce_factors exists and use_existing_outputs is set to 'true'"
    fi

    run rm -f $sce_factors &&\
                        run scran-compute-sum-factors.R\
                            --input-sce-object $sub_sce\
                            --assay-type $counts_factors_assay\
                            --output-sce-object $sce_factors

    echo "status = ${status}" #exit status
    echo "output = ${output}"

    [ "$status" -eq 0 ] 
    [ -f  "$sce_factors" ] 
}

#scran compute spike factors
@test "compute spike-in factors" {
    if [ "$use_existing_outputs" = 'true' ] && [ -f "$sce_factors_spike" ]; then
        skip "$sce_factors_spike exists and use_existing_outputs is set to 'true'"
    fi

    run rm -f $sce_factors_spike &&\
                        scran-compute-spike-factors.R\
                            --input-sce-object $sub_sce\
                            --assay-type $spike_factors_assay\
                            --output-sce-object $sce_factors_spike

    echo "status = ${status}" #exit status
    echo "output = ${output}"

    [ "$status" -eq 0 ] 
    [ -f  "$sce_factors_spike" ] 
}

#model gene variance [SCRAN 1.14]
#@test "model gene variance" {
#    if [ "$use_existing_outputs" = 'true' ] && [ -f "$GeneVar_table" ]; then
#        skip "$GeneVar_table exists and use_existing_outputs is set to 'true'"
#    fi
#
#    run rm -f $GeneVar_table &&\
#                       scran-model-gene-var.R\
#                            --input-sce-object $sub_sce\
#                            --output-geneVarSpikes-table $GeneVar_table
#
#    echo "status = ${status}" #exit status
#    echo "output = ${output}"
#
#    [ "$status" -eq 0 ] 
#    [ -f  "$GeneVar_table" ] 
#}


#model gene variance With Spikes [SCRAN 1.14]
#@test "model gene variance with spikes" {
#    if [ "$use_existing_outputs" = 'true' ] && [ -f "$GeneVarSpikes_table" ]; then
#        skip "$GeneVarSpikes_table exists and use_existing_outputs is set to 'true'"
#    fi
#
#    run rm -f $GeneVarSpikes_table &&\
#                        scran-model-gene-var-with-spikes.R\
#                            --input-sce-object $lognorm_sce\
#                            --output-geneVarSpikes-table $GeneVarSpikes_table
#
#    echo "status = ${status}" #exit status
#    echo "output = ${output}"
#
#    [ "$status" -eq 0 ] 
#    [ -f  "$GeneVarSpikes_table" ] 
#}

#trendVar - [only SCRAN 1.12, deprecated in SCRAN 1.14]
@test "Fit a mean-dependent trend to the gene-specific variances in single-cell RNA-seq data" {
    if [ "$use_existing_outputs" = 'true' ] && [ -f "$variance_trend" ]; then
        skip "$variance_trend and use_existing_outputs is set to 'true'"
    fi

    run rm -f $variance_trend &&\
                        scran-trend-var.R\
                            --input-sce-object $sub_sce\
			    --assay-type $counts_factors_assay\
			    --output-trend-var $variance_trend 
    
    echo "status = ${status}" #exit status
    echo "output = ${output}"

    [ "$status" -eq 0 ] 
    [ -f  "$variance_trend" ] 
}

#denoise PCA [requires 'technical'argument provided by modelVar, modelVarWithSpikes and fitTrendVar, all from Scran v1.14]
@test "denoise PCA" {
    if [ "$use_existing_outputs" = 'true' ] && [ -f "$sce_denoise_pca" ]; then
        skip "$sce_denoise_pca exists and use_existing_outputs is set to 'true'"
    fi

    run rm -f $sce_denoise_pca &&\
                        scran-denoise-pca.R\
                            --input-sce-object $sub_sce\
			    --assay-type $counts_factors_assay\
                            --technical $variance_trend\
                            --output-sce-object $sce_denoise_pca

    echo "status = ${status}" #exit status
    echo "output = ${output}"

    [ "$status" -eq 0 ] 
    [ -f  "$sce_denoise_pca" ] 
}

#get clustered PCs [SCRAN 1.14]
#@test "get clustered PCs" {
#    if [ "$use_existing_outputs" = 'true' ] && [ -f "$cluster_PC_sce" ]; then
#        skip "$cluster_PC_sce exists and use_existing_outputs is set to 'true'"
#    fi
#
#    run rm -f $cluster_PC_sce &&\
#                        scran-get-clustered-pcs.R\
#                            --input-sce-object $sce_denoise_pca\
#                            --output-sce-object $cluster_PC_sce
#
#    echo "status = ${status}" #exit status
#    echo "output = ${output}"
#
#    [ "$status" -eq 0 ] 
#    [ -f  "$cluster_PC_sce" ] 
#}

#buildSNNGraph
@test "build Nearest Neighbour Graph" {
    if [ "$use_existing_outputs" = 'true' ] && [ -f "$igraph_object" ]; then
        skip "$igraph_object exists and use_existing_outputs is set to 'true'"
    fi

    run rm -f $clusters_NN_sce &&\
			scran-build-snn-graph.R\
				--input-sce-object $sce_denoise_pca\
				--shared=$shared_nn_graph\
				--assay-type=$graph_assay\
				--output-igraph-object $igraph_object

    echo "status = ${status}" #exit status
    echo "output = ${output}"

    [ "$status" -eq 0 ] 
    [ -f  "$igraph_object" ] 
}
#extract clusters from igraph
@test "extract clusters from igraph" {
    if [ "$use_existing_outputs" = 'true' ] && [ -f "$sce_clusters" ]; then
        skip "$sce_clusters exists and use_existing_outputs is set to 'true'"
    fi

    run rm -f $sce_clusters &&\
			igraph_extract_clusters.R\
				--input-igraph-object $igraph_object\
                            	--input-sce-object $sce_denoise_pca\
                            	--output-sce-object $sce_clusters  

    echo "status = ${status}" #exit status
    echo "output = ${output}"

    [ "$status" -eq 0 ] 
    [ -f  "$sce_clusters" ] 
}

#Find Marker genes
@test "Find Marker genes" {
    if [ "$markers_list" = 'true' ] && [ -f "$markers_list" ]; then
        skip "$markers_list exists and use_existing_outputs is set to 'true'"
    fi

    run rm -f $markers_list &&\
                       scran-find-markers.R\
                            --input-sce-object $sce_clusters\
                            --cluster $cluster_groups\
                            --output-markers $markers_list

    echo "status = ${status}" #exit status
    echo "output = ${output}"

    [ "$status" -eq 0 ] 
    [ -f  "$markers_list" ] 
}

#Identify correlated pairs of genes
@test "Identify correlated pairs of genes" {
    if [ "$corr_gene_pairs" = 'true' ] && [ -f "$corr_gene_pairs" ]; then
        skip "$corr_gene_pairs exists and use_existing_outputs is set to 'true'"
    fi

    run rm -f $corr_gene_pairs &&\
                       scran-correlate-pairs.R\
                            --input-sce-object $sub_sce\
                            --output-pairs-df $corr_gene_pairs

    echo "status = ${status}" #exit status
    echo "output = ${output}"

#    [ "$status" -eq 0 ] 
    [ -f  "$corr_gene_pairs" ] 
}

#Identify correlated Genes
@test "Identify correlated Genes" {
    if [ "$corr_genes" = 'true' ] && [ -f "$corr_genes" ]; then
        skip "$corr_genes exists and use_existing_outputs is set to 'true'"
    fi

    run rm -f $corr_genes &&\
                       scran-correlate-genes.R\
                            --input-corr-pairs $corr_gene_pairs\
                            --output-corr-genes $corr_genes

    echo "status = ${status}" #exit status
    echo "output = ${output}"

    [ "$status" -eq 0 ] 
    [ -f  "$corr_genes" ] 
}
#Convert SCE to other formats 
@test "Convert SCE to other formats " {
    if [ "$converted_object" = 'true' ] && [ -f "$converted_object" ]; then
        skip "$converted_object exists and use_existing_outputs is set to 'true'"
    fi

    run rm -f $converted_object &&\
                        scran-convert-to.R\
                            --input-sce-object $sce_object\
                            --type $convert_to\
                            --output-converted $converted_object

    echo "status = ${status}" #exit status
    echo "output = ${output}"

    [ "$status" -eq 0 ] 
    [ -f  "$converted_object" ] 
}
