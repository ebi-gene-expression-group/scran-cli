#!/usr/bin/env bash 

#REMOVE BEFORE PUSHING
#export PATH=$(pwd):$PATH

script_name=$0 

# This is a test script designed to test that everything works in the various
# accessory scripts in this package. Parameters used have absolutely NO
# relation to best practice and this should not be taken as a sensible
# parameterisation for a workflow.

function usage {
    echo "usage: scran_post_install_tests.sh [action] [use_existing_outputs]"
    echo "  - action: what action to take, 'test' or 'clean'"
    echo "  - use_existing_outputs, 'true' or 'false'"
    exit 1
}

action=${1:-'test'}
use_existing_outputs=${2:-'false'}

if [ "$action" != 'test' ] && [ "$action" != 'clean' ]; then
    echo "Invalid action"
    usage
fi

if [ "$use_existing_outputs" != 'true' ] && [ "$use_existing_outputs" != 'false' ]; then
    echo "Invalid value ($use_existing_outputs) for 'use_existing_outputs'"
    usage
fi

test_working_dir=`pwd`/'post_install_tests'
output_dir=$test_working_dir/outputs

# Clean up if specified
if [ "$action" = 'clean' ]; then
    echo "Cleaning up $test_working_dir ..."
    rm -rf $test_working_dir
    exit 0
elif [ "$action" != 'test' ]; then
    echo "Invalid action '$action' supplied"
    exit 1
fi

# Initialise directories
mkdir -p $test_working_dir
mkdir -p $output_dir

################################################################################
# List tool outputs/inputs & parameters 
################################################################################
#get test data
export accession_code='E-MTAB-6077'
export expr_data_type='filtered'
export normalisation_method='CPM'
export data_download_dir=$test_working_dir
export get_sdrf='FALSE'
export get_marker_genes='FALSE'
#read 10X data
export input_data_dir=$test_working_dir/'10x_data'
export sce_object=$output_dir/'output_10X.rds'
#sub sample SCE
export sub_sce=$output_dir/'sub_sce.rds'
export log_transform="TRUE"
export filter_min_genes=200 
export filter_min_cells=80 #arbitrary value
export n_spikes=25 #arbitrary value
#compute Sum Factors
export counts_factors_assay='logcounts'
export sce_factors=$test_working_dir/'sce_counts_factors.rds'
#compute Spike Factors
export spike_factors_assay='logcounts'
export sce_factors_spike=$test_working_dir/'sce_spike_factors.rds'
#scran - trendVar
export variance_trend=$test_working_dir/'variance_trend.rds'
#scran - denoise PCA
export sce_denoise_pca=$test_working_dir/'sce_denoise_pca.rds'
#get clustered PCA [in SCRAN 1.14]
#export cluster_PC_sce=$test_working_dir/'cluster_PCs_sce.rds'
#biuldSNNGraph
export shared_nn_graph="TRUE"
export graph_assay="logcounts"
export dim_red_NN="PCA_sub"
export igraph_object=$test_working_dir/'igraph_object.rds'
#extract clusters from igraph
export sce_clusters=$test_working_dir/'sce_clusters.rds'
#FindMarkers
export cluster_groups="cluster"
export markers_list=$test_working_dir/'markers.rds'
#correlated pairs of genes
export corr_gene_pairs=$test_working_dir/'cor_gene_pairs.rds'
#correlated genes
export corr_genes=$test_working_dir/'cor_genes.rds'
#convert SCE to other format
export convert_to="edgeR"
export converted_object=$test_working_dir/'converted_object.rds'

################################################################################
# Test individual scripts
################################################################################

# Make the script options available to the tests so we can skip tests e.g.
# where one of a chain has completed successfullly.

export use_existing_outputs

# Derive the tests file name from the script name

tests_file="${script_name%.*}".bats

# Execute the bats tests
$tests_file
