#!/bin/bash
set -e

METADATA_FILE=Y2AVE_SingleCellDatasets.CellClusterTable.tsv
METADATA_SYNAPSE_ID=syn52252345 
DATASET_DIR=/oak/stanford/projects/igvf/Y2AVE/E2G_Predictions/inputs

# Fetch latest version of metadata file
rm -f $METADATA_FILE
synapse get $METADATA_SYNAPSE_ID

mkdir -p $DATASET_DIR

python download_datasets.py --metadata_file $METADATA_FILE --dataset_dir $DATASET_DIR
python generate_abc_biosamples.py --dataset_dir $DATASET_DIR
python qc_datasets.py --dataset_dir $DATASET_DIR
