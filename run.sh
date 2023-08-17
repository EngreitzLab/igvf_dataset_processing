#!/bin/bash

METADATA_FILE=Y2AVE_SingleCellDatasets.CellClusterTable.tsv
METADATA_SYNAPSE_ID=syn52252345 
DATASET_DIR=datasets

# Fetch latest version of metadata file
# rm $METADATA_FILE
# synapse get $METADATA_SYNAPSE_ID

mkdir -p $DATASET_DIR

python download_datasets.py --metadata_file $METADATA_FILE --dataset_dir $DATASET_DIR
