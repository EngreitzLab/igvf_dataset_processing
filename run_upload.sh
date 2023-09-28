#!/bin/bash
set -e

RESULTS_DIR=/oak/stanford/groups/engreitz/Projects/IGVF-Y2AVE/outputs/
IGVF_FOLDER_NAME=igvf
THRESHOLD=.09
SYNAPSE_UPLOAD_FOLDER=syn52408868

python convert_to_igvf_format.py --results_dir $RESULTS_DIR --igvf_folder_name $IGVF_FOLDER_NAME --threshold $THRESHOLD

python upload_predictions.py --results_dir $RESULTS_DIR --igvf_folder_name $IGVF_FOLDER_NAME --project_synapse_id $SYNAPSE_UPLOAD_FOLDER --threshold $THRESHOLD
