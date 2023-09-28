import os
import shutil
import tempfile

import click
import pandas as pd
import synapseclient
from synapseclient import File, Folder

from download_datasets import get_relevant_datasets

# Delete datasets that we've processed
# In this case, we delete datasets with num fragments < 1e6

@click.command()
@click.option("--metadata_file", type=str, default="Y2AVE_SingleCellDatasets.CellClusterTable.tsv")
@click.option("--input_dir", type=str, default="/oak/stanford/projects/igvf/Y2AVE/E2G_Predictions/inputs")
@click.option("--results_dir", type=str, default="/oak/stanford/groups/engreitz/Projects/IGVF-Y2AVE/outputs/")
@click.option("--min_num_fragments", type=int, default=1e6)
@click.option("--dataset_summary_syn_id", type=str, default="syn52576265")
def main(metadata_file, input_dir, results_dir, min_num_fragments, dataset_summary_syn_id):
    metadata_df = pd.read_csv(metadata_file, sep="\t")
    relevant_datasets = get_relevant_datasets(metadata_df, min_num_fragments)

    # This may include clusters we may not have proccessed 
    clusters_to_delete = set(metadata_df["CellClusterID"]) - set(relevant_datasets["CellClusterID"])
    print(f"Found {len(clusters_to_delete)} datasets to delete")

    for cluster in clusters_to_delete:
        for dir in [input_dir, results_dir]:
            cluster_dir = os.path.join(dir, cluster)
            try:
                shutil.rmtree(cluster_dir)
                print(f"Deleted directory: {cluster_dir}")
            except FileNotFoundError:
                # Already deleted
                pass

    print("Deleted all the directories")

    print("Deleting from synapse")
    syn = synapseclient.login()
    dataset_summary_syn = syn.get(dataset_summary_syn_id)
    # summary_df = pd.read_csv(dataset_summary_syn.path, sep="\t")
    summary_df = pd.read_csv("/home/users/atan5133/.synapseCache/311/128527311/tmpan5pxh54.tsv", sep="\t")
    for cluster in clusters_to_delete:
        if (summary_df["CellClusterID"] == cluster).any():
            cluster_folder_syn_id = summary_df[summary_df["CellClusterID"] == cluster].iloc[0]["CellClusterIDSynapseID"]
            try:
                # seems we have to get the cluster from syn in order to delete it
                cluster_syn = syn.get(cluster_folder_syn_id)
                syn.delete(cluster_syn)
                print(f"Deleted cluster from synapse: {cluster}")
            except Exception:
                # cluster folder has already been deleted
                pass

    summary_df = summary_df[~summary_df["CellClusterID"].isin(clusters_to_delete)]

    with tempfile.NamedTemporaryFile(suffix=".tsv") as temp_file:
        summary_df.to_csv(temp_file, sep="\t", index=False)
        project_folder = syn.get(dataset_summary_syn.parentId)
        summary_syn = syn.store(File(temp_file.name, name="DatasetSummary.tsv", parent=project_folder))    
    
    print(f"Uploaded new data summary file! {summary_syn.id}")

if __name__ == "__main__":
    main()