import concurrent.futures
import os
import glob
import tempfile
import click
import pandas as pd
import synapseclient
from synapseclient import Folder, File


NUM_THREADS = 20
MACS_FILE = "macs2_peaks.narrowPeak.sorted"
ENHANCERLIST_FILE = "EnhancerList.txt"
GENELIST_FILE = "GeneList.txt"
E2G_PREDICTIONS_FILE = "encode_e2g_predictions.tsv.gz"
E2G_PREDICTIONS_THRESHOLDED_FILE = "encode_e2g_predictions_threshold{threshold}.tsv.gz"

def upload_folder(syn_client, cluster, folder, parent, threshold):
    cluster_folder = Folder(cluster, parent=parent)
    cluster_folder = syn_client.store(cluster_folder)
    macs_file = File(os.path.join(folder, MACS_FILE), parent=cluster_folder)
    enh_list_file = File(os.path.join(folder, ENHANCERLIST_FILE), parent=cluster_folder)
    gene_list_file = File(os.path.join(folder, GENELIST_FILE), parent=cluster_folder)
    pred_file = File(os.path.join(folder, E2G_PREDICTIONS_FILE), parent=cluster_folder)
    pred_threshold_file = File(os.path.join(folder, E2G_PREDICTIONS_THRESHOLDED_FILE.format(threshold=threshold)), parent=cluster_folder)

    macs_syn_id = syn_client.store(macs_file).id
    enh_list_syn_id = syn_client.store(enh_list_file).id
    gene_list_syn_id = syn_client.store(gene_list_file).id
    pred_syn_id = syn_client.store(pred_file).id
    pred_threshold_syn_id = syn_client.store(pred_threshold_file).id
    return {"CellClusterID": cluster, "CellClusterIDSynapseID": cluster_folder.id, "PeaksMACS2": macs_syn_id, "EnhancerList": enh_list_syn_id, "GeneList": gene_list_syn_id, "ENCODE-rE2G": pred_syn_id, "ENCODE-rE2G_Thresholded": pred_threshold_syn_id}
    

@click.command()
@click.option("--results_dir", type=str, required=True)
@click.option("--igvf_folder_name", type=str, required=True)
@click.option("--project_synapse_id", type=str, required=True)
@click.option("--threshold", type=float, required=True)
def main(results_dir, igvf_folder_name, project_synapse_id, threshold):
    prediction_files = glob.glob(os.path.join(results_dir, '*', igvf_folder_name, E2G_PREDICTIONS_FILE))
    cell_clusters = [file.split("/")[-3] for file in prediction_files]
    syn = synapseclient.login()
    project_folder = syn.get(project_synapse_id)
    print(f"Using parallelization with {NUM_THREADS} threads")

    with concurrent.futures.ThreadPoolExecutor(max_workers=NUM_THREADS) as executor:
        futures = set()
        for cluster in cell_clusters:
            igvf_upload_folder = os.path.join(results_dir, cluster, igvf_folder_name)
            futures.add(executor.submit(upload_folder, syn, cluster, igvf_upload_folder, project_folder, threshold))

        rows = []
        for future in concurrent.futures.as_completed(futures):
            rows.append(future.result())  # raises Exceptions from threads

    with tempfile.NamedTemporaryFile(suffix=".tsv") as temp_file:
        pd.DataFrame(rows).to_csv(temp_file, sep="\t", index=False)
        syn.store(File(temp_file.name, name="Metadata.tsv", parent=project_folder))
    print("Uploaded all the datasets")

if __name__ == "__main__":
    main()
