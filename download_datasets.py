import concurrent.futures
import os
import shutil

import click
import pandas as pd
import synapseclient

NUM_THREADS = 20


def tagalign_files_exist(dir) -> bool:
    files_in_directory = os.listdir(dir)
    gz_file_exists = any(filename.endswith(".gz") for filename in files_in_directory)
    gz_tbi_file_exists = any(
        filename.endswith(".gz.tbi") for filename in files_in_directory
    )
    return gz_file_exists and gz_tbi_file_exists


def download_cell_cluster_info(syn: synapseclient.Synapse, dataset_dir, row) -> None:
    cluster = row["CellClusterID"]
    cluster_dir = os.path.join(dataset_dir, cluster)
    if os.path.exists(cluster_dir):
        if tagalign_files_exist(cluster_dir):
            return
        else:
            # redownload all files
            print(f"Incomplete download. Removing existing dir: {cluster_dir}")
            shutil.rmtree(cluster_dir)
    os.makedirs(cluster_dir)

    print(f"Downloading cell cluster: {cluster}")
    sorted_tagAlign_synapse = row["ATACtagAlignSorted"]
    tagAlign_idx_synapse = row["ATACtagAlignSortedIndex"]
    syn.get(sorted_tagAlign_synapse, downloadLocation=cluster_dir)
    if isinstance(tagAlign_idx_synapse, str) and tagAlign_idx_synapse.startswith("syn"):
        # Some rows aren't indexed yet, but we still want to process them
        syn.get(tagAlign_idx_synapse, downloadLocation=cluster_dir)


def get_relevant_datasets(metadata_df: pd.DataFrame) -> pd.DataFrame:
    return metadata_df[
        (~metadata_df["ATACtagAlignSorted"].isna())
        & (~metadata_df["ATACtagAlignSortedIndex"].isna())
        & (metadata_df["Species"] == "Human")
    ]


@click.command()
@click.option("--metadata_file")
@click.option("--dataset_dir")
def main(metadata_file, dataset_dir):
    metadata_df = pd.read_csv(metadata_file, sep="\t")

    relevant_datasets = get_relevant_datasets(metadata_df)

    syn = synapseclient.Synapse()
    syn.login()
    print(f"{len(relevant_datasets)} datasets to download")
    print(f"Using parallelization with {NUM_THREADS} threads")

    # Download datasets in parallel
    with concurrent.futures.ThreadPoolExecutor(max_workers=NUM_THREADS) as executor:
        futures = {
            executor.submit(download_cell_cluster_info, syn, dataset_dir, row)
            for _, row in relevant_datasets.iterrows()
        }

        for future in concurrent.futures.as_completed(futures):
            future.result()  # raises Exceptions from threads

        concurrent.futures.wait(futures)

    print("Downloaded all the datasets")


if __name__ == "__main__":
    main()
