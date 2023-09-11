import os
import glob
import click
import pandas as pd
import gzip
import time

MACS_FILE = "Peaks/macs2_peaks.narrowPeak.sorted"
ENHANCERLIST_FILE = "Neighborhoods/EnhancerList.txt"
GENELIST_FILE = "Neighborhoods/GeneList.txt"
E2G_PREDICTIONS_FILE = "encode_e2g_predictions.tsv.gz"
E2G_PREDICTIONS_THRESHOLDED_FILE = "encode_e2g_predictions_threshold{threshold}.tsv.gz"

MACS_COLUMNS = [
    "chr",
    "start",
    "end",
    "name",
    "score",
    "strand",
    "signalValue",
    "pValue",
    "qValue",
    "peak",
]

ELEMENT_NODE_COLUMNS = [
    "chr",
    "start",
    "end",
    "name",
    "CellClusterID",
    "BiosampleOntologyTerm",
]

GENE_NODE_COLUMNS = [
    "Symbol",
    "EnsemblID",
    "TSS",
    "CellClusterID",
    "BiosampleOntologyTerm",
]

ELEMENT_GENE_EDGE_COLUMNS = [
    "Element:chr",
    "Element:start",
    "Element:end",
    "Element:name",
    "Element:strand",
    "Gene:Symbol",
    "Gene:EnsemblID",
    "Gene:TSS",
    "Score",
    "CellClusterID",
    "BiosampleOntologyTerm",
]


def write_df_and_comments(df, comment_lines, output_file, compressed=False):
    if compressed:
        with gzip.open(output_file, mode="wt", encoding="utf-8") as gz_file:
            for comment in comment_lines:
                gz_file.write(comment + "\n")
            df.to_csv(
                gz_file, sep="\t", index=False, header=True, mode="a",
            )
    else:
        with open(output_file, "w") as file:
            for comment in comment_lines:
                file.write(comment + "\n")
            df.to_csv(
                file, sep="\t", index=False, header=True, mode="a",
            )


def convert_macs(cluster, macs_file, upload_folder):
    df = pd.read_csv(macs_file, sep="\t", header=None, names=MACS_COLUMNS)
    df["CellClusterID"] = cluster
    df["BiosampleOntologyTerm"] = ""
    additional_cols = [col for col in df.columns if col not in ELEMENT_NODE_COLUMNS]

    comment_lines = [
        "#FileType=Node:RegulatoryElement",
        "#Code=https://github.com/EngreitzLab/e2g_pipeline",
        "#Contact=Jesse Engreitz (engreitz@stanford.edu)",
        "#Genome=GRCh38",
        "#Description=MACS2 narrow peaks",
        "#BiosampleAgnostic=false",
    ]

    output_file = os.path.join(upload_folder, os.path.basename(macs_file))
    write_df_and_comments(
        df[ELEMENT_NODE_COLUMNS + additional_cols], comment_lines, output_file
    )


def convert_enhancerlist(cluster, enhancerlist_file, upload_folder):
    df = pd.read_csv(enhancerlist_file, sep="\t")
    df["CellClusterID"] = cluster
    df["BiosampleOntologyTerm"] = ""
    additional_cols = [col for col in df.columns if col not in ELEMENT_NODE_COLUMNS]

    comment_lines = [
        "#FileType=Node:RegulatoryElement",
        "#Code=https://github.com/EngreitzLab/e2g_pipeline",
        "#Contact=Jesse Engreitz (engreitz@stanford.edu)",
        "#Genome=GRCh38",
        "#Description=Candidate Enhancers",
        "#BiosampleAgnostic=false",
    ]
    output_file = os.path.join(upload_folder, os.path.basename(enhancerlist_file))
    write_df_and_comments(
        df[ELEMENT_NODE_COLUMNS + additional_cols], comment_lines, output_file
    )


def convert_genelist(cluster, genelist_file, upload_folder):
    df = pd.read_csv(genelist_file, sep="\t")
    column_name_mapping = {"symbol": "Symbol", "Ensembl_ID": "EnsemblID", "tss": "TSS"}
    df.rename(columns=column_name_mapping, inplace=True)
    df["CellClusterID"] = cluster
    df["BiosampleOntologyTerm"] = ""
    additional_cols = [col for col in df.columns if col not in GENE_NODE_COLUMNS]

    comment_lines = [
        "#FileType=Node:Gene",
        "#Code=https://github.com/EngreitzLab/e2g_pipeline",
        "#Contact=Jesse Engreitz (engreitz@stanford.edu)",
        "#Genome=GRCh38",
        "#Description=Accessibility read counts on gene bodies and gene promoter regions",
        "#BiosampleAgnostic=false",
    ]
    output_file = os.path.join(upload_folder, os.path.basename(genelist_file))
    write_df_and_comments(
        df[GENE_NODE_COLUMNS + additional_cols], comment_lines, output_file
    )


def convert_e2g_predictions(cluster, e2g_predictions_file, upload_folder, threshold=None):
    df = pd.read_csv(e2g_predictions_file, sep="\t")
    column_name_mapping = {
        "chr": "Element:chr",
        "start": "Element:start",
        "end": "Element:end",
        "name": "Element:name",
        "TargetGene": "Gene:Symbol",
        "TargetGeneEnsembl_ID": "Gene:EnsemblID",
        "TargetGeneTSS": "Gene:TSS",
        "ENCODE-E2G.Score": "Score",
        "ABC.Score": "Score.ABC",
    }
    df.rename(columns=column_name_mapping, inplace=True)
    df["Element:strand"] = "*"
    df["CellClusterID"] = cluster
    df["BiosampleOntologyTerm"] = ""
    additional_cols = [
        col for col in df.columns if col not in ELEMENT_GENE_EDGE_COLUMNS
    ]

    if threshold:
        description = f"#Description=ENCODE-rE2G predictions with scores >= {threshold}"
    else:
        description = "#Description=ENCODE-rE2G predictions"

    comment_lines = [
        "#FileType=Edge:RegulatoryElementToGene",
        "#Code=https://github.com/EngreitzLab/e2g_pipeline",
        "#Contact=Jesse Engreitz (engreitz@stanford.edu)",
        "#Genome=GRCh38",
        description,
        "#Directional=true",
        "#BiosampleAgnostic=false",
    ]
    output_file = os.path.join(upload_folder, os.path.basename(e2g_predictions_file))
    write_df_and_comments(
        df[ELEMENT_GENE_EDGE_COLUMNS + additional_cols],
        comment_lines,
        output_file,
        compressed=True,
    )


@click.command()
@click.option("--results_dir", type=str, required=True)
@click.option("--igvf_folder_name", type=str, required=True)
@click.option("--threshold", type=float, required=True)
def main(results_dir, igvf_folder_name, threshold):
    prediction_files = glob.glob(os.path.join(results_dir, '*', E2G_PREDICTIONS_FILE))
    cell_clusters = [file.split("/")[-2] for file in prediction_files]
    start = time.time()
    for i, cluster in enumerate(cell_clusters):
        print(f"Processing {cluster}")
        igvf_upload_folder = os.path.join(results_dir, cluster, igvf_folder_name)
        os.makedirs(igvf_upload_folder, exist_ok=True)

        completed_filename = os.path.join(igvf_upload_folder, ".completed")
        if os.path.exists(completed_filename):
            print(f"{cluster} has already been processed. Skipping...")
            continue

        convert_macs(
            cluster,
            os.path.join(results_dir, cluster, MACS_FILE),
            igvf_upload_folder,
        )
        convert_enhancerlist(
            cluster,
            os.path.join(results_dir, cluster, ENHANCERLIST_FILE),
            igvf_upload_folder,
        )
        convert_genelist(
            cluster,
            os.path.join(cluster, results_dir, cluster, GENELIST_FILE),
            igvf_upload_folder,
        )
        convert_e2g_predictions(
            cluster,
            os.path.join(results_dir, cluster, E2G_PREDICTIONS_FILE),
            igvf_upload_folder,
        )
        convert_e2g_predictions(
            cluster,
            os.path.join(results_dir, cluster, E2G_PREDICTIONS_THRESHOLDED_FILE.format(threshold=threshold)),
            igvf_upload_folder,
            threshold=threshold
        )
        with open(completed_filename, 'w'):
            pass

        print(f"Converted {i+1} / {len(cell_clusters)} clusters. Time taken: {int((time.time() - start)/60)} minutes")


if __name__ == "__main__":
    main()
