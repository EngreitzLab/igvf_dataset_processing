import glob
import os

import click
import pandas as pd

HiC_INFO = {
    "HiC_dir": "/oak/stanford/groups/engreitz/Projects/ABC/HiC/avg_track2",
    "HiC_type": "avg",
    "HiC_gamma": 1.024238616787792,
    "HiC_scale": 5.9594510043736655,
    "HiC_resolution": 5000,
}

CONFIG_TEMPLATE_FILE = "https://raw.githubusercontent.com/broadinstitute/ABC-Enhancer-Gene-Prediction/dev/config/config_biosamples_template.tsv"


@click.command()
@click.option("--dataset_dir", type=str, required=True)
@click.option("--config_name", type=str, default="biosamples_config.tsv")
def main(dataset_dir, config_name):
    clusters = os.listdir(dataset_dir)
    biosamples = pd.read_csv(CONFIG_TEMPLATE_FILE, sep="\t")

    tagAlign_files = glob.glob(os.path.join(dataset_dir, "*", "*tagAlign*.gz"))
    for i, tagAlign in enumerate(tagAlign_files):
        tagAlign = os.path.abspath(tagAlign)  # use full path
        cluster = tagAlign.split("/")[-2]
        biosample = {
            "biosample": cluster,
            "ATAC": tagAlign,
            "default_accessibility_feature": "ATAC",
        }
        biosample.update(HiC_INFO)
        biosamples.loc[i] = biosample

    biosamples.to_csv(config_name, sep="\t", index=False)
    print(f"Generate config: {config_name}")


if __name__ == "__main__":
    main()
