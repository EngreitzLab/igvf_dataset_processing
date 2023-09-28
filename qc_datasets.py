import glob
import gzip
import os

import click


def lines_are_sorted(line_a, line_b) -> bool:
    chr_a, start_a, end_a = line_a.split("\t")[:3]
    chr_b, start_b, end_b = line_b.split("\t")[:3]
    chr_a = int(chr_a.split("chr")[1])
    chr_b = int(chr_b.split("chr")[1])
    return chr_a <= chr_b and int(start_a) <= int(start_b)


def are_lines_sorted(lines):
    # Only look at chr1, start, end
    for i in range(len(lines) - 1):
        if not lines_are_sorted(lines[i], lines[i + 1]):
            return False
    return True


def is_sorted(filename, num_lines=1000):
    with gzip.open(filename, "rt") as file:  # 'rt' mode for reading as text
        lines = []
        for _ in range(num_lines):
            lines.append(file.readline())
        return are_lines_sorted(lines)


def has_tbi(filename):
    index_file = filename + ".tbi"
    return os.path.exists(index_file)


@click.command()
@click.option("--dataset_dir", type=str, required=True)
def main(dataset_dir):
    non_sorted_files = []
    no_tbi_files = []

    tagAlign_files = glob.glob(os.path.join(dataset_dir, "*", "*tagAlign*.gz"))
    for tagAlign in tagAlign_files:
        tagAlign = os.path.abspath(tagAlign)  # use full path
        if not is_sorted(tagAlign):
            non_sorted_files.append(tagAlign)
        if not has_tbi(tagAlign):
            no_tbi_files.append(tagAlign)

    print(f"Non sorted files: {non_sorted_files}")
    print(f"No tbi files: {no_tbi_files}")


if __name__ == "__main__":
    main()
