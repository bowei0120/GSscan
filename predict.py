###############################################################
## Copyright: Guangdong Jiyin Biotech Co. Ltd & GeneGenieDx Corp
#
## GSscan
## Description: This script is to predict GSscan score using
#   existing model.
#
## Usage:
#   python predict.py -i <bin_count_file> -w <weight_path>
###############################################################

import pandas as pd
import torch
from torch.utils.data import DataLoader
from sklearn.preprocessing import scale
import os
import argparse


from model import (
    chr_length,
    binlstm,
    evaluate,
    aggregate_chromosome,
    CNVDataset,
)


def read_bins(bin_path, hrd=None):
    """
    Read bin counts

    Params:
        bin_paths (array-like): list of bin count paths
        hrd (pd.Sereis): hrd labels for different samples

    Returns (dict):
        processed data for torch model
        pat: list of patient id
        chromosome: list of chromosome name
        x: list of bin counts array
        y: list of patient HRD labels assigned each chromosome
    """
    chroms = [str(chr) for chr in range(1, 23)]
    pat, chromosome, x, y = [], [], [], []
    bins = pd.read_csv(bin_path, sep="\t", index_col=0)

    # scale bin counts for each sample
    bins = pd.DataFrame(scale(bins.T).T, index=bins.index, columns=bins.columns)

    for pat_id in bins.index:
        pat_id = pat_id.rstrip()
        if hrd is None:
            label = 0
        else:
            if pat_id in hrd.index:
                label = 1 if hrd[pat_id] > 0.5 else 0

        # convert data for each chromosome
        for chrom in chroms:
            pat.append(pat_id)
            chromosome.append(int(chrom))
            bin_res = bins.loc[
                pat_id,
                (
                    bins.columns.str.startswith(chrom + "p")
                    + bins.columns.str.startswith(chrom + "q")
                ),
            ]
            x.append(bin_res.values)
            y.append(label)

    return {"pat": pat, "chromosome": chromosome, "x": x, "y": y}


def process_test_data(test_bin_path):
    """
    Process bin count result.

    Params:
        test_data_bin_paths (list): list of bin paths

    Returns:
        test_data_loader
    """

    test_data = CNVDataset(read_bins(test_bin_path))
    test_data_loader = DataLoader(test_data, batch_size=1, shuffle=True)

    return test_data_loader


def parse_input() -> argparse.Namespace:
    """
    Parse input argument
    """
    # argument parser
    parser = argparse.ArgumentParser(
        description="Call GSscan model to predict HRD risk score from bin count file"
    )
    parser.add_argument(
        "-i",
        "--in-path",
        required=True,
        help="Path of bin count file",
    )
    parser.add_argument(
        "-w",
        "--weight-path",
        required=True,
        help="directory where torch model and model weight is stored",
    )
    parser.add_argument(
        "-o",
        "--output-file",
        default=".",
        help="output file",
    )
    return parser.parse_args()


def main():
    args = parse_input()

    # load model
    weight_path = args.weight_path
    assert os.path.exists(weight_path), "model weight does not exist: {}".format(
        weight_path
    )
    model = binlstm(input_dim=1, hidden_dim=8)
    model.load_state_dict(torch.load(weight_path))

    # load data
    data_to_predict = process_test_data(args.in_path)

    # predict score
    data_chr = evaluate(model, data_to_predict)
    data_pat = aggregate_chromosome(data_chr)
    data_pat.to_csv(args.output_file, sep="\t")


if __name__ == "__main__":
    main()
