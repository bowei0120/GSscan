###############################################################
## Copyright: GeneGenieDx Corp 2021
## Author: Minzhe Zhang
## Version: 0.1.0
## Date of creation: 01/11/2021
#
## DACO
## Description: This script contains the torch model to predict
#   HRD risk score using CNV data.
#
## usage:
#   import model
###############################################################

# from asyncio.log import logger
import os
import json

# from sys import path
import numpy as np
import pandas as pd
import torch
from torch.utils.data import Dataset
import torch.nn as nn
from sklearn.metrics import roc_auc_score, f1_score, roc_curve, confusion_matrix
from scipy import stats


# https://en.wikipedia.org/wiki/Human_genome
chr_length = np.array(
    [
        250,
        242,
        198,
        190,
        182,
        171,
        159,
        145,
        138,
        134,
        135,
        133,
        114,
        107,
        101,
        90,
        83,
        80,
        59,
        64,
        47,
        51,
    ]
)


def aggregate_chromosome(df, include_cnv=False):
    """
    Aggregate chromosome instability score to sample level HRD risk score.

    Params:
        df (pd.DataFrame): data frame contains the ouput from
        include_cnv (bool): whether to include chromosome average cnv information

    Returns:

    """
    y_true = df.groupby(by=["pat"])["y_true"].mean()
    y_pred = df.groupby(by=["pat"])["ci_pred"].apply(
        lambda x: np.average(x, weights=chr_length)
    )
    if include_cnv:
        cnv = df.groupby(by=["pat"])["cnv"].apply(
            lambda x: np.average(x, weights=chr_length)
        )
        data = pd.concat([cnv, y_true, y_pred], axis=1)
    else:
        data = pd.concat([y_true, y_pred], axis=1)
    return data


def evaluate(model, data_loader):
    """
    Evaluate model.

    Params:
        model (torch.nn.Module): module to evaluate
        data_loader (torch.utils.data.DataLoader): validation data loader

    Returns:
        predicted output of data_loader
    """
    model.eval()
    cnv, y_true, y_pred = [], [], []
    chroms, pats = [], []
    for X, y, chrom, pat in data_loader:
        # X2 = X2.unsqueeze(-1)
        # New in bins
        X = X.unsqueeze(-1)
        y_batch = model(X.float())
        # cnv.append(chromosome_cnv(X.detach().numpy()[0]))
        y_true.append(y.detach().numpy()[0])
        y_pred.append(y_batch.detach().numpy()[0])
        chroms.append(chrom.detach().numpy()[0])
        pats.append(pat[0])
    data_chr = pd.DataFrame(
        {"pat": pats, "chr": chroms, "y_true": y_true, "ci_pred": y_pred}
    )
    return data_chr.sort_values(by=["pat", "chr"]).reset_index(drop=True)


class binlstm(nn.Module):
    """
    LSTM model to predict chromosome instability score.
    """

    def __init__(self, input_dim, hidden_dim):
        super().__init__()
        self.input_dim = input_dim
        self.hidden_dim = hidden_dim
        self.lstm = nn.LSTM(
            input_size=self.input_dim,
            hidden_size=self.hidden_dim,
            batch_first=True,
            num_layers=1,
            bidirectional=True,
        )
        self.fc1 = nn.Linear(self.hidden_dim * 2, self.hidden_dim)
        self.fc2 = nn.Linear(self.hidden_dim, 1)

    def forward(self, x):
        ouput, (hidden, cell) = self.lstm(x)
        x = hidden.reshape(-1)
        x = torch.relu(self.fc1(x))
        y = torch.sigmoid(self.fc2(x))
        return y


class CNVDataset(Dataset):
    """
    Torch dataset to store CNV segment data.
    """

    def __init__(self, data):
        # self.x1 = data["x1"]
        # self.x2 = data["x2"]
        self.x = data["x"]
        self.y = data["y"]
        self.chr = data["chromosome"]
        self.pat = data["pat"]

    def __len__(self):
        return len(self.y)

    def __getitem__(self, index):
        return (
            # self.x1[index],
            # self.x2[index],
            self.x[index],
            self.y[index],
            self.chr[index],
            self.pat[index],
        )
