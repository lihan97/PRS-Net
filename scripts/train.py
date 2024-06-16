import sys
sys.path.append("..")
import argparse
import warnings
import torch
from torch.utils.data import DataLoader
from torch import nn
from torchmetrics import AUROC
import dgl
import pandas as pd
import numpy as np
warnings.filterwarnings("ignore")

from src.dataset import Dataset
from src.utils import generate_splits, seed_worker, collate_fn, ancestry_encoding, set_random_seed
from src.model import PRSNet
from src.trainer import Trainer

def parse_args():
    parser = argparse.ArgumentParser(description="Arguments for training LiGhT")
    parser.add_argument("--data_path", type=str, required=True)
    parser.add_argument("--dataset", type=str, required=True)
    parser.add_argument("--num_workers", type=int, default=4)
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = parse_args()
    set_random_seed(22)
    ## Data loading and splits generating
    info_df = pd.read_csv(f'{args.data_path}/{args.dataset}/info.csv')
    sample_ids = info_df['sample_id'].values
    labels = torch.from_numpy(info_df['label'].values)
    ancestries = ancestry_encoding(info_df['ancestry'].values)
    splits = generate_splits(labels)
    ggi_graph = dgl.load_graphs(f'../data/ggi_graph.bin')[0][0]

    ## Device
    device = torch.device('cuda' if (torch.cuda.is_available()) else 'cpu')

    ## Validation
    for split_id, (train_ids, val_ids, test_ids) in enumerate(splits):
        train_set = Dataset(args.data_path, args.dataset, sample_ids=sample_ids[train_ids],labels=labels[train_ids], balanced_sampling=True)
        val_set = Dataset(args.data_path, args.dataset, sample_ids=sample_ids[val_ids],labels=labels[val_ids], balanced_sampling=False)
        test_set = Dataset(args.data_path, args.dataset, sample_ids=sample_ids[test_ids],labels=labels[test_ids], balanced_sampling=False)

        train_loader = DataLoader(train_set, batch_size=512, shuffle=False, num_workers=args.num_workers, worker_init_fn=seed_worker, drop_last=True, pin_memory=True, collate_fn=collate_fn)
        val_loader = DataLoader(val_set, batch_size=128, shuffle=False, num_workers=args.num_workers, worker_init_fn=seed_worker, drop_last=False, pin_memory=True, collate_fn=collate_fn)
        test_loader = DataLoader(test_set, batch_size=128, shuffle=False, num_workers=args.num_workers, worker_init_fn=seed_worker, drop_last=False, pin_memory=True, collate_fn=collate_fn)

        model = PRSNet().to(device)
        loss_fn = nn.BCEWithLogitsLoss(reduction='mean')
        optimizer = torch.optim.AdamW(model.parameters(), lr=1e-4, weight_decay=0)
        metric = AUROC(task='binary')
        trainer = Trainer(device=device)
        best_val_score, best_test_score = trainer.train_and_test(model, ggi_graph, loss_fn, optimizer, metric, train_loader, val_loader, test_loader)
        print("----------------Split {split_id} final result----------------", flush=True)
        print(f"best_val_score: {best_val_score}, best_test_score: {best_test_score}")
        
