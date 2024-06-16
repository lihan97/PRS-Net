import torch
from torch import nn
import dgl
import numpy as np
class Trainer:
    def __init__(self, device, eval_interval=100, n_steps=20000, n_early_stop=50, log_interval=20, multiple_ancestries=False):
        self.device=device
        self.eval_interval = eval_interval
        self.n_steps = n_steps
        self.n_early_stop = n_early_stop
        self.log_interval = log_interval
        self.multiple_ancestries = multiple_ancestries
    def forward_batch(self, model, ggi_graph, batch):
        feats, labels = batch
        batched_graph = dgl.batch([ggi_graph]*len(labels)).to(self.device)
        feats, labels = feats.to(self.device), labels.to(self.device)
        outputs, attn_scores = model(batched_graph, feats)
        return labels, outputs, attn_scores
    def forward_batch_ma(self, model, ggi_graph, batch):
        feats, ancestries, labels = batch
        batched_graph = dgl.batch([ggi_graph]*len(labels)).to(self.device)
        feats, ancestries, labels = feats.to(self.device), ancestries.to(self.device), labels.to(self.device)
        outputs, ph_attn_scores, anc_attn_scores = model(batched_graph, feats, ancestries)
        return labels, outputs, ph_attn_scores, anc_attn_scores
    def train_and_test(self, model, ggi_graph, loss_fn, optimizer, metric, train_loader, val_loader, test_loader):
        best_val_score, best_test_score = 0, 0
        running_loss = []
        cur_early_stop = 0

        data_iter = iter(train_loader)
        next_batch = next(data_iter)
        next_batch = [ _.cuda(non_blocking=True) for _ in next_batch ]
        print("----------------Training----------------", flush=True)
        for cur_step in range(len(train_loader)):
            ## Forward pass
            model.train()
            optimizer.zero_grad()
            batch = next_batch 
            if cur_step + 1 != len(train_loader): 
                next_batch = next(data_iter)
                next_batch = [ _.cuda(non_blocking=True) for _ in next_batch]
            if self.multiple_ancestries:
                labels, preds, ph_attn_scores, anc_attn_scores = self.forward_batch_ma(model, ggi_graph, batch)
            else:
                labels, preds, attn_scores = self.forward_batch(model, ggi_graph, batch)

            loss = loss_fn(preds, labels)

            # Backward pass and optimization
            loss.backward()
            optimizer.step()

            running_loss.append(loss.detach().cpu().numpy())
            if (cur_step+1) % self.log_interval == 0:
                print(f"[{cur_step+1}] loss: {np.mean(running_loss):.3f}", flush=True)
                running_loss = []
            if (cur_step+1) % self.eval_interval == 0:
                running_loss = []
                print("----------------Validating----------------", flush=True)
                val_score = self.evaluate(model, ggi_graph, val_loader, metric)
                if val_score > best_val_score:
                    best_val_score = val_score
                    best_test_score = self.evaluate(model, ggi_graph, test_loader, metric)
                    cur_early_stop = 0
                else:
                    cur_early_stop += 1
                    if cur_early_stop == self.n_early_stop: break
                print(f"[{cur_step+1}] val_score: {best_val_score:.3f}, test_score: {best_test_score: .3f}", flush=True)
                print("----------------Training----------------", flush=True)
            if cur_step == self.n_steps: break
        return best_val_score, best_test_score
        
    def evaluate(self, model, ggi_graph, test_loader, metric):
        with torch.no_grad():
            model.eval()
            preds_list, labels_list, ancestrys_list = [], [], []
            for batch in test_loader:
                if self.multiple_ancestries:
                    labels, preds, ph_attn_scores, anc_attn_scores = self.forward_batch_ma(model, ggi_graph, batch)
                else:
                    labels, preds, attn_scores = self.forward_batch(model, ggi_graph, batch)
                preds_list.append(preds.detach())
                labels_list.append(labels.detach())
            preds = torch.cat(preds_list).reshape(-1)
            labels = torch.cat(labels_list).reshape(-1)
            return metric(preds, labels).item()