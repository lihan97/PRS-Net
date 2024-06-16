import numpy as np
import random
import torch
import dgl

def set_random_seed(seed=22):
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    torch.backends.cudnn.deterministic = True # type: ignore
    torch.backends.cudnn.benchmark = False # type: ignore
    dgl.seed(seed)

def seed_worker(worker_id):
    worker_seed = torch.initial_seed() % 2**32
    np.random.seed(worker_seed)
    random.seed(worker_seed)

def load(file):
    if type(file) == str:
        file=open(file,"rb")
    header = file.read(128)
    if not header:
        return None
    descr = str(header[19:25], 'utf-8').replace("'","").replace(" ","")
    shape = tuple(int(num) for num in str(header[60:120], 'utf-8').replace(', }', '').replace('(', '').replace(')', '').split(','))
    datasize = np.lib.format.descr_to_dtype(descr).itemsize
    for dimension in shape:
        datasize *= dimension
    return np.ndarray(shape, dtype=descr, buffer=file.read(datasize))

def generate_splits(labels, seeds=[42,43,44,45,46,47], train_ratio=0.8, val_ratio=0.1, n_splits=3):
    splits = []
    pos_indexs = np.where(labels == 1)[0]
    neg_indexs = np.where(labels == 0)[0]
    for seed in seeds:
        np.random.seed(seed)
        np.random.shuffle(pos_indexs)
        
        shuffled_indices = np.random.permutation(len(pos_indexs))
        
        total_length = len(pos_indexs)
        
        train_size = int(total_length * train_ratio)
        val_size = int(total_length * val_ratio)
        
        # Split the shuffled indices
        train_indices = shuffled_indices[:train_size]
        val_indices = shuffled_indices[train_size:train_size + val_size]
        test_indices = shuffled_indices[train_size + val_size:]
        
        # Create train, validation, and test sets for positive class
        train_pos_data = pos_indexs[train_indices]
        val_pos_data = pos_indexs[val_indices]
        test_pos_data = pos_indexs[test_indices]
        
        # Sample negative indices for train, validation, and test sets
        train_neg_indices = np.random.choice(neg_indexs, size=int(len(neg_indexs)*0.8), replace=False)
        val_neg_indices = np.random.choice(np.setdiff1d(neg_indexs, train_neg_indices), size=int(len(neg_indexs)*0.1), replace=False)
        test_neg_indices = np.setdiff1d(np.setdiff1d(neg_indexs, train_neg_indices), val_neg_indices)
        
        # Combine positive and negative indices for each split
        train_indices = np.concatenate((train_pos_data, train_neg_indices))
        val_indices = np.concatenate((val_pos_data, val_neg_indices))
        test_indices = np.concatenate((test_pos_data, test_neg_indices))
        
        # Shuffle the combined indices
        np.random.shuffle(train_indices)
        np.random.shuffle(val_indices)
        np.random.shuffle(test_indices)
        
        splits.append((train_indices, val_indices, test_indices))
    
    return splits

def collate_fn_ma(batch):
    feats = [sample['feat'] for sample in batch]
    ancestrys = [sample['ancestry'] for sample in batch]
    labels = [sample['label'] for sample in batch]
    feats = torch.from_numpy(np.stack(feats)).to(torch.float32)
    ancestrys = torch.LongTensor(ancestrys)
    labels = torch.FloatTensor(labels).to(torch.float32).reshape(-1,1)
    return feats, ancestrys, labels

def collate_fn(batch):
    
    feats = [sample['feat'] for sample in batch]
    labels = [sample['label'] for sample in batch]
    feats = torch.from_numpy(np.stack(feats)).to(torch.float32)
    labels = torch.FloatTensor(labels).to(torch.float32).reshape(-1,1)
    return feats, labels

def ancestry_encoding(ancestries):
    encoded_ancestries = []
    for ancestry in ancestries:
        if ancestry == 'EUR':
            encoded_ancestries.append(0)
        if ancestry == 'SAS':
            encoded_ancestries.append(1)
        if ancestry == 'AFR':
            encoded_ancestries.append(2)
        else:
            encoded_ancestries.append(3)
    return encoded_ancestries