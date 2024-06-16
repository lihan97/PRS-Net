import torch
import random
from .utils import load



class Dataset(torch.utils.data.Dataset):  # type: ignore
    def __init__(self, data_path, dataset, sample_ids, labels, balanced_sampling=False, ancestries=None, multiple_ancestries=False):
        self.sample_ids = sample_ids
        self.labels = labels
        self.ancestries = ancestries
        self.data_path = data_path
        self.dataset = dataset
        self.multiple_ancestries = multiple_ancestries
        self.balanced_sampling = balanced_sampling
        if balanced_sampling == True:
            self.classes = torch.unique(self.labels)
            self.class_indices = {cls.item(): torch.where(self.labels == cls)[0] for cls in self.classes}
            self.n_samples = 10000000
        else:
            self.n_samples = len(self.labels)
    def __len__(self):
        return self.n_samples

    def __balanced_sampling__(self):
        cls = random.choice(self.classes)
        class_indices = self.class_indices[cls.item()]
        index = random.choice(class_indices)
        return index

    def __getitem__(self, index):
        if self.balanced_sampling:
            index = self.__balanced_sampling__()
        sample_id = self.sample_ids[index]
        label = self.labels[index]
        feat = load(f'{self.data_path}/{self.dataset}/feats/{sample_id}.npy')[0] # type: ignore
        if self.multiple_ancestries:
            return {"feat":feat, "ancestry":torch.FloatTensor([self.ancestries[index]]), "label":label} # type: ignore
        else:
            return {"feat":feat, "label":label}