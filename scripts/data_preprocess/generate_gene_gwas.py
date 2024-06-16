import pandas as pd
import numpy as np
import os
import sys
from multiprocessing import Pool

def process_file(file):
    try:
        gene_snp_bed = pd.read_csv(f"{ROOT_PATH}/gene_snp_bed_files_{EXTENTION}/chr{chr}/{file}", sep='\t', header=None)
        cur_gene_snps = gene_snp_bed[3].values
        gene_name = file[:-4]
        with open(f"{ROOT_PATH}/gene_gwas_{EXTENTION}/chr{chr}/{gene_name}.assoc", 'w') as f:
            f.write(' '.join(df.columns.values.astype(str)) + '\n')
            for snp in cur_gene_snps:
                f.write(' '.join(df.iloc[SNP_to_id_dict[snp]].values.astype(str)) + '\n')
    except Exception as e:
        print(e)

if __name__ == '__main__':
    ROOT_PATH = sys.argv[1]
    EXTENTION = sys.argv[2]
    df = pd.read_csv(f'{ROOT_PATH}/gwas.QC.txt',sep='\t')
    SNP_to_id_dict = {}
    snps = df['SNP'].values
    for i in range(len(df)):
        SNP_to_id_dict[snps[i]] = i

    for chr in range(1,23):
        files = os.listdir(f"{ROOT_PATH}/gene_snp_bed_files_{EXTENTION}/chr{chr}/")
        with Pool(24) as pool:
            pool.map(process_file, files)