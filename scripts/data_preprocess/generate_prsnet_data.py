import time
import os
import numpy as np
import pandas as pd
import multiprocessing
import sys
import gc


import numpy as np
import numpy.lib.format
import struct

def save(file, array):
    magic_string=b"\x93NUMPY\x01\x00v\x00"
    header=bytes(("{'descr': '"+array.dtype.descr[0][1]+"', 'fortran_order': False, 'shape': "+str(array.shape)+", }").ljust(127-len(magic_string))+"\n",'utf-8')
    if type(file) == str:
        file=open(file,"wb")
    file.write(magic_string)
    file.write(header)
    file.write(array.tobytes())

def pack(array):
    size=len(array.shape)
    return bytes(array.dtype.byteorder.replace('=','<' if sys.byteorder == 'little' else '>')+array.dtype.kind,'utf-8')+array.dtype.itemsize.to_bytes(1,byteorder='little')+struct.pack(f'<B{size}I',size,*array.shape)+array.tobytes()

def load(file):
    if type(file) == str:
        file=open(file,"rb")
    header = file.read(128)
    if not header:
        return None
    descr = str(header[19:25], 'utf-8').replace("'","").replace(" ","")
    shape = tuple(int(num) for num in str(header[60:120], 'utf-8').replace(', }', '').replace('(', '').replace(')', '').split(','))
    datasize = numpy.lib.format.descr_to_dtype(descr).itemsize
    for dimension in shape:
        datasize *= dimension
    return np.ndarray(shape, dtype=descr, buffer=file.read(datasize))

def unpack(data):
    dtype = str(data[:2],'utf-8')
    dtype += str(data[2])
    size = data[3]
    shape = struct.unpack_from(f'<{size}I', data, 4)
    datasize=data[2]
    for dimension in shape:
        datasize *= dimension
    return np.ndarray(shape, dtype=dtype, buffer=data[4+size*4:4+size*4+datasize])


DATA_PATH = sys.argv[1]
PH_PATH = sys.argv[2]
OUTPUT_PATH = sys.argv[3]
EXTENTION = sys.argv[4]
R = sys.argv[5]



def load_data(chr_gene_R_p):
    chr, gene, R, p = chr_gene_R_p
    try:
        df = pd.read_csv(f'{DATA_PATH}/gene_prs_{EXTENTION}/{R}/chr{chr}/{gene}.assoc.clumped.{p}.sscore', sep='\t')
        value = np.nan_to_num(df['SCORE1_AVG'].values.astype(np.float32))
        return value 
    except Exception as e:
        print(e)
        return np.zeros(len(samples)).astype(np.float32)

gene_list = []
genes = {}
for chr in range(1, 23):
    files = os.listdir(f'/blue/sai.zhang/lihan/GRSNet/data/ppi/gene_bed_files_10kb/chr{chr}')
    genes[f'chr{chr}'] = [file[:-4] for file in files]
    gene_list += [file[:-4] for file in files]

genes_chr = [genes[f'chr{chr}'] for chr in range(1, 23)]

samples = pd.read_csv(f'{PH_PATH}', sep=' ', header=None)[0].values

pool = multiprocessing.Pool(32)
prs = pool.map(load_data, [(chr, gene, R, p)
                            for chr in range(1, 23)
                            for gene in genes[f'chr{chr}']
                            for p in ['1e-5', '1e-4', '1e-3', '0.01', '0.05', '0.1', '0.2', '0.3', '0.4', '0.5', '1']])
pool.close()
pool.join()

prs = np.stack(prs, axis=1).astype(np.float32)
for ii in range(len(samples)):
    save(f'{OUTPUT_PATH}/{EXTENTION}_{R}/{samples[ii]}.npy', prs[ii].reshape(1,-1))