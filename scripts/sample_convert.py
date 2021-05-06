import sys
import glob
import re
import time
import pandas as pd
import numpy as np

PROCESSED="processed-data/"
RESULT="results/"

def dec2flatrgb(arr):
    b = arr & 255
    g = (arr >> 8) & 255
    r =(arr >> 16) & 255
    return np.concatenate((r,g,b), axis=None).astype(np.uint8)
    
def dec_arr2rgb_arr(arr):
    from math import sqrt,floor
    int_arr = np.int_(arr)
    no_samples = arr.shape[0]
    dim = floor(sqrt(arr.shape[1]))   # assuming square shape here
    b = (int_arr & 255)
    g = (int_arr >> 8) & 255
    r = (int_arr >> 16) & 255
    b_re = b.reshape(no_samples, dim , dim).astype(np.uint8)
    g_re = g.reshape(no_samples, dim , dim).astype(np.uint8)
    r_re = r.reshape(no_samples, dim , dim).astype(np.uint8)
    return np.stack((r_re, g_re, b_re), axis=-1)

train_samples=pd.read_csv(PROCESSED + "train_sample_labels.tsv",sep="\t")
test_samples=pd.read_csv(PROCESSED +"test_sample_labels.tsv",sep="\t")

train_data=pd.read_csv(PROCESSED + "train_data.tsv",sep="\t")
test_data=pd.read_csv(PROCESSED + "test_data.tsv",sep="\t")

# shuffling with ".sample(frac=1)" might be messing the results

#train_data = data[data['Sample'].isin(train_samples['sample_id'])]
#test_data = data[~data['Sample'].isin(train_samples['sample_id'])]
control_data = pd.read_csv(PROCESSED + "control_sample",sep="\t")
print("                        start converting rgb")

train_img = dec_arr2rgb_arr(train_data.drop(['Sample'],1))
test_img = dec_arr2rgb_arr(test_data.drop(['Sample'],1))
control_img=dec_arr2rgb_arr(control_data.drop(['Sample'],1))
print("                        converting rgb complete")


labels_train_tumor = np.array(train_samples["tumor"]).astype(np.uint8)
labels_test_tumor = np.array(test_samples["tumor"]).astype(np.uint8)

train_sample_names=np.array(train_samples['sample_id'])
test_sample_names=np.array(test_samples['sample_id'])
gene_names = np.array(list(train_data.columns)[1:])

np.savez_compressed(RESULT +"tcga_gtex_genes_data.npz",
                    x_train= train_img,
                    y_train= labels_train_tumor,
                    x_test= test_img,
                    y_test= labels_test_tumor,
                    control= control_img,
                    genes= gene_names,
                    train_samples= train_sample_names,
                    test_samples= test_sample_names)
