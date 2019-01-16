#!/usr/bin/env python

import sys
import os
import math

import numpy
from rdkit import Chem, DataStructs
from rdkit.Chem import rdMolDescriptors as rdmd

from sklearn.cluster import MiniBatchKMeans

import pandas as pd
from tqdm import tqdm
import time
import numpy as np

from scipy.spatial.distance import cdist

from docopt import docopt


class FingerprintGenerator:
    def __init__(self, fp_type, fp_bits=2048):
        """
        :param fp_type: fingerprint type
        :param fp_bits: number of fingerprint bits
        """
        self.fp_type = fp_type
        self.fp_dict = {}
        self.fp_dict['morgan2'] = [lambda m: rdmd.GetMorganFingerprintAsBitVect(m, 2, nBits=fp_bits), fp_bits]
        self.fp_dict['morgan3'] = [lambda m: rdmd.GetMorganFingerprintAsBitVect(m, 3, nBits=fp_bits), fp_bits]
        self.fp_dict['ap'] = [lambda m: rdmd.GetHashedAtomPairFingerprintAsBitVect(m, nBits=fp_bits), fp_bits]
        self.fp_dict['rdk5'] = [lambda m: Chem.RDKFingerprint(m, maxPath=5, fpSize=fp_bits, nBitsPerHash=2), fp_bits]
        if self.fp_dict.get(fp_type):
            self.fp_function = self.fp_dict[fp_type]
        else:
            print("invalid fingerprint type: %s" % fp_type)
            sys.exit(0)

    def get_fp_types(self):
        """
        Get the list of available fingerprint types
        :return: list of fingerprint type
        """
        return [x for x in self.fp_dict.keys() if x != "ref"]

    def get_names(self):
        """
        Get the names for the fingerprint bits
        :return: list of bit names
        """
        num_bits = self.fp_function[1]
        name_list = ["%s_%d" % (self.fp_type, i) for i in range(0, num_bits)]
        return name_list

    def get_fp(self, mol):
        """
        Get a fingerprint
        :param mol: input molecule
        :return: fingerprint for the molecule
        """
        return self.fp_function[0](mol)

    def get_numpy_fp(self, mol):
        """
        Get a fingerprint as a numpy array
        :param mol: input molecule
        :return: numpy array of 1 and 0 for fingerprint bits
        """
        fp = self.fp_function[0](mol)
        arr = numpy.zeros((1,), numpy.int)
        DataStructs.ConvertToNumpyArray(fp, arr)
        return arr


def molecule_supplier_from_name(input_file_name):
    """
    Get the appropriate molecule supplier based on file extension
    :param input_file_name: input file name
    :return: molecule supplier
    """
    ext = os.path.splitext(input_file_name)[-1]
    if ext == ".smi":
        suppl = Chem.SmilesMolSupplier(input_file_name, titleLine=False)
    elif ext == ".sdf":
        suppl = Chem.SDMolSupplier(input_file_name)
    elif ext == ".mae":
        suppl = Chem.MaeMolSupplier(input_file_name)
    else:
        print("%s is not a valid molecule extension" % ext)
        sys.exit(1)
    return suppl


def generate_fingerprint_df(infile_name, fp_type="morgan2", fp_bits=1024):
    """
    Read an input file and generate fingerprints to infile_name + _parquet.gz
    :param infile_name: input file name
    :param fp_type: fingerprint file
    :param fp_bits: fingerprint bits
    :return: fingerprint dataframe
    """
    fingerprint_generator = FingerprintGenerator(fp_type, fp_bits)
    suppl = molecule_supplier_from_name(infile_name)
    fp_list = []
    name_list = []
    smiles_list = []
    print(f"Generating {fp_type} fingerprints with {fp_bits} bits")
    for mol in tqdm(suppl):
        if mol:
            smiles = Chem.MolToSmiles(mol)
            fp_list.append(fingerprint_generator.get_numpy_fp(mol))
            name_list.append(mol.GetProp("_Name"))
            smiles_list.append(smiles)
    start = time.time()
    df = pd.DataFrame(np.array(fp_list), columns=fingerprint_generator.get_names())
    elapsed = time.time() - start
    df.insert(0, "SMILES", smiles_list)
    df.insert(1, "Name", name_list)
    print(f"{elapsed:.1f} sec required to generate dataframe")
    return df


def write_fingerprint_df(df, outfile_name):
    start = time.time()
    df.to_parquet(outfile_name, engine="fastparquet", compression="gzip")
    elapsed = time.time() - start
    print(f"{elapsed:.1f} sec required to write {outfile_name}")


def read_fingerprint_df(fp_file_name):
    start = time.time()
    df = pd.read_parquet(fp_file_name, engine='fastparquet')
    num_rows, num_cols = df.shape
    elapsed = time.time() - start
    print(f"Read {num_rows} rows from {fp_file_name} in {elapsed:.1f} sec, fingerprint dimension is {num_cols - 2}")
    return df


def find_cluster_centers(df,centers):
    center_set = set()
    for k,v in df.groupby("Cluster"):
        fp_list = v.values[0::,3::]
        dist_list = cdist([centers[k]],fp_list)
        min_idx = np.argmin([dist_list])
        center_set.add(v.Name.values[min_idx])
    return ["Yes" if x in center_set else "No" for x in df.Name.values]



def kmeans_cluster(df, num_clusters, outfile_name, sample_size=None):
    """
    :param df: fingerprint dataframe
    :param num_clusters: number of clusters
    :param outfile_name: output file containing molecule name and cluster id
    :param sample_size: number of molecules to use train the clustering method (with large files)
    :return: None
    """
    num_rows, num_cols = df.shape
    if num_rows > 10000:
        if sample_size:
            rows_to_sample = sample_size
        else:
            # number of samples needs to at least equal the number of clusters
            rows_to_sample = max(int(num_rows / 10),num_clusters)
        train_df = df.sample(rows_to_sample)
        print(f"Sampled {rows_to_sample} rows")
    else:
        train_df = df
    arr = np.array(train_df.values[0::, 2::], dtype=np.float16)

    start = time.time()
    km = MiniBatchKMeans(n_clusters=num_clusters, random_state=0, batch_size=3 * num_clusters)
    km.fit(arr)
    chunk_size = 500
    all_data = np.array(df.values[0::, 2::], dtype=np.bool)
    chunks = math.ceil(all_data.shape[0] / chunk_size)
    out_list = []
    # It looks like the predict method chokes if you send too much data, chunking to 500 seems to work
    cluster_id_list = []
    for row, names in tqdm(zip(np.array_split(all_data, chunks), np.array_split(df[['SMILES', 'Name']].values, chunks)),
                           total=chunks, desc="Processing chunk"):
        p = km.predict(row)
        cluster_id_list += list(p)
    elapsed = time.time() - start
    df.insert(2,"Cluster",cluster_id_list)
    center_list = find_cluster_centers(df,km.cluster_centers_)
    df.insert(3,"Center",center_list)
    out_df = df[["SMILES", "Name", "Cluster","Center"]]
    print(f"Clustered {num_rows} into {num_clusters} in {elapsed:.1f} sec")
    out_df.to_csv(outfile_name, index=False)


def main():
    command_str = """Usage:
kmeans.py all --in INPUT_FILE_NAME --clusters NUM_CLUSTERS --out OUTPUT_FILE_NAME [--fp_type FP_TYPE] [--dim FP_DIM] [--sample SAMPLE_SIZE] 
kmeans.py fp --in INPUT_FILE_NAME [--dim FP_DIM] [--fp_type FP_TYPE]
kmeans.py cluster --fp_file FP_FILE_NAME --clusters CLUSTER_FILE_NAME --out OUTPUT_FILE_NAME [--sample SAMPLE_SIZE]

Options:
--in INPUT_FILE_NAME 
--clusters NUM_CLUSTERS number of clusters to output
--out OUTPUT_FILE_NAME output csv file with SMILES, molecule name and cluster id
--dim FP_DIM number of fingerprint bits
--sample SAMPLE_SIZE number of molecules to use for training
--fp_file FP_FILE_NAME name of fingerprint file created with the "fp" option
--fp_type FP_TYPE fingerprint type, must be one of morgan2, morgan3, ap, rdkit5
"""

    cmd_input = docopt(command_str)

    infile_name = cmd_input.get("--in")
    num_clusters = cmd_input.get("--clusters")
    if num_clusters:
        num_clusters = int(num_clusters)
    fp_file_name = cmd_input.get("--fp_file")
    outfile_name = cmd_input.get("--out")
    fp_dim = cmd_input.get("--dim") or 1024
    fp_dim = int(fp_dim)
    fp_type = cmd_input.get("--fp_type") or "morgan2"
    num_sample = cmd_input.get("--sample")
    if num_sample:
        num_sample = int(num_sample)

    if cmd_input.get("all"):
        fp_df = generate_fingerprint_df(infile_name, fp_type=fp_type, fp_bits=fp_dim)
        kmeans_cluster(fp_df, num_clusters, outfile_name, sample_size=num_sample)
    elif cmd_input.get("fp"):
        fp_df = generate_fingerprint_df(infile_name, fp_type=fp_type, fp_bits=fp_dim)
        name, _ = os.path.splitext(infile_name)
        fp_file_name = name + "_parquet.gz"
        write_fingerprint_df(fp_df, fp_file_name)
    elif cmd_input.get("cluster"):
        fp_df = read_fingerprint_df(fp_file_name)
        kmeans_cluster(fp_df, num_clusters, outfile_name, sample_size=num_sample)


if __name__ == "__main__":
    main()
