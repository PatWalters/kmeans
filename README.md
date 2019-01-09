# kmeans
This script provides an implementation k-means clustering that uses the ["mini batch k-means" implementation in SciKit Learn](https://scikit-learn.org/stable/modules/generated/sklearn.cluster.MiniBatchKMeans.html) 
together with fingerprints from the [RDKit](http://rdkit.org/). 

### Installation
This script and the associated jupyter notebooks require the RDKit which can be [installed using Anaconda](https://www.rdkit.org/docs/Install.html). 

```shell
conda create -c rdkit -n my-rdkit-env rdkit
```
The script requires a few more Python libraries that can be installed with pip and conda. 

```shell
pip install sklearn tqdm docopt jupyter matplotlib seaborn
conda install fastparquet -c conda-forge
```
One of the example Jupyter notebooks also requires the "ipywidgets" package to enable interactive functionality. 
```shell
conda install -c conda-forge ipywidgets
```
I've found that when running on AWS I also have to do this in order to run Jupyter. 
```shell
pip install environment_kernels
```

### Usage

```shell
Usage:
kmeans.py all --in INPUT_FILE_NAME --clusters NUM_CLUSTERS --out OUTPUT_FILE_NAME [--fp_type FP_TYPE] [--dim FP_DIM] [--sample SAMPLE_SIZE]
kmeans.py fp --in INPUT_FILE_NAME [--dim FP_DIM] [--fp_type FP_TYPE] [--fp_bits FP_BITS]
kmeans.py cluster --fp_file FP_FILE_NAME --clusters CLUSTER_FILE_NAME --out OUTPUT_FILE_NAME [--sample SAMPLE_SIZE]
```
At the simplest level, you can just call the script with an input file, number of clusters and an ouput file. In the 
example below, we read a SMILES file with 10,000 molecules and cluster into 500 clusters. This will use the default
RDKit 1024 bit Morgan fingerprint with a radius of 2 (morgan2).  

```shell
kmeans.py --in test10K.smi all --clusters 500 --out test10K_clusters.csv
```

The script supports a few other fingerprint types, if you want to change the fingerprint type
you can use the "--fptype" flag.  Supported types are:
- morgan2 - RDKit Morgan fingerprint with a radius of 2
- morgan3 - RDKit Morgan fingerprint with a radius of 3
- ap - RDKit Atom Pair fingerprint
- rdk5 - RDKit path fingerprint with a maximum path length of 5

So if you want to use the Atom Pair fingerprints you can do this. 

```shell
kmeans.py all --in test10K.smi --clusters 500 --out test10K_clusters.csv --fp_type ap
```
You can also change the number of fingerprint bits with the "--fpbits" flag.  To change to a 512 bit
fingerprint, you could do this. 
```shell
kmeans.py all --in test10K.smi --clusters 500 --out test10K_clusters.csv --fp_bits 512
```
