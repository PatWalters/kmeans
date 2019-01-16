# kmeans
This script provides an implementation of k-means clustering that uses the ["mini batch k-means" from SciKit Learn](https://scikit-learn.org/stable/modules/generated/sklearn.cluster.MiniBatchKMeans.html) 
together with fingerprints from the [RDKit](http://rdkit.org/). 

### Installation
**Note: This script requires Python 3.6.  Seriously, Python 3.6.**

The script and the associated Jupyter notebooks require the RDKit which can be [installed using Anaconda](https://www.rdkit.org/docs/Install.html). 

```shell
conda create -c rdkit -n my-rdkit-env rdkit
```
Once you've installed the conda environment, you need to activate it. 
```shell
source activate my-rdkit-env
```
The script requires a few more Python libraries that can be installed with pip and conda. 

```shell
pip install sklearn tqdm docopt jupyter matplotlib seaborn scipy
pip install pandas -U
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
```
At the simplest level, you can just call the script with an input file, number of clusters and an output file. In the 
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
kmeans.py all --in test10K.smi --clusters 500 --out test10K_clusters.csv --dim 512
```
If you're processing a larger input file you may want to play around with some of the 
options.  As such, you may not want to regenerate the fingerprints every time. For cases 
like this, the script has an alternate workflow where you can write the fingerprints to disc
then read the fingerprints in and cluster.  The workflow is like this. 
```shell
kmeans.py fp --in test10K.smi 
kmeans.py cluster --fp_file test10K_parquet.gz --clusters 500 --out test10K_clusters.csv
```
Calling the script with the "fp" command creates the fingerprint file test10K_parquet.gz.  This 
fingerprint file is then used in the second clustering step with the "cluster" command.  All of the options discussed above
for fingerprint generation can also be used.   

One more useful trick with k-means clustering is to use a subset of the data to identify the cluster centers
then use these cluster centers to map all of the molecules onto clusters.  By default, if the dataset has more than
10,0000 molecules, the script uses 10% of the data to identify the cluster centers.  This number can be modified by
using the "--sample" flag.  For instance, if we have 25,000 molecules and we want to use 5,000 to define the
cluster centers, we can do something like this. 

```shell
kmeans.py all --in test25K.smi --cluster 500 --out test25K_clusters.csv --sample 5000
```
