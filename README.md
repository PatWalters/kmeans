# kmeans
K-means clustering

### Installation
This script and the associated jupyter notebooks require the RDKit which can be [installed using Anaconda](https://www.rdkit.org/docs/Install.html). 

```shell
conda create -c rdkit -n my-rdkit-env rdkit
```
The script requires a few more Python libraries that can be installed with pip and conda. 

```shell
pip install sklearn tqdm docopt jupyter
conda install fastparquet -c conda-forge
```
One of the example Jupyter notebooks also requires the "interact" package. 
```shell
conda install -c conda-forge ipywidgets
```


