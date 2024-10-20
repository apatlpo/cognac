# Install useful libraries for cognac:

To update with latest approach:

````
conda create -n insitu -c apatlpo -c pyviz -c conda-forge python=3.10 pynsitu jupyterlab seaborn pyTMD xrft xhistogram geojson pynmea2 rioxarray
conda activate insitu
pip install flatten_dict
# then manually install cognac
pip install freshwater
conda install -c conda-forge bambi # baysian inference
```

Download Miniconda3 (i.e. for python3) from the [conda website](https://conda.io/miniconda.html)
```
bash Miniconda3-latest-MacOSX-x86_64.sh
bash
conda update conda
conda create -n cognac -c conda-forge \
    python=3.8 xarray graphviz netCDF4 dask-jobqueue jupyterlab \
    hvplot geoviews datashader nodejs ipywidgets \
    folium cartopy gsw cmocean pyinterp pytide pyTMD \
conda activate cognac
pip install pynmea2
pip install geojsoncontour
conda install ffmpeg
pip install graphviz
pip install git+https://github.com/python-acoustics/python-acoustics.git
git clone https://github.com/apatlpo/cognac.git
cd cognac
pip install -e .
jupyter labextension install @jupyter-widgets/jupyterlab-manager \
                             jupyter-leaflet
```

Run a jupyter notebook with the following command:
```
jupyter-lab
```

In order to add the environnement to kernels available to jupyter, you need to run:
```
python -m ipykernel install --user --name cognac --display-name "COGNAC project env"
```

Uninstall library after `pip install -e .`:
- remove the egg file ( `print(distributed.__file__)` for example)
- from file `easy-install.pth`, remove the corresponding line (it should be a path to the source directory or of an egg file).

# General information about miniconda:

## Overview

Miniconda installers contain the conda package manager and Python.
Once miniconda is installed, you can use the conda command to install any other packages and create environments.

After downloading `Miniconda3-latest-Linux-x86_64.sh` or `Miniconda3-latest-MacOSX-x86_64.sh` you need to run it with: `bash Miniconda3-latest-MacOSX-x86_64.sh`

Miniconda must be used with bash. If you want to use it with csh, add in your .cshrc (not ideal solution)
```
#
#----------------------------------------------------------------
# alias Miniconda
#----------------------------------------------------------------
#
setenv PATH ${PATH}: /home/machine/username/miniconda3/bin
alias source_activate 'setenv OLDPATH ${PATH};setenv PATH /home/machine/username/miniconda3/envs/\!*/bin:${PATH}'
alias source_deactivate 'setenv PATH $OLDPATH'
```
where machine is the name of your computer and username is your username.


## Main commands:
What version, update conda
```
conda --version
conda update conda
```
Create new environment myenv
```
conda create --name myenv python
```
Switch to another environment (activate/deactivate) (or source_activate in csh)
```
conda activate myenv
```
To change your path from the current environment back to the root (or source_deactivate in csh)
```
conda deactivate
```
List all environments
```
conda info --envs
```
Delete an environment
```
conda remove --name myenv --all
```
View a list of packages and versions installed in an environmentSearch for a package
```
conda list
```
Check to see if a package is available for conda to install
```
conda search packagename
```
Install a new package
```
conda install packagename
```
Remove conda
```
rm -rf /home/machine/username/miniconda3
```
where machine is the name of your computer and username is your username.


## Install a package from Anaconda.org

For packages that are not available using conda install, we can next look on Anaconda.org. Anaconda.org is a package management service for both public and private package repositories. Anaconda.org is a Continuum Analytics product, just like Anaconda and Miniconda.

In a browser, go to http://anaconda.org. We are looking for a package named “pestc4py”
There are more than a dozen copies of petsc4py available on Anaconda.org, first select your platform, then you can sort by number of downloads by clicking the “Downloads” heading.

Select the version that has the most downloads by clicking the package name. This brings you to the Anaconda.org detail page that shows the exact command to use to download it:

Check to see that the package downloaded
```
conda list
```

## Install a package with pip

For packages that are not available from conda or Anaconda.org, we can often install the package with pip (short for “pip installs packages”).
Exporting environment

```
#conda env export > environment.yml on a machine
conda env export --from-history | cut -f 1 -d '=' > environment.yml
conda env create -f environment.yml -n $ENV_NAME on the new machine
```
