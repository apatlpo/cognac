# Prepare binder Deployment

Follow [pangeo-binder-template](https://github.com/pangeo-data/pangeo-binder-template)
and create a `binder` directory in your repository.

Add extra libraries to `binder/environment.yml` file, for example:

```
name: cognac
channels:
  - conda-forge
dependencies:
  - pangeo-notebook==2020.07.03
  - xarray
  - hvplot
  - geoviews
  - datashader
  - jupyter-panel-proxy
  - nodejs
  - cartopy
  - gsw
  - pip
  - pip:
    - pynmea2
    - graphviz
    - acoustics
    - -e ../
```

[mybinder](https://mybinder.org)

Panel binder deployment see [panel doc](https://panel.holoviz.org/user_guide/Server_Deployment.html)

pip install tweaks: [1](https://github.com/conda/conda/blob/4.7.11/tests/conda_env/support/advanced-pip/environment.yml), [2](https://remotedevdaily.com/how-to-install-packages-from-github-with-pip-and-pipenv/)

[mybinder cognac](https://mybinder.org/v2/gh/apatlpo/cognac/master?urlpath=/panel/design)
