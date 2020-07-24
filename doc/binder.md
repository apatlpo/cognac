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
  - graphviz
  - geoviews
  - datashader
  - nodejs
  - cartopy
  - gsw
  - pip
  - pip:
    - sat-search
    - pynmea2
```

[mybinder](https://mybinder.org)
