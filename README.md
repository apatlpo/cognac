# cognac
Code relevant to underwater floats, acoustic geolocationing

- acoustic/ : acoustic signal processing and modeling

- doc/ : conda and git info

- flow/ : analysis of flow numerical simulations

- geolocation/ : sandbox for geolocation methods

- instrum/ :  float dynamics and control


All scripts require python librairies that may be installed with conda according to the following instructions [here](https://github.com/apatlpo/cognac/blob/master/doc/CONDA.md)


---

In order to access jupyter working on station
```
jupyter-notebook --port=8889
```
```
ssh -N -L localhost:8889:localhost:8889
```
and go to this [address](http://localhost:8889)

