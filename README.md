# GNNPCSAFT Webapp

The GNNPCSAFT webapp is an implementation of [our project](https://github.com/wildsonbbl/gnnepcsaft/) that focuses on using Graph Neural Networks ([GNN](https://en.wikipedia.org/wiki/Graph_neural_network)) to estimate the pure-component parameters of the Equation of State [PC-SAFT](https://en.wikipedia.org/wiki/PC-SAFT). We developed this app so the scientific community can access the model's results easily.

The estimated pure-component parameters can be used to calculate thermodynamic properties and compare them with experimental data from the [ThermoML Archive](https://doi.org/10.18434/mds2-2422). More on ThermoML Archive at their [paper](https://doi.org/10.1002/jcc.26842).

You can usually find an app for the project at [SourceForge](https://sourceforge.net/projects/gnnpcsaftapp/). For this app, it will be available on releases too.

A CLI to use a model can be found at [GNNPCSAFT CLI](https://github.com/wildsonbbl/gnnepcsaftcli) and installed with [pipx](https://github.com/pypa/pipx):

```bash
pipx install gnnepcsaftcli
```
