# GNNePCSAFT App

The project focuses on using Graph Neural Networks ([GNN](https://en.wikipedia.org/wiki/Graph_neural_network)) to estimate the pure-component parameters of the Equation of State [PC-SAFT](https://en.wikipedia.org/wiki/PC-SAFT).

Currently, the model takes into account the hard-chain, dispersive, and associative terms of PC-SAFT. Future work on polar and ionic terms is being studied.

Code is being developed mainly in Pytorch ([PyG](https://pytorch-geometric.readthedocs.io/en/latest/index.html#)).

You can find a model deployed in a Desktop App at [SourceForge](https://sourceforge.net/projects/gnnepcsaft/).

A CLI to use a model can be found at [GNNePCSAFT CLI](https://github.com/wildsonbbl/gnnepcsaftcli) and installed with [pipx](https://github.com/pypa/pipx):

```bash
pipx install gnnepcsaftcli
```

---

GNN models from [GNNePCSAFT](https://github.com/wildsonbbl/gnnepcsaft/).
