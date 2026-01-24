# GNNPCSAFT Web App

The GNNPCSAFT web app is an implementation of [our project](https://github.com/wildsonbbl/gnnepcsaft/) that focuses on using Graph Neural Networks ([GNN](https://en.wikipedia.org/wiki/Graph_neural_network)) to estimate the pure-component parameters of the Equation of State [PC-SAFT](https://en.wikipedia.org/wiki/PC-SAFT). We developed this app so the scientific community can access the model's results easily.

In this app, the estimated pure-component parameters can be used to calculate thermodynamic properties and compare them with experimental data from the [ThermoML Archive](https://doi.org/10.18434/mds2-2422). More on ThermoML Archive at their [paper](https://doi.org/10.1002/jcc.26842).

On [releases](https://github.com/wildsonbbl/gnnepcsaftwebapp/releases), you find a electron app for the chat. A container image is also available on [Docker Hub](https://hub.docker.com/r/wildsonbbl/gnnpcsaftwebapp), and can be run using:

```bash
docker run -p 19770:8000 wildsonbbl/gnnpcsaftwebapp:latest
```

Other implementations with GNNPCSAFT:

- [GNNPCSAFT CLI](https://github.com/wildsonbbl/gnnepcsaftcli)
- [GNNPCSAFT APP](https://github.com/wildsonbbl/gnnpcsaftapp)
- [GNNPCSAFT MCP](https://github.com/wildsonbbl/gnnepcsaft_mcp_server)
- [GNNPCSAFT Chat](https://github.com/wildsonbbl/gnnpcsaftchat)
