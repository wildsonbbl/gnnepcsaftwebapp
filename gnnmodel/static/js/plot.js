function getplot(dataname, xlegendpos, ytitle, id) {
  var alldata = JSON.parse(dataname);

  var trace1 = {
    x: alldata["T"],
    y: alldata["TML"],
    mode: "markers",
    type: "scatter",
    name: "ThermoML Archive**",
    marker: {
      symbol: "x",
      color: "black",
    },
  };

  var trace2 = {
    x: alldata["T"],
    y: alldata["GNN"],
    mode: "lines",
    type: "scatter",
    name: "GNNePCSAFT",
  };

  var layout = {
    font: {
      family: "Times New Roman",
    },
    legend: {
      x: xlegendpos,
      y: 1,
      font: { family: "monospace", size: 10 },
    },
    paper_bgcolor: "#f8f9fa",
    plot_bgcolor: "#f8f9fa",
    margin: {
      b: 50,
      t: 50,
      l: 50,
      r: 20,
    },
    xaxis: {
      title: {
        text: "Temperature (K)",
      },
      linecolor: "black",
      ticks: "inside",
      minor: {
        ticks: "inside",
      },
      mirror: true,
      showline: true,
      showgrid: false,
    },
    yaxis: {
      title: {
        text: ytitle,
      },
      linecolor: "black",
      ticks: "inside",
      minor: {
        ticks: "inside",
      },
      mirror: true,
      showline: true,
      showgrid: false,
    },
    autosize: true,
    showlegend: true,
  };

  var plot_data = [trace1, trace2];

  Plotly.newPlot(id, plot_data, layout, {
    responsive: true,
    modeBarButtonsToRemove: ["select2d", "lasso2d"],
  });
}
