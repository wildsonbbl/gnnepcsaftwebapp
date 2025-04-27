// plots the GNN ePC-SAFT model results
// and the ThermoML archive data

function get_layout(xlegendpos, xtitle, ytitle) {
  return {
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
      b: 80,
      t: 50,
      l: 80,
      r: 20,
    },
    xaxis: {
      title: {
        text: xtitle,
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
}

function getplot(data_json, xlegendpos, ytitle, id, trace_name = "GNNePCSAFT") {
  var alldata = JSON.parse(data_json);

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
    name: trace_name,
  };

  var layout = get_layout(xlegendpos, "Temperature (K)", ytitle);

  var plot_data = [trace1, trace2];

  Plotly.newPlot(id, plot_data, layout, {
    responsive: true,
    modeBarButtonsToRemove: ["select2d", "lasso2d"],
  });
}

function get_phase_diagram(phase_diagram_data, xlegendpos, ytitle, y_pos, id) {
  var trace1 = {
    y: phase_diagram_data[y_pos],
    x: phase_diagram_data[2],
    mode: "lines",
    type: "scatter",
    name: "Liquid",
  };
  var trace2 = {
    y: phase_diagram_data[y_pos],
    x: phase_diagram_data[3],
    mode: "lines",
    type: "scatter",
    name: "Vapor",
  };

  var layout = get_layout(xlegendpos, "Density (mol / m³)", ytitle);

  var plot_data = [trace1, trace2];

  Plotly.newPlot(id, plot_data, layout, {
    responsive: true,
    modeBarButtonsToRemove: ["select2d", "lasso2d"],
  });
}
