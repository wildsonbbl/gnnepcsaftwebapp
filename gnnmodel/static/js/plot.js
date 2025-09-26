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

function get_ternary_lle_phase_diagram(ternary_lle_phase_diagram_data) {
  var trace1 = {
    a: ternary_lle_phase_diagram_data["x0"],
    b: ternary_lle_phase_diagram_data["x1"],
    c: ternary_lle_phase_diagram_data["x2"],
    mode: "markers",
    type: "scatterternary",
    name: "Liquid phase 1",
  };
  var temp_min_el = document.getElementById("id_temp_min");
  var pressure_el = document.getElementById("id_pressure");
  var temp_min = temp_min_el
    ? (temp_min_el.value || temp_min_el.textContent || "").trim()
    : "";
  var pressure = pressure_el
    ? (pressure_el.value || pressure_el.textContent || "").trim()
    : "";
  var _ternaryTitle = "LLE at " + temp_min + " K and " + pressure + " Pa";
  var trace2 = {
    a: ternary_lle_phase_diagram_data["y0"],
    b: ternary_lle_phase_diagram_data["y1"],
    c: ternary_lle_phase_diagram_data["y2"],
    mode: "markers",
    type: "scatterternary",
    name: "Liquid phase 2",
  };
  Plotly.newPlot(
    "ternary_lle_phase_diagram",
    [trace1, trace2],
    {
      title: _ternaryTitle,
      font: {
        family: "Times New Roman",
      },
      legend: {
        orientation: "h",
        font: { family: "monospace", size: 10 },
      },
      paper_bgcolor: "#f8f9fa",
      plot_bgcolor: "#f8f9fa",
      margin: {
        b: 50,
        t: 50,
        l: 50,
        r: 50,
      },
      autosize: true,
      showlegend: true,
      ternary: {
        sum: 1,
        aaxis: {
          title: { text: "A" },
          min: 0,
          linewidth: 2,
        },
        baxis: {
          title: { text: "B" },
          min: 0,
          linewidth: 2,
        },
        caxis: {
          title: { text: "C" },
          min: 0,
          linewidth: 2,
        },
      },
    },
    {
      responsive: true,
      modeBarButtonsToRemove: ["select2d", "lasso2d"],
    }
  );
}

function get_binary_lle_phase_diagram(binary_lle_phase_diagram_data) {
  var trace1 = {
    x: binary_lle_phase_diagram_data["x0"],
    y: binary_lle_phase_diagram_data["temperature"],
    mode: "lines",
    type: "scatter",
    name: "Liquid phase 1",
  };
  var trace2 = {
    x: binary_lle_phase_diagram_data["y0"],
    y: binary_lle_phase_diagram_data["temperature"],
    mode: "lines",
    type: "scatter",
    name: "Liquid phase 2",
  };
  Plotly.newPlot(
    "binary_lle_phase_diagram",
    [trace1, trace2],
    get_layout(0, "x<sub>1</sub>", "Temperature (K)"),
    {
      responsive: true,
      modeBarButtonsToRemove: ["select2d", "lasso2d"],
    }
  );
}

function get_binary_vle_phase_diagram_txy(vle_phase_diagram_data) {
  var trace1 = {
    x: vle_phase_diagram_data["x0"],
    y: vle_phase_diagram_data["temperature"],
    mode: "lines",
    type: "scatter",
    name: "Liquid phase",
  };
  var trace2 = {
    x: vle_phase_diagram_data["y0"],
    y: vle_phase_diagram_data["temperature"],
    mode: "lines",
    type: "scatter",
    name: "Vapor phase",
  };
  Plotly.newPlot(
    "vle_phase_diagram_txy",
    [trace1, trace2],
    get_layout(0, "x<sub>1</sub>", "Temperature (K)"),
    {
      responsive: true,
      modeBarButtonsToRemove: ["select2d", "lasso2d"],
    }
  );
}

function get_binary_vle_phase_diagram_xy(vle_phase_diagram_data) {
  var trace1 = {
    x: vle_phase_diagram_data["x0"],
    y: vle_phase_diagram_data["y0"],
    mode: "lines",
    type: "scatter",
    name: "VLE",
  };
  Plotly.newPlot(
    "vle_phase_diagram_xy",
    [trace1],
    get_layout(0, "x<sub>1</sub>", "y<sub>1</sub>"),
    {
      responsive: true,
      modeBarButtonsToRemove: ["select2d", "lasso2d"],
    }
  );
}
