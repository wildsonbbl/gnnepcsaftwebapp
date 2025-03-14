let element = document.querySelector("#molplot");
let config = { backgroundColor: "#f8f9fa" };
let v = $3Dmol.createViewer(element, config);

function loadmol(data) {
  v.addModel(data, "sdf");
  v.setStyle(
    {},
    {
      stick: { color: "spectrum", radius: 0.2 },
      sphere: { color: "spectrum", scale: 0.333 },
    }
  );
  v.render();
  v.zoomTo();
}

function reset() {
  v.removeAllLabels();
  v.setStyle(
    {},
    {
      stick: { color: "spectrum", radius: 0.2 },
      sphere: { color: "spectrum", scale: 0.333 },
    }
  );
  v.render();
}

function addlabel() {
  v.addPropertyLabels(
    "atom",
    { not: { elem: "H" } },
    {
      fontColor: " black ",
      font: " sans-serif ",
      fontSize: 28,
      showBackground: false,
      alignment: " center ",
    }
  );
  v.render();
}
