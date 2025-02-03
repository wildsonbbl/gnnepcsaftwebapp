let submitQuery = document.getElementById("submitQuery");

submitQuery.addEventListener("submit", (event) => {
  event.preventDefault();
  let query = document.getElementById("id_query");
  const inchip = getinchi(query.value);
  const smilesp = getsmiles(query.value);
  smilesp
    .then((valor) => {
      if (!valor) {
        inchip
          .then((valor) => {
            if (!valor) {
              document.getElementById("output").innerHTML =
                "InChI/SMILES not available in local database";
              document.getElementById("output").className =
                "text-center alert alert-danger";
              document.getElementById("m").innerHTML = "";
              document.getElementById("sigma").innerHTML = "";
              document.getElementById("e").innerHTML = "";
              document.getElementById("k_ab").innerHTML = "";
              document.getElementById("e_ab").innerHTML = "";
              document.getElementById("mu").innerHTML = "";
              document.getElementById("na").innerHTML = "";
              document.getElementById("nb").innerHTML = "";
              document.getElementById("inchi").innerHTML = "";
              document.getElementById("smiles").innerHTML = "";
            } else {
              document.getElementById("output").innerHTML = "";
              document.getElementById("output").className = "";

              for (prop in valor) {
                document.getElementById(prop).innerHTML = valor[prop];
              }
            }
          })
          .catch((err) => {
            console.log(err);
          });
      } else {
        document.getElementById("output").innerHTML = "";
        document.getElementById("output").className = "";

        for (prop in valor) {
          document.getElementById(prop).innerHTML = valor[prop];
        }
      }
    })
    .catch((err) => {
      console.log(err);
    });
});
