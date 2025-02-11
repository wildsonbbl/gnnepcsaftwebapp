let submitQueryPred = document.getElementById("submitQueryPred");

submitQueryPred.addEventListener("submit", (event) => {
  let loadingElem = document.getElementById("loadingIndicator");
  if (loadingElem) {
    loadingElem.style.display = "flex";
  }
});
