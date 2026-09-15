document.addEventListener("DOMContentLoaded", function () {
  const copyableParameters = document.querySelectorAll(".copyable-parameter");

  async function copyParameter(parameterCell) {
    const value = parameterCell.dataset.copyValue;
    try {
      if (navigator.clipboard && window.isSecureContext) {
        await navigator.clipboard.writeText(value);
      } else {
        throw new Error("Modern clipboard API unavailable");
      }
    } catch (error) {
      const copyInput = document.createElement("textarea");
      copyInput.value = value;
      copyInput.style.position = "fixed";
      copyInput.style.opacity = "0";
      document.body.appendChild(copyInput);
      copyInput.focus();
      copyInput.select();
      document.execCommand("copy");
      copyInput.remove();
    }
    parameterCell.title = "Copied";
    window.setTimeout(() => {
      parameterCell.title = "Copy parameter";
    }, 1200);
  }

  copyableParameters.forEach((parameterCell) => {
    parameterCell.addEventListener("click", () => copyParameter(parameterCell));
    parameterCell.addEventListener("keydown", (event) => {
      if (event.key === "Enter" || event.key === " ") {
        event.preventDefault();
        copyParameter(parameterCell);
      }
    });
  });

  const pressureButtons = document.querySelectorAll(".pressure-btn");
  const pressureInput = document.getElementById("id_pressure");
  pressureButtons.forEach((button) => {
    button.addEventListener("click", function () {
      if (pressureInput) {
        pressureInput.value = this.getAttribute("data-pressure") * 1000;
        pressureInput.dispatchEvent(new Event("change"));
      }
    });
  });

  const mixturePressureInput = document.getElementById("id_pressure");
  const tempMinInput = document.getElementById("id_temp_min");
  const moleFractionsInput = document.getElementById("id_mole_fractions");
  const kijInput = document.getElementById("id_kijs");
  const rhoBinaryButtons = document.querySelectorAll(".rho-px-data-b-btn");
  const rhoTernaryButtons = document.querySelectorAll(".rho-px-data-t-btn");
  const lleBinaryButtons = document.querySelectorAll(".lle-p-data-b-btn");
  const lleTernaryButtons = document.querySelectorAll(".lle-pt-data-t-btn");
  const vlePressureButtons = document.querySelectorAll(".vle-p-data-b-btn");
  const vleTemperatureButtons = document.querySelectorAll(".vle-t-data-b-btn");
  const vleCompositionButtons = document.querySelectorAll(".vle-x-data-b-btn");
  const vlePressureTemperatureButtons =
    document.querySelectorAll(".vle-pt-data-t-btn");
  const vleTemperatureCompositionButtons =
    document.querySelectorAll(".vle-tx-data-t-btn");
  const experimentalDataButtons = [
    ...rhoBinaryButtons,
    ...rhoTernaryButtons,
    ...lleBinaryButtons,
    ...lleTernaryButtons,
    ...vlePressureButtons,
    ...vleTemperatureButtons,
    ...vleCompositionButtons,
    ...vlePressureTemperatureButtons,
    ...vleTemperatureCompositionButtons,
  ];

  function resetButtonStyle(button) {
    button.classList.remove("btn-success");
    button.classList.add("btn-outline-primary");
  }

  function fillFormAndHighlight(
    button,
    pressureValue,
    tempValue,
    moleFractions = null,
  ) {
    experimentalDataButtons.forEach(resetButtonStyle);

    if (pressureValue !== null && mixturePressureInput) {
      mixturePressureInput.value = pressureValue * 1000;
      mixturePressureInput.dispatchEvent(new Event("change"));
    }
    if (tempValue !== null && tempMinInput) {
      tempMinInput.value = tempValue;
      tempMinInput.dispatchEvent(new Event("change"));
    }
    if (moleFractions !== null && moleFractionsInput) {
      moleFractionsInput.value = moleFractions
        .map((value) => value.toFixed(4))
        .join(" ");
      moleFractionsInput.dispatchEvent(new Event("input"));
      moleFractionsInput.dispatchEvent(new Event("change"));
    }

    button.classList.add("btn-success");
    button.classList.remove("btn-outline-primary");
  }

  rhoBinaryButtons.forEach((button) => {
    button.addEventListener("click", function (event) {
      event.preventDefault();
      const pressure = parseFloat(this.getAttribute("data-pressure"));
      const fraction = parseFloat(this.getAttribute("data-x"));
      fillFormAndHighlight(this, pressure, null, [fraction, 1 - fraction]);
    });
  });
  rhoTernaryButtons.forEach((button) => {
    button.addEventListener("click", function (event) {
      event.preventDefault();
      const pressure = parseFloat(this.getAttribute("data-pressure"));
      const firstFraction = parseFloat(this.getAttribute("data-x1"));
      const secondFraction = parseFloat(this.getAttribute("data-x2"));
      fillFormAndHighlight(this, pressure, null, [
        firstFraction,
        secondFraction,
        1 - firstFraction - secondFraction,
      ]);
    });
  });
  lleBinaryButtons.forEach((button) => {
    button.addEventListener("click", function (event) {
      event.preventDefault();
      fillFormAndHighlight(
        this,
        parseFloat(this.getAttribute("data-pressure")),
        null,
      );
    });
  });
  lleTernaryButtons.forEach((button) => {
    button.addEventListener("click", function (event) {
      event.preventDefault();
      fillFormAndHighlight(
        this,
        parseFloat(this.getAttribute("data-pressure")),
        parseFloat(this.getAttribute("data-tk")),
      );
    });
  });
  vlePressureButtons.forEach((button) => {
    button.addEventListener("click", function (event) {
      event.preventDefault();
      fillFormAndHighlight(
        this,
        parseFloat(this.getAttribute("data-pressure")),
        null,
      );
    });
  });
  vleTemperatureButtons.forEach((button) => {
    button.addEventListener("click", function (event) {
      event.preventDefault();
      fillFormAndHighlight(
        this,
        null,
        parseFloat(this.getAttribute("data-tk")),
      );
    });
  });
  vleCompositionButtons.forEach((button) => {
    button.addEventListener("click", function (event) {
      event.preventDefault();
      const fraction = parseFloat(this.getAttribute("data-x"));
      fillFormAndHighlight(this, null, null, [fraction, 1 - fraction]);
    });
  });
  vlePressureTemperatureButtons.forEach((button) => {
    button.addEventListener("click", function (event) {
      event.preventDefault();
      fillFormAndHighlight(
        this,
        parseFloat(this.getAttribute("data-pressure")),
        parseFloat(this.getAttribute("data-tk")),
      );
    });
  });
  vleTemperatureCompositionButtons.forEach((button) => {
    button.addEventListener("click", function (event) {
      event.preventDefault();
      const ratio = parseFloat(this.getAttribute("data-ratio"));
      fillFormAndHighlight(
        this,
        null,
        parseFloat(this.getAttribute("data-tk")),
        [0, ratio, 1 - ratio],
      );
    });
  });

  if (kijInput) {
    const mixtureConfig = document.getElementById("mixture-config");
    const kijValue = mixtureConfig?.dataset.initialKijValue;
    if (kijValue) {
      kijInput.value = parseFloat(kijValue).toFixed(4);
      kijInput.dispatchEvent(new Event("input"));
      kijInput.dispatchEvent(new Event("change"));
    }
  }
});
