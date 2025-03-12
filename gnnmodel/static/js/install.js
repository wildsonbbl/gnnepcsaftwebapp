let deferredPrompt;
window.addEventListener("beforeinstallprompt", (e) => {
  deferredPrompt = e;
  if (document.getElementById("installApp") !== null) {
    return;
  }
  const button = document.createElement("button");
  button.type = "button";
  button.id = "installApp";
  button.ariaLabel = "Install App Button";
  button.innerHTML = "Install";
  button.className = "btn btn-outline-success text-light";
  button.addEventListener("click", async () => {
    if (deferredPrompt !== null) {
      deferredPrompt.prompt();
      const { outcome } = await deferredPrompt.userChoice;
      if (outcome === "accepted") {
        deferredPrompt = null;
      }
    }
  });
  const installButton = document.getElementById("installButton");
  installButton.appendChild(button);
});
