let deferredPrompt;
window.addEventListener("beforeinstallprompt", (e) => {
  deferredPrompt = e;
  const installApp = document.getElementById("installApp");
  installApp.innerHTML = "Install";
  installApp.className = "btn btn-outline-success text-light";
});

const installApp = document.getElementById("installApp");
installApp.addEventListener("click", async () => {
  if (deferredPrompt !== null) {
    deferredPrompt.prompt();
    const { outcome } = await deferredPrompt.userChoice;
    if (outcome === "accepted") {
      deferredPrompt = null;
    }
  }
});
