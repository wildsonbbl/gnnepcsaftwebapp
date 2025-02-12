const { app, BrowserWindow, shell } = require("electron");
const { spawn } = require("child_process");
const controller = new AbortController();
const { signal } = controller;

if (require("electron-squirrel-startup")) app.quit();

const createWindow = async () => {
  const djangoBackend = startDjangoServer();
  await waitForDjangoServer(djangoBackend);

  const win = new BrowserWindow({
    width: 600,
    height: 600,
    titleBarStyle: "default",
  });

  win.menuBarVisible = false;

  win.loadURL("http://localhost:19770");

  win.webContents.setWindowOpenHandler(({ url }) => {
    shell.openExternal(url);
    return { action: "deny" };
  });
};

app.whenReady().then(() => {
  createWindow();

  app.on("activate", () => {
    if (BrowserWindow.getAllWindows().length === 0) {
      createWindow();
    }
  });
});

app.on("window-all-closed", () => {
  if (process.platform !== "darwin") {
    app.quit();
    controller.abort();
  }
});

const startDjangoServer = () => {
  const djangoBackend = spawn(
    `.\\gnnepcsaftwebapp\\gnnepcsaftwebapp.exe`,
    ["runserver", "--noreload", "localhost:19770"],
    { signal }
  );

  djangoBackend.stdout.on("data", (data) => {
    console.log(`stdout:\n${data}`);
  });
  djangoBackend.stderr.on("data", (data) => {
    console.log(`stderr:\n${data}`);
  });
  djangoBackend.on("error", (error) => {
    console.log(`error:\n${error.message}`);
  });
  djangoBackend.on("close", (code) => {
    console.log(`child process exited with code ${code}`);
    app.quit();
  });
  djangoBackend.on("message", (message) => {
    console.log(`message:\n${message}`);
  });
  return djangoBackend;
};

const waitForDjangoServer = (djangoBackend) => {
  return new Promise((resolve, reject) => {
    const onData = (data) => {
      const text = data.toString();
      console.log("Waiting for Django server...\n", text);
      if (text.includes("Starting development server")) {
        // Once we detect the server is running, remove this listener.
        djangoBackend.stdout.off("data", onData);
        resolve();
      }
    };
    djangoBackend.stdout.on("data", onData);
    // Optionally add a timeout or error handling here if needed.
  });
};
