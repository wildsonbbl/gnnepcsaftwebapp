const { app, BrowserWindow, shell } = require("electron");
const { spawn } = require("child_process");
const controller = new AbortController();
const { signal } = controller;
const path = require("path");
const log = require("electron-log/main");
const fetch = require("node-fetch");
// Optional, initialize the logger for any renderer process
log.initialize();

if (require("electron-squirrel-startup")) app.quit();

const createWindow = async () => {
  const djangoBackend = startDjangoServer();

  const win = new BrowserWindow({
    width: 768,
    height: 600,
    titleBarStyle: "default",
  });

  win.menuBarVisible = false;

  win.loadFile(path.join(__dirname, "index.html")); //from loading.io

  await waitForDjangoServer();

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
  let appPath;
  if (process.platform === "win32") {
    appPath = path.join(
      process.resourcesPath,
      "gnnepcsaftwebapp/gnnepcsaftwebapp.exe"
    );
  } else {
    appPath = path.join(
      process.resourcesPath,
      "gnnepcsaftwebapp/gnnepcsaftwebapp"
    );
  }
  const djangoBackend = spawn(
    appPath,
    ["runserver", "--noreload", "--skip-checks", "localhost:19770"],
    { signal }
  );

  djangoBackend.stdout.on("data", (data) => {
    log.info(`stdout:\n${data}`);
  });
  djangoBackend.stderr.on("data", (data) => {
    log.info(`stderr:\n${data}`);
  });
  djangoBackend.on("error", (error) => {
    log.error(`error:\n${error.message}`);
  });
  djangoBackend.on("close", (code) => {
    log.info(`child process exited with code ${code}`);
    app.quit();
  });
  djangoBackend.on("message", (message) => {
    log.info(`message:\n${message}`);
  });
  return djangoBackend;
};

const waitForDjangoServer = () => {
  return new Promise((resolve, reject) => {
    const interval = setInterval(() => {
      fetch("http://localhost:19770")
        .then((response) => {
          if (response.status === 200) {
            clearInterval(interval);
            log.info("Django server is running");
            resolve();
          }
        })
        .catch((error) => {
          log.info("Django server is not running");
          log.error(error.message);
        });
    }, 1000);
  });
};
