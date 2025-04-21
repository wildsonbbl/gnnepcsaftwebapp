const { app, BrowserWindow, shell, Menu } = require("electron");
const { spawn } = require("child_process");
const controller = new AbortController();
const { signal } = controller;
const path = require("path");
const log = require("electron-log/main");
const fetch = require("node-fetch");
const fs = require("fs");

// Optional, initialize the logger for any renderer process
log.initialize();

if (require("electron-squirrel-startup")) app.quit();
Menu.setApplicationMenu(null); // Hide the menu bar

const createWindow = async () => {
  // Set up user data directory for database
  const userDataPath = app.getPath("userData");
  const dbDir = path.join(userDataPath, "db");

  // Create directory if it doesn't exist
  if (!fs.existsSync(dbDir)) {
    fs.mkdirSync(dbDir, { recursive: true });
  }

  // Check if we need to copy the initial database
  copyDB(dbDir, "gnnepcsaft.db");
  copyDB(dbDir, "gnnepcsaft.chat.db");

  // Start Django with the user database path
  const djangoBackend = startDjangoServer(dbDir);

  const win = new BrowserWindow({
    width: 1280,
    height: 760,
    minWidth: 600,
    minHeight: 680,
    titleBarStyle: "default",
    title: "GNNePCSAFT",
  });

  win.loadFile(path.join(__dirname, "index.html"));

  await waitForDjangoServer();

  win.loadURL("http://localhost:19770");
  win.maximize();

  win.webContents.setWindowOpenHandler(({ url }) => {
    if (url.startsWith("http://localhost:19770")) {
      return {
        action: "allow",
        overrideBrowserWindowOptions: {
          fullscreen: false,
          width: 600,
          height: 680,
          minWidth: 600,
          minHeight: 680,
          title: "GNNePCSAFT",
        },
      };
    }
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

const startDjangoServer = (dbDir) => {
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

  // Pass the database path as an environment variable
  const dbPath = path.join(dbDir, "gnnepcsaft.db");
  const dbChatPath = path.join(dbDir, "gnnepcsaft.chat.db");
  const env = {
    ...process.env,
    GNNEPCSAFT_DB_PATH: dbPath,
    GNNEPCSAFT_DB_CHAT_PATH: dbChatPath,
  };

  const djangoBackend = spawn(
    appPath,
    ["runserver", "--noreload", "--skip-checks", "localhost:19770"],
    { signal, env }
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

function copyDB(dbDir, dbName) {
  // Path to the database in user data directory
  const userDbPath = path.join(dbDir, dbName);

  log.info(`User database path: ${userDbPath}`);

  // Check if the database already exists in user data directory
  if (!fs.existsSync(userDbPath)) {
    const resourceDbPath = path.join(
      process.resourcesPath,
      "gnnepcsaftwebapp/_internal",
      dbName
    );

    // Copy the database if it exists in resources
    if (fs.existsSync(resourceDbPath)) {
      fs.copyFileSync(resourceDbPath, userDbPath);
      log.info(`Copied initial database to: ${userDbPath}`);
    }
  }
}
