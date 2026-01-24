const { FusesPlugin } = require("@electron-forge/plugin-fuses");
const { FuseV1Options, FuseVersion } = require("@electron/fuses");

module.exports = {
  packagerConfig: {
    asar: true,
    name: "gnnpcsaftwebapp",
    icon: "./storelogo_scale_400_hf3_icon.ico",
    extraResource: ["../app_pkg/dist/gnnpcsaftwebapp"],
  },
  rebuildConfig: {},
  makers: [
    {
      name: "@electron-forge/maker-squirrel",
      config: {
        setupIcon: "./storelogo_scale_400_hf3_icon.ico",
        iconUrl:
          "https://raw.githubusercontent.com/wildsonbbl/gnnepcsaftwebapp/refs/heads/main/electron/storelogo_scale_400_hf3_icon.ico",
      },
    },
    {
      name: "@electron-forge/maker-deb",
      config: {
        options: {
          maintainer: "Wildson Lima",
          homepage: "https://github.com/wildsonbbl/gnnepcsaftwebapp",
          icon: "./storelogo_scale_400_hf3_icon.ico",
        },
      },
    },
    {
      name: "@electron-forge/maker-dmg",
      config: {
        format: "ULFO",
        icon: "./storelogo_scale_400_hf3_icon.ico",
      },
    },
  ],
  publishers: [
    {
      name: "@electron-forge/publisher-github",
      config: {
        repository: {
          owner: "wildsonbbl", //"github-user-name",
          name: "gnnepcsaftwebapp", //"github-repo-name",
        },
        prerelease: false,
        draft: true,
        generateReleaseNotes: true,
      },
    },
  ],
  plugins: [
    {
      name: "@electron-forge/plugin-auto-unpack-natives",
      config: {},
    },
    // Fuses are used to enable/disable various Electron functionality
    // at package time, before code signing the application
    new FusesPlugin({
      version: FuseVersion.V1,
      [FuseV1Options.RunAsNode]: false,
      [FuseV1Options.EnableCookieEncryption]: true,
      [FuseV1Options.EnableNodeOptionsEnvironmentVariable]: false,
      [FuseV1Options.EnableNodeCliInspectArguments]: false,
      [FuseV1Options.EnableEmbeddedAsarIntegrityValidation]: true,
      [FuseV1Options.OnlyLoadAppFromAsar]: true,
    }),
  ],
};
