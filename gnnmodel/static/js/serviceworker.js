// Service Worker for PWA
// This service worker is used to cache files for offline access in a PWA application.

var staticCacheName = "pwa-v" + new Date().getTime();
var filesToCache = [
  "/offline/",
  "/static/fontawesomefree/css/brands.css",
  "/static/fontawesomefree/css/fontawesome.css",
  "/static/fontawesomefree/css/solid.css",
  "/static/css/styles.css",
  "/static/css/chat.css",
  "/static/fontawesomefree/webfonts/fa-brands-400.woff2",
  "/static/fontawesomefree/webfonts/fa-brands-400.ttf",
  "/static/fontawesomefree/webfonts/fa-solid-900.woff2",
  "/static/fontawesomefree/webfonts/fa-solid-900.ttf",
  "/static/js/myidb.js",
  "/static/js/queries.js",
  "/static/js/install.js",
  "/static/images/icons/android/android-launchericon-48-48.png",
  "/static/images/icons/android/android-launchericon-72-72.png",
  "/static/images/icons/android/android-launchericon-96-96.png",
  "/static/images/icons/android/android-launchericon-144-144.png",
  "/static/images/icons/android/android-launchericon-192-192.png",
  "/static/images/icons/android/android-launchericon-512-512.png",
  "/static/images/icons/windows11/SplashScreen.scale-100.png",
  "/static/images/icons/windows11/SplashScreen.scale-125.png",
  "/static/images/icons/windows11/SplashScreen.scale-150.png",
  "/static/images/icons/windows11/SplashScreen.scale-200.png",
  "/static/images/icons/windows11/SplashScreen.scale-400.png",
  "/static/images/icons/ios/180.png",
  "/static/manifest.json",
];

// Cache on install
self.addEventListener("install", (event) => {
  this.skipWaiting();
  event.waitUntil(
    caches.open(staticCacheName).then((cache) => {
      return cache.addAll(filesToCache);
    })
  );
});

// Clear cache on activate
self.addEventListener("activate", (event) => {
  event.waitUntil(
    caches.keys().then((cacheNames) => {
      return Promise.all(
        cacheNames
          .filter((cacheName) => cacheName.startsWith("pwa-"))
          .filter((cacheName) => cacheName !== staticCacheName)
          .map((cacheName) => caches.delete(cacheName))
      );
    })
  );
});

// Serve from Cache
self.addEventListener("fetch", (event) => {
  event.respondWith(
    caches
      .match(event.request)
      .then((response) => {
        return response || fetch(event.request);
      })
      .catch(() => {
        return caches.match("/offline/");
      })
  );
});
