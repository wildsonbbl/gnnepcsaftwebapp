// Base Service Worker implementation.  To use your own Service Worker, set the PWA_SERVICE_WORKER_PATH variable in settings.py

var staticCacheName = "django-pwa-v" + new Date().getTime();
var filesToCache = [
  "/offline/",
  "/static/css/django-pwa-app.css",
  "/static/fontawesomefree/css/brands.css",
  "/static/fontawesomefree/css/fontawesome.css",
  "/static/fontawesomefree/webfonts/fa-brands-400.woff2",
  "/static/fontawesomefree/webfonts/fa-brands-400.ttf",
  "/static/js/myidb.js",
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
          .filter((cacheName) => cacheName.startsWith("django-pwa-"))
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
