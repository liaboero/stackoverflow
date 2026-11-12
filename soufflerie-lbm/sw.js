/* Service worker : application utilisable hors ligne.
   Ressources locales : « stale-while-revalidate » (réponse immédiate depuis le
   cache, mise à jour en arrière-plan pour la prochaine ouverture). */
const VERSION = 'soufflerie-v5';
const FONTS = 'soufflerie-fonts';
const ASSETS = [
  './',
  './index.html',
  './css/app.css',
  './js/app.js',
  './js/geometry.js',
  './js/lbm.js',
  './js/lbm-gpu.js',
  './js/render-gpu.js',
  './js/shapes.js',
  './js/sim.js',
  './js/stats.js',
  './js/render.js',
  './js/chart.js',
  './js/panel.js',
  './js/wall.js',
  './js/vlm.js',
  './js/wing-tab.js',
  './js/lbm3d.js',
  './js/bodies3d.js',
  './js/render3d.js',
  './js/lab3d.js',
  './manifest.webmanifest',
  './icons/icon.svg',
  './icons/icon-192.png',
  './icons/icon-512.png',
  './icons/apple-touch-icon.png',
];

self.addEventListener('install', (event) => {
  event.waitUntil(caches.open(VERSION).then((c) => c.addAll(ASSETS)).then(() => self.skipWaiting()));
});

self.addEventListener('activate', (event) => {
  event.waitUntil(
    caches
      .keys()
      .then((keys) => Promise.all(keys.filter((k) => k !== VERSION && k !== FONTS).map((k) => caches.delete(k))))
      .then(() => self.clients.claim()),
  );
});

function staleWhileRevalidate(event, cacheName) {
  return caches.open(cacheName).then((cache) =>
    cache.match(event.request, { ignoreSearch: true }).then((cached) => {
      const network = fetch(event.request)
        .then((res) => {
          if (res && (res.ok || res.type === 'opaque')) cache.put(event.request, res.clone());
          return res;
        })
        .catch(() => cached);
      if (cached) {
        event.waitUntil(network);
        return cached;
      }
      return network;
    }),
  );
}

self.addEventListener('fetch', (event) => {
  const req = event.request;
  if (req.method !== 'GET') return;
  const url = new URL(req.url);
  if (url.origin === self.location.origin) {
    event.respondWith(staleWhileRevalidate(event, VERSION));
  } else if (url.hostname === 'fonts.googleapis.com' || url.hostname === 'fonts.gstatic.com') {
    event.respondWith(staleWhileRevalidate(event, FONTS));
  }
});
