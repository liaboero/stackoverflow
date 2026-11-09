/**
 * Construit une version « fichier unique » de l'application (CSS et JS en
 * ligne), utilisable hors ligne en ouvrant simplement le fichier.
 *
 *   node tools/build.mjs                    → dist/soufflerie-de-poche.html
 *   node tools/build.mjs --fragment <out>   → fragment sans <html>/<head>/<body>
 *                                             (hébergement qui fournit le squelette)
 *
 * Mini-empaqueteur adapté au style du projet : imports nommés sur une ligne,
 * exports par `export function|class|const|let`.
 */
import { readFileSync, writeFileSync, mkdirSync } from 'node:fs';
import { dirname, join } from 'node:path';
import { fileURLToPath } from 'node:url';

const root = join(dirname(fileURLToPath(import.meta.url)), '..');
const read = (p) => readFileSync(join(root, p), 'utf8');

// Ordre de dépendance.
const MODULES = ['geometry', 'lbm', 'shapes', 'stats', 'sim', 'lbm-gpu', 'render', 'render-gpu', 'chart', 'app'];
const ident = (name) => name.replace(/-/g, '_');

function bundle() {
  const parts = ['window.__SOUFFLERIE_INLINE__ = true;'];
  for (const name of MODULES) {
    let src = read(`js/${name}.js`);
    const exportsList = [];
    src = src.replace(/^export (function|class|const|let) ([A-Za-z_$][\w$]*)/gm, (_, kind, id) => {
      exportsList.push(id);
      return `${kind} ${id}`;
    });
    // Réexportations : export { a, b };
    src = src.replace(/^export \{([^}]+)\};$/gm, (_, names) => {
      exportsList.push(...names.split(',').map((n) => n.trim()).filter(Boolean));
      return '';
    });
    src = src.replace(/^import \{([^}]+)\} from '\.\/([\w-]+)\.js';$/gm, (_, names, mod) => {
      if (!MODULES.includes(mod)) throw new Error(`Module inconnu : ${mod}`);
      return `const {${names}} = __mod_${ident(mod)};`;
    });
    if (/^\s*(import|export)\s/m.test(src)) throw new Error(`Import/export non pris en charge dans ${name}.js`);
    parts.push(`const __mod_${ident(name)} = (() => {\n${src}\nreturn { ${exportsList.join(', ')} };\n})();`);
  }
  return parts.join('\n');
}

function between(html, a, b) {
  const i = html.indexOf(a);
  const j = html.indexOf(b);
  if (i < 0 || j < 0) throw new Error(`Repères ${a} / ${b} introuvables`);
  return html.slice(i + a.length, j).trim();
}

const html = read('index.html');
const head = between(html, '<!-- HEAD START -->', '<!-- HEAD END -->');
const body = between(html, '<!-- BODY START -->', '<!-- BODY END -->');
const css = read('css/app.css');
const js = bundle().replace(/<\/script/gi, '<\\/script');
const icon = `data:image/svg+xml,${encodeURIComponent(read('icons/icon.svg').trim())}`;

const args = process.argv.slice(2);
if (args[0] === '--fragment') {
  const out = args[1];
  if (!out) throw new Error('Chemin de sortie manquant');
  const title = head.match(/<title>[^<]*<\/title>/)[0];
  const fonts = head.match(/<link rel="stylesheet" href="https:\/\/fonts[^>]+>/)[0];
  const frag = `${title}
<meta name="claude:viewport-fit-sides" content="cover">
${fonts}
<style>
${css}
</style>
${body}
<script>
${js}
</script>
`;
  mkdirSync(dirname(out), { recursive: true });
  writeFileSync(out, frag);
  console.log(`Fragment écrit : ${out} (${(frag.length / 1024).toFixed(0)} ko)`);
} else {
  const doc = `<!doctype html>
<html lang="fr">
<head>
<meta charset="utf-8">
${head}
<link rel="icon" href="${icon}">
<style>
${css}
</style>
</head>
<body>
${body}
<script>
${js}
</script>
</body>
</html>
`;
  const out = join(root, 'dist', 'soufflerie-de-poche.html');
  mkdirSync(dirname(out), { recursive: true });
  writeFileSync(out, doc);
  console.log(`Version autonome écrite : ${out} (${(doc.length / 1024).toFixed(0)} ko)`);
}
