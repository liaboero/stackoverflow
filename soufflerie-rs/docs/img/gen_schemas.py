"""Schémas et graphiques de la documentation (docs/img/*.svg).

python docs/img/gen_schemas.py  (depuis soufflerie-rs/)

SVG autonomes, thème clair et sombre (prefers-color-scheme), palette de
référence : bleu, orange, aqua. Les données des graphiques viennent du
README (sections « Validation 3D » et « Décrochage »).
"""
import math
import os

OUT = os.path.dirname(os.path.abspath(__file__))

STYLE = """<style>
.bg{fill:#fcfcfb}
.t{fill:#0b0b0b;font:14px system-ui,-apple-system,'Segoe UI',sans-serif}
.tb{fill:#0b0b0b;font:600 15px system-ui,-apple-system,'Segoe UI',sans-serif}
.t2{fill:#52514e;font:12px system-ui,-apple-system,'Segoe UI',sans-serif}
.tm{fill:#52514e;font:italic 13px Georgia,'Times New Roman',serif}
.ln{stroke:#52514e;fill:none}
.lnf{fill:#52514e}
.mu{stroke:#898781;fill:none}
.grid{stroke:#e1e0d9;fill:none}
.body{fill:#dedcd5;stroke:#52514e}
.b{stroke:#2a78d6;fill:none}.bf{fill:#2a78d6}.bl{fill:#2a78d6;fill-opacity:.14}
.o{stroke:#eb6834;fill:none}.of{fill:#eb6834}.ol{fill:#eb6834;fill-opacity:.16}
.a{stroke:#1baf7a;fill:none}.af{fill:#1baf7a}.al{fill:#1baf7a;fill-opacity:.18}
.sf{fill:#fcfcfb}
@media (prefers-color-scheme:dark){
.bg{fill:#1a1a19}.sf{fill:#1a1a19}
.t,.tb{fill:#f0efec}.t2,.tm{fill:#c3c2b7}
.ln{stroke:#c3c2b7}.lnf{fill:#c3c2b7}.grid{stroke:#2c2c2a}
.body{fill:#3f3e3b;stroke:#c3c2b7}
.b{stroke:#3987e5}.bf{fill:#3987e5}.bl{fill:#3987e5}
.o{stroke:#d95926}.of{fill:#d95926}.ol{fill:#d95926}
.a{stroke:#199e70}.af{fill:#199e70}.al{fill:#199e70}
}
</style>"""


def markers():
    m = ""
    for c in ["ln", "b", "o", "a"]:
        cls = {"ln": "lnf", "b": "bf", "o": "of", "a": "af"}[c]
        m += (f'<marker id="m{c}" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" '
              f'markerHeight="7" orient="auto-start-reverse"><path d="M0,1 L10,5 L0,9 z" class="{cls}"/></marker>')
    m += ('<pattern id="hach" width="8" height="8" patternUnits="userSpaceOnUse" patternTransform="rotate(45)">'
          '<line x1="0" y1="0" x2="0" y2="8" class="ln" stroke-width="1.2"/></pattern>')
    return f"<defs>{m}</defs>"


def svg(name, w, h, body, title):
    s = (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {h}" width="{w}" height="{h}" '
         f'role="img" aria-label="{title}"><title>{title}</title>{STYLE}{markers()}'
         f'<rect class="bg" x="0" y="0" width="{w}" height="{h}" rx="10"/>{body}</svg>\n')
    with open(os.path.join(OUT, name), "w", encoding="utf-8") as f:
        f.write(s)


def T(x, y, s, cls="t", anchor="start", extra=""):
    s = s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
    return f'<text x="{x:.1f}" y="{y:.1f}" class="{cls}" text-anchor="{anchor}" {extra}>{s}</text>'


def L(x1, y1, x2, y2, cls="ln", w=1.5, mk=None, dash=None):
    m = f' marker-end="url(#m{mk})"' if mk else ""
    d = f' stroke-dasharray="{dash}"' if dash else ""
    return f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" class="{cls}" stroke-width="{w}"{m}{d}/>'


def C(x, y, r, cls, extra=""):
    return f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{r}" class="{cls}" {extra}/>'


def P(pts, cls, w=1.5, close=False, extra=""):
    d = "M" + " L".join(f"{x:.1f},{y:.1f}" for x, y in pts) + (" Z" if close else "")
    return f'<path d="{d}" class="{cls}" stroke-width="{w}" {extra}/>'


# ---------------------------------------------------------------- D3Q19
def d3q19():
    cx, cy, s = 190, 195, 100
    ca, sa = math.cos(math.radians(35)), math.sin(math.radians(35))

    def pr(x, y, z):
        # perspective cavalière : y en profondeur, vers le haut à droite
        return cx + s * (x + 0.5 * y * ca), cy - s * (z + 0.5 * y * sa)

    b = ""
    # cube [-1,1]^3
    corners = [(x, y, z) for x in (-1, 1) for y in (-1, 1) for z in (-1, 1)]
    for p in corners:
        for q in corners:
            if p < q and sum(a != c for a, c in zip(p, q)) == 1:
                b += L(*pr(*p), *pr(*q), cls="grid", w=1.2)
    vel = [(x, y, z) for x in (-1, 0, 1) for y in (-1, 0, 1) for z in (-1, 0, 1)
           if 0 < abs(x) + abs(y) + abs(z) <= 2]
    # arêtes d'abord (derrière), faces ensuite
    for v in sorted(vel, key=lambda v: -(abs(v[0]) + abs(v[1]) + abs(v[2]))):
        n = abs(v[0]) + abs(v[1]) + abs(v[2])
        c = "b" if n == 1 else "o"
        x0, y0 = pr(0, 0, 0)
        x1, y1 = pr(*v)
        k = 0.93
        b += L(x0, y0, x0 + k * (x1 - x0), y0 + k * (y1 - y0), cls=c, w=2 if n == 1 else 1.6, mk=c)
    b += C(*pr(0, 0, 0), 6, "af")
    # axes
    ox, oy = 40, 350
    for v, lab in [((0.45, 0, 0), "x"), ((0, 0.45, 0), "y"), ((0, 0, 0.45), "z")]:
        x1, y1 = pr(*v)
        x0, y0 = pr(0, 0, 0)
        b += L(ox, oy, ox + (x1 - x0), oy + (y1 - y0), w=1.3, mk="ln")
        b += T(ox + 1.25 * (x1 - x0) - 4, oy + 1.25 * (y1 - y0) + 4, lab, "t2")
    # légende
    lx = 400
    b += T(lx, 70, "Réseau D3Q19", "tb")
    rows = [("af", "repos", "1 vitesse", "w₀ = 1/3"),
            ("bf", "faces", "6 vitesses, |c| = 1", "w = 1/18"),
            ("of", "arêtes", "12 vitesses, |c| = √2", "w = 1/36")]
    for k, (cls, a, bb, w) in enumerate(rows):
        y = 108 + 38 * k
        b += f'<rect x="{lx}" y="{y-11}" width="14" height="14" rx="3" class="{cls}"/>'
        b += T(lx + 24, y, f"{a} : {bb}", "t") + T(lx + 24, y + 17, w, "t2")
    b += T(lx, 240, "c_s² = 1/3,   τ = ½ + 3ν", "tm")
    b += T(lx, 268, "Pas de vitesse de coin (±1, ±1, ±1) :", "t2")
    b += T(lx, 285, "c_α³ = c_α, d'où le défaut cubique", "t2")
    b += T(lx, 302, "corrigé en option (MACH3_BIT).", "t2")
    svg("d3q19.svg", 680, 380, b, "Les 19 vitesses du réseau D3Q19 et leurs poids")


# ---------------------------------------------------------------- Bouzidi
def bouzidi():
    b = ""
    dx = 120
    rows = [(0.3, 95, "q < ½"), (0.7, 245, "q ≥ ½")]
    for q, y, lab in rows:
        x0 = 120
        xs = [x0, x0 + dx, x0 + 2 * dx]
        xw = xs[1] + q * dx
        # solide
        b += f'<rect x="{xw:.1f}" y="{y-45}" width="{xs[2]+70-xw:.1f}" height="90" fill="url(#hach)" opacity=".55"/>'
        b += L(xw, y - 50, xw, y + 50, w=2.5)
        b += T(xw + 6, y - 52, "paroi", "t2")
        b += L(xs[0] - 30, y, xs[2] + 40, y, cls="grid", w=1.2)
        for k, x in enumerate(xs):
            if k < 2:
                b += C(x, y, 8, "bf")
            else:
                b += C(x, y, 8, "sf", 'stroke-width="2" stroke="currentColor" style="stroke:#898781"')
        b += T(xs[0], y + 30, "x − c_i", "t2", "middle")
        b += T(xs[1], y + 30, "x", "t2", "middle")
        b += T(xs[2], y + 30, "nœud solide", "t2", "middle")
        # q
        b += L(xs[1], y - 26, xw, y - 26, w=1, mk="ln")
        b += L(xw, y - 26, xs[1], y - 26, w=1, mk="ln")
        b += T((xs[1] + xw) / 2, y - 32, "q", "tm", "middle")
        # trajet
        b += L(xs[1] + 4, y + 12, xw - 3, y + 12, cls="o", w=2, mk="o")
        b += L(xw - 3, y + 18, xs[1] + 4, y + 18, cls="b", w=2, mk="b")
        b += T(36, y + 5, lab, "tb")
        if q < 0.5:
            f = "f_ī(x) = 2q f_i*(x) + (1 − 2q) f_i*(x − c_i) + 6w_i c_ī·u_w"
            n = "interpolation avant le rebond, entre x − c_i et x"
        else:
            f = "f_ī(x) = [f_i*(x) + (2q − 1) f_ī*(x) + 6w_i c_ī·u_w] / 2q"
            n = "interpolation après le rebond, entre x et x + c_i"
        b += T(450, y - 6, f, "tm")
        b += T(450, y + 16, n, "t2")
    b += T(450, 330, "u_w : vitesse de glissement de la loi de paroi,", "t2")
    b += T(450, 347, "projetée dans le plan de la face touchée.", "t2")
    b += L(120, 340, 160, 340, cls="o", w=2, mk="o") + T(168, 344, "f_i* sortant (après collision)", "t2")
    b += L(120, 362, 160, 362, cls="b", w=2, mk="b") + T(168, 366, "f_ī reconstruit", "t2")
    svg("bouzidi.svg", 900, 380, b, "Rebond interpolé de Bouzidi selon la fraction q du lien")


# ---------------------------------------------------------------- couche pariétale
def couche():
    X0, Y0, W, H = 90, 330, 380, 280  # u de 0 à 1.1, y de 0 à 4 mailles

    def px(u):
        return X0 + u / 1.1 * W

    def py(y):
        return Y0 - y / 4 * H

    b = ""
    b += f'<rect x="{X0-30}" y="{Y0}" width="{W+40}" height="22" fill="url(#hach)" opacity=".55"/>'
    b += L(X0 - 30, Y0, X0 + W + 10, Y0, w=2.5)
    for k in range(1, 5):
        b += L(X0, py(k), X0 + W, py(k), cls="grid", w=1)
        b += T(X0 - 8, py(k) + 4, f"{k}", "t2", "end")
    b += T(X0 - 8, Y0 + 4, "0", "t2", "end")
    b += L(X0, Y0, X0, py(4.2), w=1.3, mk="ln") + T(X0 - 8, py(4.25) - 4, "y (mailles)", "t2", "end")
    b += L(X0, Y0, px(1.12), Y0, w=0)
    b += T(px(1.08), Y0 + 40, "u", "tm", "middle")
    # profil réel (loi log, couche visqueuse très mince)
    pts = [(px(0), py(0))]
    for i in range(60):
        y = 0.01 * (4 / 0.01) ** (i / 59)
        u = min(1.0, 0.9 + 0.125 * math.log(y / 2))
        pts.append((px(max(u, 0)), py(y)))
    b += P(pts, "mu", 2, extra='stroke-dasharray="5 4"')
    b += T(px(0.97), py(3.2), "profil réel", "t2", "end")
    # couche pariétale linéaire [0, h]
    us, slope = 0.62, 0.12
    b += P([(px(us), py(0)), (px(us + slope * 2), py(2))], "b", 2.5)
    yn = 0.55
    un = us + slope * yn
    uh = 0.9
    b += P([(px(us + slope * 2), py(2)), (px(uh), py(2))], "b", 1, extra='stroke-dasharray="2 3"')
    b += f'<rect x="{px(us)-6:.1f}" y="{py(0)-6:.1f}" width="12" height="12" rx="2" class="bf"/>'
    b += C(px(un), py(yn), 6, "bf")
    b += C(px(uh), py(2), 7, "of")
    # annotations
    ax = 600
    b += T(px(uh) + 14, py(2) - 10, "u_h, relevée à h = 2", "t")
    b += T(px(un) - 12, py(yn) + 5, "nœud pariétal (y_n)", "t", "end")
    b += T(px(us) - 12, py(0) - 8, "glissement u_s", "t", "end")
    b += T(px(us + slope) - 10, py(1) + 4, "couche linéaire [0, h]", "t2", "end")
    b += T(ax, 70, "Loi de paroi (Spalding)", "tb")
    lines = [
        ("1.", "u_τ de Spalding en (h, u_h)"),
        ("2.", "couche [0, h] : ν_w = κ u_τ h,"),
        ("", "pente τ_p / (ν + ν_w)"),
        ("3.", "glissement du nœud :"),
        ("", "u_s = min(u_n, u_h) − pente · y_n,"),
        ("", "plancher u_h − pente · h"),
        ("4.", "relaxation : u_w += 0,1 (u_s − u_w)"),
        ("5.", "Bouzidi en paroi mobile u_w"),
    ]
    y = 104
    for a, s in lines:
        b += T(ax, y, a, "t2") + T(ax + 22, y, s, "t2")
        y += 22
    b += T(ax, y + 12, "La vitesse extrapolée à la paroi ne", "t2")
    b += T(ax, y + 29, "dépend plus de la place du nœud.", "t2")
    svg("couche_parietale.svg", 900, 380, b, "Couche pariétale de la loi de paroi à h = 2 mailles")


# ---------------------------------------------------------------- profil NACA
def naca(n=60):
    up, lo = [], []
    for i in range(n + 1):
        x = (1 - math.cos(math.pi * i / n)) / 2
        t = 0.6 * (0.2969 * math.sqrt(x) - 0.126 * x - 0.3516 * x**2 + 0.2843 * x**3 - 0.1036 * x**4)
        up.append((x, t))
        lo.append((x, -t))
    return up[::-1] + lo[1:]


def raffinement():
    b = ""
    # panneau gauche : régions imbriquées
    ox, oy, s = 150, 200, 190
    a = math.radians(6)

    def tr(x, y):
        xr = (x - 0.25) * math.cos(a) + y * math.sin(a)
        yr = -(x - 0.25) * math.sin(a) + y * math.cos(a)
        return ox + s * xr, oy - s * yr

    prof = [tr(x, y) for x, y in naca()]
    b += f'<rect x="20" y="40" width="420" height="320" class="grid" stroke-width="1.5"/>'
    b += T(30, 60, "veine (L)", "t2")
    b += f'<rect x="55" y="90" width="355" height="220" rx="2" class="b" stroke-width="1.3" stroke-dasharray="6 4"/>'
    b += T(62, 106, "boîte ×2", "t2")
    b += f'<rect x="80" y="125" width="270" height="150" rx="2" class="b" stroke-width="1.6"/>'
    b += T(87, 141, "boîte ×4", "t2")
    d = "M" + " L".join(f"{x:.1f},{y:.1f}" for x, y in prof) + " Z"
    xt, yt = tr(1.0, 0)
    b += f'<path d="M{xt:.1f},{yt:.1f} L{xt+70:.1f},{yt+4:.1f}" class="o" stroke-width="16" stroke-opacity=".35" stroke-linecap="round"/>'
    b += f'<path d="{d}" class="o" stroke-width="30" stroke-opacity=".35" stroke-linejoin="round"/>'
    b += f'<path d="{d}" class="o" stroke-width="12" stroke-opacity=".5" stroke-linejoin="round"/>'
    b += f'<path d="{d}" class="body" stroke-width="1.5"/>'
    b += T(250, 330, "couches creuses ×8, ×16", "t", "middle")
    b += L(250, 318, 220, 222, cls="o", w=1, mk="o")
    b += T(xt + 40, yt - 22, "sillage balayé", "t2", "middle")
    b += T(230, 30, "Grilles imbriquées autour de la maquette", "tb", "middle")
    # panneau droit : raccord
    X0, Y0, g = 500, 80, 44
    b += T(690, 30, "Raccord parente → fine (2D, schéma)", "tb", "middle")
    nx, ny = 9, 6
    xf0 = 3  # début de la grille fine (mailles parentes)
    for i in range(nx):
        for j in range(ny):
            x, y = X0 + i * g, Y0 + j * g
            cls = "af" if i >= xf0 + 2 and 1 <= j <= ny - 2 else "lnf"
            b += C(x, y, 5.5 if cls == "af" else 4, cls)
    # fine
    for i in range(2 * xf0, 2 * (nx - 1) + 1):
        for j in range(0, 2 * (ny - 1) + 1):
            x, y = X0 + i * g / 2, Y0 + j * g / 2
            ghost = i == 2 * xf0 or j == 0 or j == 2 * (ny - 1)
            if ghost:
                b += C(x, y, 3.6, "sf", 'style="stroke:#eb6834" stroke-width="1.8"')
            else:
                b += C(x, y, 2.4, "bf")
    b += f'<rect x="{X0 + (xf0 + 2) * g - 14}" y="{Y0 + g - 14}" width="{(nx - 1 - xf0 - 2) * g + 28}" height="{(ny - 3) * g + 28}" rx="6" class="a" stroke-width="1.5" stroke-dasharray="5 3"/>'
    ly = Y0 + ny * g + 4
    b += C(X0, ly, 4, "lnf") + T(X0 + 12, ly + 4, "nœud de la grille parente", "t2")
    b += C(X0, ly + 22, 2.6, "bf") + T(X0 + 12, ly + 26, "maille de la grille fine (Δx/2, Δt/2)", "t2")
    b += C(X0, ly + 44, 3.6, "sf", 'style="stroke:#eb6834" stroke-width="1.8"') + T(X0 + 12, ly + 48, "fantôme : interpolation cubique en espace, linéaire en temps", "t2")
    b += C(X0, ly + 66, 5.5, "af") + T(X0 + 12, ly + 70, "recopie fine → parente (à MARGIN = 2 du bord)", "t2")
    b += T(X0, ly + 96, "Π⁽¹⁾_fine = τ_fine / (2 τ_parente) · Π⁽¹⁾_parente", "tm")
    svg("raffinement.svg", 900, 450, b, "Grilles imbriquées et raccord entre une grille et sa grille fine")


# ---------------------------------------------------------------- blocs creux
def blocs():
    b = ""
    X0, Y0, bs = 40, 60, 34
    nbx, nby = 16, 8
    pts = [(X0 + 40 + 470 * x, Y0 + 210 - 520 * y) for x, y in naca(80) if y >= 0]
    D = 30

    def dist(px_, py_):
        return min(math.hypot(px_ - x, py_ - y) for x, y in pts)

    kept = 0
    for i in range(nbx):
        for j in range(nby):
            x, y = X0 + i * bs, Y0 + j * bs
            dc = dist(x + bs / 2, y + bs / 2)
            keep = dc < D + bs * 0.71 and y + bs / 2 < Y0 + 214
            cls = "ol" if keep else "bg"
            kept += keep
            b += f'<rect x="{x}" y="{y}" width="{bs}" height="{bs}" class="{cls}"/>'
            b += f'<rect x="{x}" y="{y}" width="{bs}" height="{bs}" class="grid" stroke-width="1"/>'
    # extrados (surface) et couche D
    d = "M" + " L".join(f"{x:.1f},{y:.1f}" for x, y in pts)
    b += f'<path d="{d}" class="b" stroke-width="{2*D}" stroke-opacity=".18" stroke-linecap="round"/>'
    b += f'<path d="{d}" class="ln" stroke-width="2.5"/>'
    body = d + f" L{pts[-1][0]:.1f},{Y0 + 214} L{pts[0][0]:.1f},{Y0 + 214} Z"
    b += f'<path d="{body}" class="body" stroke-width="0" fill-opacity=".8"/>'
    b += f'<rect x="{X0}" y="{Y0 + 214}" width="{nbx*bs}" height="{(nby)*bs - 214}" fill="url(#hach)" opacity=".4"/>'
    b += T(X0 + nbx * bs / 2, 36, "Couche creuse : blocs de 8³ mailles", "tb", "middle")
    lx = X0 + nbx * bs + 30
    b += f'<rect x="{lx}" y="80" width="16" height="16" class="ol"/><rect x="{lx}" y="80" width="16" height="16" class="grid"/>'
    b += T(lx + 24, 93, "bloc gardé (maille active", "t2") + T(lx + 24, 109, "ou fantôme à remplir)", "t2")
    b += f'<rect x="{lx}" y="130" width="16" height="16" class="bg"/><rect x="{lx}" y="130" width="16" height="16" class="grid"/>'
    b += T(lx + 24, 143, "bloc absent : toute lecture", "t2") + T(lx + 24, 159, "va au bloc vide n° 0", "t2")
    b += f'<rect x="{lx}" y="180" width="16" height="16" class="bl" fill-opacity=".18"/>'
    b += T(lx + 24, 193, "couche d'épaisseur D", "t2")
    b += T(lx, 232, "rang = slot · 512 + local", "tm")
    b += T(lx, 256, "id_of(x, y, z) : table des blocs", "t2")
    b += T(lx, 276, "pos_of(rang) : coordonnées du bloc", "t2")
    b += T(lx, 306, "TR-1208, couche ×8 : 2,6 M mailles", "t2")
    b += T(lx, 322, "gardées au lieu de 16 M (boîte)", "t2")
    svg("blocs_creux.svg", 860, 360, b, "Stockage creux par blocs de 8 au cube d'une couche pariétale")


# ---------------------------------------------------------------- efforts
def efforts():
    b = ""
    titles = ["1. Échange de quantité de mouvement", "2. Intégration sur la surface", "3. Bilan sur boîtes"]
    for k, t in enumerate(titles):
        b += T(20 + 300 * k + 135, 34, t, "tb", "middle")
    # 1. escalier et liens
    x0, y0, g = 30, 90, 30
    surf = lambda x: y0 + 150 - 0.35 * (x - x0)
    b += f'<path d="M{x0},{surf(x0):.1f} L{x0+250},{surf(x0+250):.1f} L{x0+250},{y0+200} L{x0},{y0+200} Z" fill="url(#hach)" opacity=".5"/>'
    b += L(x0, surf(x0), x0 + 250, surf(x0 + 250), w=2.5)
    for i in range(9):
        for j in range(6):
            x, y = x0 + 5 + i * g, y0 + j * g
            if y < surf(x) - 2:
                b += C(x, y, 3.2, "bf")
                if y + g > surf(x):
                    for dxl in (-1, 0, 1):
                        xe, ye = x + dxl * g, y + g
                        if ye > surf(xe):
                            t = (surf(x) - y) / (g * (1 + 0.35 * dxl))
                            xi, yi = x + t * dxl * g, y + t * g
                            b += L(x, y, xi, yi, cls="o", w=1.6)
                            b += C(xi, yi, 2.6, "of")
            else:
                b += C(x, y, 3, "sf", 'style="stroke:#898781" stroke-width="1.2"')
    b += T(x0, 330, "F = Σ c_j (f_j* + f_ī) − u_w (f_j* − f_ī)", "tm")
    b += T(x0, 350, "liens coupés (orange) ; surestime la", "t2")
    b += T(x0, 366, "traînée aux marches de l'escalier", "t2")
    # 2. surface
    x0 = 330
    cxs, cys, R = x0 + 135, 255, 140
    arc = [(cxs + R * math.cos(t), cys - R * math.sin(t)) for t in [math.pi * (0.2 + 0.6 * i / 40) for i in range(41)]]
    b += P(arc + [(arc[-1][0], 300), (arc[0][0], 300)], "body", 0, close=True, extra='fill-opacity=".7"')
    b += P(arc, "ln", 2.5)
    for k in range(5):
        tt = math.pi * (0.28 + 0.11 * k)
        x, y = cxs + R * math.cos(tt), cys - R * math.sin(tt)
        nx_, ny_ = math.cos(tt), -math.sin(tt)
        b += L(x, y, x + 34 * nx_, y + 34 * ny_, cls="mu", w=1, dash="2 2")
        b += C(x + 13 * nx_, y + 13 * ny_, 3, "af") + C(x + 26 * nx_, y + 26 * ny_, 3, "af")
        b += L(x + 26 * nx_, y + 26 * ny_, x + 3 * nx_, y + 3 * ny_, cls="b", w=1.8, mk="b")
        b += L(x, y, x + 24 * ny_, y - 24 * nx_, cls="o", w=1.8, mk="o")
    b += T(x0, y0 + 230 + 90 - 90, "", "t2")
    b += T(x0 + 5, 330, "−p_p n dA + τ_p t dA par triangle", "tm")
    b += T(x0 + 5, 350, "p_p extrapolée depuis ½ et 1 maille (aqua),", "t2")
    b += T(x0 + 5, 366, "pression en bleu, frottement en orange", "t2")
    # 3. boîte
    x0 = 630
    b += P([(x0 + 70 + 140 * x, 190 - 90 * y) for x, y in naca()], "body", 1.5, close=True)
    b += f'<rect x="{x0+30}" y="105" width="225" height="160" class="b" stroke-width="2" stroke-dasharray="7 4"/>'
    for yy in (135, 185, 235):
        b += L(x0 - 5, yy, x0 + 28, yy, cls="o", w=1.6, mk="o")
        b += L(x0 + 257, yy, x0 + 290, yy, cls="o", w=1.6, mk="o")
    b += T(x0 + 5, 330, "F = −∮ (ρ u (u·n) + p n) dA", "tm")
    b += T(x0 + 5, 350, "indépendant du traitement de la paroi ;", "t2")
    b += T(x0 + 5, 366, "Cx extrapolé à volume nul", "t2")
    svg("efforts.svg", 930, 380, b, "Les trois mesures des efforts")


# ---------------------------------------------------------------- graphiques
def axes(X0, Y0, W, H, xt, yt, fx, fy, xl, yl):
    b = ""
    for v, lab in yt:
        y = fy(v)
        b += L(X0, y, X0 + W, y, cls="grid", w=1)
        b += T(X0 - 8, y + 4, lab, "t2", "end")
    for v, lab in xt:
        x = fx(v)
        b += L(x, Y0, x, Y0 + 5, cls="mu", w=1)
        b += T(x, Y0 + 20, lab, "t2", "middle")
    b += L(X0, Y0, X0 + W, Y0, cls="mu", w=1)
    b += T(X0 + W, Y0 + 42, xl, "t2", "end")
    b += T(X0, Y0 - H - 14, yl, "t2", "start")
    return b


def series(pts, c, fx, fy, dash=None):
    pp = [(fx(x), fy(y)) for x, y in pts]
    b = P(pp, c, 2, extra=f'stroke-dasharray="{dash}"' if dash else "")
    for x, y in pp:
        b += C(x, y, 4.5, c + "f", 'style="stroke:#fcfcfb" stroke-width="1.5"')
    return b


def sphere():
    X0, Y0, W, H = 80, 320, 520, 250
    lx = lambda r: math.log10(r)
    fx = lambda r: X0 + (lx(r) - 1.2) / (3.1 - 1.2) * W
    fy = lambda c: Y0 - (math.log10(c) - math.log10(0.4)) / (math.log10(3.2) - math.log10(0.4)) * H
    re = [20, 50, 100, 200, 300, 500, 1000]
    d12 = [2.848, 1.665, 1.168, 0.856, 0.739, 0.671, 0.580]
    d16 = [2.838, 1.648, 1.145, 0.823, 0.695, 0.617, 0.533]
    ref = [2.610, 1.539, 1.094, 0.810, 0.691, 0.576, 0.466]
    b = T(20, 30, "Sphère : coefficient de traînée selon Re", "tb")
    b += axes(X0, Y0, W, H, [(r, str(r)) for r in [20, 50, 100, 200, 500, 1000]],
              [(v, f"{v:g}".replace(".", ",")) for v in [0.5, 1, 2, 3]], fx, fy, "Re", "Cx (échelles log)")
    b += P([(fx(r), fy(c)) for r, c in zip(re, ref)], "ln", 2, extra='stroke-dasharray="6 4"')
    b += series(list(zip(re, d12)), "o", fx, fy)
    b += series(list(zip(re, d16)), "b", fx, fy)
    b += T(fx(1000) + 10, fy(0.580) - 4, "D = 12", "t2")
    b += T(fx(1000) + 10, fy(0.533) + 9, "D = 16", "t2")
    b += T(fx(1000) + 10, fy(0.466) + 14, "référence", "t2")
    lgx = 330
    b += L(lgx, 60, lgx + 24, 60, cls="b", w=2) + T(lgx + 30, 64, "D = 16 mailles", "t2")
    b += L(lgx, 80, lgx + 24, 80, cls="o", w=2) + T(lgx + 30, 84, "D = 12 mailles", "t2")
    b += L(lgx, 100, lgx + 24, 100, w=2, dash="6 4") + T(lgx + 30, 104, "Schiller–Naumann, Clift", "t2")
    svg("sphere_cx.svg", 680, 370, b, "Traînée de la sphère selon le nombre de Reynolds, soufflerie et référence")


def decrochage():
    X0, Y0, W, H = 70, 360, 560, 290
    fx = lambda a: X0 + a / 20 * W
    fy = lambda c: Y0 - (c + 0.1) / 2.0 * H
    al = [0, 4, 8, 10, 12, 14, 15, 16, 17, 18, 19]
    s64 = [-0.013, 0.458, 0.971, 1.115, 1.311, 1.534, 1.491, 1.621, 1.641, 1.706, 1.693]
    s128 = [-0.013, 0.458, 0.899, 1.112, 1.317, 1.517, 1.591, 1.710, 1.790, 1.714, 1.765]
    tmr = [(10, 1.091), (15, 1.546), (16, 1.62), (17, 1.67), (18, 1.69), (19, 1.60)]
    lad = [(0, 0, 0), (4, .43, .44), (8, .87, .89), (10, 1.07, 1.08), (12, 1.26, 1.27), (14, 1.44, 1.45),
           (15, 1.51, 1.53), (16, 1.57, 1.58), (17, 1.61, 1.63), (18, 1.0, 1.2)]
    b = T(20, 30, "NACA 0012, Re = 6·10⁶ : portance jusqu'au décrochage", "tb")
    b += axes(X0, Y0, W, H, [(a, f"{a}°") for a in range(0, 21, 2)],
              [(v, f"{v:g}".replace(".", ",")) for v in [0, 0.5, 1, 1.5]], fx, fy, "incidence α", "Cz")
    # Ladson : barres d'étendue
    for a, lo, hi in lad:
        x = fx(a)
        b += f'<rect x="{x-5:.1f}" y="{fy(hi)-2:.1f}" width="10" height="{fy(lo)-fy(hi)+4:.1f}" rx="2" class="lnf" fill-opacity=".55"/>'
    b += series(tmr, "a", fx, fy)
    b += series(list(zip(al, s64)), "o", fx, fy)
    b += series(list(zip(al, s128)), "b", fx, fy)
    b += T(fx(18) - 8, fy(1.1) + 4, "Ladson : chute à 18°", "t2", "end")
    lgx, lgy = 90, 70
    b += L(lgx, lgy, lgx + 24, lgy, cls="b", w=2) + T(lgx + 30, lgy + 4, "soufflerie 3D, 128 mailles par corde", "t2")
    b += L(lgx, lgy + 20, lgx + 24, lgy + 20, cls="o", w=2) + T(lgx + 30, lgy + 24, "soufflerie 3D, 64 mailles par corde", "t2")
    b += L(lgx, lgy + 40, lgx + 24, lgy + 40, cls="a", w=2) + T(lgx + 30, lgy + 44, "NASA TMR, Spalart–Allmaras", "t2")
    b += f'<rect x="{lgx+7}" y="{lgy+52}" width="10" height="16" rx="2" class="lnf" fill-opacity=".55"/>' + T(lgx + 30, lgy + 64, "mesure de Ladson (étendue des grains)", "t2")
    svg("decrochage.svg", 680, 410, b, "Portance du NACA 0012 jusqu'au décrochage : soufflerie, TMR et Ladson")


if __name__ == "__main__":
    d3q19()
    bouzidi()
    couche()
    raffinement()
    blocs()
    efforts()
    sphere()
    decrochage()
    print("ok")
