"""Figura 3.3: árbol completo sin nombres (A) y ampliación de un clado (B).
Entrada: arbol_nj.nwk (árbol de MEGA-CC, ver mega/). Salida: figura_3_3_arbol.png
"""
import collections
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.patches import Rectangle, ConnectionPatch
from Bio import Phylo

plt.rcParams['font.family'] = 'Arial'
tree = Phylo.read('arbol_nj.nwk', 'newick')
tree.root_at_midpoint()
tree.ladderize()


def especie(name):
    p = name.split('|')
    return p[1].replace('_', ' ') if len(p) > 1 else name


def genero(name):
    p = name.split('|')
    return p[2].strip('[]') if len(p) > 2 else '?'


leaves = tree.get_terminals()
print('hojas', len(leaves))
depth = tree.depths()
ypos = {c: i for i, c in enumerate(leaves)}


def y_of(clade):
    if clade in ypos:
        return ypos[clade]
    ys = [y_of(c) for c in clade.clades]
    ypos[clade] = (min(ys) + max(ys)) / 2
    return ypos[clade]


y_of(tree.root)


def segmentos(root):
    segs = []
    for cl in root.find_clades(order='preorder'):
        if cl.clades:
            ys = [ypos[c] for c in cl.clades]
            segs.append([(depth[cl], min(ys)), (depth[cl], max(ys))])
            for c in cl.clades:
                segs.append([(depth[cl], ypos[c]), (depth[c], ypos[c])])
    return segs


GRUPOS = [('Penicillium', '#1f77b4'), ('Aspergillus', '#2ca02c'), ('Fusarium', '#d62728'),
          ('Cladosporium', '#9467bd')]
COLOR = dict(GRUPOS)

# ---- clado ampliado: el de Trichosporon y sus parientes ----
tri = [l for l in leaves if genero(l.name) in ('Trichosporon', 'Cutaneotrichosporon')]
clado = tree.common_ancestor(tri)
sub = clado.get_terminals()
print('clado ampliado:', len(sub), collections.Counter(genero(l.name) for l in sub))

fig = plt.figure(figsize=(7.2, 9.6), dpi=300)
ax = fig.add_axes([0.02, 0.06, 0.36, 0.90])
az = fig.add_axes([0.47, 0.17, 0.52, 0.66])

ax.add_collection(LineCollection(segmentos(tree.root), colors='#333333', linewidths=0.25))
xmax = max(depth[l] for l in leaves)
for l in leaves:
    g = genero(l.name)
    ax.plot([xmax * 1.03, xmax * 1.10], [ypos[l]] * 2, color=COLOR.get(g, '#bbbbbb'), lw=0.6, solid_capstyle='butt')
ax.set_xlim(0, xmax * 1.12); ax.set_ylim(len(leaves) + 2, -2); ax.axis('off')
# barra de escala
ax.plot([0.02 * xmax, 0.02 * xmax + 0.1], [len(leaves) + 1] * 2, color='black', lw=0.8)
ax.text(0.02 * xmax + 0.05, len(leaves) + 1 - 8, '0,1 diferencias por posición', ha='left', va='bottom', fontsize=6)

# recuadro sobre el clado
y0 = min(ypos[l] for l in sub) - 2; y1 = max(ypos[l] for l in sub) + 2
x0 = depth[clado] - 0.01
rect = Rectangle((x0, y0), xmax * 1.11 - x0, y1 - y0, fill=False, ec='#e07b00', lw=1.0)
ax.add_patch(rect)

# ampliación
ys = {l: i for i, l in enumerate(sub)}
ypos_sub = dict(ys)


def y_sub(c):
    if c in ypos_sub:
        return ypos_sub[c]
    v = [y_sub(k) for k in c.clades]
    ypos_sub[c] = (min(v) + max(v)) / 2
    return ypos_sub[c]


y_sub(clado)
segs = []
for cl in clado.find_clades(order='preorder'):
    if cl.clades:
        v = [ypos_sub[c] for c in cl.clades]
        segs.append([(depth[cl], min(v)), (depth[cl], max(v))])
        for c in cl.clades:
            segs.append([(depth[cl], ypos_sub[c]), (depth[c], ypos_sub[c])])
az.add_collection(LineCollection(segs, colors='#333333', linewidths=0.6))
xm = max(depth[l] for l in sub)
for l in sub:
    az.text(xm + 0.004, ypos_sub[l], especie(l.name), fontsize=7, va='center', style='italic')
az.set_xlim(depth[clado] - 0.005, xm + (xm - depth[clado]) * 1.25)
az.set_ylim(len(sub), -1); az.axis('off')
for s in ('top', 'bottom', 'left', 'right'):
    az.spines[s].set_visible(False)
frame = Rectangle((0, 0), 1, 1, transform=az.transAxes, fill=False, ec='#e07b00', lw=1.0, clip_on=False)
az.add_patch(frame)

# líneas de zoom
for (ya, yb) in ((y0, 1.0), (y1, 0.0)):
    con = ConnectionPatch(xyA=(xmax * 1.11, ya), coordsA=ax.transData, xyB=(0, yb), coordsB=az.transAxes,
                          color='#e07b00', lw=0.8)
    fig.add_artist(con)

# leyenda
hand = [plt.Line2D([], [], color=c, lw=4, label=g) for g, c in GRUPOS] + \
       [plt.Line2D([], [], color='#bbbbbb', lw=4, label='Otros géneros')]
leg = fig.legend(handles=hand, loc='lower center', ncol=5, fontsize=6.5, frameon=False,
           prop={'style': 'italic', 'size': 6.5}, bbox_to_anchor=(0.5, 0.005), title='Género (barra de color)', title_fontsize=7)
fig.text(0.02, 0.975, 'A. Árbol completo: 1.136 secuencias ITS', fontsize=8.5, weight='bold')
fig.text(0.47, 0.845, 'B. Ampliación del clado recuadrado (30 secuencias)', fontsize=8.5, weight='bold')
fig.savefig('figura_3_3_arbol.png', dpi=300, facecolor='white')
print('ok')
