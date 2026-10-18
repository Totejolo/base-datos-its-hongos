"""Figura 1.1: esquema de la unidad de ADN ribosómico de hongos y posición de los cebadores.
Salida: figura_1_1_esquema_its.png
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

plt.rcParams['font.family'] = 'Arial'
fig, ax = plt.subplots(figsize=(10, 3.6), dpi=300)
ax.set_xlim(0, 100); ax.set_ylim(0, 40); ax.axis('off')

Y, H = 18, 5          # altura de la línea de la unidad
GENE = '#9e9e9e'      # genes ribosómicos (conservados)
SPACER = '#ffffff'    # espaciadores (variables)

# Espaciador intergénico (IGS) a ambos lados
ax.plot([1, 10], [Y + H / 2] * 2, color='black', lw=1.5)
ax.plot([90, 99], [Y + H / 2] * 2, color='black', lw=1.5)
ax.text(5.5, Y - 3.5, 'IGS', ha='center', fontsize=9)
ax.text(94.5, Y - 3.5, 'IGS', ha='center', fontsize=9)

bloques = [  # (inicio, fin, etiqueta, relleno)
    (10, 32, '18S', GENE),
    (32, 44, 'ITS1', SPACER),
    (44, 51, '5.8S', GENE),
    (51, 63, 'ITS2', SPACER),
    (63, 90, '26S/28S', GENE),
]
for x0, x1, lab, col in bloques:
    ax.add_patch(FancyBboxPatch((x0, Y), x1 - x0, H, boxstyle='square,pad=0', fc=col, ec='black', lw=1.2))
    ax.text((x0 + x1) / 2, Y + H / 2, lab, ha='center', va='center', fontsize=10,
            fontweight='bold' if col == GENE else 'normal')

# Región ITS completa
ax.annotate('', xy=(32, Y + H + 3), xytext=(63, Y + H + 3),
            arrowprops=dict(arrowstyle='-', lw=1))
for x in (32, 63):
    ax.plot([x, x], [Y + H + 2.2, Y + H + 3.8], color='black', lw=1)
ax.text(47.5, Y + H + 4.6, 'Región ITS (≈ 400-800 pb)', ha='center', fontsize=9.5)


def primer(x, sentido, nombre, nivel, color):
    """Flecha de primer: sentido +1 (forward, hacia la derecha) o -1 (reverse)."""
    y = Y - 4 - 4.2 * nivel
    largo = 4.5
    ax.annotate('', xy=(x + sentido * largo, y), xytext=(x, y),
                arrowprops=dict(arrowstyle='-|>', lw=1.6, color=color, mutation_scale=11))
    ax.text(x + sentido * largo / 2, y - 2.6, nombre, ha='center', fontsize=8.5, color=color)


NEGRO, AZUL = '#333333', '#1f5fa8'
primer(28.5, +1, 'ITS1', 0, NEGRO)          # extremo 3' de 18S
primer(44.5, +1, 'ITS3', 0, NEGRO)          # dentro de 5.8S
primer(44.5, +1, 'LEVNGS F', 1, AZUL)       # misma zona que ITS3
primer(70.0, -1, 'ITS4', 0, NEGRO)          # extremo 5' de 26S/28S
primer(67.0, -1, 'LEVNGS R', 1, AZUL)       # extremo 5' de 26S/28S, antes de ITS4

# Amplicón de los primers diseñados
ax.plot([44.5, 67.0], [3.0, 3.0], color=AZUL, lw=1.2)
for x in (44.5, 67.0):
    ax.plot([x, x], [2.2, 3.8], color=AZUL, lw=1.2)
ax.text(55.75, 0.4, 'Amplicón LEVNGS F/R (≈ 330 pb)', ha='center', fontsize=8.5, color=AZUL)

# Leyenda
ax.add_patch(FancyBboxPatch((1, 34.5), 2.5, 1.8, boxstyle='square,pad=0', fc=GENE, ec='black', lw=0.8))
ax.text(4.2, 35.4, 'Genes de ARNr (conservados)', va='center', fontsize=8)
ax.add_patch(FancyBboxPatch((30, 34.5), 2.5, 1.8, boxstyle='square,pad=0', fc=SPACER, ec='black', lw=0.8))
ax.text(33.2, 35.4, 'Espaciadores (variables)', va='center', fontsize=8)

fig.tight_layout(pad=0.2)
fig.savefig('figura_1_1_esquema_its.png', dpi=300, facecolor='white')
print('ok')
