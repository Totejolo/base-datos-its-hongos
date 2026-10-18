"""Figura 3.2: rotula la captura del alineamiento de MEGA y deja resaltado solo el tramo del cebador directo.
Entrada: captura_alineamiento_mega.png (captura de MEGA). Salida: figura_3_2_alineamiento.png
"""
from PIL import Image, ImageDraw, ImageFont
import numpy as np

im = Image.open('captura_alineamiento_mega.png').convert('RGB')
a = np.asarray(im).copy()
W, H = im.size
X0, X1 = 598, 998            # columnas 7 a 26 del recorte: GCATCGATGAAGAACGCAGC (posición 492)
amarillo = (a[:, :, 0] > 200) & (a[:, :, 1] > 200) & (a[:, :, 2] < 120)
fuera = np.zeros_like(amarillo)
fuera[:, 482:X0] = True
fuera[:, X1:] = True
a[amarillo & fuera] = (255, 255, 255)
base = Image.fromarray(a)

S = 2
TOP = 120
c = Image.new('RGB', (W * S, (H + TOP) * S), 'white')
c.paste(base.resize((W * S, H * S), Image.LANCZOS), (0, TOP * S))
d = ImageDraw.Draw(c)
f = ImageFont.truetype('arial.ttf', 17 * S)
fb = ImageFont.truetype('arialbd.ttf', 17 * S)
y = (TOP - 14) * S
d.line((X0 * S, y, X1 * S, y), fill='black', width=2 * S)
d.line((X0 * S, y, X0 * S, y + 10 * S), fill='black', width=2 * S)
d.line((X1 * S, y, X1 * S, y + 10 * S), fill='black', width=2 * S)
mid = (X0 + X1) / 2 * S
t1 = 'Cebador directo candidato (posición 492)'
t2 = 'Conservación: 96,12 % de las 1.136 secuencias'
d.text((mid - d.textlength(t1, font=fb) / 2, (TOP - 92) * S), t1, font=fb, fill='black')
d.text((mid - d.textlength(t2, font=f) / 2, (TOP - 62) * S), t2, font=f, fill='black')
t3 = 'Secuencias de Aspergillus (MEGA)'
d.text((8 * S, (TOP - 40) * S), t3, font=f, fill='#444444')
c.save('figura_3_2_alineamiento.png')
print(c.size)
