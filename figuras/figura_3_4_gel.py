"""Figura 3.4: foto del gel en escala de grises (canal verde) con el contraste ajustado por igual en toda la imagen,
y rotulado de carriles y bandas. Entrada: foto_gel.png. Salida: figura_3_4_gel.png
"""
from PIL import Image, ImageDraw, ImageFont
import numpy as np
a = np.asarray(Image.open('foto_gel.png').convert('RGB')).astype(float)[:, :, 1]
lo, hi = np.percentile(a, 40), np.percentile(a, 99.9)
g = Image.fromarray((np.clip((a - lo) / (hi - lo), 0, 1) * 255).astype('uint8'))
S = 2
g = g.resize((g.width * S, g.height * S), Image.LANCZOS)
W, H = (690 + 190) * S, (480 + 62 + 40) * S
c = Image.new('RGB', (W, H), 'white')
c.paste(g.convert('RGB'), (0, 62 * S))
d = ImageDraw.Draw(c)
f = ImageFont.truetype('arial.ttf', 20 * S); fb = ImageFont.truetype('arialbd.ttf', 20 * S)
for lab, x in zip(['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', 'M'],
                  [78, 122, 163, 207, 249, 290, 332, 376, 418, 462, 506, 547]):
    font = fb if lab == 'M' else f
    d.text((x * S - d.textlength(lab, font=font) / 2, 16 * S), lab, font=font, fill='black')


def flecha(y, textos):
    d.line((744 * S, y * S, 704 * S, y * S), fill='black', width=2 * S)
    d.polygon([(700 * S, y * S), (712 * S, (y - 6) * S), (712 * S, (y + 6) * S)], fill='black')
    for k, t in enumerate(textos):
        d.text((752 * S, (y - 11 + 24 * k - 12 * (len(textos) - 1)) * S), t, font=f, fill='black')


flecha(345, ['≈ 500 pb'])
flecha(441, ['Dímeros de', 'cebadores'])
ley = 'Carriles 1-11: muestras  ·  M: marcador de peso molecular  ·  hibridación a 55 °C'
d.text((8 * S, (480 + 62 + 8) * S), ley, font=ImageFont.truetype('arial.ttf', 17 * S), fill='black')
c.save('figura_3_4_gel.png'); print(c.size)
