"""Cebadores inversos en sentido 5'->3' y tamaño del fragmento amplificado, calculados sobre el alineamiento múltiple.

Para cada candidato a cebador inverso muestra el tramo de la hebra alineada y su complementaria inversa, que es la
secuencia real del cebador. Después cuenta, en cada secuencia que contiene los dos sitios de unión, las bases que hay
entre el inicio del cebador directo y el final del inverso.

Uso: python 18_calcular_tamano_amplicon.py alineamiento.fas
"""
import collections
import statistics
import sys

from Bio import SeqIO
from Bio.Seq import Seq

ALINEAMIENTO = sys.argv[1] if len(sys.argv) > 1 else 'micobiotapostmega_CORREGIDO_ALINEADO_PEP.fas'
DIRECTO = 494                       # primera columna del cebador directo (ATCGATGAAGAACGCAGC, 18 nt)
INVERSO_INICIO, INVERSO_LARGO = 994, 19   # tramo cuya complementaria inversa es el cebador inverso (19 nt)
COLAS = 33 + 34                     # nucleótidos de las colas de los adaptadores de Illumina (directo + inverso)
CANDIDATOS_INVERSOS = [(993, 20), (994, 20), (970, 20)]   # posición y longitud (Tabla 3.5)

secuencias = [str(r.seq).upper() for r in SeqIO.parse(ALINEAMIENTO, 'fasta')]
print(f'Secuencias en el alineamiento: {len(secuencias)}; columnas: {len(secuencias[0])}')


def consenso(inicio, largo):
    """Base más frecuente (sin contar huecos) en cada columna del tramo; las columnas empiezan en 1."""
    tramo = ''
    for c in range(inicio - 1, inicio - 1 + largo):
        cuenta = collections.Counter(s[c] for s in secuencias if s[c] != '-')
        tramo += cuenta.most_common(1)[0][0]
    return tramo


print('\nCandidatos a cebador inverso')
for pos, largo in CANDIDATOS_INVERSOS:
    tramo = consenso(pos, largo)
    print(f'  posición {pos}: hebra alineada {tramo} -> cebador 5\'-{Seq(tramo).reverse_complement()}-3\'')

print('\nPareja de cebadores')
print(f'  directo 5\'-{consenso(DIRECTO, 18)}-3\'')
print(f'  inverso 5\'-{Seq(consenso(INVERSO_INICIO, INVERSO_LARGO)).reverse_complement()}-3\'')

fin = INVERSO_INICIO + INVERSO_LARGO - 1
largos = []
for s in secuencias:
    if s[DIRECTO - 1] == '-' or s[fin - 1] == '-':
        continue
    largos.append(sum(1 for b in s[DIRECTO - 1:fin] if b != '-'))
largos.sort()
q = statistics.quantiles(largos, n=20)
print(f'\nSecuencias con los dos sitios de unión: {len(largos)}')
print(f'Fragmento amplificado: mediana {statistics.median(largos):.0f} pb; el 90 % entre {q[0]:.0f} y {q[-1]:.0f} pb '
      f'(mínimo {largos[0]}, máximo {largos[-1]})')
print(f'Con las colas de los adaptadores (+{COLAS} nt): mediana {statistics.median(largos) + COLAS:.0f} pb; '
      f'el 90 % entre {q[0] + COLAS:.0f} y {q[-1] + COLAS:.0f} pb')
