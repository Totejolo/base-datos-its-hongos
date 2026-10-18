# ordenar_fasta_por_genero.py
from Bio import SeqIO
from collections import defaultdict
import re

# 1. Leer el .fas
fasta_file = "ITS_13009_ORDENADO.fasta"  # ← TU ARCHIVO
records = list(SeqIO.parse(fasta_file, "fasta"))

# 2. Extraer género y especie del ID: >NR_..._Genero_especie
def extraer_genero_especie(header):
    match = re.search(r'_([A-Za-z]+)_([a-z]+)', header)
    if match:
        genero = match.group(1)
        especie = match.group(2)
        return genero, especie
    return "Unknown", "unknown"

# 3. Agrupar por género
grupos = defaultdict(list)
for rec in records:
    genero, especie = extraer_genero_especie(rec.id)
    grupos[genero].append((especie, rec))

# 4. Ordenar: géneros alfabéticamente → especies alfabéticamente
generos_ordenados = sorted(grupos.keys())

# 5. Escribir nuevo .fas ordenado
output_file = "ITS_13009_ORDENADO_POR_GENERO.fas"
with open(output_file, "w") as f:
    for genero in generos_ordenados:
        # Ordenar especies dentro del género
        secuencias = sorted(grupos[genero], key=lambda x: x[0])
        for especie, rec in secuencias:
            f.write(f">{rec.id} [{genero} {especie}]\n")
            f.write(str(rec.seq) + "\n")

print(f"¡Listo! Archivo guardado como: {output_file}")
print(f"Total secuencias: {len(records)} | Géneros: {len(generos_ordenados)}")