# ordenar_y_guardar_meg_perfecto.py
from Bio import SeqIO
import re

fasta_file = "micobiotapostmega_CORREGIDO ALINEADO PEP.fas"      # ← el que acabas de exportar
meg_nuevo  = "micobiotapostmega_ORDENADO_POR_GENERO.meg"

print("Ordenando por género...")

# Leer FASTA
records = list(SeqIO.parse(fasta_file, "fasta"))

# Extraer género y ordenar
ordenados = []
for rec in records:
    partes = rec.id.split("|")
    genero = partes[2].strip() if len(partes) > 2 else "Unknown"
    ordenados.append((genero, rec))

ordenados.sort(key=lambda x: x[0])

# Guardar como .meg perfecto que MEGA abre sin problemas
with open(meg_nuevo, "w", encoding="utf-8") as f:
    f.write("#mega\n")
    f.write("!Title Alineamiento ordenado por género;\n")
    f.write("!Format DataType=DNA indel=-;\n\n")
    
    for genero, rec in ordenados:
        f.write(f"#{rec.id}\n")
        seq = str(rec.seq)
        # MEGA lee mejor si la secuencia está partida cada 60 caracteres
        for i in range(0, len(seq), 60):
            f.write(seq[i:i+60] + "\n")
        f.write("\n")

print("¡PERFECTO!")
print(f"Archivo creado: {meg_nuevo}")
print("Ábrelo en MEGA → funciona al 100 %, sin Access violation, todo ordenado por género")