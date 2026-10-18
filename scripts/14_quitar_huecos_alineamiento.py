from Bio import SeqIO
from Bio.Seq import Seq

# ================== CONFIGURACIÓN ==================
input_file  = "micobiotapostmega_ORDENADO_POR_GENERO.fas"
output_file = "micobiotapostmega_SIN_GAPS.fas"
# ==================================================

print("Cargando archivo alineado...")
records = list(SeqIO.parse(input_file, "fasta"))

print(f"Eliminando gaps de {len(records)} secuencias...")

for rec in records:
    # Forma compatible con TODAS las versiones de Biopython
    rec.seq = Seq(str(rec.seq).replace("-", ""))

# Guardar el nuevo FASTA sin gaps
SeqIO.write(records, output_file, "fasta")

print(f"\n✅ ¡Listo!")
print(f"Archivo creado: **{output_file}**")
print(f"Procesadas {len(records)} secuencias sin ningún gap.")
print("El archivo está en la misma carpeta que el script.")