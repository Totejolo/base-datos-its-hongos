from Bio import Entrez
import pandas as pd
import time
import sys
import re

# ==================== CONFIGURACIÓN ====================
Entrez.email = "tu_correo@ejemplo.com"  # correo con el que NCBI identifica las consultas
fasta_actual = "micobiotapostmega_ORDENADO_POR_GENERO.fas"
csv_faltantes = "GENEROS_FALTANTES_BLAST_NCBI.csv"
fasta_final = "micobiotapostmega_BASE_COMPLETA_FINAL.fas"  # Nuevo: FASTA combinado
meg_final = "micobiotapostmega_BASE_COMPLETA_FINAL.meg"    # Opcional: MEG
# =======================================================

# 1. Cargar CSV
df = pd.read_csv(csv_faltantes)
if 'Accession' not in df.columns:
    print("Error: El CSV no tiene la columna 'Accession'")
    sys.exit(1)

df['Secuencia'] = ''
df['Nombre_Cientifico'] = ''

print("Iniciando descarga de secuencias y nombres científicos desde NCBI...\n")

# 2. Fetch secuencia + nombre científico (si ya lo tienes, puedes skippear manualmente)
for idx, row in df.iterrows():
    acc = str(row['Accession']).strip()
    try:
        # Descargar secuencia (FASTA)
        handle = Entrez.efetch(db="nucleotide", id=acc, rettype="fasta", retmode="text")
        fasta_text = handle.read()
        seq = ''.join(fasta_text.splitlines()[1:]).upper()
        
        # Limpiar secuencia (solo ACGTN-)
        seq = re.sub(r'[^ACGTN-]', '', seq)

        # Descargar nombre científico (GenBank)
        handle_gb = Entrez.efetch(db="nucleotide", id=acc, rettype="gb", retmode="text")
        gb_text = handle_gb.read()
        
        # Extraer ORGANISM (Genus species)
        organism = "Unknown"
        for line in gb_text.splitlines():
            if line.strip().startswith("ORGANISM"):
                organism = line.split("ORGANISM")[1].strip()
                break

        df.at[idx, 'Secuencia'] = seq
        df.at[idx, 'Nombre_Cientifico'] = organism

        print(f"✓ {acc} → {organism}  ({len(seq)} bp)")

    except Exception as e:
        df.at[idx, 'Secuencia'] = "ERROR"
        df.at[idx, 'Nombre_Cientifico'] = f"ERROR: {e}"
        print(f"✗ Error con {acc}: {e}")

    time.sleep(0.4)   # Delay seguro

# Guardar CSV con la nueva información
try:
    df.to_csv("GENEROS_FALTANTES_con_secuencias.csv", index=False)
except PermissionError:
    print("ERROR: Cierra 'GENEROS_FALTANTES_con_secuencias.csv' y vuelve a ejecutar.")
    sys.exit(1)

# 3. Leer FASTA actual
actuales = []
with open(fasta_actual, "r") as f:
    header = ""
    seq = ""
    for line in f:
        line = line.strip()
        if line.startswith(">"):
            if header and seq:
                actuales.append((header, seq))
            header = line[1:]
            seq = ""
        elif line:
            seq += line
    if header and seq:
        actuales.append((header, seq))

# 4. Generar FASTA combinado (PRINCIPAL)
try:
    with open(fasta_final, "w", encoding="utf-8") as f:
        # Secuencias originales
        for header, seq in actuales:
            f.write(f">{header}\n")
            for i in range(0, len(seq), 60):
                f.write(seq[i:i+60] + "\n")
        
        # Secuencias nuevas
        añadidas = 0
        for _, row in df.iterrows():
            if row['Secuencia'] == "ERROR" or len(row['Secuencia']) < 50:
                continue

            header_nuevo = f"{row['Accession']}|{row['Nombre_Cientifico']}"
            header_nuevo = re.sub(r'[^\w| .-]', '_', header_nuevo)  # Limpieza

            f.write(f">{header_nuevo}\n")
            seq = row['Secuencia']
            for i in range(0, len(seq), 60):
                f.write(seq[i:i+60] + "\n")
            añadidas += 1
except Exception as e:
    print(f"Error generando FASTA: {e}")
    sys.exit(1)

print(f"\n¡FASTA GENERADO! Archivo: {fasta_final}")
print(f"Secuencias añadidas: {añadidas}")

# 5. Opcional: Generar .meg también (por si quieres comparar)
try:
    with open(meg_final, "w", encoding="utf-8") as f:
        f.write("#mega\n")
        f.write("!Title= Base completa + secuencias añadidas;\n")
        f.write("!Format DataType=DNA indel=-;\n\n")

        # Secuencias originales
        for header, seq in actuales:
            f.write(f"#{header}\n")
            for i in range(0, len(seq), 60):
                f.write(seq[i:i+60] + "\n")
            f.write("\n")

        # Secuencias nuevas
        for _, row in df.iterrows():
            if row['Secuencia'] == "ERROR" or len(row['Secuencia']) < 50:
                continue

            header_nuevo = f"{row['Accession']}|{row['Nombre_Cientifico']}"
            header_nuevo = re.sub(r'[^\w| .-]', '_', header_nuevo)

            f.write(f"#{header_nuevo}\n")
            seq = row['Secuencia']
            for i in range(0, len(seq), 60):
                f.write(seq[i:i+60] + "\n")
            f.write("\n")
    print(f"¡MEG también generado para comparar: {meg_final}")
except Exception as e:
    print(f"Error generando MEG opcional: {e}")

print("Ahora, abre el FASTA en MEGA: File > Open > selecciona .fas > luego File > Save As > MEGA format (.meg).")
print("Si MEGA da error, prueba el FASTA en otro software como BioEdit o UGENE para convertir.")