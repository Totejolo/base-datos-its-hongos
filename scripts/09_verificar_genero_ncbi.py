# detectar_generos_faltantes_fasta.py
from Bio import Entrez, SeqIO
import pandas as pd
import time

# CAMBIA ESTO
Entrez.email = "tu_correo@ejemplo.com"  # correo con el que NCBI identifica las consultas
fasta_file = "micobiotapostmega_ORDENADO_POR_GENERO.fas"  # ← tu archivo

# Géneros clínicos esperados en micobiota intestinal / clínica
GENEROS_ESPERADOS = {
    'Candida', 'Saccharomyces', 'Malassezia', 'Trichosporon', 'Aspergillus',
    'Penicillium', 'Pichia', 'Debaryomyces', 'Fusarium', 'Cladosporium',
    'Rhodotorula', 'Cryptococcus', 'Meyerozyma', 'Wickerhamomyces',
    'Cyberlindnera', 'Nakaseomyces', 'Yarrowia', 'Kazachstania'
}

print("Leyendo FASTA y extrayendo accessions...")
records = list(SeqIO.parse(fasta_file, "fasta"))

data = []
for rec in records:
    header = rec.id
    partes = header.split("|")
    accession = partes[0] if len(partes) > 0 else "Unknown"
    genero_header = partes[2].strip() if len(partes) > 2 else "Unknown"
    data.append({
        'Accession': accession,
        'Header_completo': header,
        'Genero_header': genero_header,
        'Secuencia': str(rec.seq)
    })

df = pd.DataFrame(data)
print(f"Total secuencias: {len(df)}")

# BLAST automático a NCBI con Max target sequences = 5000
print("Haciendo consulta a NCBI (puede tardar 3-8 minutos)...")
generos_ncbi = []

for i, acc in enumerate(df['Accession']):
    if acc == "Unknown":
        generos_ncbi.append("Unknown")
        continue
    try:
        print(f"   {i+1}/{len(df)} → {acc}", end="\r")
        handle = Entrez.efetch(db="nucleotide", id=acc, rettype="gb", retmode="text")
        record = SeqIO.read(handle, "genbank")
        handle.close()
        organismo = record.annotations.get("organism", "Unknown")
        genero_ncbi = organismo.split()[0]  # primera palabra = género
        generos_ncbi.append(genero_ncbi)
    except Exception as e:
        print(f"\nError con {acc}: {e}")
        generos_ncbi.append("Error_BLAST")
    time.sleep(0.3)  # ser amable con NCBI

df['Genero_NCBI'] = generos_ncbi

# Detectar los que FALTAN (no están en lista clínica)
df['Es_faltante'] = ~df['Genero_NCBI'].isin(GENEROS_ESPERADOS)

# Guardar resultados
faltantes = df[df['Es_faltante'] == True].copy()
faltantes = faltantes[['Accession', 'Genero_header', 'Genero_NCBI', 'Header_completo']]
faltantes.to_csv("GENEROS_FALTANTES_BLAST_NCBI.csv", index=False)

todos = df[['Accession', 'Genero_header', 'Genero_NCBI']].copy()
todos.to_csv("TODOS_LOS_GENEROS_BLAST_NCBI.csv", index=False)

# RESUMEN BONITO
print("\n" + "="*60)
print("RESUMEN FINAL")
print("="*60)
print(f"Total secuencias analizadas      : {len(df)}")
print(f"Géneros clínicos esperados       : {len(GENEROS_ESPERADOS)}")
print(f"Secuencias con género FALTANTE   : {len(faltantes)}")
print("\nGéneros raros / inesperados más comunes:")
print(faltantes['Genero_NCBI'].value_counts().head(15))
print("\nArchivos generados:")
print("   → GENEROS_FALTANTES_BLAST_NCBI.csv  (solo los raros)")
print("   → TODOS_LOS_GENEROS_BLAST_NCBI.csv  (todos con taxonomía real)")
print("="*60)