# filtrar_gastrointestinal.py (CORREGIDO PARA ;)
import pandas as pd
from Bio import Entrez
import time

# 1. TU EMAIL
Entrez.email = "tu_correo@ejemplo.com"  # correo con el que NCBI identifica las consultas

# 2. CARGAR CSV CON SEPARADOR ;
df = pd.read_csv('hongos_gastrointestinales.csv', sep=';', encoding='utf-8')

# 3. LIMPIAR NOMBRES DE COLUMNAS
df.columns = df.columns.str.strip().str.lower()

# 4. VERIFICAR
print("Columnas:", df.columns.tolist())
if 'accession' not in df.columns:
    print("ERROR: No hay 'accession'")
    exit()

# 5. DESCARGAR METADATOS
def get_metadata(acc):
    try:
        handle = Entrez.efetch(db="nucleotide", id=acc, rettype="gb", retmode="text")
        record = handle.read()
        source = ""
        host = ""
        for line in record.split('\n'):
            if '/isolation_source="' in line:
                source = line.split('"')[1]
            if '/host="' in line:
                host = line.split('"')[1]
        return source, host
    except Exception as e:
        print(f"Error {acc}: {e}")
        return "", ""
    time.sleep(0.3)

# 6. AÑADIR
sources = []
hosts = []
for acc in df['accession']:
    s, h = get_metadata(acc)
    sources.append(s)
    hosts.append(h)
    print(f"{acc} → {s} | {h}")

df['isolation_source'] = sources
df['host'] = hosts

# 7. GUARDAR CON COMA
df.to_csv('hongos_gastrointestinales_FINAL.csv', index=False, sep=',', encoding='utf-8')
print("\n¡LISTO! Archivo: hongos_gastrointestinales_FINAL.csv")