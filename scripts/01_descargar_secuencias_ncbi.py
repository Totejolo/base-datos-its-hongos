import os
import pandas as pd
import psycopg2
from Bio import Entrez, SeqIO
import re
import sys
import time

# Configuración
Entrez.email = "tu_correo@ejemplo.com"  # correo con el que NCBI identifica las consultas
DB_HOST = "localhost"
DB_PORT = "5432"
DB_NAME = "tfg_Its"
DB_USER = "postgres"
DB_PASS = os.environ.get("PGPASSWORD", "")  # contraseña de PostgreSQL, en la variable de entorno PGPASSWORD
CSV_INPUT = r"C:\Users\hp\Desktop\TFG\secuencias_completas_sin_igs.csv"
CSV_OUTPUT = r"C:\Users\hp\Desktop\TFG\secuencias_its_exportado.csv"

def extraer_regiones(record):
    regions = {}
    for feature in record.features:
        if "product" in feature.qualifiers:
            product = feature.qualifiers["product"][0].lower()
            location = str(feature.location)
            if any(x in product for x in ["its1", "internal transcribed spacer 1"]): regions["ITS1"] = location
            elif "5.8s" in product or "5.8s rrna" in product: regions["5.8S"] = location
            elif any(x in product for x in ["its2", "internal transcribed spacer 2"]): regions["ITS2"] = location
            elif "18s" in product or "18s rrna" in product: regions["18S"] = location
            elif any(x in product for x in ["26s", "28s", "large subunit rrna"]): regions["26S/28S"] = location
    return regions

def parse_coordenada(location):
    start_match = re.search(r'(?:<)?(\d+)', location)
    end_match = re.search(r'(\d+)(?:>)?', location)
    start = int(start_match.group(1)) - 1 if start_match else None  # Base 0
    end = int(end_match.group(1)) if end_match else None
    return start, end

# Leer accession del CSV
df = pd.read_csv(CSV_INPUT, delimiter=';')
accessions = df['accession'].tolist()
print(f"Total de accessions en el CSV: {len(accessions)}")

# Conexión a la base de datos
try:
    conn = psycopg2.connect(host=DB_HOST, port=DB_PORT, database=DB_NAME, user=DB_USER, password=DB_PASS)
    cur = conn.cursor()
    print("Conexión a BD OK")
except Exception as e:
    print(f"Error BD: {e}")
    sys.exit(1)

# Crear tabla
cur.execute("""
    CREATE TABLE IF NOT EXISTS secuencias_its (
        id SERIAL PRIMARY KEY,
        especie TEXT,
        region TEXT,
        secuencia TEXT,
        longitud INT,
        chr_id TEXT,
        start INT,
        end_position INT,
        accession TEXT,
        description TEXT,
        regions_detected JSONB,
        CONSTRAINT unique_accession_region UNIQUE (accession, region)
    );
""")
conn.commit()

# Procesar solo los accessions del CSV
inserted = 0
skipped = 0
batch_size = 10
for start in range(0, len(accessions), batch_size):
    batch = accessions[start:start + batch_size]
    print(f"Lote {start//batch_size + 1}: {len(batch)} accessions")
    for seq_id in batch:
        try:
            handle = Entrez.efetch(db="nucleotide", id=seq_id, rettype="gb", retmode="text")
            record = SeqIO.read(handle, "genbank")
            handle.close()
            
            regions = extraer_regiones(record)
            if not regions:
                skipped += 1
                print(f"Saltada (sin regiones): {seq_id}")
                continue
            
            if "basidiomycota" in record.description.lower() or "agaricomycetes" in record.description.lower():
                skipped += 1
                print(f"Saltada (seta probable): {seq_id}")
                continue
            
            especie = record.annotations.get("organism", "Unknown")
            full_sequence = re.sub(r'[^ACGT]', 'N', str(record.seq).upper())
            
            for region_name in regions:
                start, end = parse_coordenada(regions[region_name])
                if start is None or end is None or start >= end:
                    print(f"Coordenadas inválidas para {region_name} en {seq_id}")
                    continue
                start = max(0, start)
                end = min(len(full_sequence), end)
                subsequence = full_sequence[start:end]
                longitud = len(subsequence)
                
                cur.execute("""
                    INSERT INTO secuencias_its (especie, region, secuencia, longitud, start, end_position, accession, description, regions_detected)
                    VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s)
                    ON CONFLICT (accession, region) DO NOTHING
                """, (especie, region_name, subsequence, longitud, start + 1, end, seq_id, record.description, json.dumps(regions)))
                inserted += 1
            
            conn.commit()
            time.sleep(0.4)
        except Exception as e:
            skipped += 1
            print(f"Error {seq_id}: {e}")
            conn.rollback()
    
    print(f"Insertadas: {inserted}, Saltadas: {skipped}")

# Exportar a CSV
cur.execute("""
    COPY (
        SELECT accession, especie, region, secuencia, longitud, start, end_position, regions_detected
        FROM secuencias_its
        WHERE regions_detected ? 'ITS1'
        ORDER BY accession, region
    ) TO %s WITH (FORMAT CSV, DELIMITER ';', HEADER)
""", (CSV_OUTPUT,))
conn.commit()

cur.close()
conn.close()
print(f"Terminado! Insertadas: {inserted} | Saltadas: {skipped}")
print(f"Datos exportados a {CSV_OUTPUT}")