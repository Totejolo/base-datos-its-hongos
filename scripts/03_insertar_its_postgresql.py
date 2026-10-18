# coding:UTF8
import os
from Bio import SeqIO
import psycopg2
import sys
try:
    conn= psycopg2.connect(
        host="localhost",
        port="5432",
        database="tfg_Its",
        user="postgres",
        password=os.environ.get("PGPASSWORD", ""),
        client_encoding= "UTF8"
    )
    cur=conn.cursor()
except Exception as e:
    print(f"Error al conectar a Postfres : {e}")
    sys.exit(1)
fasta_file = "its_extraidas.fasta"
try:
    for record in SeqIO.parse(fasta_file, "fasta"):
        especie = "Saccharomyces cerevisiae"
        region = record.id.split('_')[2]  # ITS1 o ITS2
        secuencia = str(record.seq)
        longitud = len(secuencia)
        chr_id = "NC_001144.5"  # chrXII
        start = 451574 if region == "ITS1" else 452226
        end_position = 452074 if region == "ITS1" else 452456

        cur.execute("""
            INSERT INTO secuencias_its (especie, region, secuencia, longitud, chr_id, start, end_position)
            VALUES (%s, %s, %s, %s, %s, %s, %s)
        """, (especie, region, secuencia, longitud, chr_id, start, end_position))
except Exception as e:
    print(f"Error al procesar FASTA o insertar: {e}")
    sys.exit(1)

conn.commit()
cur.close()
conn.close()

print("ITS insertadas en la base de datos. Verifica con: SELECT * FROM secuencias_its;")