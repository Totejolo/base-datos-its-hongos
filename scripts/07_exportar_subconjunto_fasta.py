import pandas as pd

# Ruta al archivo CSV (corrigiendo el nombre con espacios)
csv_path = 'micobiota intestinal final.csv'
fas_path = 'micobiota_full.fas'

# Leer el CSV
df = pd.read_csv(csv_path)

# Sin filtro: usa todas las secuencias (todas son de micobiota intestinal)
df_filtered = df  # O quita esta línea y usa df directamente abajo

# Crear el archivo FASTA
with open(fas_path, 'w') as fas_file:
    for _, row in df_filtered.iterrows():
        header = f">{row['accession']}|{row['especie'].replace(' ', '_')}|{row['genero']}|{row['longitud_secuencia']}bp\n"
        sequence = row['secuencia'] + '\n'
        fas_file.write(header + sequence)

print(f"Archivo FASTA generado: {fas_path}")
print(f"Secuencias exportadas: {len(df_filtered)}")