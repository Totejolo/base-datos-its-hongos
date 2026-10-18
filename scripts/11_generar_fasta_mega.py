# generar_fasta_ordenado.py
import pandas as pd

# TU CSV
csv_path = 'secuencias_ordenadas_MEGA.csv'
fasta_path = 'ITS_13009_ORDENADO.fasta'

df = pd.read_csv(csv_path, sep=';', encoding='utf-8')
df = df.sort_values(by='especie')  # Orden alfabético por especie

with open(fasta_path, 'w', encoding='utf-8') as f:
    for _, row in df.iterrows():
        acc = row['accession']
        especie = row['especie'].replace(' ', '_').replace('[', '').replace(']', '')
        seq = row['secuencia'].upper()
        f.write(f">{acc}_{especie}\n")
        for i in range(0, len(seq), 80):
            f.write(seq[i:i+80] + '\n')

print("FASTA ORDENADO LISTO →", fasta_path)