import pandas as pd

# Cargar CSV
df = pd.read_csv('todas_secuencias_13009.csv')

# Extraer género (primera palabra, limpiar corchetes)
df['genero'] = df['especie'].str.split().str[0]
df['genero'] = df['genero'].str.replace(r'[\[\]]', '', regex=True)

# Géneros únicos
generos = df['genero'].dropna().unique()
total_generos = len(generos)
print(f"Total géneros distintos: {total_generos}")

# Conteo por género (top 20)
conteo = df['genero'].value_counts().head(20)
print("\nTop 10000 géneros:")
print(conteo)

# Guardar lista completa
pd.DataFrame(generos, columns=['genero']).to_csv('generos_1203.csv', index=False)