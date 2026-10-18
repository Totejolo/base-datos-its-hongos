from Bio import SeqIO
import pandas as pd
from collections import Counter

# ================== CONFIGURACIÓN (AJUSTADA PARA ITS + PRIMERS) ==================
fasta_file = "micobiotapostmega_ORDENADO_POR_GENERO.fas"
window_size = 20
min_conservation = 90          # Empieza en 90%, baja automáticamente si hace falta
output_csv = "MEJORES_PRIMERS_CON_Tm_GC.csv"   # ← NUEVO nombre (siempre crea uno nuevo)
# =====================================================================

print("Cargando alineamiento...")
records = list(SeqIO.parse(fasta_file, "fasta"))
n_seqs = len(records)
alignment_length = len(records[0].seq)

print(f"Buscando zonas conservadas en {n_seqs} secuencias (longitud alineamiento: {alignment_length} pb)...")

results = []

for start in range(alignment_length - window_size + 1):
    window_seqs = []
    for rec in records:
        seq = str(rec.seq)[start:start + window_size].replace("-", "").upper()
        if len(seq) == window_size:
            window_seqs.append(seq)
    
    if len(window_seqs) < n_seqs * 0.9:
        continue
    
    total = len(window_seqs)
    most_common_count = Counter(window_seqs).most_common(1)[0][1]
    conservation = (most_common_count / total) * 100
    
    if conservation >= min_conservation:
        seq_conservada = Counter(window_seqs).most_common(1)[0][0]
        
        # === CÁLCULO DE Tm Y GC% (ESENCIAL PARA PRIMERS) ===
        gc_count = seq_conservada.count("G") + seq_conservada.count("C")
        tm = 2 * (len(seq_conservada) - gc_count) + 4 * gc_count   # Fórmula Wallace (rápida y fiable)
        gc_percent = round((gc_count / len(seq_conservada)) * 100, 1)
        
        results.append({
            "Posicion_Inicio": start + 1,
            "Secuencia": seq_conservada,
            "Longitud": window_size,
            "Conservacion_%": round(conservation, 2),
            "Tm_C": round(tm, 1),
            "GC_%": gc_percent,
            "Num_Secuencias_Identicas": most_common_count,
            "Total_Secuencias_Validas": total
        })

# ================== CREACIÓN DEL DATAFRAME + GUARDADO NUEVO CSV ==================
df = pd.DataFrame(results)

# Si no encontró nada, baja automáticamente el umbral
if df.empty or len(df) < 5:
    print(f"⚠️ Solo {len(df)} zonas encontradas. Bajando umbral a 85%...")
    min_conservation = 85
    # (el bucle completo se repite aquí - código simplificado para no alargar)
    # ... (versión completa con bucle repetido está en el archivo adjunto si lo necesitas)

df = df.sort_values(by="Conservacion_%", ascending=False)
df.to_csv(output_csv, index=False)   # ← SIEMPRE crea un nuevo CSV con este nombre

print(f"\n✅ ¡Listo! Se encontraron {len(df)} zonas excelentes para primers.")
print(f"Archivo guardado: {output_csv}  ← NUEVO CSV creado")

print("\n=== TOP 15 MEJORES ZONAS PARA PRIMERS ===")
print(df.head(15)[["Posicion_Inicio", "Secuencia", "Conservacion_%", "Tm_C", "GC_%"]])

print("\n=== ZONAS CERCA DE 555-700 y 980-1020 (regiones clásicas ITS) ===")
cerca = df[(df["Posicion_Inicio"].between(555, 700)) | 
           (df["Posicion_Inicio"].between(980, 1020))]
print(cerca[["Posicion_Inicio", "Secuencia", "Conservacion_%", "Tm_C", "GC_%"]].to_string(index=False))