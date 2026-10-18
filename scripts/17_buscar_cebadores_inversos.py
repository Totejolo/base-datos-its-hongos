from Bio import SeqIO
import pandas as pd
from collections import Counter

fasta_file = "micobiotapostmega_ORDENADO_POR_GENERO.fas"
window_size = 20

print("🔍 Buscando zonas conservadas en la PARTE FINAL (995+)...")

records = list(SeqIO.parse(fasta_file, "fasta"))
alignment_length = len(records[0].seq)

results = []
min_conservation = 80

for start in range(950, alignment_length - window_size + 1):
    window_seqs = []
    for rec in records:
        seq = str(rec.seq)[start:start + window_size].replace("-", "").upper()
        if len(seq) == window_size:
            window_seqs.append(seq)
    
    if len(window_seqs) < len(records) * 0.85:
        continue
    
    total = len(window_seqs)
    most_common_count = Counter(window_seqs).most_common(1)[0][1]
    conservation = (most_common_count / total) * 100
    
    if conservation >= min_conservation:
        seq_conservada = Counter(window_seqs).most_common(1)[0][0]
        gc_count = seq_conservada.count("G") + seq_conservada.count("C")
        tm = 2 * (len(seq_conservada) - gc_count) + 4 * gc_count
        
        results.append({
            "Posicion_Inicio": start + 1,
            "Secuencia": seq_conservada,
            "Conservacion_%": round(conservation, 2),
            "Tm_C": round(tm, 1),
            "GC_%": round((gc_count / len(seq_conservada)) * 100, 1)
        })

df = pd.DataFrame(results)

# === SEGURIDAD: si no encontró nada, baja el umbral automáticamente ===
if df.empty:
    print(f"⚠️ No se encontraron zonas con ≥{min_conservation}%. Bajando umbral a 75%...")
    min_conservation = 75
    # (el bucle se repite aquí - versión simplificada)
    results = []
    for start in range(950, alignment_length - window_size + 1):
        window_seqs = [str(rec.seq)[start:start+window_size].replace("-","").upper() 
                       for rec in records if len(str(rec.seq)[start:start+window_size].replace("-","")) == window_size]
        if len(window_seqs) < len(records)*0.85: continue
        total = len(window_seqs)
        most = Counter(window_seqs).most_common(1)[0]
        cons = (most[1]/total)*100
        if cons >= min_conservation:
            seq = most[0]
            gc = seq.count("G")+seq.count("C")
            tm = 2*(len(seq)-gc) + 4*gc
            results.append({"Posicion_Inicio": start+1, "Secuencia": seq,
                            "Conservacion_%": round(cons,2), "Tm_C": round(tm,1),
                            "GC_%": round((gc/len(seq))*100,1)})
    df = pd.DataFrame(results)

df = df.sort_values(by="Conservacion_%", ascending=False)
df.to_csv("ZONAS_FINALES_REVERSAS_PRIMERS.csv", index=False)

print(f"\n✅ Archivo creado: **ZONAS_FINALES_REVERSAS_PRIMERS.csv**")
print(f"Se encontraron {len(df)} zonas en la parte final.")

print("\n=== TOP 10 MEJORES ZONAS REVERSAS ===")
print(df.head(10)[["Posicion_Inicio", "Secuencia", "Conservacion_%", "Tm_C", "GC_%"]].to_string(index=False))

print("\n=== RECOMENDACIÓN MÁS FUERTE PARA TU TFG ===")
print("Usa el primer reverse CLÁSICO (funciona en >98% de hongos):")
print("ITS4 → 5'-TCCTCCGCTTATTGATATGC-3'   (Tm ≈ 58°C | GC 45%)")
print("Combínalo con tu mejor forward de posición 492:")
print("Forward → GCATCGATGAAGAACGCAGC (96.12% conservación)")