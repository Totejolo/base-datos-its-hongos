import pandas as pd
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

# RUTAS
input_fas = 'micobiotapostmega.fas'
output_fas = 'micobiotapostmega_CORREGIDO.fas'

# Motivo conservado (forward)
FORWARD_MOTIF = "TCCGTAGGTGAACCTGCGGAAGGATCA"

# Leer FASTA
records = list(SeqIO.parse(input_fas, "fasta"))

corrected_records = []
for rec in records:
    seq = str(rec.seq).upper()
    rev_seq = str(Seq(seq).reverse_complement()).upper()

    # Comparar con motivo
    if seq.startswith(FORWARD_MOTIF):
        corrected_seq = seq
        orientation = "forward"
    elif rev_seq.startswith(FORWARD_MOTIF):
        corrected_seq = rev_seq
        orientation = "REVERSE_COMPLEMENT_CORREGIDO"
    else:
        # Si no es claro, usar score
        fwd_score = sum(1 for i in range(len(FORWARD_MOTIF)) if i < len(seq) and seq[i] == FORWARD_MOTIF[i])
        rev_score = sum(1 for i in range(len(FORWARD_MOTIF)) if i < len(rev_seq) and rev_seq[i] == FORWARD_MOTIF[i])
        if rev_score > fwd_score:
            corrected_seq = rev_seq
            orientation = "REVERSE_COMPLEMENT_CORREGIDO_SCORE"
        else:
            corrected_seq = seq
            orientation = "forward"

    # NUEVO: Evitar caracteres especiales
    new_rec = SeqRecord(
        Seq(corrected_seq),
        id=rec.id,
        description=f"{rec.description.split('|', 1)[1] if '|' in rec.description else ''} | {orientation}"
    )
    corrected_records.append(new_rec)

# GUARDAR CON UTF-8 (compatible con Windows)
with open(output_fas, 'w', encoding='utf-8') as f:
    SeqIO.write(corrected_records, f, "fasta")

print(f"Archivo corregido: {output_fas}")
print(f"Secuencias procesadas: {len(corrected_records)}")
print("Todas en orientación FORWARD - Listo para MEGA!")