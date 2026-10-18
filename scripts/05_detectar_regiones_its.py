"""
Detecta las regiones ITS1, 5.8S e ITS2 en cada secuencia de hongos_its_definitiva
y actualiza la tabla añadiendo columnas con las posiciones y subsecuencias de cada región.

Método: búsqueda de motivos conservados (primers universales ITS1/ITS4 y 5.8S)
"""

import os
import psycopg2
import re

# ─────────────────────────────────────────
#  CONFIGURACIÓN
# ─────────────────────────────────────────
DB_HOST = "localhost"
DB_PORT = "5432"
DB_NAME = "tfg_Its"
DB_USER = "postgres"
DB_PASS = os.environ.get("PGPASSWORD", "")  # contraseña de PostgreSQL, en la variable de entorno PGPASSWORD
# ─────────────────────────────────────────


# ─────────────────────────────────────────
#  MOTIVOS CONSERVADOS PARA DETECCIÓN
#  Basados en primers universales ITS
# ─────────────────────────────────────────

# Fin del 18S / inicio de ITS1
MOTIVO_18S_FIN = [
    "TCCGTAGGTGAACCTGCGG",
    "TCCGTAGGTGAACCTGCG",
    "GTAGGTGAACCTGCGG",
]

# Inicio del 5.8S (= fin de ITS1)
MOTIVO_58S_INI = [
    "TCGATGAAGAACGCAGCGAAAT",
    "TCGATGAAGAACGCAGCG",
    "ATGAAGAACGCAGCGAAAT",
    "GAAGAACGCAGCGAAAT",
]

# Fin del 5.8S (= inicio de ITS2)
MOTIVO_58S_FIN = [
    "TGCGATAACTAATGTGAATTGCAG",
    "TGCGATACGTAATGTGAATTGCAG",
    "CGATAACTAATGTGAATTGCAG",
    "TGTGAATTGCAGAATTCAGTGAAT",
    "TGAATTGCAGAATTCAGTGAAT",
]

# Inicio del 28S (= fin de ITS2)
MOTIVO_28S_INI = [
    "TCGGATCAGGTAGGGATACCCGCT",
    "GGATCAGGTAGGGATACCCGCT",
    "TCGGATCAGGT",
    "AGGTAGGGATACCCGCT",
]


def buscar_motivo(seq, motivos, desde=0):
    """Busca el primer motivo encontrado en la secuencia desde una posición dada."""
    for motivo in motivos:
        match = re.search(motivo, seq[desde:], re.IGNORECASE)
        if match:
            return desde + match.start(), desde + match.end()
    return None, None


def detectar_regiones(seq):
    """
    Detecta posiciones de ITS1, 5.8S e ITS2.
    Devuelve dict con posiciones (0-based) y subsecuencias.
    """
    seq = seq.upper()
    r = {}

    # 1. Fin del 18S → inicio de ITS1
    _, fin_18s = buscar_motivo(seq, MOTIVO_18S_FIN)
    its1_start = fin_18s if fin_18s else 0

    # 2. Inicio del 5.8S → fin de ITS1
    ini_58s, fin_58s_ini = buscar_motivo(seq, MOTIVO_58S_INI, its1_start)
    its1_end = ini_58s if ini_58s else int(len(seq) * 0.30)

    r['ITS1_inicio'] = its1_start
    r['ITS1_fin']    = its1_end
    r['ITS1_seq']    = seq[its1_start:its1_end]
    r['5.8S_inicio'] = its1_end

    # 3. Fin del 5.8S → inicio de ITS2
    buscar_desde = r['5.8S_inicio']
    ini_its2, _ = buscar_motivo(seq, MOTIVO_58S_FIN, buscar_desde)
    fin_58s = ini_its2 if ini_its2 else r['5.8S_inicio'] + 160

    r['5.8S_fin']    = fin_58s
    r['5.8S_seq']    = seq[r['5.8S_inicio']:fin_58s]
    r['ITS2_inicio'] = fin_58s

    # 4. Inicio del 28S → fin de ITS2
    ini_28s, _ = buscar_motivo(seq, MOTIVO_28S_INI, r['ITS2_inicio'])
    its2_end = ini_28s if ini_28s else len(seq)

    r['ITS2_fin'] = its2_end
    r['ITS2_seq'] = seq[r['ITS2_inicio']:its2_end]

    return r


def main():
    print("=" * 55)
    print("  DETECCIÓN DE REGIONES ITS — TFG HONGOS")
    print("=" * 55)

    print(f"\nConectando a PostgreSQL ({DB_HOST}/{DB_NAME})...")
    try:
        conn = psycopg2.connect(
            host=DB_HOST, port=DB_PORT,
            database=DB_NAME, user=DB_USER, password=DB_PASS
        )
        cur = conn.cursor()
        print("  → Conexión exitosa\n")
    except Exception as e:
        print(f"ERROR de conexión: {e}")
        return

    # Añadir columnas si no existen
    print("Añadiendo columnas de regiones...")
    columnas = [
        ("its1_inicio", "INTEGER"), ("its1_fin", "INTEGER"), ("its1_seq", "TEXT"),
        ("s58_inicio",  "INTEGER"), ("s58_fin",  "INTEGER"), ("s58_seq",  "TEXT"),
        ("its2_inicio", "INTEGER"), ("its2_fin", "INTEGER"), ("its2_seq", "TEXT"),
    ]
    for col, tipo in columnas:
        cur.execute(f"ALTER TABLE public.hongos_its_definitiva ADD COLUMN IF NOT EXISTS {col} {tipo};")
    conn.commit()
    print("  → Columnas añadidas\n")

    # Procesar secuencias
    cur.execute("SELECT id, accession, especie, secuencia FROM public.hongos_its_definitiva ORDER BY id;")
    filas = cur.fetchall()
    print(f"Procesando {len(filas)} secuencias...\n")

    ok = parcial = sin_motivos = 0

    for id_, accession, especie, secuencia in filas:
        r = detectar_regiones(secuencia)

        tiene_its1 = len(r.get('ITS1_seq', '')) > 50
        tiene_58s  = len(r.get('5.8S_seq', '')) > 50
        tiene_its2 = len(r.get('ITS2_seq', '')) > 50

        if tiene_its1 and tiene_58s and tiene_its2:
            ok += 1
        elif tiene_its1 or tiene_its2:
            parcial += 1
        else:
            sin_motivos += 1
            print(f"  ⚠ Sin motivos claros: {accession} ({especie})")

        cur.execute("""
            UPDATE public.hongos_its_definitiva SET
                its1_inicio = %s, its1_fin = %s, its1_seq = %s,
                s58_inicio  = %s, s58_fin  = %s, s58_seq  = %s,
                its2_inicio = %s, its2_fin = %s, its2_seq = %s
            WHERE id = %s;
        """, (
            r.get('ITS1_inicio'), r.get('ITS1_fin'), r.get('ITS1_seq'),
            r.get('5.8S_inicio'), r.get('5.8S_fin'), r.get('5.8S_seq'),
            r.get('ITS2_inicio'), r.get('ITS2_fin'), r.get('ITS2_seq'),
            id_
        ))

    conn.commit()
    cur.close()
    conn.close()

    print("\n" + "=" * 55)
    print("  RESULTADO FINAL")
    print("=" * 55)
    print(f"  Detección completa (ITS1+5.8S+ITS2) : {ok}")
    print(f"  Detección parcial                    : {parcial}")
    print(f"  Sin motivos detectados               : {sin_motivos}")
    print("=" * 55)
    print("\n¡Listo! Revisa en pgAdmin las columnas its1_seq, s58_seq, its2_seq")


if __name__ == "__main__":
    main()
