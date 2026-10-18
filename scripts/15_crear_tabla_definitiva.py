"""
Script para crear la tabla definitiva de hongos gastrointestinales
e importar las 1136 especies desde el FASTA alineado de MEGA.

Formato del header: >NR_077146.1|Aspergillus_caesiellus|Aspergillus|598bp_...|forward
"""

import os
import psycopg
import re

# ─────────────────────────────────────────
#  CONFIGURACIÓN — ajusta solo esto
# ─────────────────────────────────────────
FASTA_FILE = r"C:\Users\Usuario\Desktop\TFG (2)\TFG\micobiotapostmega_ORDENADO_POR_GENERO - copia.fas"

DB_HOST = "localhost"
DB_PORT = "5432"
DB_NAME = "TFG"
DB_USER = "postgres"
DB_PASS = os.environ.get("PGPASSWORD", "")  # contraseña de PostgreSQL, en la variable de entorno PGPASSWORD
# ─────────────────────────────────────────


def limpiar_secuencia(seq: str) -> str:
    """Elimina guiones de alineamiento y espacios, deja solo ACGTN."""
    seq = seq.replace("-", "").replace(" ", "").replace("\n", "").replace("\r", "")
    seq = seq.upper()
    seq = re.sub(r"[^ACGTN]", "N", seq)
    return seq


def parsear_header(header: str):
    """
    Extrae accession, especie y genero del header.
    Formato: >NR_077146.1|Aspergillus_caesiellus|Aspergillus|598bp_...|forward
    """
    header = header.lstrip(">").strip()
    partes = header.split("|")

    accession = partes[0].strip() if len(partes) > 0 else "desconocido"
    especie   = partes[1].replace("_", " ").strip() if len(partes) > 1 else "desconocida"
    genero    = partes[2].strip() if len(partes) > 2 else especie.split()[0]

    return accession, especie, genero


def leer_fasta(filepath: str):
    """Lee el FASTA alineado y devuelve lista de (accession, especie, genero, secuencia)."""
    registros = []
    accession = especie = genero = None
    seq_buffer = []

    with open(filepath, "r", encoding="latin-1") as f:
        for linea in f:
            linea = linea.rstrip()
            if linea.startswith(">"):
                # Guardar el anterior
                if accession and seq_buffer:
                    seq_limpia = limpiar_secuencia("".join(seq_buffer))
                    if len(seq_limpia) > 0:
                        registros.append((accession, especie, genero, seq_limpia))
                # Parsear nuevo header
                accession, especie, genero = parsear_header(linea)
                seq_buffer = []
            else:
                seq_buffer.append(linea)

        # Último registro
        if accession and seq_buffer:
            seq_limpia = limpiar_secuencia("".join(seq_buffer))
            if len(seq_limpia) > 0:
                registros.append((accession, especie, genero, seq_limpia))

    return registros


def crear_tabla(cur):
    """Elimina tablas antiguas innecesarias y crea la tabla definitiva."""

    print("Eliminando tablas intermedias antiguas...")
    tablas_a_eliminar = [
        "secuencias_ordenadas",
        "hongos_gastrointestinales",
        "micobiota_intestinal_full",
        "micobiota_219_perfecta",
        "micobiota_219_familia",
        "micobiota_364",
        "micobiota_364_full",
        "micobiota_364_tfg",
        "micobiota_animales",
        "micobiota_animales_full",
        "secuencias_its_completas",
        "secuencias_its_con_posiciones",
        "secuencias_its_definitiva",      # por si existe de antes
    ]
    for tabla in tablas_a_eliminar:
        cur.execute(f'DROP TABLE IF EXISTS public."{tabla}" CASCADE;')
        print(f"  → {tabla} eliminada (si existía)")

    cur.execute("DROP VIEW IF EXISTS public.secuencias_completas_7000 CASCADE;")
    cur.execute("DROP VIEW IF EXISTS public.todas_secuencias_13009 CASCADE;")
    print("  → Vistas eliminadas")

    print("\nCreando tabla definitiva: hongos_its_definitiva ...")
    cur.execute("""
        CREATE TABLE IF NOT EXISTS public.hongos_its_definitiva (
            id          SERIAL PRIMARY KEY,
            accession   TEXT UNIQUE NOT NULL,
            especie     TEXT NOT NULL,
            genero      TEXT,
            secuencia   TEXT NOT NULL,
            longitud    INTEGER GENERATED ALWAYS AS (length(secuencia)) STORED
        );
    """)

    # Índices para búsqueda rápida (importante para el pipeline NGS)
    cur.execute("CREATE INDEX IF NOT EXISTS idx_especie ON public.hongos_its_definitiva (especie);")
    cur.execute("CREATE INDEX IF NOT EXISTS idx_genero  ON public.hongos_its_definitiva (genero);")
    print("  → Tabla e índices creados correctamente")


def insertar_registros(cur, registros):
    """Inserta todos los registros en la tabla definitiva."""
    insertados = 0
    duplicados = 0

    for accession, especie, genero, secuencia in registros:
        try:
            cur.execute("""
                INSERT INTO public.hongos_its_definitiva (accession, especie, genero, secuencia)
                VALUES (%s, %s, %s, %s)
                ON CONFLICT (accession) DO NOTHING;
            """, (accession, especie, genero, secuencia))
            if cur.rowcount > 0:
                insertados += 1
            else:
                duplicados += 1
        except Exception as e:
            print(f"  ⚠ Error insertando {accession}: {e}")

    return insertados, duplicados


def main():
    print("=" * 55)
    print("  IMPORTACIÓN TABLA DEFINITIVA — TFG HONGOS ITS")
    print("=" * 55)

    # 1. Leer FASTA
    print(f"\nLeyendo FASTA: {FASTA_FILE}")
    registros = leer_fasta(FASTA_FILE)
    print(f"  → {len(registros)} secuencias encontradas en el FASTA")

    if len(registros) == 0:
        print("ERROR: No se encontraron secuencias. Revisa la ruta del FASTA.")
        return

    # 2. Conectar a BD
    print(f"\nConectando a PostgreSQL ({DB_HOST}/{DB_NAME})...")
    try:
        conn = psycopg.connect(
        host=DB_HOST, port=DB_PORT,
        dbname=DB_NAME, user=DB_USER, password=DB_PASS
        )
        conn.autocommit = False
        cur = conn.cursor()
        print("  → Conexión exitosa")
    except Exception as e:
        print(f"ERROR de conexión: {e}")
        return

    # 3. Crear tabla definitiva
    try:
        crear_tabla(cur)
        conn.commit()
    except Exception as e:
        conn.rollback()
        print(f"ERROR al crear tabla: {e}")
        conn.close()
        return

    # 4. Insertar registros
    print(f"\nInsertando {len(registros)} registros...")
    try:
        insertados, duplicados = insertar_registros(cur, registros)
        conn.commit()
    except Exception as e:
        conn.rollback()
        print(f"ERROR al insertar: {e}")
        conn.close()
        return

    # 5. Verificación final
    cur.execute("SELECT COUNT(*) FROM public.hongos_its_definitiva;")
    total = cur.fetchone()[0]

    cur.execute("SELECT COUNT(DISTINCT genero) FROM public.hongos_its_definitiva;")
    n_generos = cur.fetchone()[0]

    cur.execute("SELECT genero, COUNT(*) FROM public.hongos_its_definitiva GROUP BY genero ORDER BY COUNT(*) DESC LIMIT 5;")
    top_generos = cur.fetchall()

    print("\n" + "=" * 55)
    print("  RESULTADO FINAL")
    print("=" * 55)
    print(f"  Registros insertados : {insertados}")
    print(f"  Duplicados omitidos  : {duplicados}")
    print(f"  Total en BD          : {total}")
    print(f"  Géneros distintos    : {n_generos}")
    print(f"\n  Top 5 géneros:")
    for genero, count in top_generos:
        print(f"    {genero:30s} → {count} especies")
    print("=" * 55)

    cur.close()
    conn.close()
    print("\n¡Importación completada! Revisa pgAdmin.")


if __name__ == "__main__":
    main()
