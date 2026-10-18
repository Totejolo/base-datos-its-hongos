import os
import re
from Bio import SeqIO
import gffutils
import psycopg2
import sys
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML

# Configuración de rutas para Fusarium
fasta_file = r"C:\Users\hp\Desktop\TFG\ncbi_dataset fusarium\data\GCF_013085055.1\GCF_013085055.1_ASM1308505v1_genomic.fna"
gff_file = r"C:\Users\hp\Desktop\TFG\ncbi_dataset fusarium\data\GCF_013085055.1\genomic.gff"
output_fasta = r"C:\Users\hp\Desktop\TFG\fusarium_its.fasta"
especie = "Fusarium oxysporum"

# Configuración de la base de datos
DB_HOST = "localhost"
DB_PORT = "5432"
DB_NAME = "tfg_Its"
DB_USER = "postgres"
DB_PASS = os.environ.get("PGPASSWORD", "")  # contraseña de PostgreSQL, en la variable de entorno PGPASSWORD

# Paso 1: Cargar genoma FASTA
print("Paso 1: Cargando genoma FASTA...")
try:
    genome = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
    print(f"Genoma cargado: {len(genome)} contigs. Contigs: {list(genome.keys())}")
except Exception as e:
    print(f"Error al cargar FASTA: {e}")
    sys.exit(1)

# Paso 2: Cargar GFF y buscar rRNA
print("Paso 2: Cargando GFF...")
try:
    db = gffutils.create_db(str(gff_file), dbfn='temp_fusarium.db', force=True, keep_order=True, merge_strategy='merge')
except Exception as e:
    print(f"Error al cargar GFF: {e}")
    sys.exit(1)

# Paso 3: Buscar regiones rRNA
print("Paso 3: Buscando regiones rRNA...")
rrna_features = []
for feature in db.all_features():
    rRNA_name = feature.attributes.get('Name', [''])[0].lower()
    rRNA_product = feature.attributes.get('product', [''])[0].lower()
    if feature.featuretype in ['rRNA', 'region', 'gene'] and any(x in rRNA_name + rRNA_product for x in ['18s', '5.8s', '28s', 'ribosomal rna']):
        print(f"rRNA encontrada: {feature.seqid}:{feature.start}-{feature.end}, Type={feature.featuretype}, Product={rRNA_product}")
        rrna_features.append({
            'seqid': feature.seqid,
            'start': feature.start - 1,  # 0-based
            'end': feature.end,
            'name': rRNA_product or rRNA_name
        })

# Paso 4: Extraer ITS desde GFF
its_sequences = []
if rrna_features:
    print("Paso 4: Extrayendo ITS desde GFF...")
    rrna_features.sort(key=lambda x: x['start'])
    for i in range(len(rrna_features) - 1):
        curr = rrna_features[i]
        next_feat = rrna_features[i + 1]
        if curr['seqid'] != next_feat['seqid']:
            continue
        if '18s' in curr['name'] and '5.8s' in next_feat['name']:
            its_start = curr['end'] + 1
            its_end = next_feat['start']
            its_length = its_end - its_start + 1
            if 100 <= its_length <= 500:
                try:
                    its_seq = genome[curr['seqid']].seq[its_start-1:its_end]
                    its_seq_clean = re.sub(r'[^ACGT]', 'N', str(its_seq).upper())
                    if len(its_seq_clean) > 0:
                        its_sequences.append({
                            'region': 'ITS1',
                            'seq': its_seq_clean,
                            'length': len(its_seq_clean),
                            'chr_id': curr['seqid'],
                            'start': its_start,
                            'end': its_end
                        })
                        print(f"ITS1 extraída: {curr['seqid']}:{its_start}-{its_end} ({its_length} bp)")
                except KeyError:
                    print(f"Error: Contig {curr['seqid']} no en FASTA")
        elif '5.8s' in curr['name'] and '28s' in next_feat['name']:
            its_start = curr['end'] + 1
            its_end = next_feat['start']
            its_length = its_end - its_start + 1
            if 100 <= its_length <= 500:
                try:
                    its_seq = genome[curr['seqid']].seq[its_start-1:its_end]
                    its_seq_clean = re.sub(r'[^ACGT]', 'N', str(its_seq).upper())
                    if len(its_seq_clean) > 0:
                        its_sequences.append({
                            'region': 'ITS2',
                            'seq': its_seq_clean,
                            'length': len(its_seq_clean),
                            'chr_id': curr['seqid'],
                            'start': its_start,
                            'end': its_end
                        })
                        print(f"ITS2 extraída: {curr['seqid']}:{its_start}-{its_end} ({its_length} bp)")
                except KeyError:
                    print(f"Error: Contig {curr['seqid']} no en FASTA")
else:
    print("Paso 4: No se encontraron rRNA en GFF. Usando BLAST local...")

    # Paso 5: BLAST local - CORREGIDO PARA ITS1 E ITS2
    print("Paso 5: Ejecutando BLAST local...")
    try:
        blastn_cline = NcbiblastnCommandline(
            cmd="blastn",
            query=fasta_file,
            db=r"C:\Users\hp\Desktop\TFG\its_db",
            out=r"C:\Users\hp\Desktop\TFG\blast_result.xml",
            outfmt=5,  # Formato XML
            evalue=1e-5  # Relajado para capturar más hits
        )
        stdout, stderr = blastn_cline()
        if stderr:
            print(f"Error en BLAST local: {stderr}")
            sys.exit(1)

        # Parsear resultados y segmentar ITS1 e ITS2
        print("Parseando resultados BLAST...")
        blast_records = NCBIXML.parse(open(r"C:\Users\hp\Desktop\TFG\blast_result.xml"))
        for record in blast_records:
            # Extraer contig_id del query (genoma), no del subject
            query_def = record.query_def if hasattr(record, 'query_def') else record.query
            contig_id = query_def.split()[0]  # Primer palabra del query_def, ej. NC_072840.1
            
            for alignment in record.alignments:
                for hsp in alignment.hsps:
                    if hsp.expect < 1e-5:  # Usar el evalue relajado
                        # Coordenadas del query (genoma)
                        its_start = hsp.query_start
                        its_end = hsp.query_end
                        its_length = its_end - its_start + 1
                        
                        if 100 <= its_length <= 800:  # Rango más amplio para ITS+5.8S
                            try:
                                # Extraer secuencia del genoma (query)
                                its_seq = genome[contig_id].seq[its_start-1:its_end]
                                its_seq_clean = re.sub(r'[^ACGT]', 'N', str(its_seq).upper())
                                
                                if len(its_seq_clean) > 0:
                                    total_length = len(its_seq_clean)
                                    print(f"Hit encontrado: {contig_id}:{its_start}-{its_end} ({total_length} bp), e-value: {hsp.expect}")
                                    
                                    # Segmentar ITS1 e ITS2 si el hit es largo
                                    if total_length >= 400:  # ITS1 + 5.8S + ITS2
                                        # División aproximada: ITS1 (primer 40%), 5.8S (20%), ITS2 (40%)
                                        its1_end = int(total_length * 0.4)
                                        its2_start = int(total_length * 0.6)
                                        its1_length = its1_end
                                        its2_length = total_length - its2_start
                                        
                                        its1_seq = its_seq_clean[:its1_end]
                                        its2_seq = its_seq_clean[its2_start:]
                                        
                                        # ITS1
                                        its_sequences.append({
                                            'region': 'ITS1',
                                            'seq': its1_seq,
                                            'length': len(its1_seq),
                                            'chr_id': contig_id,
                                            'start': its_start,
                                            'end': its_start + len(its1_seq) - 1
                                        })
                                        print(f"ITS1 segmentada: {contig_id}:{its_start}-{its_start + len(its1_seq) - 1} ({len(its1_seq)} bp)")
                                        
                                        # ITS2
                                        its_sequences.append({
                                            'region': 'ITS2',
                                            'seq': its2_seq,
                                            'length': len(its2_seq),
                                            'chr_id': contig_id,
                                            'start': its_start + its2_start,
                                            'end': its_end
                                        })
                                        print(f"ITS2 segmentada: {contig_id}:{its_start + its2_start}-{its_end} ({len(its2_seq)} bp)")
                                        
                                    else:
                                        # Si es corto, guardar como ITS genérica
                                        its_sequences.append({
                                            'region': 'ITS',
                                            'seq': its_seq_clean,
                                            'length': len(its_seq_clean),
                                            'chr_id': contig_id,
                                            'start': its_start,
                                            'end': its_end
                                        })
                                        print(f"ITS encontrada (corta): {contig_id}:{its_start}-{its_end} ({len(its_seq_clean)} bp)")
                                        
                            except KeyError as e:
                                print(f"Error: Contig {contig_id} no en FASTA. Disponibles: {list(genome.keys())}")
                            except Exception as e:
                                print(f"Error al extraer secuencia: {e}")
                                
    except Exception as e:
        print(f"Error en BLAST local: {e}")
        sys.exit(1)

# Paso 6: Guardar en FASTA
if its_sequences:
    with open(output_fasta, "w", encoding='utf-8') as out:
        for i, data in enumerate(its_sequences, 1):
            header = f">{especie}_{data['region']}_{i} chr_id={data['chr_id']} start={data['start']} end={data['end']}"
            print(f"Guardando header: {header}")
            out.write(f"{header}\n{data['seq']}\n")
    print(f"Paso 6: Secuencias ITS guardadas en {output_fasta} ({len(its_sequences)} secuencias)")
else:
    print("Error: No se extrajeron ITS ni con GFF ni con BLAST. Revisa los archivos de entrada.")
    sys.exit(1)

# Paso 7: Conectar a PostgreSQL e insertar/actualizar
print("Paso 7: Conectando a PostgreSQL e insertando/actualizando...")
try:
    conn = psycopg2.connect(
        host=DB_HOST, port=DB_PORT, database=DB_NAME,
        user=DB_USER, password=DB_PASS
    )
    cur = conn.cursor()
    updated = 0
    inserted = 0
    for data in its_sequences:
        # Buscar duplicados por especie, región, chr_id, start y end
        cur.execute("""
            SELECT id FROM secuencias_its 
            WHERE especie = %s AND region = %s AND chr_id = %s AND start = %s
        """, (especie, data['region'], data['chr_id'], data['start']))
        existing_id = cur.fetchone()
        if existing_id:
            cur.execute("""
                UPDATE secuencias_its 
                SET end_position = %s, longitud = %s, secuencia = %s
                WHERE id = %s
            """, (data['end'], data['length'], data['seq'], existing_id[0]))
            updated += 1
        else:
            cur.execute("""
                INSERT INTO secuencias_its (especie, region, secuencia, longitud, chr_id, start, end_position)
                VALUES (%s, %s, %s, %s, %s, %s, %s)
            """, (especie, data['region'], data['seq'], data['length'], data['chr_id'], data['start'], data['end']))
            inserted += 1
    conn.commit()
    print(f"Actualizadas {updated} secuencias y insertadas {inserted} nuevas en BD.")

    # Verificar
    cur.execute("SELECT COUNT(*) FROM secuencias_its WHERE especie = %s;", (especie,))
    count = cur.fetchone()[0]
    print(f"Verificación: {count} entradas para {especie}.")

    # Mostrar las filas
    cur.execute("""
        SELECT id, especie, region, longitud, chr_id, start, end_position 
        FROM secuencias_its 
        WHERE especie = %s 
        ORDER BY region, start
    """, (especie,))
    rows = cur.fetchall()
    for row in rows:
        print(f"ID {row[0]}: {row[1]} {row[2]} ({row[3]} bp) en {row[4]}:{row[5]}-{row[6]}")

    cur.close()
    conn.close()
except Exception as e:
    print(f"Error en BD: {e}")
    if 'conn' in locals():
        conn.close()
    sys.exit(1)

print("¡Procesamiento completado! Revisa pgAdmin.")