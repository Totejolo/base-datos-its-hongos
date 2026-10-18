# Base de datos de secuencias ITS de hongos para análisis metagenómico

Código del Trabajo Fin de Grado *Elaboración de una Base de Datos de secuencias de la región nuclear ribosómica para
análisis metagenómico en ADN de hongos* (Grado en Bioquímica, Universidad de Extremadura).

El flujo de trabajo descarga secuencias ITS de hongos de NCBI, las organiza en una base de datos PostgreSQL, localiza las
regiones ITS1, 5.8S e ITS2, selecciona los géneros de interés para la micobiota intestinal, prepara los ficheros para el
alineamiento en MEGA y busca zonas conservadas en las que diseñar cebadores para amplificar la región ITS2.

## Scripts, en el orden del flujo de trabajo

| Paso | Script | Qué hace |
|---|---|---|
| 00 | `scripts/00_comprobar_conexion_postgresql.py` | Comprueba la conexión con la base de datos PostgreSQL. |
| 01 | `scripts/01_descargar_secuencias_ncbi.py` | Descarga masiva de las secuencias ITS de NCBI Nucleotide (GenBank y FASTA) en lotes de 500, extrae las posiciones de las regiones anotadas y las inserta en PostgreSQL. |
| 02 | `scripts/02_extraer_its_genomas_gff.py` | Extrae las regiones ITS de genomas de referencia a partir de su anotación GFF3 (con BLAST como alternativa) y las inserta en la base de datos. |
| 03 | `scripts/03_insertar_its_postgresql.py` | Inserta en PostgreSQL las subsecuencias ITS extraídas. |
| 04 | `scripts/04_listar_generos.py` | Obtiene la lista de géneros presentes en las 13.009 secuencias descargadas. |
| 05 | `scripts/05_detectar_regiones_its.py` | Localiza ITS1, 5.8S e ITS2 mediante motivos conservados y guarda sus posiciones y secuencias. |
| 06 | `scripts/06_corregir_orientacion_secuencias.py` | Detecta las secuencias depositadas en orientación inversa y las sustituye por su complementaria inversa. |
| 07 | `scripts/07_exportar_subconjunto_fasta.py` | Exporta a formato FASTA las secuencias seleccionadas. |
| 08 | `scripts/08_descargar_hospedador_y_origen.py` | Descarga de NCBI el origen de aislamiento y el hospedador de las secuencias seleccionadas. |
| 09 | `scripts/09_verificar_genero_ncbi.py` | Comprueba el género de cada secuencia frente a NCBI Taxonomy e identifica los géneros que no están en la lista de referencia. |
| 10 | `scripts/10_descargar_secuencias_revisadas.py` | Vuelve a descargar las secuencias revisadas y genera el FASTA combinado. |
| 11 | `scripts/11_generar_fasta_mega.py` | Genera los ficheros FASTA que se importan en MEGA. |
| 12 | `scripts/12_ordenar_secuencias_por_especie.py` | Ordena las secuencias por especie y género antes del alineamiento. |
| 13 | `scripts/13_ordenar_alineamiento_por_genero.py` | Reordena el alineamiento por género y especie. |
| 14 | `scripts/14_quitar_huecos_alineamiento.py` | Elimina los huecos (gaps) del alineamiento. |
| 15 | `scripts/15_crear_tabla_definitiva.py` | Crea la tabla definitiva en PostgreSQL e importa las secuencias desde el FASTA alineado. |
| 16 | `scripts/16_buscar_cebadores_directos.py` | Analiza la conservación de cada posición del alineamiento, calcula Tm y %GC y selecciona los candidatos a cebador directo. |
| 17 | `scripts/17_buscar_cebadores_inversos.py` | Hace el mismo análisis en la región final del alineamiento para los candidatos a cebador inverso. |
| 18 | `scripts/18_calcular_tamano_amplicon.py` | Escribe los candidatos a cebador inverso en sentido 5'→3' (complementaria inversa) y calcula el tamaño del fragmento que amplifica la pareja de cebadores en cada secuencia del alineamiento. |

## Consultas SQL

| Fichero | Qué hace |
|---|---|
| `sql/01_contar_generos.sql` | Cuenta cuántas secuencias hay de cada género en el conjunto descargado. |
| `sql/02_seleccion_inicial_generos.sql` | Primera selección por nombre de diez géneros de levaduras y mohos. |

## Requisitos

- Python 3.10 o posterior y las librerías de `requirements.txt` (`pip install -r requirements.txt`).
- PostgreSQL con una base de datos local. La contraseña se lee de la variable de entorno `PGPASSWORD`.
- BLAST+ instalado en el sistema (solo para `02_extraer_its_genomas_gff.py`).
- MEGA para el alineamiento múltiple (MUSCLE), la matriz de distancias y el árbol.

Antes de ejecutar los scripts que consultan NCBI, sustituye `tu_correo@ejemplo.com` por tu correo, como piden las normas
de uso de las E-utilities de NCBI.

## Datos

Los datos no se incluyen en el repositorio. Las secuencias se obtienen de NCBI Nucleotide con la búsqueda
`ITS ribosomal fungi`; el script `01_descargar_secuencias_ncbi.py` automatiza la descarga.

## Árbol de distancias (MEGA-CC)

El árbol de la figura 3.3 se reconstruye con MEGA-CC 12 y el fichero de opciones `mega/arbol_nj_distancia_p.mao`
(unión de vecinos, distancia p, eliminación por pares de los huecos):

```
megacc -a mega/arbol_nj_distancia_p.mao -d alineamiento.fas -o arbol_nj
```

Se usa la distancia p porque la distancia de Kimura de 2 parámetros no puede calcularse entre algunos pares de
secuencias muy diferentes del alineamiento.

## Figuras

La carpeta `figuras/` contiene los scripts que generan las figuras 1.1, 3.2, 3.3 y 3.4 del trabajo. Necesitan
matplotlib y Pillow, además de la fuente Arial.

| Script | Figura |
|---|---|
| `figuras/figura_1_1_esquema_its.py` | Esquema de la unidad de ADN ribosómico y posición de los cebadores |
| `figuras/figura_3_2_alineamiento.py` | Fragmento del alineamiento con el cebador directo rotulado |
| `figuras/figura_3_3_arbol.py` | Árbol completo y ampliación de un clado |
| `figuras/figura_3_4_gel.py` | Electroforesis de los productos de PCR |

## Licencia

MIT (ver `LICENSE`).

## Nombres anteriores de los scripts

| Nombre anterior | Nombre actual |
|---|---|
| `test_conexion.py` | `00_comprobar_conexion_postgresql.py` |
| `blast.py` | `01_descargar_secuencias_ncbi.py` |
| `procesar_its_generico.py` | `02_extraer_its_genomas_gff.py` |
| `insertar_its.py` | `03_insertar_its_postgresql.py` |
| `generos.py` | `04_listar_generos.py` |
| `detectar_regiones_its.py` | `05_detectar_regiones_its.py` |
| `micobiota_complementarias.py` | `06_corregir_orientacion_secuencias.py` |
| `micobiota.py` | `07_exportar_subconjunto_fasta.py` |
| `filtrar_gastrointestinal.py` | `08_descargar_hospedador_y_origen.py` |
| `geenros_faltantes.py` | `09_verificar_genero_ncbi.py` |
| `geenrar_meg_definitivo.py` | `10_descargar_secuencias_revisadas.py` |
| `crudo_generar.py` | `11_generar_fasta_mega.py` |
| `automatizar_hongos_its.py` | `12_ordenar_secuencias_por_especie.py` |
| `agrupar_los_generos_para_mega.py` | `13_ordenar_alineamiento_por_genero.py` |
| `quitar_gaps.py` | `14_quitar_huecos_alineamiento.py` |
| `importar_tabla_definitiva.py` | `15_crear_tabla_definitiva.py` |
| `primer.py` | `16_buscar_cebadores_directos.py` |
| `primers_995.py` | `17_buscar_cebadores_inversos.py` |
| `numero de generos.sql` | `01_contar_generos.sql` |
| `gastrointetsinal.sql` | `02_seleccion_inicial_generos.sql` |

## Autor

Jorge Sánchez Germán. Trabajo Fin de Grado, Universidad de Extremadura.
