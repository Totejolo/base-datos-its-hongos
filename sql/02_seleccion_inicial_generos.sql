-- 3. CREAR TABLA GASTROINTESTINAL (SIN METADATOS AÚN)
CREATE TABLE hongos_gastrointestinales AS
SELECT *
FROM secuencias_its_con_posiciones
WHERE 
    LOWER(especie) LIKE ANY(ARRAY[
        '%candida%', '%saccharomyces%', '%malassezia%', '%debaryomyces%',
        '%pichia%', '%yarrowia%', '%rhodotorula%', '%cryptococcus%',
        '%penicillium%', '%aspergillus%'
    ]);

-- 4. VER CUÁNTAS SALIERON
SELECT COUNT(*) FROM hongos_gastrointestinales;

-- 5. EXPORTAR A CSV
COPY hongos_gastrointestinales 
TO 'ruta/a/tu/carpeta/hongos_gastrointestinales.csv' 
WITH (FORMAT CSV, HEADER TRUE, DELIMITER ',', ENCODING 'UTF8');