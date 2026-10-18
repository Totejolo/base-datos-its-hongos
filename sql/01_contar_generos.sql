-- TODOS LOS GÉNEROS DISTINTOS + CONTEO → LISTO PARA CSV
SELECT
    TRIM(BOTH '[]' FROM SPLIT_PART(especie, ' ', 1)) AS genero,
    COUNT(*) AS num_especies
FROM secuencias_its_con_posiciones
WHERE especie IS NOT NULL 
  AND especie != ''
  AND SPLIT_PART(especie, ' ', 1) != ''
GROUP BY genero
ORDER BY genero ASC;