# -*- coding: utf-8 -*-
import os
import psycopg2
import sys
try:
    conn = psycopg2.connect(
        host="localhost",
        port="5432",
        database="tfg_Its",
        user="postgres",
        password=os.environ.get("PGPASSWORD", ""),
        client_encoding="UTF8"
    )
    print("Conexion exitosa")
    conn.close()
except Exception as e:
    print(f"Error al conectar: {e}")
    sys.exit(1)