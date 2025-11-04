"""
Tarea 1 - Análisis funcional de genes
------------------------------------------------------------
Autor: Patricia Rodríguez Lidueña

Este script realiza un análisis funcional de los genes COX4I2, ND1 y ATP6
usando las librerías MyGene y GSEApy.

"""

# ========================================================
# 1. Importación de librerías 
# =========================================================
import argparse
import pandas as pd
import mygene
import gseapy as gp


# =========================================================
# 2. Funciones auxiliares
# =========================================================

def leer_genes(ruta_input):
    """
    Lee el archivo de genes de entrada y devuelve una lista limpia de símbolos.
    """
    with open(ruta_input) as f:
        content = f.read().strip()
    genes = [g.strip() for g in content.replace("\n", ",").split(",") if g.strip()]
    return genes


def obtener_info_mygene(genes):
    """
    Consulta la API MyGene.info para obtener información básica de los genes:
    - Identificadores Ensembl y UniProt
    - Anotaciones de procesos biológicos (GO)
    """
    mg = mygene.MyGeneInfo()
    info = mg.querymany(
        genes,
        scopes="symbol",
        fields="ensembl.gene,uniprot,go.BP.name",
        species="human",
    )
    return pd.DataFrame(info)


def analisis_funcional_enrichr(genes, ruta_salida):
    """
    Realiza el análisis funcional mediante Enrichr (GSEApy).
    Busca términos enriquecidos en bases de datos de procesos biológicos (GO)
    y rutas metabólicas (KEGG).
    """
    enr = gp.enrichr(
        gene_list=genes,
        gene_sets=["GO_Biological_Process_2023", "KEGG_2021_Human"],
        organism="Human",
        outdir="results",   # Carpeta donde se guardan ls resultados 
        cutoff=0.05,
    )

    # Guardar resultados en formato CSV
    enr.results.to_csv(ruta_salida, index=False)

    # Mostrar un resumen con los términos más enriquecidos
    print("\nAnálisis funcional completado correctamente.")
    print("Principales términos enriquecidos:")
    print(enr.results[["Term", "Adjusted P-value"]].head(5))


# =========================================================
# 3. Configuración de argumentos CLI (para ejecurae en la terminal)
# =========================================================

def parse_args():
    """
    Define los argumentos que se pasan al script desde la línea de comandos.
    """
    parser = argparse.ArgumentParser(
        description="Análisis funcional de genes usando MyGene y Enrichr (GSEApy)."
    )
    parser.add_argument("input_file", help="Ruta al archivo de entrada con los genes (ej. data/genes_input.txt).")
    parser.add_argument("output_file", help="Ruta al archivo de salida CSV donde se guardarán los resultados.")
    return parser.parse_args()


# =========================================================
# 4. Función principal
# =========================================================

def main():
    args = parse_args()

    # 1️ Leer genes desde el archivo
    genes = leer_genes(args.input_file)
    print(f"Genes cargados: {genes}")

    # 2️ Obtener información de los genes desde MyGene
    info_df = obtener_info_mygene(genes)
    print("\nInformación básica obtenida de MyGene.info:")
    print(info_df[["query", "ensembl", "uniprot", "notfound"]])

    # 3️ Realizar análisis funcional con Enrichr
    analisis_funcional_enrichr(genes, args.output_file)

    print(f"\nResultados guardados en: {args.output_file}")


# =========================================================
# 5. Punto de entrada del script
# =========================================================

if __name__ == "__main__":
    main()
