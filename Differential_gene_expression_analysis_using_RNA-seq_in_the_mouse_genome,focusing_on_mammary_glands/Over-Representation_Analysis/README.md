# Análisis de Sobre-representación en Genes Diferencialmente Expresados (DEG)

## Introducción

En un análisis de expresión diferencial, identificamos genes que presentan cambios significativos en su expresión entre condiciones experimentales. Sin embargo, el listado de DEG por sí solo no proporciona información sobre los procesos biológicos afectados. Para interpretar estos resultados en un contexto funcional, se realiza un análisis de sobre-representación (ORA, por sus siglas en inglés: Over-Representation Analysis).

## ¿Por qué realizar un análisis de sobre-representación?

El ORA permite responder preguntas clave, como:
- ¿Cuáles son las rutas metabólicas o procesos biológicos enriquecidos en los DEG?
- ¿Cuáles son las categorías funcionales más representadas en los genes regulados al alza o a la baja?

Este análisis se basa en comparar el conjunto de DEG con bases de datos de referencia como Gene Ontology, KEGG o Reactome para identificar si algún conjunto de genes está significativamente representado en la lista de genes estudiados.

## Pasos del análisis de sobre-representación

1. **Preprocesamiento de datos:**
   - Obtención de DEG.

2. **Filtrado de genes:**
   - Eliminación de genes no significativos según valores de p ajustados, en este caso **FDR > 0.05**.
   - Clasificación en genes regulados al alza y a la baja.

3. **Definición del fondo de genes:**
   - Se usa el conjunto total de genes medidos en el experimento como referencia para evitar sesgos.

4. **Análisis de enriquecimiento:**
   - Se emplean herramientas como `clusterProfiler` en R para realizar la prueba de sobre-representación.
   - Se comparan los DEG con bases de datos de vías y procesos biológicos, en este caso con GO.

5. **Filtrado de resultados:**
   - Se seleccionan solo los términos significativamente enriquecidos según un umbral de p-ajustada (< 0.05) y un mínimo de genes por vía.

6. **Visualización de resultados:**
   - Se generan gráficos como mapas de calor, gráficos de dot plot o redes de interacciones.

## Implementación en R

Para realizar este análisis en R, se utilizó el paquete `clusterProfiler`, basado en el tutorial de BiostatSquid [(enlace al tutorial)](https://biostatsquid.com/pathway-enrichment-analysis-tutorial-clusterprofiler/).

El código completo utilizado para este análisis se encuentra en el siguiente script: 
- [clusterProfiler.R](./Recursos/clusterProfiler.R)

En este archivo se incluyen todos los pasos descritos anteriormente para realizar el análisis de sobre-representación en un conjunto de DEG.

Quiero agradecer a Laura, la persona detrás de biostatsquid, por su dedicación al crear estos tutoriales. Su esfuerzo por enseñar y divulgar conceptos complejos de forma sencilla es de gran ayuda para bioinformáticos recién graduados.

(Dentro de nada, actualizaré este readme para visualizar los [resultados obtenidos](./Recursos/GO.RDS_resclusterp.csv) en el script.) 
---


