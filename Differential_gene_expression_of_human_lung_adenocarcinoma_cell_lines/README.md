# Introducción

Este proyecto tiene como objetivo profundizar en los aspectos básicos de **RNA-seq** y las buenas prácticas para el análisis transcriptómico utilizando la plataforma Galaxy. A través de este análisis, se busca **desarrollar conocimientos avanzados sobre el análisis de expresión diferencial**, abordando temas clave como:

- La creación de gráficos para evaluar el estado de las muestras antes del análisis.
- La visualización de los efectos del filtrado de genes de baja expresión mediante gráficos específicos.
- El uso de gráficos como **MDS**, **PCA** y **barplots** para explorar la distribución del p-valor y otros parámetros relevantes.
- Aplicación de las herramientas de **EdgeR** como `plotBCV` y gráficos relacionados con la binomial negativa, examinando y explicando cómo estas afectan a los datos.

Además, se busca experimentar con el método de **pseudomapping** utilizando la herramienta **Salmon** como parte del flujo de trabajo. También se incluye un análisis básico de expresión diferencial (DGE) empleando **R**.

El proyecto utiliza como referencia los datos del estudio [*Long and short-read transcriptome profiling of human lung cancer cell lines*](https://pmc.ncbi.nlm.nih.gov/articles/PMC7545141/), proporcionando un contexto práctico para ilustrar cada aspecto abordado.

## Información de la muestra

Se basa en datos transcriptómicos generados en el marco del estudio titulado [*Long and short-read transcriptome profiling of human lung cancer cell lines*](https://pmc.ncbi.nlm.nih.gov/articles/PMC7545141/) de Dong et al. (2020). Este estudio examinó los perfiles transcriptómicos de las líneas celulares **HCC827** y **H1975**. Estas son líneas celulares presentan diferentes mutaciones en el gen **EGFR**.

Los datos utilizados incluyen lecturas _paired end_ de **RNA-Seq** generadas con tecnología **Illumina™**. Estos datos están disponibles en el repositorio SRA:

- **H1975**:
  - [SRR14286057](https://www.ncbi.nlm.nih.gov/sra/?term=SRR14286057)
  - [SRR14286058](https://www.ncbi.nlm.nih.gov/sra/?term=SRR14286058)
  
- **HCC_827**:
  - [SRR14286065](https://www.ncbi.nlm.nih.gov/sra/?term=SRR14286065)
  - [SRR14286066](https://www.ncbi.nlm.nih.gov/sra/?term=SRR14286066)

# Metodología

## Análisis de Control de Calidad

Para no alargar este README, se informa que las lecturas presentan una buena calidad. Los valores de GC están dentro de lo esperado para datos de RNA-seq, los duplicados son aceptables, el contenido de adaptadores es despreciable, y el número de secuencias es adecuado para un análisis de expresión diferencial. 

Las imágenes relacionadas con el control de calidad se pueden visualizar en el siguiente enlace: [Control de Calidad](./Differential_gene_expression_of_human_lung_adenocarcinoma_cell_lines/Control_Calidad)

## Pseudomapping

Se utilizó la configuración predeterminada de la herramienta **Salmon**, ajustando el parámetro **Kmer length** a un valor de 31. Esta decisión se basó en las recomendaciones de la [documentación oficial de Salmon](https://salmon.readthedocs.io/en/latest/salmon.html#preparing-transcriptome-indices-mapping-based-mode), donde indican _We find that a k of 31 seems to work well for reads of 75bp or longer_, en este caso las lecturas tienen 80bp de longitud, así que cumplen esta condición.

El análisis se llevó a cabo utilizando el transcriptoma de **GENCODE v47**, que se encuentra disponible en el siguiente enlace: [GENCODE v47](https://www.gencodegenes.org/human/).

Una vez generados los resultados con **Salmon**, estos se procesaron en **R** para llevar a cabo análisis adicionales, como la exploración y el análisis diferencial de expresión génica.

## Exploración de los datos

Para la exploración de los datos, se definieron los nombres de las muestras y sus grupos experimentales. Se localizaron los archivos de cuantificación de Salmon y se generó un mapa de transcritos a genes basado en GENCODE v47. Los datos se importaron con tximport para obtener una matriz de conteos, que se limpió y renombró. Finalmente, se anotaron los genes utilizando org.Hs.eg.db, eliminando duplicados y genes sin anotación válida. Los datos están preparados para el análisis DGE.

### Gráfica preliminar

Para explorar las características iniciales de los datos RNA-seq antes de aplicar análisis diferenciales o correcciones metodológicas, se hizo uso de una gráfica de tipo violin. En el eje **y** se encuentran los valores de expresión en **log2 Counts per Million (logCPM)**, mientras que el eje **x** muestra las distintas muestras: **H1975_1**, **H1975_2**, **HCC827_1** y **HCC827_2**.

Cada violín ilustra la densidad de los datos de expresión para cada muestra, proporcionando una idea de cómo se distribuyen los genes. La línea negra en el centro de cada violín señala la **mediana** de la expresión, que sirve como referencia para los niveles centrales de cada muestra.

<p align="center">
  <img src="./Differential_gene_expression_of_human_lung_adenocarcinoma_cell_lines/Recursos/LogCPMnofiltnonTMM.png">
</p>

Se observa una acumulación significativa de genes con baja expresión por debajo de la mediana. Estos genes suelen introducir ruido en el análisis y no aportan información biológicamente relevante. El paquete **edgeR** tiene una opción denominada `filterbyExpr` para eliminar genes de baja expresión y asegurar que solo se incluyan aquellos que están biológicamente activos y relevantes para el análisis de expresión diferencial.

### Gráfica tras aplicar `filterByExpr`

Tras aplicar el filtro `filterByExpr`, se observan cambios significativos en la distribución de los datos RNA-seq, reflejados en la gráfica:

<p align="center">
  <img src="./Differential_gene_expression_of_human_lung_adenocarcinoma_cell_lines/Recursos/LogCPMfilterednonTMM.png">
</p>

Los genes con baja expresión, que generan ruido en los datos y afectan la interpretación, han sido eliminados. Esto se traduce en unas distribuciones más homogéneas y limpias. Tras el filtrado de genes, la mediana está ajustada a valores de expresión logCPM donde los genes ya no son ruidosos, proporcionando una referencia adecuada para comparaciones entre muestras.
Aunque los datos han sido filtrados, es necesario proceder con la normalización por TMM para corregir diferencias técnicas entre muestras, asegurando que las variaciones observadas sean de origen biológico.

### Gráfica tras aplicar TMM

Después de aplicar la normalización por TMM (Trimmed Mean of M-values), los datos RNA-seq muestran una distribución ajustada y comparable entre las muestras, reflejada en la gráfica:

<p align="center">
  <img src="./Differential_gene_expression_of_human_lung_adenocarcinoma_cell_lines/Recursos/LogCPMfilteredTMM.png">
</p>

Las distribuciones de logCPM son ahora más consistentes entre las muestras, lo que indica que se han corregido las diferencias técnicas, como la profundidad de secuenciación o la composición de las librerías.
Las líneas negras que representan las medianas de las distribuciones están ahora alineadas entre las muestras, mostrando que la normalización ha equilibrado los niveles de expresión para permitir comparaciones confiables.

Gracias a TMM, las variaciones no biológicas entre las muestras han sido mitigadas, garantizando que cualquier diferencia observada sea mayoritariamente de origen biológico.


El plot MDS (Multidimensional Scaling Plot) es una herramienta para visualizar las relaciones entre las muestras basándose en sus perfiles de expresión génica.










  
