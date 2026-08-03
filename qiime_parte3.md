# Guía paso a paso de QIIME2 — Parte 3

*Entrenamiento de un clasificador de taxonomía (SILVA, V3-V4) y asignación de taxonomía*

> Antes de poder asignar taxonomía a `rep-seqs.qza`, se necesita un **clasificador entrenado específicamente para la región y los cebadores usados** (341F/806R, V3-V4). Un clasificador entrenado sobre el gen 16S completo, o sobre otra región (p. ej. V4 con 515F/806R), da resultados menos precisos que uno entrenado solo sobre la región realmente amplificada.
>
> Como el entorno conda es `qiime2-amplicon-2024.10`, el objetivo correcto es **SILVA 138.2 SSURef_NR99**.
>
> La descarga y el entrenamiento se hacen **una sola vez** — el clasificador resultante (`SILVA_138.2_V3V4_341F_806R_CLASSIFIER.qza`) se reutiliza para cualquier muestra futura con los mismos cebadores, no hay que repetir este proceso por proyecto. Todos los archivos intermedios y el clasificador viven en el subdirectorio `classifier_SILVA/`, separado de los outputs de tus muestras.

---

## Índice

- [11. Descarga y prepara SILVA](#11-descarga-y-prepara-silva)
- [12. Reverse-transcribe](#12-reverse-transcribe)
- [13. Extraer región V3-V4 con nuestros primers](#13-extraer-región-v3-v4-con-nuestros-primers)
- [14. Entrenar el clasificador](#14-entrenar-el-clasificador-el-más-pesado)
- [15. Asignar taxonomía a tus ASVs](#15-asignar-taxonomía-a-tus-asvs)
- [16. Filtrado y limpieza de la tabla](#16-filtrado-y-limpieza-de-la-tabla)
- [17. Árbol filogenético y diversidad (ASVs) — pendiente](#17-árbol-filogenético-y-diversidad-asvs--pendiente)
- [18. Convertir ASVs a OTUs — pendiente](#18-convertir-asvs-a-otus--pendiente)
- [19. Diversidad para OTUs — pendiente](#19-diversidad-para-otus--pendiente)
- [Recursos utilizados](#recursos-utilizados)

---

## 11. Descarga y prepara SILVA

En el directorio de scripts, crea el script de descarga:

```bash
nano 07_get_silva_data_v1.sh
```

Contenido del script:

```bash
#!/bin/bash
#SBATCH -p cu
#SBATCH -n 4
#SBATCH -N 1
#SBATCH -D /Tayra-Share/home/Universidades/UCaldas/Users/eariasa/Rhogeessa
#SBATCH -J silva_download
#SBATCH --mail-user=edwin.2051826349@ucaldas.edu.co
#SBATCH --mail-type=END,FAIL
#SBATCH --mem=4G
#SBATCH -o /Tayra-Share/home/Universidades/UCaldas/Users/eariasa/Rhogeessa/analysis_outputs/logs/silva_download_%j.out
#SBATCH -e /Tayra-Share/home/Universidades/UCaldas/Users/eariasa/Rhogeessa/analysis_outputs/logs/silva_download_%j.err
#SBATCH --time=01:00:00

source /Tayra-Share/miniconda/etc/profile.d/conda.sh
conda activate qiime2-amplicon-2024.10

OUTPUT_DIR=/Tayra-Share/home/Universidades/UCaldas/Users/eariasa/Rhogeessa/analysis_outputs/classifier_SILVA

mkdir -p $OUTPUT_DIR

qiime rescript get-silva-data \
  --p-version 138.2 \
  --o-silva-sequences $OUTPUT_DIR/SILVA_138.2_SSURef_NR99_RNAseq.qza \
  --o-silva-taxonomy $OUTPUT_DIR/taxmap_slv_ssu_ref_nr_138.2.qza
```

### Explicación de los flags

- **`qiime rescript get-silva-data`**
  Acción del plugin q2-rescript que descarga y empaqueta directamente en formato QIIME2 la base de datos SILVA (secuencias + taxonomía).

- **`--p-version 138.2`**
  La versión de SILVA a descargar. Debe coincidir con la que soporta tu entorno conda (`qiime2-amplicon-2024.10` → 138.2). Al no especificar `--p-target`, se usa el valor por defecto del plugin (`SSURef_NR99`, la referencia no-redundante al 99% — la variante estándar para 16S bacteriano).

- **`--o-silva-sequences`**
  Salida en formato `FeatureData[RNASequence]` — nota que son secuencias de **ARN**, no de ADN todavía, por eso el siguiente paso las transcribe.

- **`--o-silva-taxonomy`**
  El mapa de taxonomía correspondiente a cada secuencia (necesario más adelante para entrenar el clasificador).

---

## 12. Reverse-transcribe

En el directorio de scripts, crea el script de transcripción:

```bash
nano 08_reverse_transcribe_v1.sh
```

Contenido del script:

```bash
#!/bin/bash
#SBATCH -p cu
#SBATCH -n 4
#SBATCH -N 1
#SBATCH -D /Tayra-Share/home/Universidades/UCaldas/Users/eariasa/Rhogeessa
#SBATCH -J silva_reverse_transcribe
#SBATCH --mail-user=edwin.2051826349@ucaldas.edu.co
#SBATCH --mail-type=END,FAIL
#SBATCH --mem=2G
#SBATCH -o /Tayra-Share/home/Universidades/UCaldas/Users/eariasa/Rhogeessa/analysis_outputs/logs/silva_reverse_transcribe_%j.out
#SBATCH -e /Tayra-Share/home/Universidades/UCaldas/Users/eariasa/Rhogeessa/analysis_outputs/logs/silva_reverse_transcribe_%j.err
#SBATCH --time=00:30:00

source /Tayra-Share/miniconda/etc/profile.d/conda.sh
conda activate qiime2-amplicon-2024.10

CLASSIFIER_DIR=/Tayra-Share/home/Universidades/UCaldas/Users/eariasa/Rhogeessa/analysis_outputs/classifier_SILVA

qiime rescript reverse-transcribe \
  --i-rna-sequences $CLASSIFIER_DIR/SILVA_138.2_SSURef_NR99_RNAseq.qza \
  --o-dna-sequences $CLASSIFIER_DIR/SILVA_138.2_SSURef_NR99_DNA.qza
```

### Explicación de los flags

- **`qiime rescript reverse-transcribe`**
  Convierte las bases U a T en cada secuencia.

- **`--i-rna-sequences`**
  El artefacto `FeatureData[RNASequence]` descargado en el paso 11.

- **`--o-dna-sequences`**
  El mismo conjunto de secuencias, ahora como `FeatureData[Sequence]` (ADN), listo para que `extract-reads` pueda buscar los cebadores en él.

---

## 13. Extraer región V3-V4 con nuestros primers

En el directorio de scripts, crea el script de extracción de región:

```bash
nano 09_extract_v3v4_v1.sh
```

Contenido del script:

```bash
#!/bin/bash
#SBATCH -p cu
#SBATCH --cpus-per-task=12
#SBATCH -N 1
#SBATCH -D /Tayra-Share/home/Universidades/UCaldas/Users/eariasa/Rhogeessa
#SBATCH -J 10_extract_v3v4
#SBATCH --mail-user=edwin.2051826349@ucaldas.edu.co
#SBATCH --mail-type=END,FAIL
#SBATCH --mem=80G
#SBATCH -o /Tayra-Share/home/Universidades/UCaldas/Users/eariasa/Rhogeessa/analysis_outputs/logs/extract_v3v4_%j.out
#SBATCH -e /Tayra-Share/home/Universidades/UCaldas/Users/eariasa/Rhogeessa/analysis_outputs/logs/extract_v3v4_%j.err
#SBATCH --time=03:00:00

source /Tayra-Share/miniconda/etc/profile.d/conda.sh
conda activate qiime2-amplicon-2024.10

CLASSIFIER_DIR=/Tayra-Share/home/Universidades/UCaldas/Users/eariasa/Rhogeessa/analysis_outputs/classifier_SILVA

qiime feature-classifier extract-reads \
  --i-sequences $CLASSIFIER_DIR/SILVA_138.2_SSURef_NR99_DNA.qza \
  --p-f-primer CCTAYGGGDBGCWGCAG \
  --p-r-primer GAMTACNVGGGTHTCTAATCC \
  --p-min-length 100 \
  --p-max-length 550 \
  --p-n-jobs 12 \
  --o-reads $CLASSIFIER_DIR/SILVA_138.2_V3V4_341F_806R_CUSTOM.qza
```

### Explicación de los flags

- **`qiime feature-classifier extract-reads`**
  Busca la posición de nuestros cebadores dentro de las secuencias de referencia completas y recorta solo la región amplificada entre ellos — así el clasificador se entrena únicamente sobre la parte del gen que realmente se secuenció.

- **`--p-f-primer CCTAYGGGDBGCWGCAG`**
  Cebador forward 341F — el mismo usado en `cutadapt` en la parte 2. Debe ser idéntico para que la región extraída coincida con las lecturas reales.

- **`--p-r-primer GAMTACNVGGGTHTCTAATCC`**
  Cebador reverse 806R, igual al usado en `cutadapt`.

- **`--p-min-length 100`**
  Longitud mínima (en pb) que debe tener una secuencia ya extraída para conservarse. Nuestros primers tienen bases IUPAC ambiguas (Y, D, B, W en el forward; N, V, H en el reverse), que a veces hacen match en un lugar equivocado del genoma de referencia y generan fragmentos "extraídos" muy cortos que no son el amplicón V3-V4 real, sino ruido. Este filtro descarta esa basura antes de que llegue al entrenamiento.

- **`--p-max-length 550`**
  Lo mismo que `--p-min-length`, pero por el lado superior — descarta fragmentos anormalmente largos (posibles matches espurios o secuencias mal anotadas en SILVA). El rango 100-550 pb es generoso alrededor del amplicón V3-V4 esperado.

- **`--o-reads`**
  El artefacto de salida: solo el fragmento V3-V4 de cada secuencia de referencia, listo para entrenar el clasificador. El nombre del archivo (`V3V4_341F_806R_CUSTOM`) deja claro que es un clasificador hecho a la medida de nuestros cebadores, no uno genérico descargado ya entrenado.

- **`--p-n-jobs`**
  Número de procesos en paralelo que `extract-reads` lanza internamente para buscar los cebadores en las secuencias de referencia. Por defecto (`--p-n-jobs 1`) el comando es de un solo hilo, sin importar cuántas CPUs le pidas a Slurm — este flag es el que realmente le dice al programa que reparta el trabajo entre varios núcleos. Debe ir acompañado de una reserva de CPUs equivalente en el script (`--cpus-per-task=12` en este caso).
  > **Nota:** recomiendo asignar más recursos de cómputo; tiempo de ejecución real: ~25 min.

---

## 14. Entrenar el clasificador (el más pesado)

En el directorio de scripts, crea el script de entrenamiento:

```bash
nano 10_silva_train_classifier.sh
```

Contenido del script:

```bash
#!/bin/bash
#SBATCH -p cu
#SBATCH --cpus-per-task=4
#SBATCH -N 1
#SBATCH -D /Tayra-Share/home/Universidades/UCaldas/Users/eariasa/Rhogeessa
#SBATCH -J silva_train_classifier
#SBATCH --mail-user=edwin.2051826349@ucaldas.edu.co
#SBATCH --mail-type=END,FAIL
#SBATCH --mem=64G
#SBATCH -o /Tayra-Share/home/Universidades/UCaldas/Users/eariasa/Rhogeessa/analysis_outputs/logs/silva_train_classifier_%j.out
#SBATCH -e /Tayra-Share/home/Universidades/UCaldas/Users/eariasa/Rhogeessa/analysis_outputs/logs/silva_train_classifier_%j.err
#SBATCH --time=08:00:00

source /Tayra-Share/miniconda/etc/profile.d/conda.sh
conda activate qiime2-amplicon-2024.10

CLASSIFIER_DIR=/Tayra-Share/home/Universidades/UCaldas/Users/eariasa/Rhogeessa/analysis_outputs/classifier_SILVA

qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads $CLASSIFIER_DIR/SILVA_138.2_V3V4_341F_806R_CUSTOM.qza \
  --i-reference-taxonomy $CLASSIFIER_DIR/taxmap_slv_ssu_ref_nr_138.2.qza \
  --o-classifier $CLASSIFIER_DIR/SILVA_138.2_V3V4_341F_806R_CLASSIFIER.qza
```

### Explicación de los flags

- **`qiime feature-classifier fit-classifier-naive-bayes`**
  Entrena un clasificador Naive Bayes multinomial usando las secuencias de referencia recortadas (paso 13) y su taxonomía asociada, para luego poder predecir la taxonomía de nuestros ASVs (`rep-seqs.qza`).

- **`--i-reference-reads`**
  Las secuencias de referencia ya recortadas a la región V3-V4 (salida del paso 13).

- **`--i-reference-taxonomy`**
  El mapa de taxonomía de SILVA descargado en el paso 11 (no cambia entre pasos, es el mismo para toda la base).

- **`--o-classifier`**
  El clasificador entrenado y listo, en formato `.qza`. Este es el archivo que se usa en el paso 15 contra `rep-seqs.qza` para asignar taxonomía real a los ASVs.

---

## 15. Asignar taxonomía a tus ASVs

En el directorio de scripts, crea el script de clasificación:

```bash
nano 11_qiime2_classify.sh
```

Contenido del script:

```bash
#!/bin/bash
#SBATCH -p cu
#SBATCH --cpus-per-task=8
#SBATCH -N 1
#SBATCH -D /Tayra-Share/home/Universidades/UCaldas/Users/eariasa/Rhogeessa
#SBATCH -J 11_quiime2_classify
#SBATCH --mail-user=edwin.2051826349@ucaldas.edu.co
#SBATCH --mail-type=END,FAIL
#SBATCH --mem=40G
#SBATCH -o /Tayra-Share/home/Universidades/UCaldas/Users/eariasa/Rhogeessa/analysis_outputs/logs/qiime2_classify_%j.out
#SBATCH -e /Tayra-Share/home/Universidades/UCaldas/Users/eariasa/Rhogeessa/analysis_outputs/logs/qiime2_classify_%j.err
#SBATCH --time=02:00:00

source /Tayra-Share/miniconda/etc/profile.d/conda.sh
conda activate qiime2-amplicon-2024.10

OUTPUT_DIR=/Tayra-Share/home/Universidades/UCaldas/Users/eariasa/Rhogeessa/analysis_outputs
CLASSIFIER_DIR=$OUTPUT_DIR/classifier_SILVA

qiime feature-classifier classify-sklearn \
  --i-classifier $CLASSIFIER_DIR/SILVA_138.2_V3V4_341F_806R_CLASSIFIER.qza \
  --i-reads $OUTPUT_DIR/rep-seqs.qza \
  --p-n-jobs 8 \
  --o-classification $OUTPUT_DIR/taxonomy.qza

qiime metadata tabulate \
  --m-input-file $OUTPUT_DIR/taxonomy.qza \
  --o-visualization $OUTPUT_DIR/taxonomy.qzv

qiime taxa barplot \
  --i-table $OUTPUT_DIR/table.qza \
  --i-taxonomy $OUTPUT_DIR/taxonomy.qza \
  --o-visualization $OUTPUT_DIR/barplot.qzv
```

### Explicación de los flags

- **`qiime feature-classifier classify-sklearn`**
  Aplica el clasificador Naive Bayes entrenado (paso 14) sobre cada ASV de `rep-seqs.qza`, prediciendo su taxonomía (dominio → especie, según qué tan lejos llegue la resolución de SILVA para cada secuencia).

- **`--i-classifier`**
  El clasificador entrenado a la medida de tus cebadores (`V3V4_341F_806R`), no un clasificador genérico.

- **`--i-reads $OUTPUT_DIR/rep-seqs.qza`**
  Tus propias secuencias representativas (salida de DADA2 en la parte 2, paso 9), no las de referencia de SILVA.

- **`--o-classification $OUTPUT_DIR/taxonomy.qza`**
  Tabla de asignación taxonómica por ASV: para cada `Feature ID` de `rep-seqs.qza`, la cadena taxonómica asignada y su nivel de confianza.

- **`qiime metadata tabulate` / `--o-visualization taxonomy.qzv`**
  Convierte `taxonomy.qza` en una tabla visualizable en [view.qiime2.org](https://view.qiime2.org/), igual que se hizo con `denoising-stats.qzv` en la parte 2 (paso 10).

---

## 16. Filtrado y limpieza de la tabla

En el directorio de scripts, crea el script de filtrado y limpieza:

```bash
nano 12_filter_and_cleanup_v1.sh
```

Contenido del script:

```bash
#!/bin/bash
#SBATCH -p cu
#SBATCH --cpus-per-task=4
#SBATCH -N 1
#SBATCH -D /Tayra-Share/home/Universidades/UCaldas/Users/eariasa/Rhogeessa
#SBATCH -J 12_filter_cleanup
#SBATCH --mail-user=edwin.2051826349@ucaldas.edu.co
#SBATCH --mail-type=END,FAIL
#SBATCH --mem=40G
#SBATCH -o /Tayra-Share/home/Universidades/UCaldas/Users/eariasa/Rhogeessa/analysis_outputs/logs/12_filter_cleanup_%j.out
#SBATCH -e /Tayra-Share/home/Universidades/UCaldas/Users/eariasa/Rhogeessa/analysis_outputs/logs/12_filter_cleanup_%j.err
#SBATCH --time=01:00:00

source /Tayra-Share/miniconda/etc/profile.d/conda.sh
conda activate qiime2-amplicon-2024.10

OUTPUT_DIR=/Tayra-Share/home/Universidades/UCaldas/Users/eariasa/Rhogeessa/analysis_outputs

# 1. Filtrar solo bacterias
qiime taxa filter-table \
  --i-table $OUTPUT_DIR/table.qza \
  --i-taxonomy $OUTPUT_DIR/taxonomy.qza \
  --p-include "d__Bacteria" \
  --p-exclude "Mitochondria,Chloroplast" \
  --o-filtered-table $OUTPUT_DIR/table-bacteria-clean.qza

# 2. Filtrar secuencias representativas
qiime feature-table filter-seqs \
  --i-data $OUTPUT_DIR/rep-seqs.qza \
  --i-table $OUTPUT_DIR/table-bacteria-clean.qza \
  --o-filtered-data $OUTPUT_DIR/rep-seqs-bacteria-clean.qza

# 3. Resumen tabla limpia
qiime feature-table summarize \
  --i-table $OUTPUT_DIR/table-bacteria-clean.qza \
  --o-visualization $OUTPUT_DIR/table-bacteria-clean.qzv

# 4. Secuencias limpias tabuladas
qiime feature-table tabulate-seqs \
  --i-data $OUTPUT_DIR/rep-seqs-bacteria-clean.qza \
  --o-visualization $OUTPUT_DIR/rep-seqs-bacteria-clean.qzv

# 5. Barplot limpio (solo bacterias)
qiime taxa barplot \
  --i-table $OUTPUT_DIR/table-bacteria-clean.qza \
  --i-taxonomy $OUTPUT_DIR/taxonomy.qza \
  --o-visualization $OUTPUT_DIR/barplot-bacteria-clean.qzv
```

### Explicación de los flags

- **`qiime taxa filter-table`**
  Filtra la tabla de features conservando solo lo que coincide con `--p-include` y descartando lo que coincide con `--p-exclude` — aquí, se queda con lecturas de bacterias y descarta contaminación de origen mitocondrial o de cloroplasto.

- **`qiime feature-table filter-seqs`**
  Aplica el mismo filtro, pero sobre las secuencias representativas (`rep-seqs.qza`), para que la tabla y las secuencias queden sincronizadas tras la limpieza.

- **`qiime feature-table summarize` / `qiime feature-table tabulate-seqs`**
  Generan visualizaciones (`.qzv`) de la tabla y las secuencias ya filtradas, para revisar cuántas features y lecturas sobrevivieron a la limpieza.

- **`qiime taxa barplot`**
  Genera el barplot de composición taxonómica, ya usando solo la tabla limpia de bacterias.

---

## 17. Árbol filogenético y diversidad (ASVs) — pendiente

> **Estado: no ejecutado.** Estos pasos requieren varias muestras para tener sentido (un `--p-sampling-depth` y una comparación de grupos no aplican con una sola muestra). Se deja aquí como plantilla, en el mismo formato de script que el resto de la guía, para correr tal cual cuando el proyecto tenga más muestras.

En el directorio de scripts, crea el script de árbol y diversidad:

```bash
nano 13_phylogeny_diversity_v1.sh
```

Contenido del script (plantilla):

```bash
#!/bin/bash
#SBATCH -p cu
#SBATCH --cpus-per-task=8
#SBATCH -N 1
#SBATCH -D /Tayra-Share/home/Universidades/UCaldas/Users/eariasa/Rhogeessa
#SBATCH -J 13_phylogeny_diversity
#SBATCH --mail-user=edwin.2051826349@ucaldas.edu.co
#SBATCH --mail-type=END,FAIL
#SBATCH --mem=40G
#SBATCH -o /Tayra-Share/home/Universidades/UCaldas/Users/eariasa/Rhogeessa/analysis_outputs/logs/13_phylogeny_diversity_%j.out
#SBATCH -e /Tayra-Share/home/Universidades/UCaldas/Users/eariasa/Rhogeessa/analysis_outputs/logs/13_phylogeny_diversity_%j.err
#SBATCH --time=03:00:00

source /Tayra-Share/miniconda/etc/profile.d/conda.sh
conda activate qiime2-amplicon-2024.10

OUTPUT_DIR=/Tayra-Share/home/Universidades/UCaldas/Users/eariasa/Rhogeessa/analysis_outputs

# 1. Alinear y construir el árbol filogenético
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences $OUTPUT_DIR/rep-seqs-bacteria-clean.qza \
  --o-alignment $OUTPUT_DIR/aligned-rep-seqs.qza \
  --o-masked-alignment $OUTPUT_DIR/masked-aligned-rep-seqs.qza \
  --o-tree $OUTPUT_DIR/unrooted-tree.qza \
  --o-rooted-tree $OUTPUT_DIR/rooted-tree.qza

# 2. Métricas de diversidad (alfa y beta) a partir del árbol
qiime diversity core-metrics-phylogenetic \
  --i-phylogeny $OUTPUT_DIR/rooted-tree.qza \
  --i-table $OUTPUT_DIR/table-bacteria-clean.qza \
  --p-sampling-depth 22000 \
  --m-metadata-file $OUTPUT_DIR/seqs_and_metadata/metadata2.txt \
  --output-dir $OUTPUT_DIR/core-metrics-results

# 3. Significancia de diversidad alfa — Faith's PD
qiime diversity alpha-group-significance \
  --i-alpha-diversity $OUTPUT_DIR/core-metrics-results/faith_pd_vector.qza \
  --m-metadata-file $OUTPUT_DIR/seqs_and_metadata/metadata2.txt \
  --o-visualization $OUTPUT_DIR/core-metrics-results/faith-pd-group-significance.qzv

# 4. Significancia de diversidad alfa — Features observadas
qiime diversity alpha-group-significance \
  --i-alpha-diversity $OUTPUT_DIR/core-metrics-results/observed_features_vector.qza \
  --m-metadata-file $OUTPUT_DIR/seqs_and_metadata/metadata2.txt \
  --o-visualization $OUTPUT_DIR/core-metrics-results/observed_features-group-significance.qzv
```

### Explicación de los flags

- **`qiime phylogeny align-to-tree-mafft-fasttree`**
  Pipeline de 4 pasos en uno: alinea las secuencias representativas con MAFFT, enmascara las posiciones muy variables del alineamiento, y construye un árbol filogenético con FastTree (primero sin raíz, luego enraizado). Este árbol es lo que le da el componente "filogenético" a las métricas de diversidad (p. ej. Faith's PD, UniFrac).

- **`--i-sequences`**
  Las secuencias representativas ya filtradas a solo bacterias (salida del paso 16).

- **`--o-alignment` / `--o-masked-alignment`**
  El alineamiento múltiple crudo y su versión enmascarada (columnas de baja calidad o muy variables removidas), respectivamente. Se guardan como artefactos intermedios, normalmente no se inspeccionan directamente.

- **`--o-tree` / `--o-rooted-tree`**
  El árbol sin raíz (topología simple) y el árbol enraizado (necesario para métricas como Faith's PD y UniFrac, que requieren un punto de referencia común).

- **`qiime diversity core-metrics-phylogenetic`**
  Calcula de una sola vez un conjunto estándar de métricas de diversidad alfa (dentro de cada muestra) y beta (entre muestras), tanto filogenéticas como no filogenéticas.

- **`--i-phylogeny`**
  El árbol enraizado generado arriba.

- **`--i-table`**
  La tabla de features ya filtrada a bacterias (salida del paso 16), con el conteo de cada ASV por muestra.

- **`--p-sampling-depth 22000`**
  Profundidad de rarefacción: el número de lecturas al que se submuestrea cada muestra antes de calcular diversidad, para que todas las muestras se comparen en igualdad de condiciones. Se elige normalmente mirando la distribución de profundidad por muestra en `table-bacteria-clean.qzv` (paso 16) y tomando un valor que conserve la mayoría de las muestras sin descartar demasiadas lecturas. El valor de 22000 conviene revisarlo contra tus propias muestras antes de correrlo.

- **`--m-metadata-file`**
  Archivo de metadatos por muestra (grupos, tratamientos, etc.) usado para las comparaciones de significancia más adelante. El archivo `metadata2.txt` aún no existe en este pipeline y debe crearse antes de correr el script.

- **`--output-dir core-metrics-results`**
  Carpeta donde `core-metrics-phylogenetic` deja todos sus artefactos de salida (vectores de diversidad alfa, matrices de distancia beta, PCoA, etc.) de una sola vez.

- **`qiime diversity alpha-group-significance`**
  Prueba si las diferencias en una métrica de diversidad alfa (Faith's PD, o número de features observadas) entre los grupos definidos en los metadatos son estadísticamente significativas, y genera un `.qzv` con los resultados y gráficos de caja por grupo.

---

## 18. Convertir ASVs a OTUs — pendiente

> **Estado: no ejecutado.** Queda como script de referencia para cuando el proyecto lo requiera.

En el directorio de scripts, crea el script de clustering:

```bash
nano 14_cluster_otus_v1.sh
```

Contenido del script (plantilla):

```bash
#!/bin/bash
#SBATCH -p cu
#SBATCH --cpus-per-task=4
#SBATCH -N 1
#SBATCH -D /Tayra-Share/home/Universidades/UCaldas/Users/eariasa/Rhogeessa
#SBATCH -J 14_cluster_otus
#SBATCH --mail-user=edwin.2051826349@ucaldas.edu.co
#SBATCH --mail-type=END,FAIL
#SBATCH --mem=40G
#SBATCH -o /Tayra-Share/home/Universidades/UCaldas/Users/eariasa/Rhogeessa/analysis_outputs/logs/14_cluster_otus_%j.out
#SBATCH -e /Tayra-Share/home/Universidades/UCaldas/Users/eariasa/Rhogeessa/analysis_outputs/logs/14_cluster_otus_%j.err
#SBATCH --time=02:00:00

source /Tayra-Share/miniconda/etc/profile.d/conda.sh
conda activate qiime2-amplicon-2024.10

OUTPUT_DIR=/Tayra-Share/home/Universidades/UCaldas/Users/eariasa/Rhogeessa/analysis_outputs

# 1. Agrupar ASVs en OTUs de novo al 97% de identidad
qiime vsearch cluster-features-de-novo \
  --i-sequences $OUTPUT_DIR/rep-seqs-bacteria-clean.qza \
  --i-table $OUTPUT_DIR/table-bacteria-clean.qza \
  --p-perc-identity 0.97 \
  --p-threads 4 \
  --output-dir $OUTPUT_DIR/cluster_denovo_dada

# 2. Resumen de la tabla de OTUs
qiime feature-table summarize \
  --i-table $OUTPUT_DIR/cluster_denovo_dada/clustered_table.qza \
  --o-visualization $OUTPUT_DIR/cluster_denovo_dada/clustered_table.qzv

# 3. Barplot taxonómico sobre la tabla de OTUs
qiime taxa barplot \
  --i-table $OUTPUT_DIR/cluster_denovo_dada/clustered_table.qza \
  --i-taxonomy $OUTPUT_DIR/taxonomy.qza \
  --m-metadata-file $OUTPUT_DIR/seqs_and_metadata/metadata.txt \
  --o-visualization $OUTPUT_DIR/cluster_denovo_dada/barplot-otus.qzv
```

### Explicación de los flags

- **`qiime vsearch cluster-features-de-novo`**
  Agrupa los ASVs en OTUs mediante clustering de novo (sin referencia externa): ASVs cuya similitud de secuencia supera el umbral dado se colapsan en un solo OTU. Es el paso inverso, en cierto sentido, a la resolución de base única que da DADA2.

- **`--i-sequences` / `--i-table`**
  Las secuencias representativas y la tabla de features ya filtradas a bacterias (salidas del paso 16) — los mismos insumos que usa el árbol filogenético en el paso 17.

- **`--p-perc-identity 0.97`**
  Umbral de identidad de secuencia (97%) para agrupar dos ASVs en el mismo OTU. Es el valor convencional en estudios de 16S, equivalente aproximadamente a diferenciación a nivel de especie.

- **`--p-threads 4`**
  Número de hilos para el clustering; debe ir acompañado de una reserva equivalente de CPUs en el script (`--cpus-per-task=4`).

- **`--output-dir cluster_denovo_dada`**
  Carpeta donde `vsearch` deja `clustered_table.qza` (tabla de OTUs) y `clustered_sequences.qza` (secuencia representativa de cada OTU), entre otros artefactos.

- **`qiime taxa barplot`** (aplicado aquí a la tabla de OTUs)
  Mismo uso que en el paso 15, pero reemplazando la tabla de ASVs por la tabla de OTUs recién generada, para comparar la composición taxonómica bajo ambos esquemas.

---

## 19. Diversidad para OTUs — pendiente

> **Estado: no ejecutado.** Repite los mismos pasos del punto 17, pero sobre la tabla y las secuencias ya agrupadas en OTUs (paso 18) en lugar de los ASVs originales.

En el directorio de scripts, crea el script de diversidad para OTUs:

```bash
nano 15_diversity_otus_v1.sh
```

Contenido del script (plantilla):

```bash
#!/bin/bash
#SBATCH -p cu
#SBATCH --cpus-per-task=8
#SBATCH -N 1
#SBATCH -D /Tayra-Share/home/Universidades/UCaldas/Users/eariasa/Rhogeessa
#SBATCH -J 15_diversity_otus
#SBATCH --mail-user=edwin.2051826349@ucaldas.edu.co
#SBATCH --mail-type=END,FAIL
#SBATCH --mem=40G
#SBATCH -o /Tayra-Share/home/Universidades/UCaldas/Users/eariasa/Rhogeessa/analysis_outputs/logs/15_diversity_otus_%j.out
#SBATCH -e /Tayra-Share/home/Universidades/UCaldas/Users/eariasa/Rhogeessa/analysis_outputs/logs/15_diversity_otus_%j.err
#SBATCH --time=03:00:00

source /Tayra-Share/miniconda/etc/profile.d/conda.sh
conda activate qiime2-amplicon-2024.10

OUTPUT_DIR=/Tayra-Share/home/Universidades/UCaldas/Users/eariasa/Rhogeessa/analysis_outputs
CLASSIFIER_DIR=$OUTPUT_DIR/classifier_SILVA
CLUSTER_DIR=$OUTPUT_DIR/cluster_denovo_dada

# 1. Tabular las secuencias representativas de los OTUs
qiime feature-table tabulate-seqs \
  --i-data $CLUSTER_DIR/clustered_sequences.qza \
  --o-visualization $CLUSTER_DIR/clustered_sequences.qzv

# 2. Clasificar taxonómicamente los OTUs
qiime feature-classifier classify-sklearn \
  --i-classifier $CLASSIFIER_DIR/SILVA_138.2_V3V4_341F_806R_CLASSIFIER.qza \
  --i-reads $CLUSTER_DIR/clustered_sequences.qza \
  --o-classification $CLUSTER_DIR/taxonomyOTUS.qza

# 3. Visualizar la taxonomía de los OTUs
qiime metadata tabulate \
  --m-input-file $CLUSTER_DIR/taxonomyOTUS.qza \
  --o-visualization $CLUSTER_DIR/taxonomyOTUS.qzv

# 4. Árbol filogenético de los OTUs
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences $CLUSTER_DIR/clustered_sequences.qza \
  --o-alignment $CLUSTER_DIR/aligned-rep-seqs.qza \
  --o-masked-alignment $CLUSTER_DIR/masked-aligned-rep-seqs.qza \
  --o-tree $CLUSTER_DIR/unrootedOTUS-tree.qza \
  --o-rooted-tree $CLUSTER_DIR/rootedOTUS-tree.qza

# 5. Métricas de diversidad sobre los OTUs
qiime diversity core-metrics-phylogenetic \
  --i-phylogeny $CLUSTER_DIR/rootedOTUS-tree.qza \
  --i-table $CLUSTER_DIR/clustered_table.qza \
  --p-sampling-depth 22000 \
  --m-metadata-file $OUTPUT_DIR/seqs_and_metadata/metadata2.txt \
  --output-dir $CLUSTER_DIR/OTUScore-metrics-results

# 6. Significancia de diversidad alfa — Faith's PD (OTUs)
qiime diversity alpha-group-significance \
  --i-alpha-diversity $CLUSTER_DIR/OTUScore-metrics-results/faith_pd_vector.qza \
  --m-metadata-file $OUTPUT_DIR/seqs_and_metadata/metadata2.txt \
  --o-visualization $CLUSTER_DIR/OTUScore-metrics-results/faith-pd-group-significance.qzv

# 7. Significancia de diversidad alfa — Features observadas (OTUs)
qiime diversity alpha-group-significance \
  --i-alpha-diversity $CLUSTER_DIR/OTUScore-metrics-results/observed_features_vector.qza \
  --m-metadata-file $OUTPUT_DIR/seqs_and_metadata/metadata2.txt \
  --o-visualization $CLUSTER_DIR/OTUScore-metrics-results/observed_features-group-significance.qzv
```

### Explicación de los flags

- **`qiime feature-table tabulate-seqs`**
  Genera una tabla visualizable con cada secuencia representativa de OTU y su longitud — útil para revisar rápidamente cuántos OTUs quedaron y su distribución de tamaños, igual que se hizo con `rep-seqs-bacteria-clean.qzv` en el paso 16.

- **`qiime feature-classifier classify-sklearn`** (aplicado aquí a OTUs)
  Mismo uso que en el paso 15, pero clasificando `clustered_sequences.qza` (una secuencia por OTU) en vez de `rep-seqs.qza` (una por ASV). El resultado es una tabla de taxonomía paralela, específica para el esquema de OTUs.

- **Pasos de árbol y diversidad (`align-to-tree-mafft-fasttree`, `core-metrics-phylogenetic`, `alpha-group-significance`)**
  Funcionalmente idénticos a los del paso 17, aplicados sobre `clustered_sequences.qza` y `clustered_table.qza` en lugar de las secuencias y tabla de ASVs. El objetivo es poder comparar directamente los resultados de diversidad bajo el esquema de ASVs (paso 17) contra el esquema de OTUs (este paso).

---

## Recursos utilizados

- [Training a classifier with EMP primers — UWEC HPC docs](https://docs.hpc.uwec.edu/classes/qiime2-testing/#training-a-classifier-with-emp-primers)
- [Tutorial Moving Pictures QIIME 2 en español, paso a paso — Microbioma Lab](https://microbioma-lab.com/tutorial-moving-pictures-qiime-2-en-espanol-paso-a-paso/)
- [RESCRIPt — bokulich-lab (GitHub)](https://github.com/bokulich-lab/RESCRIPt)
- [QIIME 2 Amplicon Docs](https://amplicon-docs.qiime2.org/en/stable/)
- [HCGS Metabarcoding Tutorials — Joseph7e (GitHub)](https://github.com/Joseph7e/HCGS_Metabarcoding_Tutorials)
- [Amplicon Data Analysis Using QIIME2 — ELIXIR Slovenia (PDF)](https://elixir.mf.uni-lj.si/pluginfile.php/4319/mod_resource/content/6/1.AmpliconDataAnalysisUsingQIIME2.pdf)