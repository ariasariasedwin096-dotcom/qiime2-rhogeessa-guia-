# Guía paso a paso de QIIME2 — Parte 2

*Control de calidad, recorte de cebadores con cutadapt y denoising con DADA2*

---

## Índice

- [7. Recorta los primers con cutadapt](#7-recorta-los-primers-con-cutadapt)
- [8. Verifica la calidad tras el recorte](#8-verifica-la-calidad-tras-el-recorte)
- [9. Denoising con DADA2](#9-denoising-con-dada2)
- [10. Genera la tabla de estadísticas de denoising](#10-genera-la-tabla-de-estadísticas-de-denoising)

---

## 7. Recorta los primers con cutadapt

En el directorio de scripts, crea el script de recorte de cebadores:

```bash
nano 03_cutadapt_trim_primers.sh
```

Contenido del script:

```bash
#!/bin/bash
#SBATCH -p cu
#SBATCH -n 4
#SBATCH -N 1
#SBATCH -D /Tayra-Share/home/Universidades/UCaldas/Users/eariasa/Rhogeessa
#SBATCH -J qiime2_cutadapt
#SBATCH --mail-user=edwin.2051826349@ucaldas.edu.co
#SBATCH --mail-type=END,FAIL
#SBATCH --mem=40G
#SBATCH -o /Tayra-Share/home/Universidades/UCaldas/Users/eariasa/Rhogeessa/analysis_outputs/logs/qiime2_cutadapt_%j.out
#SBATCH -e /Tayra-Share/home/Universidades/UCaldas/Users/eariasa/Rhogeessa/analysis_outputs/logs/qiime2_cutadapt_%j.err
#SBATCH --time=02:00:00

source /Tayra-Share/miniconda/etc/profile.d/conda.sh
conda activate qiime2-amplicon-2024.10

OUTPUT_DIR=/Tayra-Share/home/Universidades/UCaldas/Users/eariasa/Rhogeessa/analysis_outputs

qiime cutadapt trim-paired \
  --i-demultiplexed-sequences $OUTPUT_DIR/demux-paired-end.qza \
  --p-front-f CCTAYGGGDBGCWGCAG \
  --p-front-r GAMTACNVGGGTHTCTAATCC \
  --o-trimmed-sequences $OUTPUT_DIR/demux-trimmed.qza
```

#### Explicación de los flags

- **`qiime cutadapt trim-paired`**
  Subcomando del plugin q2-cutadapt para recortar secuencias (típicamente adaptadores) de lecturas pareadas, buscando la secuencia real del primer en cada read en vez de asumir un recorte ciego de N bases fijas.

- **`--i-demultiplexed-sequences $OUTPUT_DIR/demux-paired-end.qza`**
  El artefacto de entrada: nuestras lecturas ya importadas y organizadas por muestra (el resultado del paso 5), todavía con los primers puestos.

- **`--p-front-f CCTAYGGGDBGCWGCAG`**
  La secuencia del primer forward (341F) a buscar y recortar del extremo 5' de cada lectura forward.

- **`--p-front-r GAMTACNVGGGTHTCTAATCC`**
  Lo mismo, pero para el primer reverse (806R), aplicado al extremo 5' de cada lectura reverse.

- **`--o-trimmed-sequences $OUTPUT_DIR/demux-trimmed.qza`**
  El artefacto de salida: lecturas pareadas, ya sin los primers, listas para el control de calidad y el denoising.

---

## 8. Verifica la calidad tras el recorte

Ya que cutadapt cambió la longitud y el contenido de las lecturas, necesitamos un nuevo resumen de calidad — el gráfico anterior (paso 6) ya no sirve para elegir los valores de solapamiento, porque incluía los primers.

En el directorio de scripts, crea el script de control de calidad:

```bash
nano 04_quality_control_trimmed_v1.sh
```

Contenido del script:

```bash
#!/bin/bash
#SBATCH -p cu
#SBATCH -n 4
#SBATCH -N 1
#SBATCH -D /Tayra-Share/home/Universidades/UCaldas/Users/eariasa/Rhogeessa
#SBATCH -J qiime2_qc_trimmed
#SBATCH --mail-user=edwin.2051826349@ucaldas.edu.co
#SBATCH --mail-type=END,FAIL
#SBATCH --mem=40G
#SBATCH -o /Tayra-Share/home/Universidades/UCaldas/Users/eariasa/Rhogeessa/analysis_outputs/logs/qiime2_qc_trimmed_%j.out
#SBATCH -e /Tayra-Share/home/Universidades/UCaldas/Users/eariasa/Rhogeessa/analysis_outputs/logs/qiime2_qc_trimmed_%j.err
#SBATCH --time=02:00:00

source /Tayra-Share/miniconda/etc/profile.d/conda.sh
conda activate qiime2-amplicon-2024.10

OUTPUT_DIR=/Tayra-Share/home/Universidades/UCaldas/Users/eariasa/Rhogeessa/analysis_outputs

qiime demux summarize \
  --i-data $OUTPUT_DIR/demux-trimmed.qza \
  --o-visualization $OUTPUT_DIR/demux-trimmed.qzv
```

> `qiime demux summarize` es la misma acción usada en el paso 6, aplicada esta vez sobre `demux-trimmed.qza`.

---

## 9. Denoising con DADA2

Con la calidad post-recorte revisada, definimos los puntos de solapamiento y corremos el denoising.

**Cálculo del overlap:**

- Forward sin primer: 284 pb → `trunc-len-f 267`
- Reverse sin primer: 214 pb → `trunc-len-r 193`
- Overlap = (267 + 193) − 430 = 30 pb

> **Nota sobre el valor de 430 pb:** es una estimación aproximada del largo del amplicón V3-V4 (341F/806R) tomada de la literatura, que reporta un rango de ~440-465 pb incluyendo primers. No es un dato medido directamente de nuestros datos — conviene confirmarlo empíricamente revisando la tasa de fusión (% merged) en `denoising-stats.qza` una vez corrido este script, o la distribución de longitudes de `rep-seqs.qza` con `qiime feature-table tabulate-seqs`.

En el directorio de scripts, crea el script de denoising:

```bash
nano 05_dada2_denoise.sh
```

Contenido del script:

```bash
#!/bin/bash
#SBATCH -p cu
#SBATCH -n 4
#SBATCH -N 1
#SBATCH -D /Tayra-Share/home/Universidades/UCaldas/Users/eariasa/Rhogeessa
#SBATCH -J qiime2_dada2
#SBATCH --mail-user=edwin.2051826349@ucaldas.edu.co
#SBATCH --mail-type=END,FAIL
#SBATCH --mem=40G
#SBATCH -o /Tayra-Share/home/Universidades/UCaldas/Users/eariasa/Rhogeessa/analysis_outputs/logs/qiime2_dada2_%j.out
#SBATCH -e /Tayra-Share/home/Universidades/UCaldas/Users/eariasa/Rhogeessa/analysis_outputs/logs/qiime2_dada2_%j.err
#SBATCH --time=04:00:00

source /Tayra-Share/miniconda/etc/profile.d/conda.sh
conda activate qiime2-amplicon-2024.10

OUTPUT_DIR=/Tayra-Share/home/Universidades/UCaldas/Users/eariasa/Rhogeessa/analysis_outputs

qiime dada2 denoise-paired \
  --i-demultiplexed-seqs $OUTPUT_DIR/demux-trimmed.qza \
  --p-trunc-len-f 267 \
  --p-trunc-len-r 193 \
  --o-table $OUTPUT_DIR/table.qza \
  --o-representative-sequences $OUTPUT_DIR/rep-seqs.qza \
  --o-denoising-stats $OUTPUT_DIR/denoising-stats.qza
```

#### Explicación de los flags

- **`qiime dada2 denoise-paired`**
  Acción del plugin q2-dada2 que corrige errores de secuenciación (denoising), fusiona los pares forward/reverse y agrupa las lecturas en ASVs (variantes de secuencia de amplicón) — el equivalente moderno a OTUs, pero con resolución de una sola base.

- **`--i-demultiplexed-seqs $OUTPUT_DIR/demux-trimmed.qza`**
  El artefacto de entrada: nuestras lecturas ya sin primers (salida del paso 7).

- **`--p-trunc-len-f 267` / `--p-trunc-len-r 193`**
  Posición en la que se truncan las lecturas forward y reverse respectivamente, elegida a partir del `.qzv` del paso 8. Todo lo que esté después de esa posición se descarta; lecturas más cortas que este valor se descartan por completo. Como en este script no hay `trim-left`, estos valores representan directamente el largo final de la lectura.

- **`--o-table $OUTPUT_DIR/table.qza`**
  La tabla de features: cuántas veces aparece cada ASV en cada muestra. Es el insumo principal para los análisis de diversidad.

- **`--o-representative-sequences $OUTPUT_DIR/rep-seqs.qza`**
  Las secuencias representativas de cada ASV (una secuencia por variante detectada) — se usan después para asignación taxonómica y para, por ejemplo, construir un árbol filogenético.

- **`--o-denoising-stats $OUTPUT_DIR/denoising-stats.qza`**
  Estadísticas por muestra de cada etapa del proceso: cuántas lecturas entraron, cuántas pasaron el filtrado, cuántas se fusionaron (merge) exitosamente, y cuántas quedaron tras remover quimeras. Este es el archivo que debes revisar primero después de correr el script — una tasa de fusión (merged) baja confirma si el overlap de 30 pb calculado arriba fue suficiente en la práctica o si hay que volver a ajustar los valores de truncamiento.

---

## 10. Genera la tabla de estadísticas de denoising

Ejecuta el siguiente script para obtener una tabla `.qzv` a partir de la salida del denoising. Esto te permitirá ver cuántas lecturas se conservaron durante el paso de filtrado de DADA2.

En el directorio de scripts, crea el script de estadísticas:

```bash
nano 06_denoising_stats_v1.sh
```

Contenido del script:

```bash
#!/bin/bash
#SBATCH -p cu
#SBATCH -n 4
#SBATCH -N 1
#SBATCH -D /Tayra-Share/home/Universidades/UCaldas/Users/eariasa/Rhogeessa
#SBATCH -J qiime2_stats
#SBATCH --mail-user=edwin.2051826349@ucaldas.edu.co
#SBATCH --mail-type=END,FAIL
#SBATCH --mem=40G
#SBATCH -o /Tayra-Share/home/Universidades/UCaldas/Users/eariasa/Rhogeessa/analysis_outputs/logs/qiime2_stats_%j.out
#SBATCH -e /Tayra-Share/home/Universidades/UCaldas/Users/eariasa/Rhogeessa/analysis_outputs/logs/qiime2_stats_%j.err
#SBATCH --time=00:30:00

source /Tayra-Share/miniconda/etc/profile.d/conda.sh
conda activate qiime2-amplicon-2024.10

OUTPUT_DIR=/Tayra-Share/home/Universidades/UCaldas/Users/eariasa/Rhogeessa/analysis_outputs

qiime metadata tabulate \
  --m-input-file $OUTPUT_DIR/denoising-stats.qza \
  --o-visualization $OUTPUT_DIR/denoising-stats.qzv
```

Puedes visualizar el archivo `denoising-stats.qzv` cargándolo en la herramienta de visualización de QIIME2: [view.qiime2.org](https://view.qiime2.org/)

---

*Continúa en la Parte 3: entrenamiento de un clasificador de taxonomía con SILVA.*