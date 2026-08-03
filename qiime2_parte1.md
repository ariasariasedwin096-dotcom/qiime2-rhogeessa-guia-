# Guía paso a paso de QIIME2 — Proyecto Rhogeessa

*Desde la importación hasta la clasificación taxonómica*

> **Grupo de investigación:** GEBIOME
> **Servidor:** Tayra — Centro de Bioinformática y Biología Computacional de Colombia (BIOS)
> **Autor:** E. Antonio Arias-Arias
> **Contacto para acceder a las muestras usadas en esta guía:** Diego Ceballos — `Diego.ceballos05837@ucaldas.edu.co`

Esta guía está pensada para integrantes de GEBIOME con usuario en el servidor Tayra. Con mucho cariño, comparto este contenido para que sirva de referencia al equipo.

---

## Índice

- [Parte 1: Importación y visualización](#parte-1-importación-y-visualización)
  - [1. Subir las secuencias al servidor](#1-subir-las-secuencias-al-servidor)
  - [2. Organiza tus archivos](#2-organiza-tus-archivos)
  - [3. Crea un archivo manifest](#3-crea-un-archivo-manifest-necesario-para-qiime2)
  - [4. Crea un directorio de scripts ordenado](#4-crea-un-directorio-de-scripts-ordenado)
  - [5. Importa las secuencias a QIIME2](#5-importa-las-secuencias-a-qiime2)
  - [6. Visualiza la calidad](#6-visualiza-la-calidad)

---

## Parte 1: Importación y visualización

### 1. Subir las secuencias al servidor

Para esto utilizamos el comando `scp`, que sigue la siguiente estructura:

```bash
scp archivo1 archivo2 archivo3 usuario@servidor:/ruta/carpeta_destino/
```

**Ejemplo:**

```bash
scp /home/e-antonio-arias-arias/Rhogeessa/zr27477_19V3V4_R1.fastq.gz \
    /home/e-antonio-arias-arias/Rhogeessa/zr27477_19V3V4_R2.fastq.gz \
    eariasa@10.0.80.100:/Tayra-Share/home/Universidades/UCaldas/Users/eariasa/Rhogeessa/
```

**Estado de la ejecución:**

```
Rhogeessa/
├── zr27477_19V3V4_R1.fastq.gz  # Forward reads
└── zr27477_19V3V4_R2.fastq.gz  # Reverse reads
```

> Como no tenemos el archivo de barcodes, esto significa que nuestros datos **YA ESTÁN DEMULTIPLEXADOS**, por lo que iremos directo a QIIME2.

---

### 2. Organiza tus archivos

```bash
cd ~/Rhogeessa
mkdir -p qiime2_analysis/reads

# Mueve tus archivos
mv zr27477_19V3V4_R1.fastq.gz qiime2_analysis/reads/
mv zr27477_19V3V4_R2.fastq.gz qiime2_analysis/reads/

cd qiime2_analysis/reads
```

---

### 3. Crea un archivo manifest (necesario para QIIME2)

```bash
# Primero asegúrate de estar en reads
cd qiime2_analysis/reads

# Luego crea el manifest file
printf "sample-id\tforward-absolute-filepath\treverse-absolute-filepath\n" > manifest.tsv
printf "Rhogeessa\t/Tayra-Share/home/Universidades/UCaldas/Users/eariasa/Rhogeessa/qiime2_analysis/reads/zr27477_19V3V4_R1.fastq.gz\t/Tayra-Share/home/Universidades/UCaldas/Users/eariasa/Rhogeessa/qiime2_analysis/reads/zr27477_19V3V4_R2.fastq.gz\n" >> manifest.tsv
```

---

### 4. Crea un directorio de scripts ordenado

Desde tu carpeta principal de análisis (`Rhogeessa`):

```bash
# Crea carpetas para scripts y outputs
mkdir -p scripts
mkdir -p analysis_outputs

# Verifica la estructura
tree -L 2 ~/Rhogeessa
```

---

### 5. Importa las secuencias a QIIME2

En el directorio de scripts, crea el primer script de importación:

```bash
nano 01_import_sequences_v1.sh
```

Contenido del script:

```bash
#!/bin/bash
#SBATCH -p cu
#SBATCH -n 4
#SBATCH -N 1
#SBATCH -D /Tayra-Share/home/Universidades/UCaldas/Users/eariasa/Rhogeessa
#SBATCH -J qiime2_import
#SBATCH --mail-user=edwin.2051826349@ucaldas.edu.co
#SBATCH --mail-type=END,FAIL
#SBATCH --mem=40G
#SBATCH -o /Tayra-Share/home/Universidades/UCaldas/Users/eariasa/Rhogeessa/analysis_outputs/logs/qiime2_import_%j.out
#SBATCH -e /Tayra-Share/home/Universidades/UCaldas/Users/eariasa/Rhogeessa/analysis_outputs/logs/qiime2_import_%j.err
#SBATCH --time=02:00:00

source /Tayra-Share/miniconda/etc/profile.d/conda.sh
conda activate qiime2-amplicon-2024.10

DATA_DIR=/Tayra-Share/home/Universidades/UCaldas/Users/eariasa/Rhogeessa/qiime2_analysis/reads
OUTPUT_DIR=/Tayra-Share/home/Universidades/UCaldas/Users/eariasa/Rhogeessa/analysis_outputs

mkdir -p $OUTPUT_DIR/logs

qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path $DATA_DIR/manifest.tsv \
  --input-format PairedEndFastqManifestPhred33V2 \
  --output-path $OUTPUT_DIR/demux-paired-end.qza
```

#### Explicación de los flags

- **`qiime tools import`**
  Subcomando genérico de QIIME2 para traer datos externos (no generados por QIIME2) al formato de artefacto `.qza`.

- **`--type 'SampleData[PairedEndSequencesWithQuality]'`**
  Le dice a QIIME2 qué tipo semántico de dato es:
  - `SampleData[...]` = datos organizados por muestra
  - `PairedEndSequencesWithQuality` = lecturas pareadas (forward + reverse) con información de calidad (los scores Phred del FASTQ)

- **`--input-path $DATA_DIR/manifest.tsv`**
  La ruta al archivo que QIIME2 va a leer. Ojo: aquí no apunta directamente a los FASTQ, sino al archivo manifiesto (`manifest.tsv`), que es una tabla que mapea cada `sample-id` a las rutas de sus archivos forward/reverse. QIIME2 lee ese manifiesto y va a buscar los FASTQ a las rutas que ahí se indican.

- **`--input-format PairedEndFastqManifestPhred33V2`**
  Especifica el formato exacto en que está estructurado el archivo de `--input-path`, para que QIIME2 sepa cómo parsearlo.

- **`--output-path $OUTPUT_DIR/demux-paired-end.qza`**
  Donde se guarda el resultado: un `.qza`, que es un zip + metadatos.

---

### 6. Visualiza la calidad

En el directorio de scripts, crea el segundo script de visualización:

```bash
nano 02_quality_control_v1.sh
```

Contenido del script:

```bash
#!/bin/bash
#SBATCH -p cu
#SBATCH -n 4
#SBATCH -N 1
#SBATCH -D /Tayra-Share/home/Universidades/UCaldas/Users/eariasa/Rhogeessa
#SBATCH -J qiime2_qc
#SBATCH --mail-user=edwin.2051826349@ucaldas.edu.co
#SBATCH --mail-type=END,FAIL
#SBATCH --mem=40G
#SBATCH -o /Tayra-Share/home/Universidades/UCaldas/Users/eariasa/Rhogeessa/analysis_outputs/logs/qiime2_qc_%j.out
#SBATCH -e /Tayra-Share/home/Universidades/UCaldas/Users/eariasa/Rhogeessa/analysis_outputs/logs/qiime2_qc_%j.err
#SBATCH --time=02:00:00

source /Tayra-Share/miniconda/etc/profile.d/conda.sh
conda activate qiime2-amplicon-2024.10

OUTPUT_DIR=/Tayra-Share/home/Universidades/UCaldas/Users/eariasa/Rhogeessa/analysis_outputs

qiime demux summarize \
  --i-data $OUTPUT_DIR/demux-paired-end.qza \
  --o-visualization $OUTPUT_DIR/demux-paired-end.qzv
```

Puedes visualizar el archivo `demux-paired-end.qzv` cargándolo en la herramienta de visualización de QIIME2: [view.qiime2.org](https://view.qiime2.org/)

---

*Continúa en la Parte 2: control de calidad, recorte de cebadores con cutadapt y denoising con DADA2.*