# Guía paso a paso de QIIME2 — Proyecto Rhogeessa

Guía de análisis de secuencias 16S (región V3-V4, cebadores 341F/806R) con QIIME2, documentada para el grupo de investigación **GEBIOME**. Cubre desde la importación de las secuencias hasta la clasificación taxonómica, corriendo en el servidor **Tayra** del Centro de Bioinformática y Biología Computacional de Colombia (BIOS) vía SLURM.

## Contenido

| Parte | Descripción | Archivo |
|---|---|---|
| 1 | Importación de secuencias a QIIME2 y visualización de calidad inicial | [`qiime2_parte1.md`](./qiime2_parte1.md) |
| 2 | Recorte de cebadores con cutadapt y denoising con DADA2 | [`qiime2_parte2.md`](./qiime2_parte2.md) |
| 3 | Entrenamiento del clasificador SILVA, asignación de taxonomía y scripts de referencia (árbol filogenético, diversidad, OTUs) | [`qiime2_parte3.md`](./qiime2_parte3.md) |

## Requisitos

- Usuario en el clúster Tayra (BIOS)
- Entorno conda `qiime2-amplicon-2024.10`
- Datos demultiplexados en formato FASTQ pareado (forward/reverse)

## Autor

E. Antonio Arias-Arias
[![DOI](https://zenodo.org/badge/1322161926.svg)](https://doi.org/10.5281/zenodo.22803027)
