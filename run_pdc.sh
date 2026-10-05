#!/bin/bash
cd /home/thor/Downloads/pseudomonas-aeruginosa
DATASETS=/home/thor/miniconda3/envs/ncbi_datasets/bin/datasets
LOG=ncbi_genomes_pdc_only/dehydrated_download.log
ZIP=ncbi_genomes_pdc_only/all_835_dehydrated.zip
DIR=ncbi_genomes_pdc_only/all
mkdir -p ncbi_genomes_pdc_only
echo "=== Started $(date) ===" >> "$LOG"
rm -f "$ZIP"
"$DATASETS" download genome accession \
  --inputfile P_aeruginosa_PDC_only/PDC_assembly_accessions.txt \
  --assembly-source genbank \
  --include genome \
  --dehydrated \
  --filename "$ZIP" >> "$LOG" 2>&1
echo "download exit: $? $(date)" >> "$LOG"
ls -lh "$ZIP" >> "$LOG"
rm -rf "$DIR"
mkdir -p "$DIR"
unzip -q "$ZIP" -d "$DIR" >> "$LOG" 2>&1
echo "unzip exit: $? $(date)" >> "$LOG"
"$DATASETS" rehydrate --directory "$DIR" --max-workers 5 >> "$LOG" 2>&1
echo "rehydrate exit: $? $(date)" >> "$LOG"
echo "fna count: $(find "$DIR" -name '*.fna' | wc -l) / 835" >> "$LOG"
echo "=== Done $(date) ===" >> "$LOG"
