# PFE Benchmark — VG Giraffe Pangenome Variant Calling

Part of the thesis: *Benchmarking Variant Calling Paradigms on HG002 GRCh38 (30x)*
USTHB | CERIST, Algiers — 2025/2026

## What this pipeline does

Aligns reads to the **HPRC v1.1 pangenome** (47 human haplotypes + GRCh38 + CHM13)
using **vg Giraffe v1.70.0** and calls variants with two strategies:
Pipeline A (VG+BCFtools):
FASTQ → vg giraffe (--output-format BAM) → reheader → BCFtools call → VCF
Pipeline B (VG+vgcall, graph-native):
FASTQ → vg giraffe (--output-format gam) | vg pack → vg call → VCF

## Results (HG002 GRCh38, GIAB v4.2.1, chr1-22, PASS filter)

| Pipeline | SNP F1 | INDEL F1 | Time | Machine |
|----------|--------|----------|------|---------|
| VG+BCFtools t8 | 99.06% | **91.17%** | ~32h | ServerA (64 GB) |
| VG+vgcall t12 | 96.17% | **79.02%** | --- | ServerB (128 GB) |
| *BCFtools t8 (reference)* | *99.05%* | *94.12%* | --- | *ServerA* |

**Key finding**: Pangenome alignment paired with a linear caller (BCFtools)
**degrades** INDEL accuracy by 3% vs BWA+BCFtools. VG+vgcall produces 164k
false positive INDELs due to a protobuf overflow bug in vg v1.70.0 at
full-genome scale (issue [#4721](https://github.com/vgteam/vg/issues/4721)).

**Critical note on indexes**: HPRC v1.1 indexes distributed with the pangenome
were built with vg v1.46 and are **incompatible** with vg v1.70.0
(SnarlDistanceIndex format v2 → v3). Indexes must be rebuilt:

```bash
# Requires 128 GB RAM, ~2.5h
vg index -s -t 16 -j grch38.dist grch38.giraffe.gbz          # 30 min
vg minimizer -t 16 -d grch38.dist -z grch38.shortread.zipcodes \
    -o grch38.shortread.withzip.min grch38.giraffe.gbz         # 2h
```

The rebuilt indexes (vg v1.70.0 compatible) are available at: **[ZENODO DOI TBD]**

## Requirements

| Software | Version |
|----------|---------|
| Nextflow | 25.10.2 |
| vg | v1.70.0 "Zebedassi" |
| BCFtools | 1.23.1 |
| SAMtools | 1.23.1 |
| Apptainer | 1.2.5 |

```bash
# Pull vg container
apptainer pull vg_1.70.0.sif docker://quay.io/vgteam/vg:v1.70.0
```

## Input data

| File | Size | Source |
|------|------|--------|
| HG002 HiSeq 30x FASTQ | 83 GB | GIAB FTP |
| HPRC v1.1 GBZ | 3.1 GB | [HPRC S3](s3://human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-grch38/) |
| Minimizer index (rebuilt) | 36 GB | Zenodo (see above) |
| Distance index (rebuilt) | 8.1 GB | Zenodo (see above) |
| GIAB v4.2.1 truth | 1 GB | GIAB FTP |

## Usage — Pipeline A (VG + BCFtools)

```bash
nextflow -C pipelines/grch38_vg_t8.config \
    run pipelines/grch38_vg_giraffe_fixed.nf \
    -work-dir work/grch38_vg_t8 \
    2>&1 | tee run_vg_bcf.log
```

## Usage — Pipeline B (VG + vg call, requires ≥128 GB RAM)

```bash
# Requires the pack file from Pipeline A (or run VG_ALIGN_PACK process)
nextflow -C pipelines/grch38_vg_call_ursda.config \
    run pipelines/grch38_vg_call_ursda.nf \
    --pack /path/to/HG002_vg.pack \
    -work-dir work/grch38_vg_call \
    2>&1 | tee run_vg_call.log
```

## Known issues

| Issue | Description | Workaround |
|-------|-------------|------------|
| Index incompatibility | vg v1.46 indexes crash vg v1.70 | Rebuild indexes (see above) |
| CHM13 coordinate space | Without `--ref-paths`, vg outputs CHM13 coords | Always pass `--ref-paths ref_paths_grch38.txt` |
| BAM contig naming | vg outputs `GRCh38#0#chr1`, BCFtools needs `chr1` | `samtools reheader` step included |
| vg call protobuf overflow | >2 GB message on large chromosomes | Requires 128 GB RAM, sequential per-chr |
| hap.py segfault on MNPs | vg call produces 83k MNPs incompatible with hap.py | `bcftools view -v snps,indels` before hap.py |
