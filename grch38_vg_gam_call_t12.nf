#!/usr/bin/env nextflow
nextflow.enable.dsl=2

VG_SIF    = "/home/etud2026/PFE_VC/tools/containers/vg_1.70.0.sif"
BCFTOOLS  = "/home/etud2026/miniconda3/envs/pfe_benchmark/bin/bcftools"
TABIX     = "/home/etud2026/miniconda3/envs/pfe_benchmark/bin/tabix"
HAPPY_SIF = "/home/etud2026/PFE_VC/tools/containers/hap.py.sif"
GBZ       = "/home/etud2026/PFE_VC/data_grch38/vg_indexes/grch38.giraffe.gbz"
MIN_IDX   = "/home/etud2026/PFE_VC/data_grch38/vg_indexes/grch38.shortread.withzip.min"
DIST_IDX  = "/home/etud2026/PFE_VC/data_grch38/vg_indexes/grch38.dist"
ZIP_IDX   = "/home/etud2026/PFE_VC/data_grch38/vg_indexes/grch38.shortread.zipcodes"

params.fastq_r1  = "/home/etud2026/PFE_VC/data_grch38/reads/HG002_HiSeq30x_subsampled_R1.fastq.gz"
params.fastq_r2  = "/home/etud2026/PFE_VC/data_grch38/reads/HG002_HiSeq30x_subsampled_R2.fastq.gz"
params.ref       = "/home/etud2026/PFE_VC/data_grch38/reference/GRCh38_no_alt.fa"
params.truth_vcf = "/home/etud2026/PFE_VC/data_grch38/truth/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
params.truth_bed = "/home/etud2026/PFE_VC/data_grch38/truth/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"
params.outdir    = "/home/etud2026/PFE_VC/results/grch38_t12/vg_giraffe_vgcall"
params.tool_tag  = "vg_giraffe_vgcall_t12"
params.threads   = 12

// ── Step 1: Align (GAM) + Pack streamed in one step 
process VG_ALIGN_PACK {
    cpus   12
    memory '70 GB'
    time   '120h'
    tag    "vg_align_pack_t12"

    publishDir "${params.outdir}", mode: 'copy',
               pattern: 'vg_*.txt'

    input:
    path r1; path r2

    output:
    path("HG002_vg_t12.pack"),     emit: pack
    path("vg_align_timing.txt"),   emit: timing

    script:
    """
    echo "=== VG Giraffe → GAM → Pack started: \$(date) ==="

    /usr/bin/time -v -o vg_align_timing.txt \
    apptainer exec --bind /home/etud2026:/home/etud2026 \
        ${VG_SIF} vg giraffe \
            -Z ${GBZ} \
            -m ${MIN_IDX} \
            -d ${DIST_IDX} \
            -z ${ZIP_IDX} \
            -f ${r1} -f ${r2} \
            --threads ${params.threads} \
            --output-format gam \
    2>vg_giraffe_stderr.txt \
    | apptainer exec --bind /home/etud2026:/home/etud2026 \
        ${VG_SIF} vg pack \
            -x ${GBZ} \
            -g - \
            -Q 5 \
            -t ${params.threads} \
            -o HG002_vg_t12.pack \
    2>vg_pack_stderr.txt

    echo "=== Giraffe stderr ===" && cat vg_giraffe_stderr.txt
    echo "=== Pack stderr ===" && cat vg_pack_stderr.txt
    echo "=== Pack size: \$(du -sh HG002_vg_t12.pack) ==="
    echo "=== Align+Pack done: \$(date) ==="
    cat vg_align_timing.txt | grep -E "wall clock|Maximum|Percent"
    """
}

// ── Step 2: vg call per chromosome (22 parallel, max 3 concurrent) 
process VG_CALL_CHR {
    errorStrategy 'ignore'
    cpus   4
    memory '16 GB'
    time   '8h'
    tag    "vgcall_${chr}"

    input:
    tuple val(chr), path(pack)

    output:
    tuple val(chr), path("${chr}.vcf.gz"), path("${chr}.vcf.gz.tbi"), emit: vcf

    script:
    """
    echo "=== vg call ${chr} started: \$(date) ==="

    # Write to file first — avoids broken pipe if vg call produces empty/malformed output
    apptainer exec --bind /home/etud2026:/home/etud2026 \
        ${VG_SIF} vg call \
            ${GBZ} \
            -k ${pack} \
            -p GRCh38#0#${chr} \
            -t ${task.cpus} \
            --ploidy 2 \
            -s HG002 \
        > raw_${chr}.vcf \
        2>vg_call_${chr}.log
    echo "=== vg call exit: \$? lines: \$(wc -l < raw_${chr}.vcf) ==="
    cat vg_call_${chr}.log

    sed 's/GRCh38#0#chr/chr/g' raw_${chr}.vcf \
        | ${BCFTOOLS} sort -Oz -o ${chr}.vcf.gz
    rm -f raw_${chr}.vcf
    ${TABIX} -p vcf ${chr}.vcf.gz

    VARS=\$(${BCFTOOLS} stats ${chr}.vcf.gz | grep "^SN.*records" | awk '{print \$NF}')
    echo "=== ${chr}: \${VARS} variants at \$(date) ==="
    """
}

// ── Step 3: Merge all chromosomes 
process MERGE_VCFS {
    cpus 4; memory '8 GB'; time '2h'; tag "merge"

    publishDir "${params.outdir}", mode: 'copy',
               pattern: '*.{vcf.gz,vcf.gz.tbi}'

    input:
    path vcfs; path tbis

    output:
    path("HG002_${params.tool_tag}.vcf.gz"),     emit: vcf
    path("HG002_${params.tool_tag}.vcf.gz.tbi"), emit: tbi

    script:
    """
    for c in \$(seq 1 22); do echo "chr\${c}.vcf.gz"; done > order.txt
    ${BCFTOOLS} concat -f order.txt -Oz --threads ${task.cpus} \
        -o HG002_${params.tool_tag}.vcf.gz
    ${TABIX} -p vcf HG002_${params.tool_tag}.vcf.gz
    TOTAL=\$(${BCFTOOLS} stats HG002_${params.tool_tag}.vcf.gz | grep "^SN.*records" | awk '{print \$NF}')
    echo "=== Total variants: \${TOTAL} ==="
    """
}

// ── Step 4: hap.py evaluation 
process HAPPY_EVAL {
    cpus 4; memory '14 GB'; time '4h'; tag "happy"

    publishDir "${params.outdir}/happy", mode: 'copy'

    input:
    path ref; path fai
    path vcf; path tbi
    path truth_vcf; path truth_bed

    output:
    path("${params.tool_tag}_happy.*")

    script:
    """
    apptainer exec --bind /home/etud2026:/home/etud2026 \
        ${HAPPY_SIF} /usr/local/bin/hap.py \
            ${truth_vcf} ${vcf} \
            -r ${ref} -f ${truth_bed} \
            -o ${params.tool_tag}_happy \
            --threads ${task.cpus}
    """
}

// ── Workflow 
workflow {
    r1        = file(params.fastq_r1)
    r2        = file(params.fastq_r2)
    ref_fa    = file(params.ref)
    ref_fai   = file("${params.ref}.fai")
    truth_vcf = file(params.truth_vcf)
    truth_bed = file(params.truth_bed)

    // 1. Align + Pack
    VG_ALIGN_PACK(r1, r2)

    // 2. Call 22 chromosomes in parallel
    chr_ch = Channel.from(1..22)
        .combine(VG_ALIGN_PACK.out.pack)
        .map { c, pack -> tuple("chr${c}", pack) }

    called = VG_CALL_CHR(chr_ch)

    // 3. Merge
    vcfs   = called.vcf.map { it[1] }.collect()
    tbis   = called.vcf.map { it[2] }.collect()
    merged = MERGE_VCFS(vcfs, tbis)

    // 4. Evaluate
    HAPPY_EVAL(ref_fa, ref_fai,
               merged.vcf, merged.tbi,
               truth_vcf, truth_bed)
}
