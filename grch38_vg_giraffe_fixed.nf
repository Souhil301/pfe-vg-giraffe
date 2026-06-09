#!/usr/bin/env nextflow
nextflow.enable.dsl=2

SAMTOOLS  = "/home/etud2026/miniconda3/envs/pfe_benchmark/bin/samtools"
BCFTOOLS  = "/home/etud2026/miniconda3/envs/pfe_benchmark/bin/bcftools"
TABIX     = "/home/etud2026/miniconda3/envs/pfe_benchmark/bin/tabix"
VG_SIF    = "/home/etud2026/PFE_VC/tools/containers/vg_1.70.0.sif"
HAPPY_SIF = "/home/etud2026/PFE_VC/tools/containers/hap.py.sif"

// Index paths used as absolute (avoids Nextflow staging zipcodes issues)
VG_GBZ    = "/home/etud2026/PFE_VC/data_grch38/vg_indexes/grch38.giraffe.gbz"
VG_MIN    = "/home/etud2026/PFE_VC/data_grch38/vg_indexes/grch38.shortread.withzip.min"
VG_DIST   = "/home/etud2026/PFE_VC/data_grch38/vg_indexes/grch38.dist"
VG_PATHS  = "/home/etud2026/PFE_VC/data_grch38/vg_indexes/ref_paths_grch38.txt"

params.ref       = "/home/etud2026/PFE_VC/data_grch38/reference/GRCh38_no_alt.fa"
params.fastq_r1  = "/home/etud2026/PFE_VC/data_grch38/reads/HG002_HiSeq30x_subsampled_R1.fastq.gz"
params.fastq_r2  = "/home/etud2026/PFE_VC/data_grch38/reads/HG002_HiSeq30x_subsampled_R2.fastq.gz"
params.truth_vcf = "/home/etud2026/PFE_VC/data_grch38/truth/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
params.truth_bed = "/home/etud2026/PFE_VC/data_grch38/truth/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"
params.sample_id = "HG002"
params.tool_tag  = "vg_giraffe_t8"
params.outdir    = "/home/etud2026/PFE_VC/results/grch38_t8/vg_giraffe"
params.threads   = 8

process VG_ALIGN {
    cpus   { params.threads as int }
    memory '50 GB'
    time   '96h'
    tag    "vg_align_t${params.threads}"

    publishDir "${params.outdir}/mapping", mode: 'copy',
               pattern: '*.{bam,bai}'

    input:
    path r1; path r2

    output:
    path("${params.sample_id}_${params.tool_tag}.bam"),     emit: bam
    path("${params.sample_id}_${params.tool_tag}.bam.bai"), emit: bai

    script:
    def T     = params.threads as int
    def sortT = Math.max(1, T.intdiv(2))
    """
    echo "=== VG Giraffe v1.70.0 started: \$(date) ==="

    /usr/bin/time -v -o vg_align_timing.txt \
    apptainer exec --bind /home/etud2026:/home/etud2026 \
        ${VG_SIF} vg giraffe \
            -Z ${VG_GBZ} \
            -m ${VG_MIN} \
            -d ${VG_DIST} \
            --ref-paths ${VG_PATHS} \
            -f ${r1} -f ${r2} \
            --threads ${T} \
            --output-format BAM \
    2>vg_stderr.txt \
    | ${SAMTOOLS} sort \
        -@ ${sortT} -m 2G \
        -o tmp_sorted.bam -

    echo "=== vg stderr ===" && cat vg_stderr.txt

    # Rename GRCh38#0#chr1 -> chr1 in header (same count/order, just renamed)
    ${SAMTOOLS} view -H tmp_sorted.bam \
        | sed 's/SN:GRCh38#0#chr/SN:chr/g' \
        > new_header.sam

    ${SAMTOOLS} reheader new_header.sam tmp_sorted.bam \
        > ${params.sample_id}_${params.tool_tag}.bam

    rm -f tmp_sorted.bam new_header.sam

    ${SAMTOOLS} index -@ ${T} \
        ${params.sample_id}_${params.tool_tag}.bam

    echo "=== Header check (must show chr1) ==="
    ${SAMTOOLS} view -H ${params.sample_id}_${params.tool_tag}.bam \
        | grep "^@SQ" | head -3
    echo "=== Mapped reads on chr1 (must be > 0) ==="
    ${SAMTOOLS} view -c -F 4 \
        ${params.sample_id}_${params.tool_tag}.bam
    echo "=== Done: \$(date) ==="
    """
}

process CALL_CHR {
    cpus   2
    memory '8 GB'
    time   '4h'
    tag    "call_${chr}"

    input:
    tuple val(chr), path(ref), path(fai), path(bam), path(bai)

    output:
    tuple val(chr), path("${chr}.vcf.gz"), path("${chr}.vcf.gz.tbi"), emit: vcf

    script:
    """
    ${BCFTOOLS} mpileup \
        -f ${ref} -r ${chr} -q 20 -Q 20 \
        -Ou ${bam} \
    | ${BCFTOOLS} call -mv -Oz -o ${chr}.vcf.gz
    ${TABIX} -p vcf ${chr}.vcf.gz
    """
}

process MERGE_VCFS {
    cpus   4
    memory '8 GB'
    time   '2h'
    tag    "merge"

    publishDir "${params.outdir}", mode: 'copy',
               pattern: '*.{vcf.gz,vcf.gz.tbi}'

    input:
    path vcfs
    path tbis

    output:
    path("${params.sample_id}_${params.tool_tag}.vcf.gz"),     emit: vcf
    path("${params.sample_id}_${params.tool_tag}.vcf.gz.tbi"), emit: tbi

    script:
    """
    for c in \$(seq 1 22); do echo "chr\${c}.vcf.gz"; done > order.txt
    ${BCFTOOLS} concat -f order.txt -Oz \
        --threads ${task.cpus} \
        -o ${params.sample_id}_${params.tool_tag}.vcf.gz
    ${TABIX} -p vcf ${params.sample_id}_${params.tool_tag}.vcf.gz
    """
}

process HAPPY_EVAL {
    cpus   4
    memory '14 GB'
    time   '4h'
    tag    "happy"

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

workflow {
    ref_fa    = file(params.ref)
    ref_fai   = file("${params.ref}.fai")
    r1        = file(params.fastq_r1)
    r2        = file(params.fastq_r2)
    truth_vcf = file(params.truth_vcf)
    truth_bed = file(params.truth_bed)

    VG_ALIGN(r1, r2)

    // .combine() ensures CALL_CHR waits for VG_ALIGN to finish
    chr_ch = Channel.from(1..22)
        .combine(VG_ALIGN.out.bam)
        .combine(VG_ALIGN.out.bai)
        .map { c, bam, bai ->
            tuple("chr${c}", ref_fa, ref_fai, bam, bai)
        }

    called = CALL_CHR(chr_ch)
    vcfs   = called.vcf.map { it[1] }.collect()
    tbis   = called.vcf.map { it[2] }.collect()
    merged = MERGE_VCFS(vcfs, tbis)

    HAPPY_EVAL(ref_fa, ref_fai,
               merged.vcf, merged.tbi,
               truth_vcf, truth_bed)
}