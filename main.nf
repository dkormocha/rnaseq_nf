#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// ═══════════════════════════════════════════════════════════
//  PARAMETERS
// ═══════════════════════════════════════════════════════════
params.reads       = "data/*_{1,2}.fastq.gz"
params.genome      = "ref/genome.fa"
params.gtf         = "ref/annotation.gtf"
params.star_index  = "ref/star_index"
params.outdir      = "results"
params.threads     = 8

// ═══════════════════════════════════════════════════════════
//  PROCESS 1: FastQC on raw reads
// ═══════════════════════════════════════════════════════════
process FASTQC_RAW {
    tag "$sample_id raw"
    publishDir "${params.outdir}/fastqc/raw", mode: 'copy'
    conda 'bioconda::fastqc=0.11.9'

    input:
    tuple val(sample_id), path(reads)

    output:
    path "*.{html,zip}", emit: qc_files

    script:
    """
    fastqc --threads ${params.threads} --outdir . ${reads}
    """
}

// ═══════════════════════════════════════════════════════════
//  PROCESS 2: Trim adapters + low-quality bases
// ═══════════════════════════════════════════════════════════
process TRIM_READS {
    tag "$sample_id"
    publishDir "${params.outdir}/trimmed", mode: 'copy'
    conda 'bioconda::trimmomatic=0.39'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id),
          path("${sample_id}_R{1,2}_trimmed.fastq.gz"),
          emit: trimmed_reads
    path "trimmomatic.log", emit: log

    script:
    """
    trimmomatic PE -threads ${params.threads} \\
        ${reads[0]} ${reads[1]} \\
        ${sample_id}_R1_trimmed.fastq.gz /dev/null \\
        ${sample_id}_R2_trimmed.fastq.gz /dev/null \\
        ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 \\
        LEADING:3 TRAILING:3 \\
        SLIDINGWINDOW:4:15 MINLEN:36 \\
        2> trimmomatic.log
    """
}

// ═══════════════════════════════════════════════════════════
//  PROCESS 3: FastQC on trimmed reads
// ═══════════════════════════════════════════════════════════
process FASTQC_TRIMMED {
    tag "$sample_id trimmed"
    publishDir "${params.outdir}/fastqc/trimmed", mode: 'copy'
    conda 'bioconda::fastqc=0.11.9'

    input:
    tuple val(sample_id), path(reads)

    output:
    path "*.{html,zip}", emit: qc_files

    script:
    """
    fastqc --threads ${params.threads} --outdir . ${reads}
    """
}

// ═══════════════════════════════════════════════════════════
//  PROCESS 4: STAR alignment
// ═══════════════════════════════════════════════════════════
process STAR_ALIGN {
    tag "$sample_id"
    publishDir "${params.outdir}/star", mode: 'copy'
    conda 'bioconda::star=2.7.10a'

    input:
    tuple val(sample_id), path(reads)
    path star_index    // passed as a single value channel, shared across samples

    output:
    tuple val(sample_id),
          path("${sample_id}*.bam"),
          emit: bam
    path "${sample_id}Log.final.out", emit: log

    script:
    """
    STAR \\
        --runThreadN ${params.threads} \\
        --genomeDir ${star_index} \\
        --readFilesIn ${reads[0]} ${reads[1]} \\
        --readFilesCommand zcat \\
        --outSAMtype BAM SortedByCoordinate \\
        --outSAMattributes NH HI AS NM \\
        --outFileNamePrefix ${sample_id} \\
        --quantMode GeneCounts \\
        --outFilterMultimapNmax 20 \\
        --alignSJoverhangMin 8 \\
        --alignSJDBoverhangMin 1
    """
}

// ═══════════════════════════════════════════════════════════
//  PROCESS 5: Sort + index BAM
// ═══════════════════════════════════════════════════════════
process SAMTOOLS_SORT_INDEX {
    tag "$sample_id"
    publishDir "${params.outdir}/bam", mode: 'copy'
    conda 'bioconda::samtools=1.17'

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id),
          path("${sample_id}.sorted.bam"),
          path("${sample_id}.sorted.bam.bai"),
          emit: sorted_bam
    path "${sample_id}.flagstat.txt", emit: flagstat

    script:
    """
    samtools sort -@ ${params.threads} -o ${sample_id}.sorted.bam ${bam}
    samtools index ${sample_id}.sorted.bam
    samtools flagstat ${sample_id}.sorted.bam > ${sample_id}.flagstat.txt
    """
}

// ═══════════════════════════════════════════════════════════
//  PROCESS 6: featureCounts — gene-level counts
// ═══════════════════════════════════════════════════════════
process FEATURECOUNTS {
    tag "all samples"
    publishDir "${params.outdir}/counts", mode: 'copy'
    conda 'bioconda::subread=2.0.3'

    input:
    path bams          // collected list of ALL sorted BAMs
    path gtf

    output:
    path "counts.txt",         emit: counts
    path "counts.txt.summary", emit: summary

    script:
    """
    featureCounts \\
        -T ${params.threads} \\
        -t exon \\
        -g gene_id \\
        -s 2 \\               # stranded: 2 = reverse-stranded (common for TruSeq)
        -p \\                 # paired-end
        -a ${gtf} \\
        -o counts.txt \\
        ${bams.join(' ')}
    """
}

// ═══════════════════════════════════════════════════════════
//  PROCESS 7: MultiQC — aggregate report
// ═══════════════════════════════════════════════════════════
process MULTIQC {
    publishDir "${params.outdir}/multiqc", mode: 'copy'
    conda 'bioconda::multiqc=1.14'

    input:
    path qc_files     // all QC + log files, collected

    output:
    path "multiqc_report.html"
    path "multiqc_data"

    script:
    """
    multiqc . --outdir .
    """
}

// ═══════════════════════════════════════════════════════════
//  WORKFLOW — wire it all together
// ═══════════════════════════════════════════════════════════
workflow {

    // ── 1. Create input channel from paired FASTQ files ────
    reads_ch = Channel
        .fromFilePairs(params.reads, checkIfExists: true)

    // ── 2. QC raw reads ────────────────────────────────────
    FASTQC_RAW(reads_ch)

    // ── 3. Trim ────────────────────────────────────────────
    TRIM_READS(reads_ch)

    // ── 4. QC trimmed reads ────────────────────────────────
    FASTQC_TRIMMED(TRIM_READS.out.trimmed_reads)

    // ── 5. Align with STAR ─────────────────────────────────
    star_index_ch = Channel.value(file(params.star_index))
    STAR_ALIGN(TRIM_READS.out.trimmed_reads, star_index_ch)

    // ── 6. Sort and index BAM ──────────────────────────────
    SAMTOOLS_SORT_INDEX(STAR_ALIGN.out.bam)

    // ── 7. Count reads per gene ────────────────────────────
    // .collect() waits for ALL bam files before running featureCounts
    all_bams = SAMTOOLS_SORT_INDEX.out.sorted_bam
        .map { id, bam, bai -> bam }
        .collect()
    FEATURECOUNTS(all_bams, file(params.gtf))

    // ── 8. Collect ALL QC files and run MultiQC ────────────
    all_qc = FASTQC_RAW.out.qc_files
        .mix(FASTQC_TRIMMED.out.qc_files)
        .mix(TRIM_READS.out.log)
        .mix(STAR_ALIGN.out.log)
        .mix(FEATURECOUNTS.out.summary)
        .collect()
    MULTIQC(all_qc)
}