#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
========================================================================================
    RNA-seq Analysis Pipeline
========================================================================================
    Quality Control -> Trimming -> Pseudoalignment (Kallisto) -> Differential Expression
----------------------------------------------------------------------------------------
*/

// Print pipeline info
log.info """\
    RNA-SEQ ANALYSIS PIPELINE
    ===================================
    input       : ${params.input}
    outdir      : ${params.outdir}
    transcriptome : ${params.transcriptome}
    gtf         : ${params.gtf}
    """
    .stripIndent()

/*
========================================================================================
    PROCESSES
========================================================================================
*/

// Process 1: Quality Control with FastQC
process FASTQC {
    tag "$sample_id"
    container 'docker://biocontainers/fastqc:v0.11.9_cv8'

    input:
    tuple val(sample_id), path(reads)

    output:
    path "*.{html,zip}"

    script:
    """
    fastqc -t ${task.cpus} ${reads}
    """
}

// Process 2: Read Trimming with fastp
process TRIM_READS {
    tag "$sample_id"
    container 'biocontainers/fastp:v0.20.1_cv1'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_trimmed_R{1,2}.fastq.gz"), emit: trimmed_reads
    path "${sample_id}_fastp.{html,json}", emit: report

    script:
    """
    fastp \
        -i ${reads[0]} \
        -I ${reads[1]} \
        -o ${sample_id}_trimmed_R1.fastq.gz \
        -O ${sample_id}_trimmed_R2.fastq.gz \
        --detect_adapter_for_pe \
        --qualified_quality_phred ${params.min_quality} \
        --length_required ${params.min_length} \
        --thread ${task.cpus} \
        --compression 6 \
        --html ${sample_id}_fastp.html \
        --json ${sample_id}_fastp.json
    """
}

// Process 3a: Create tx2gene mapping
process CREATE_TX2GENE {
    container 'docker://rnakato/rumball'
    publishDir "${params.outdir}/reference", mode: 'copy'

    input:
    path transcriptome

    output:
    path "tx2gene.tsv", emit: tx2gene

    script:
    """
    Rscript ${projectDir}/create_tx2gene.R ${transcriptome} tx2gene.tsv
    """
}

// Process 3b: Kallisto Index (run once)
process KALLISTO_INDEX {
    container 'docker://zlskidmore/kallisto'

    input:
    path transcriptome

    output:
    path "transcriptome.idx", emit: index

    script:
    """
    kallisto index \
        -i transcriptome.idx \
        ${transcriptome}
    """
}

// Process 4: Kallisto Quantification
process KALLISTO_QUANT {
    tag "$sample_id"
    container 'docker://zlskidmore/kallisto'

    input:
    tuple val(sample_id), path(reads)
    path kallisto_index

    output:
    tuple val(sample_id), path("${sample_id}"), emit: quant_dir
    path "${sample_id}/abundance.tsv", emit: abundance
    path "${sample_id}/run_info.json", emit: run_info

    script:
    """
    kallisto quant \
        -i ${kallisto_index} \
        -o ${sample_id} \
        -t ${task.cpus} \
        -b ${params.bootstrap_samples} \
        ${reads[0]} ${reads[1]}
    """
}

// Process 7: MultiQC Report
process MULTIQC {
    publishDir "${params.outdir}/multiqc", mode: 'copy', overwrite: true
    container 'docker://ewels/multiqc:v1.19'

    input:
    path 'qc_data/*/*', stageAs: 'qc_data/?/*'

    output:
    path "multiqc_report.html"
    path "multiqc_data"

    script:
    """
    multiqc qc_data/
    """
}

// Process 5: Merge Kallisto Results
process MERGE_KALLISTO {
    publishDir "${params.outdir}/final_counts", mode: 'copy', overwrite: true
    container 'docker://rnakato/rumball'

    input:
    path quant_dirs
    path tx2gene

    output:
    path "kallisto_count_matrix.txt", emit: count_matrix
    path "kallisto_tpm_matrix.txt", emit: tpm_matrix

    script:
    """
    Rscript ${projectDir}/merge_kallisto.R ${tx2gene}
    """
}

// Process 6: Differential Expression Analysis
process DIFFERENTIAL_EXPRESSION {
    publishDir "${params.outdir}/differential_expression", mode: 'copy', overwrite: true
    container 'docker://rnakato/rumball'

    input:
    path count_matrix
    path metadata

    output:
    path "differential_expression_results.csv", emit: de_results
    path "normalized_counts.csv", emit: normalized_counts
    path "volcano_plot.png", emit: volcano_plot
    path "upregulated_genes.csv", emit: up_genes, optional: true
    path "downregulated_genes.csv", emit: down_genes, optional: true

    script:
    """
    Rscript ${projectDir}/differential_expression.R \
        ${count_matrix} \
        ${metadata} \
        .
    """
}

/*
========================================================================================
    WORKFLOW
========================================================================================
*/

workflow {
    // Read input samplesheet
    Channel
        .fromPath(params.input)
        .splitCsv(header:true)
        .map { row -> 
            def sample_id = row.Sample
            def r1 = file(row.R1)
            def r2 = file(row.R2)
            return tuple(sample_id, [r1, r2])
        }
        .set { samples_ch }

    // Main processing chain: Trim -> Kallisto -> Merge -> DE
    // Step 1: Trim reads
    TRIM_READS(samples_ch)

    // Step 2a: Create tx2gene mapping (if not provided)
    if (params.tx2gene) {
        tx2gene_ch = Channel.fromPath(params.tx2gene)
    } else {
        transcriptome_for_tx2gene = Channel.fromPath(params.transcriptome)
        CREATE_TX2GENE(transcriptome_for_tx2gene)
        tx2gene_ch = CREATE_TX2GENE.out.tx2gene
    }

    // Step 2b: Build Kallisto index (if not provided)
    if (params.kallisto_index) {
        kallisto_index_ch = Channel.fromPath(params.kallisto_index)
    } else {
        transcriptome = Channel.fromPath(params.transcriptome)
        KALLISTO_INDEX(transcriptome)
        kallisto_index_ch = KALLISTO_INDEX.out.index
    }

    // Step 3: Kallisto quantification on ALL trimmed samples
    // The index needs to be combined with each sample
    KALLISTO_QUANT(
        TRIM_READS.out.trimmed_reads,
        kallisto_index_ch.first()
    )

    // Step 4: Wait for ALL Kallisto quantifications, then merge
    // Extract just the directories from the tuple and collect all of them
    quant_dirs_only = KALLISTO_QUANT.out.quant_dir
        .map { sample_id, dir -> dir }
        .toList()
    
    MERGE_KALLISTO(quant_dirs_only, tx2gene_ch)

    // Step 5: Differential Expression and Gene Ontology Analysis
    metadata_ch = Channel.fromPath(params.input)
    DIFFERENTIAL_EXPRESSION(MERGE_KALLISTO.out.count_matrix, metadata_ch)

    // Step 6: Quality Control on raw reads (runs in parallel for all samples)
    FASTQC(samples_ch)

    // Step 7: Generate MultiQC report (waits for ALL QC processes to complete)
    // Collect all QC outputs and flatten into a single list before running MultiQC once
    all_qc_files = FASTQC.out.collect()
        .concat(TRIM_READS.out.report.collect())
        .concat(KALLISTO_QUANT.out.run_info.collect())
        .flatten()
        .collect()
    
    MULTIQC(all_qc_files)
}

workflow.onComplete {
    log.info """\
        Pipeline completed!
        ------------------
        Status    : ${workflow.success ? 'SUCCESS' : 'FAILED'}
        Output    : ${params.outdir}
        """
        .stripIndent()
}
