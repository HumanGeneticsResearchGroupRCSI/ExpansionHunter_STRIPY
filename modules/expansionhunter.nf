/*
    EXPANSIONHUNTER MODULE
    Repeat expansion calling with ExpansionHunter v5.0.0
    Single-threaded execution (1 CPU)

    Input BAM: mini-BAM extracted by SAMTOOLS_EXTRACT_REGIONS
               (only repeat-locus-covering reads, ~50-100MB)
    Mode:      seeking - required to produce realigned BAM for REViewer
*/

process EXPANSIONHUNTER {
    tag   "${meta.id}"
    label 'process_medium'

    container "${params.containers.expansionhunter}"

    // Publish VCF and JSON - the final deliverables from EH
    // The unsorted realigned BAM is NOT published here:
    //   - it cannot be indexed (unsorted) → no scientific use
    //   - it is passed via channel to SAMTOOLS_INDEX_REALIGNED
    //   - the sorted+indexed version is published by that process
    publishDir "${params.outdir}/expansionhunter/${meta.id}",
               mode: params.publish_dir_mode,
               pattern: "*.vcf"

    publishDir "${params.outdir}/expansionhunter/${meta.id}",
               mode: params.publish_dir_mode,
               pattern: "*.json"

    publishDir "${params.outdir}/expansionhunter/${meta.id}",
               mode: params.publish_dir_mode,
               pattern: "*_realigned.bam"

    publishDir "${params.outdir}/expansionhunter/${meta.id}",
               mode: params.publish_dir_mode,
               pattern: "versions.yml"

    input:
    tuple val(meta), path(bam), path(bai)
    path  reference
    path  reference_index
    path  variant_catalog

    output:
    tuple val(meta), path("${meta.id}.vcf"),               emit: vcf
    tuple val(meta), path("${meta.id}.json"),              emit: json
    tuple val(meta), path("${meta.id}_realigned.bam"),     emit: realigned_bam, optional: true
    path  "versions.yml",                                  emit: versions

    script:
    def sex = (meta.sex == '1') ? 'male' : 'female'
    // NOTE: --region-extension-length is an INTEGER (bp), not a file path.
    // Region targeting is handled upstream by SAMTOOLS_EXTRACT_REGIONS.
    // params.regions_bed is for samtools, not passed to EH.
    """
    ExpansionHunter \\
        --reads             ${bam} \\
        --reference         ${reference} \\
        --variant-catalog   ${variant_catalog} \\
        --output-prefix     ${meta.id} \\
        --sex               ${sex} \\
        --analysis-mode     ${params.analysis_mode}

    # NOTE: realigned BAM indexing is handled downstream by SAMTOOLS_INDEX_REALIGNED
    # using the samtools container - expansionhunter.sif does not include samtools

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ExpansionHunter: \$(ExpansionHunter --version 2>&1 | grep -oP 'v[0-9]+\\.[0-9]+\\.[0-9]+' || echo 'v5.0.0')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.vcf
    touch ${meta.id}.json
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ExpansionHunter: v5.0.0
        samtools: 1.17.0
    END_VERSIONS
    """
}
