/*
========================================================================================
    SAMTOOLS MODULES
========================================================================================
    SAMTOOLS_INDEX            - index BAM files if .bai is missing
    SAMTOOLS_EXTRACT_REGIONS  - extract repeat-locus reads into mini-BAM
                                for hybrid seeking mode (WES optimisation)
----------------------------------------------------------------------------------------
*/

// ── SAMTOOLS_INDEX ────────────────────────────────────────────────────────────

process SAMTOOLS_INDEX {
    tag   "${meta.id}"
    label 'process_low'

    container "${params.containers.samtools}"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path(bam), path("${bam}.bai"), emit: bam_bai
    path  "versions.yml",                           emit: versions

    script:
    """
    samtools index -@ ${task.cpus} ${bam}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -n1 | sed 's/samtools //')
    END_VERSIONS
    """
}


// ── SAMTOOLS_EXTRACT_REGIONS ──────────────────────────────────────────────────

process SAMTOOLS_EXTRACT_REGIONS {
    tag   "${meta.id}"
    label 'process_low'

    container "${params.containers.samtools}"

    // Do not publishDir - mini-BAM is an intermediate file
    // Realigned BAM from EH is the useful output

    input:
    tuple val(meta), path(bam), path(bai)
    path  regions_bed   // BED file of padded repeat loci (generated from catalog)

    output:
    tuple val(meta), path("${meta.id}.regions.bam"), path("${meta.id}.regions.bam.bai"), emit: mini_bam
    path  "versions.yml",                                                                  emit: versions

    script:
    """
    # Extract reads overlapping repeat loci (with padding) into a mini-BAM
    # This reduces the BAM from ~5-15GB to ~50-100MB for WES data
    # EH seeking mode on this mini-BAM is fast - no index-jumping overhead
    samtools view \\
        -b \\
        -L ${regions_bed} \\
        -o ${meta.id}.regions.bam \\
        --write-index \\
        ${bam}

    # If samtools version does not support --write-index, index separately
    # The || true prevents failure if index was already written
    [ -f "${meta.id}.regions.bam.bai" ] || samtools index ${meta.id}.regions.bam

    # Sanity check - mini-BAM must have reads
    N_READS=\$(samtools view -c ${meta.id}.regions.bam)
    if [ "\${N_READS}" -eq "0" ]; then
        echo "[ERROR] Mini-BAM has 0 reads for ${meta.id}" >&2
        echo "[ERROR] Check that regions_bed coordinates match BAM chromosome naming" >&2
        exit 1
    fi
    echo "[INFO] Mini-BAM reads: \${N_READS} for ${meta.id}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -n1 | sed 's/samtools //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.regions.bam
    touch ${meta.id}.regions.bam.bai
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: 1.17.0
    END_VERSIONS
    """
}


// ── SAMTOOLS_INDEX_REALIGNED ──────────────────────────────────────────────────
// Index the realigned BAM produced by ExpansionHunter seeking mode.
// Runs in samtools container - expansionhunter.sif does not include samtools.
// Optional: only runs when realigned BAM exists (seeking mode).

process SAMTOOLS_INDEX_REALIGNED {
    tag   "${meta.id}"
    label 'process_low'

    container "${params.containers.samtools}"

    publishDir "${params.outdir}/expansionhunter/${meta.id}",
               mode: params.publish_dir_mode,
               pattern: "*.sorted.bam"

    publishDir "${params.outdir}/expansionhunter/${meta.id}",
               mode: params.publish_dir_mode,
               pattern: "*.sorted.bam.bai"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${meta.id}_realigned.sorted.bam"),
                     path("${meta.id}_realigned.sorted.bam.bai"), emit: bam_bai
    path  "versions.yml",                                         emit: versions

    script:
    // ExpansionHunter realigned BAM is NOT coordinate-sorted (graph-alignment order).
    // Must sort before indexing - samtools index will fail on unsorted input.
    // TMPDIR is set to node-local scratch by SLURM on RCSI HPC.
    """
    samtools sort \\
        -@ ${task.cpus} \\
        -T \${TMPDIR:-/tmp}/${meta.id}_sort \\
        -o ${meta.id}_realigned.sorted.bam \\
        ${bam}

    samtools index -@ ${task.cpus} ${meta.id}_realigned.sorted.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -n1 | sed 's/samtools //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_realigned.sorted.bam
    touch ${meta.id}_realigned.sorted.bam.bai
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: 1.23.0
    END_VERSIONS
    """
}
