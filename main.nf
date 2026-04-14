#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
========================================================================================
    EXPANSIONHUNTER REPEAT EXPANSION DETECTION PIPELINE v1.1
========================================================================================
    PED-driven repeat expansion calling with ExpansionHunter v5.0.0
    Hybrid seeking mode: samtools region extraction → EH seeking → REViewer-ready

    Author: Deepak Bharti
    Version: 1.1.0
----------------------------------------------------------------------------------------
*/

log.info """
╔══════════════════════════════════════════════════════════════╗
║      EXPANSIONHUNTER REPEAT EXPANSION PIPELINE v1.1          ║
╚══════════════════════════════════════════════════════════════╝

INPUT FILES
  PED file        : ${params.ped ?: 'NOT SET'}
  BAM directory   : ${params.bam_dir ?: 'NOT SET'}
  BAM suffix      : ${params.bam_suffix}

REFERENCE
  Genome          : ${params.genome}
  FASTA           : ${params.fasta ?: 'NOT SET'}

EXPANSIONHUNTER
  Variant catalog : ${params.variant_catalog ?: 'NOT SET'}
  Analysis mode   : ${params.analysis_mode}
  Regions BED     : ${params.regions_bed ?: 'NOT SET (run generate_regions_bed.py first)'}

STRIPY POST-PROCESSING
  Locus ref       : ${params.stripy_locus_ref}
  HPO file        : ${params.hpo_file}
  HPO filter      : ${params.hpo_filter_file ?: 'disabled'}

OUTPUT
  Output dir      : ${params.outdir}
  Batch name      : ${params.batch_name}

RESOURCES
  Max CPUs        : ${params.max_cpus}
  Max memory      : ${params.max_memory}
  Queue size      : ${params.slurm_queue_size}
──────────────────────────────────────────────────────────────
""".stripIndent()

// ============================================================================
// PARAMETER VALIDATION
// ============================================================================

def required_params = ['ped', 'bam_dir', 'fasta', 'variant_catalog', 'regions_bed']
required_params.each { p ->
    if (!params[p]) {
        log.error "[PIPELINE ERROR] Missing required parameter: --${p}"
        System.exit(1)
    }
}

[
    [params.ped,             'PED file'],
    [params.bam_dir,         'BAM directory'],
    [params.fasta,           'Reference FASTA'],
    [params.variant_catalog, 'Variant catalog'],
    [params.regions_bed,     'Repeat regions BED'],
    [params.hpo_file,        'HPO file'],
    [params.stripy_locus_ref,'Stripy locus reference'],
].each { path, label ->
    if (path && !file(path).exists()) {
        log.error "[PIPELINE ERROR] ${label} not found: ${path}"
        System.exit(1)
    }
}

// ============================================================================
// INCLUDE MODULES
// ============================================================================

include { SAMTOOLS_INDEX          } from './modules/samtools'
include { SAMTOOLS_EXTRACT_REGIONS} from './modules/samtools'
include { SAMTOOLS_INDEX_REALIGNED} from './modules/samtools'
include { EXPANSIONHUNTER         } from './modules/expansionhunter'
include { STRIPY_POSTPROCESS      } from './subworkflows/local/stripy_postprocess'

// ============================================================================
// HELPER FUNCTIONS
// ============================================================================

def parse_ped(ped_file) {
    def samples = []
    file(ped_file).eachLine { line ->
        line = line.trim()
        if (!line || line.startsWith('#')) return

        def cols = line.tokenize()
        if (cols.size() < 6) {
            log.warn "[PED] Skipping malformed line (need ≥6 columns): ${line}"
            return
        }

        samples << [
            family_id : cols[0],
            sample_id : cols[1],
            father_id : cols[2],
            mother_id : cols[3],
            sex       : cols[4],
            phenotype : cols[5]
        ]
    }

    if (samples.isEmpty()) {
        log.error "[PIPELINE ERROR] No valid samples found in PED file: ${ped_file}"
        System.exit(1)
    }

    log.info "[PED] Loaded ${samples.size()} samples from ${ped_file}"
    return samples
}

def resolve_bam_files(sample_id, bam_dir, bam_suffix) {
    def bam = file("${bam_dir}/${sample_id}${bam_suffix}")

    if (!bam.exists()) {
        log.warn "[BAM] Not found for sample '${sample_id}': ${bam}"
        return null
    }

    def bai = file("${bam}.bai")
    if (!bai.exists()) {
        def stem = bam_suffix.replaceAll(/\.bam$/, '')
        bai = file("${bam_dir}/${sample_id}${stem}.bai")
    }
    if (!bai.exists()) {
        log.warn "[BAI] Index not found for ${bam.name} - will create with samtools index"
        bai = null
    }

    return [bam, bai]
}

// ============================================================================
// MAIN WORKFLOW
// ============================================================================

workflow {

    // ── Parse PED file ────────────────────────────────────────────────────
    def ped_samples = parse_ped(params.ped)

    // ── Match samples to BAM files ────────────────────────────────────────
    def matched_samples = ped_samples.collect { s ->
        def result = resolve_bam_files(s.sample_id, params.bam_dir, params.bam_suffix)
        if (!result) return null

        def (bam, bai) = result
        def meta = [
            id        : s.sample_id,
            family_id : s.family_id,
            sex       : s.sex,
            phenotype : s.phenotype
        ]
        return [meta, bam, bai]
    }.findAll { it != null }

    if (matched_samples.isEmpty()) {
        log.error "[PIPELINE ERROR] No BAM files matched. Check --bam_dir and --bam_suffix"
        System.exit(1)
    }

    log.info "[PIPELINE] ${matched_samples.size()} / ${ped_samples.size()} samples matched"

    // ── Create channels ───────────────────────────────────────────────────
    def (samples_with_bai, samples_without_bai) = matched_samples
        .split { meta, bam, bai -> bai != null }

    ch_to_index = Channel.fromList(samples_without_bai ?: [])
        .map { meta, bam, bai -> [meta, bam] }

    ch_indexed = Channel.fromList(samples_with_bai ?: [])

    // Reference files
    ch_fasta     = Channel.value(file(params.fasta))
    ch_fasta_fai = params.fasta_fai ?
        Channel.value(file(params.fasta_fai)) :
        Channel.value(file("${params.fasta}.fai"))
    ch_catalog   = Channel.value(file(params.variant_catalog))

    // Repeat regions BED - single file broadcast to all extract jobs
    ch_regions_bed = Channel.value(file(params.regions_bed))

    // ── Index BAM files if needed ─────────────────────────────────────────
    SAMTOOLS_INDEX(ch_to_index)

    // ── Combine all samples (indexed + already-indexed) ───────────────────
    ch_all_samples = ch_indexed.mix(SAMTOOLS_INDEX.out.bam_bai)

    // ── Extract repeat-locus reads into mini-BAM ──────────────────────────
    // Reduces BAM from ~5-15GB to ~50-100MB before EH seeking mode
    SAMTOOLS_EXTRACT_REGIONS(
        ch_all_samples,
        ch_regions_bed
    )

    // ── Run ExpansionHunter in seeking mode on mini-BAM ───────────────────
    // seeking mode produces realigned BAM required for REViewer
    EXPANSIONHUNTER(
        SAMTOOLS_EXTRACT_REGIONS.out.mini_bam,
        ch_fasta,
        ch_fasta_fai,
        ch_catalog
    )

    // ── Index realigned BAM using samtools container ───────────────────────
    // expansionhunter.sif does not include samtools - index separately
    // filter() ensures this only runs for samples that produced a realigned BAM
    SAMTOOLS_INDEX_REALIGNED(
        EXPANSIONHUNTER.out.realigned_bam
    )

    // ── Generate completion manifest ──────────────────────────────────────
    // collectFile is a Nextflow operator - no process/container needed
    // Writes one row per sample as they complete, header seeded once
    EXPANSIONHUNTER.out.json
        .map { meta, json ->
            "${meta.id}\t${meta.family_id}\t${meta.sex}\t${meta.phenotype}\t${json}"
        }
        .collectFile(
            name:     'completed_samples.tsv',
            storeDir: "${params.outdir}/pipeline_info",
            newLine:  true,
            seed:     "sample_id\tfamily_id\tsex\tphenotype\tjson_path"
        )

    // ── STRipy post-processing ────────────────────────────────────────────
    STRIPY_POSTPROCESS(
        EXPANSIONHUNTER.out.vcf,
        ch_catalog,
        ch_fasta,
        file(params.ped)
    )
}

// ============================================================================
// WORKFLOW COMPLETION
// ============================================================================

workflow.onComplete {
    def status = workflow.success ? 'SUCCESS ✔' : 'FAILED ✘'

    log.info """
──────────────────────────────────────────────────────────────
  Pipeline Status : ${status}
  Duration        : ${workflow.duration}
  Output dir      : ${params.outdir}
  Exit status     : ${workflow.exitStatus}
──────────────────────────────────────────────────────────────
""".stripIndent()
}

workflow.onError {
    log.error "[PIPELINE ERROR] ${workflow.errorReport}"
}
