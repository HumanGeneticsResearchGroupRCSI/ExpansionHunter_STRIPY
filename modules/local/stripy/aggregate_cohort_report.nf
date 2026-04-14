/*
========================================================================================
    STRIPY - AGGREGATE COHORT REPORT
========================================================================================
    Collect all per-sample TSV reports and generate per-family Excel workbooks.
    One Excel per family: Summary sheet + one sheet per family member.
    Runs once after all per-sample reports are complete.
----------------------------------------------------------------------------------------
*/

process STRIPY_AGGREGATE_COHORT_REPORT {
    tag        "cohort_report"
    label      'process_aggregate'

    container  params.containers.python

    publishDir "${params.outdir}/cohort",
               mode: 'copy'

    input:
    path tsv_files   // collected list of all *.stripy_report.tsv files
    path ped         // PED file for family structure

    output:
    path "*.xlsx",        emit: cohort_excel
    path "versions.yml",  emit: versions

    script:
    """
    # Collect all TSVs into a working directory
    mkdir -p tsv_input
    for f in ${tsv_files}; do
        cp \$f tsv_input/
    done

    python3 ${projectDir}/scripts/aggregate_cohort_report.py \
        --tsv-dir   tsv_input/ \
        --ped       ${ped} \
        --outdir    .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | sed 's/Python //')
        openpyxl: \$(python3 -c "import openpyxl; print(openpyxl.__version__)")
    END_VERSIONS
    """

    stub:
    """
    touch cohort.stripy_report.xlsx
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11.0
        openpyxl: 3.1.0
    END_VERSIONS
    """
}
