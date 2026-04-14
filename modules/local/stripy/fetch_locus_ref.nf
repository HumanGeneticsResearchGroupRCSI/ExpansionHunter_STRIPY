/*
========================================================================================
    STRIPY - FETCH LOCUS REFERENCE
========================================================================================
    Fetch /locus metadata from the Stripy API for all loci in the variant catalog.
    Runs ONCE per pipeline run - result is cached and reused across all samples.

    Behaviour:
      - If params.stripy_locus_ref already exists on disk → copy to work dir (fast path)
      - If not → run fetch_stripy_locus_ref.py to download from API → save permanently

    Output:
      stripy_locus_reference.json   - emitted to downstream STRIPY_BUILD_SAMPLE_REPORT
      versions.yml
----------------------------------------------------------------------------------------
*/

process STRIPY_FETCH_LOCUS_REF {
    tag   "stripy_locus_ref"
    label 'process_low'

    container params.containers.python

    // Not using publishDir - file is written directly to params.stripy_locus_ref
    // by the script block below (permanent resource location)

    input:
    path variant_catalog

    output:
    path "stripy_locus_reference.json", emit: locus_ref
    path "versions.yml",                emit: versions

    script:
    def permanent_path = params.stripy_locus_ref
    """
    if [ -f "${permanent_path}" ]; then
        echo "[INFO] Locus reference already exists: ${permanent_path}"
        echo "[INFO] Copying to work directory for pipeline use"
        cp "${permanent_path}" stripy_locus_reference.json
    else
        echo "[INFO] Locus reference not found - fetching from Stripy API"
        python3 ${projectDir}/scripts/fetch_stripy_locus_ref.py \\
            --catalog ${variant_catalog} \\
            --output  stripy_locus_reference.json

        # Persist to permanent resource path for future runs (avoids re-downloading)
        mkdir -p \$(dirname "${permanent_path}")
        cp stripy_locus_reference.json "${permanent_path}"
        echo "[INFO] Saved to: ${permanent_path}"
    fi

    # Validate output is non-empty JSON array
    python3 -c "
import json, sys
with open('stripy_locus_reference.json') as f:
    data = json.load(f)
if not isinstance(data, list) or len(data) == 0:
    print('[ERROR] stripy_locus_reference.json is empty or not a list', file=sys.stderr)
    sys.exit(1)
print(f'[INFO] Locus reference validated: {len(data)} loci')
"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | sed 's/Python //')
        stripy_locus_api: "https://api.stripy.org/locus"
    END_VERSIONS
    """

    stub:
    """
    echo '[]' > stripy_locus_reference.json
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.14.2
        stripy_locus_api: "https://api.stripy.org/locus"
    END_VERSIONS
    """
}
