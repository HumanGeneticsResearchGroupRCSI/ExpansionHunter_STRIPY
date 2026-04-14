/*
========================================================================================
    STRIPY - ANNOTATE VCF
========================================================================================
*/

process STRIPY_ANNOTATE_VCF {
    tag   "${meta.id}"
    label 'process_low'

    container params.containers.python

    publishDir "${params.outdir}/stripy/annotated_vcf",
               mode: params.publish_dir_mode,
               pattern: "*.Reannotated.vcf"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("${meta.id}.Reannotated.vcf"), emit: annotated_vcf
    path "versions.yml",                                  emit: versions

    script:
    """
    python3 - <<'PYEOF'
import json
import sys
import time
import urllib.error
import urllib.request
import uuid
from pathlib import Path

vcf_path    = "${vcf}"
sample_id   = "${meta.id}"
out_path    = f"{sample_id}.Reannotated.vcf"
api_url     = "https://api.stripy.org/annotateVCF"
max_retries = 3
retry_wait  = 10
timeout     = 120

def post_vcf(vcf_path: str, url: str, timeout: int) -> bytes:
    # POST a VCF file to the Stripy /annotateVCF endpoint as multipart/form-data.
    # Returns response bytes.
    boundary = uuid.uuid4().hex
    vcf_bytes = Path(vcf_path).read_bytes()
    vcf_name  = Path(vcf_path).name

    body = (
        f"--{boundary}\\r\\n"
        f'Content-Disposition: form-data; name="file"; filename="{vcf_name}"\\r\\n'
        f"Content-Type: text/plain\\r\\n"
        f"\\r\\n"
    ).encode() + vcf_bytes + f"\\r\\n--{boundary}--\\r\\n".encode()

    req = urllib.request.Request(
        url,
        data=body,
        headers={"Content-Type": f"multipart/form-data; boundary={boundary}"},
        method="POST",
    )
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return resp.read()

# Retry loop
last_error = None
for attempt in range(1, max_retries + 1):
    try:
        print(f"[INFO] Attempt {attempt}/{max_retries}: POSTing {vcf_path} to {api_url}")
        response_bytes = post_vcf(vcf_path, api_url, timeout)
        response_text  = response_bytes.decode("utf-8", errors="replace")
        break
    except urllib.error.HTTPError as e:
        last_error = f"HTTP {e.code}: {e.reason}"
        print(f"[WARN] Attempt {attempt} failed: {last_error}", file=sys.stderr)
    except Exception as e:
        last_error = str(e)
        print(f"[WARN] Attempt {attempt} failed: {last_error}", file=sys.stderr)

    if attempt < max_retries:
        print(f"[INFO] Retrying in {retry_wait}s...", file=sys.stderr)
        time.sleep(retry_wait)
else:
    print(f"[ERROR] All {max_retries} attempts failed for {sample_id}: {last_error}", file=sys.stderr)
    sys.exit(1)

# Validate response is a VCF
if "#CHROM" not in response_text:
    print(f"[ERROR] Stripy /annotateVCF returned invalid response for {sample_id}", file=sys.stderr)
    sys.exit(1)

# Count annotated loci
n_lines      = sum(1 for l in response_text.splitlines() if not l.startswith('#'))
n_annotated  = sum(1 for l in response_text.splitlines() if not l.startswith('#') and 'DISID=' in l)

Path(out_path).write_text(response_text)
print(f"[INFO] Written  : {out_path}")
print(f"[INFO] VCF rows : {n_lines}  annotated: {n_annotated}")
PYEOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | sed 's/Python //')
        stripy_annotateVCF_api: "https://api.stripy.org/annotateVCF"
    END_VERSIONS
    """
}
