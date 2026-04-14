#!/usr/bin/env python3
"""
fetch_stripy_locus_ref.py
=========================
One-time bootstrap: fetch /locus metadata from Stripy API for all loci
in the ExpansionHunter variant catalog and save as a local JSON reference.

Called by the FETCH_STRIPY_LOCUS_REF Nextflow module.
Logic:
  - If output file already exists → print message and exit 0 (skip download)
  - If not exists → fetch all loci and save

Usage
-----
    python3 fetch_stripy_locus_ref.py \\
        --catalog  /path/to/variant_catalog.json \\
        --output   /home/data/ONE/resources/stripy_locus_reference.json

Version : 1.0.0
"""

import argparse
import json
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

STRIPY_API     = "https://api.stripy.org"
API_TIMEOUT    = 15
API_RETRY      = 3
API_RETRY_WAIT = 5
POLITE_DELAY   = 0.25   # seconds between requests - avoid rate limiting

# VCF REPID → Stripy API locus ID normalisation
LOCUS_ID_MAP = {
    'RFC1:AAGGG': 'RFC1',
    'RFC1:ACAGG': 'RFC1',
}


def normalise_locus_id(locus_id: str) -> str:
    return LOCUS_ID_MAP.get(locus_id, locus_id)


def fetch_locus(locus_id: str) -> dict:
    """
    GET /locus/{locus_id} with retry.
    Returns parsed JSON or {"_error": "<reason>", "_locus_id": locus_id}.
    """
    api_id = normalise_locus_id(locus_id)
    url    = f"{STRIPY_API}/locus/{api_id}"

    for attempt in range(1, API_RETRY + 1):
        try:
            with urllib.request.urlopen(url, timeout=API_TIMEOUT) as resp:
                data = json.loads(resp.read())
                data['_source_repid'] = locus_id   # keep original VCF ID
                return data
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return {"_error": "NOT_FOUND", "_locus_id": locus_id}
            if attempt < API_RETRY:
                time.sleep(API_RETRY_WAIT)
                continue
            return {"_error": f"HTTP_{e.code}", "_locus_id": locus_id}
        except Exception as exc:
            if attempt < API_RETRY:
                time.sleep(API_RETRY_WAIT)
                continue
            return {"_error": str(exc), "_locus_id": locus_id}

    return {"_error": "max_retries_exceeded", "_locus_id": locus_id}


def parse_args():
    p = argparse.ArgumentParser(
        description="Fetch Stripy /locus metadata for all catalog loci"
    )
    p.add_argument('--catalog', required=True,
                   help="ExpansionHunter variant catalog JSON")
    p.add_argument('--output',  required=True,
                   help="Output path for stripy_locus_reference.json")
    return p.parse_args()


def main():
    args = parse_args()
    out_path = Path(args.output)

    # ── Check if already exists ───────────────────────────────────────────
    if out_path.exists():
        print(f"[INFO] Locus reference already exists: {out_path}")
        print(f"[INFO] Skipping download - delete file to force refresh")
        sys.exit(0)

    # ── Load catalog ──────────────────────────────────────────────────────
    catalog_path = Path(args.catalog)
    if not catalog_path.exists():
        print(f"[ERROR] Catalog not found: {catalog_path}", file=sys.stderr)
        sys.exit(1)

    with open(catalog_path) as f:
        catalog = json.load(f)

    # Extract unique locus IDs - deduplicate RFC1 variants
    seen     = set()
    locus_ids = []
    for entry in catalog:
        if not isinstance(entry, dict) or 'LocusId' not in entry:
            continue
        raw_id = entry['LocusId']
        api_id = normalise_locus_id(raw_id)
        if api_id not in seen:
            seen.add(api_id)
            locus_ids.append(raw_id)    # keep original for display

    print(f"[INFO] Catalog          : {catalog_path}")
    print(f"[INFO] Unique loci      : {len(locus_ids)}")
    print(f"[INFO] Output           : {out_path}")
    print(f"[INFO] Starting fetch...")
    print()

    # ── Fetch all loci ────────────────────────────────────────────────────
    results   = []
    ok        = not_found = errors = 0
    log_lines = []

    for i, locus_id in enumerate(locus_ids, 1):
        api_id = normalise_locus_id(locus_id)
        print(f"  [{i:>3}/{len(locus_ids)}] {locus_id:<25} (API: {api_id})",
              end='', flush=True)

        data = fetch_locus(locus_id)

        if '_error' not in data:
            status = "OK"
            results.append(data)
            ok += 1
        else:
            status = data['_error']
            if status == "NOT_FOUND":
                not_found += 1
            else:
                errors += 1

        print(f"  {status}")
        log_lines.append(f"{locus_id}\t{api_id}\t{status}")
        time.sleep(POLITE_DELAY)

    # ── Write output ──────────────────────────────────────────────────────
    out_path.parent.mkdir(parents=True, exist_ok=True)

    with open(out_path, 'w') as f:
        json.dump(results, f, indent=2)

    log_path = out_path.with_suffix('.log')
    with open(log_path, 'w') as f:
        f.write("LocusId\tApiId\tStatus\n")
        f.write('\n'.join(log_lines) + '\n')

    # ── Summary ───────────────────────────────────────────────────────────
    print()
    print(f"{'='*60}")
    print(f"  OK              : {ok}")
    print(f"  Not in DB       : {not_found}  "
          f"← these loci will have '.' for metadata fields")
    print(f"  Errors          : {errors}")
    print(f"  Written         : {out_path}")
    print(f"  Log             : {log_path}")
    print(f"{'='*60}")

    if errors > 0:
        print(f"\n[WARN] {errors} loci failed to fetch - check log file")
        sys.exit(1)

    print(f"\n[INFO] Done - reference file ready for pipeline use")


if __name__ == '__main__':
    main()
