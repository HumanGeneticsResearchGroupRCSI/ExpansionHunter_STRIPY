#!/usr/bin/env python3
"""
generate_regions_bed.py
=======================
Extract repeat loci genomic regions from an ExpansionHunter variant catalog
and write a padded, merged BED file for use with samtools view -b -L.

The resulting BED file is used to extract only repeat-locus-covering reads
from a full WES BAM before running EH in seeking mode (hybrid approach).

Usage
-----
    python3 generate_regions_bed.py \\
        --catalog  variant_catalog_v2.json \\
        --output   repeat_regions.bed \\
        --padding  1000

Output
------
    BED file (0-based, half-open) with padded, merged, sorted intervals.
    One line per interval, sorted by chrom then start.
    Duplicate/overlapping intervals are merged.

Notes
-----
    - Only canonical chromosomes (chr1-22, chrX, chrY, chrM) are included
    - ReferenceRegion format: 'chrN:start-end' (0-based) or list of same
    - RFC1 appears as RFC1:AAGGG and RFC1:ACAGG - same coords → deduplicated
    - Padding is applied symmetrically; start clamped to 0

Version: 1.0.0
"""

import argparse
import json
import sys
from pathlib import Path

CANONICAL_CHROMS = set(
    [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY', 'chrM']
)


def parse_reference_region(region) -> list:
    """
    Parse ReferenceRegion field - handles both string and list formats.
    Returns list of (chrom, start, end) tuples (0-based, half-open).

    Examples:
        'chr6:16327635-16327722'         → [('chr6', 16327635, 16327722)]
        ['chr4:39348424-39348485']       → [('chr4', 39348424, 39348485)]
        'chrX:25013529-25013565'         → [('chrX', 25013529, 25013565)]
    """
    regions = [region] if isinstance(region, str) else region
    result  = []
    for r in regions:
        r = r.strip()
        if not r:
            continue
        try:
            chrom, coords = r.rsplit(':', 1)
            start_str, end_str = coords.split('-')
            start = int(start_str)
            end   = int(end_str)
            result.append((chrom, start, end))
        except (ValueError, AttributeError):
            print(f"[WARN] Cannot parse ReferenceRegion: {r!r} - skipping",
                  file=sys.stderr)
    return result


def merge_intervals(intervals: list) -> list:
    """
    Sort and merge overlapping/adjacent intervals.
    Input/output: list of (chrom, start, end) tuples.
    """
    if not intervals:
        return []

    # Sort by chrom then start
    sorted_ivs = sorted(intervals, key=lambda x: (x[0], x[1]))
    merged = [sorted_ivs[0]]

    for chrom, start, end in sorted_ivs[1:]:
        last_chrom, last_start, last_end = merged[-1]
        if chrom == last_chrom and start <= last_end:
            # Overlapping or adjacent - extend
            merged[-1] = (last_chrom, last_start, max(last_end, end))
        else:
            merged.append((chrom, start, end))

    return merged


def main():
    p = argparse.ArgumentParser(
        description="Generate padded BED file from EH variant catalog"
    )
    p.add_argument('--catalog', required=True,
                   help="ExpansionHunter variant catalog JSON")
    p.add_argument('--output',  required=True,
                   help="Output BED file path")
    p.add_argument('--padding', type=int, default=1000,
                   help="Padding (bp) to add on each side of each locus (default: 1000)")
    args = p.parse_args()

    catalog_path = Path(args.catalog)
    if not catalog_path.exists():
        print(f"[ERROR] Catalog not found: {catalog_path}", file=sys.stderr)
        sys.exit(1)

    with open(catalog_path) as f:
        catalog = json.load(f)

    # Normalise: catalog may be list or dict
    if isinstance(catalog, dict):
        entries = list(catalog.values())
    elif isinstance(catalog, list):
        entries = catalog
    else:
        print("[ERROR] Unexpected catalog format (not list or dict)", file=sys.stderr)
        sys.exit(1)

    # Collect all padded intervals
    raw_intervals = []
    skipped       = 0

    for entry in entries:
        locus_id = entry.get('LocusId', 'unknown')
        ref_region = entry.get('ReferenceRegion')

        if not ref_region:
            print(f"[WARN] No ReferenceRegion for locus {locus_id} - skipping",
                  file=sys.stderr)
            skipped += 1
            continue

        regions = parse_reference_region(ref_region)

        for chrom, start, end in regions:
            # Skip non-canonical contigs
            if chrom not in CANONICAL_CHROMS:
                print(f"[WARN] Skipping non-canonical contig: {chrom} ({locus_id})",
                      file=sys.stderr)
                skipped += 1
                continue

            # Apply padding - clamp start to 0
            padded_start = max(0, start - args.padding)
            padded_end   = end + args.padding

            raw_intervals.append((chrom, padded_start, padded_end))

    # Deduplicate before merging (handles RFC1:AAGGG / RFC1:ACAGG etc.)
    unique_intervals = list(set(raw_intervals))

    # Sort and merge overlapping intervals
    merged = merge_intervals(unique_intervals)

    # Write BED file
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    with open(out_path, 'w') as f:
        for chrom, start, end in merged:
            f.write(f"{chrom}\t{start}\t{end}\n")

    # Summary
    print(f"[INFO] Catalog entries   : {len(entries)}")
    print(f"[INFO] Raw intervals     : {len(raw_intervals)}")
    print(f"[INFO] After dedup/merge : {len(merged)}")
    print(f"[INFO] Skipped           : {skipped}")
    print(f"[INFO] Padding           : {args.padding} bp")
    print(f"[INFO] Written to        : {out_path}")


if __name__ == '__main__':
    main()
