#!/usr/bin/env python3
"""
build_sample_report.py
======================
Generate a per-sample STR annotation TSV report by combining:

  1. STRipy-annotated ExpansionHunter VCF  (*.Reannotated.vcf)
  2. Pre-fetched Stripy /locus metadata    (stripy_locus_reference.json)
  3. Local HPO gene-to-phenotype file      (genes_to_phenotype.txt)
  4. PED file                              (family relationships for DeNovo)
  5. All sample VCFs directory             (parent REPCN lookup for DeNovo)
  6. Live Stripy /compare API              (called loci only)

Output: per-sample TSV - one row per disease per locus.

Usage
-----
    python3 build_sample_report.py \\
        --vcf        LH0014A.Reannotated.vcf \\
        --sample-id  LH0014A \\
        --locus-ref  stripy_locus_reference.json \\
        --hpo-file   genes_to_phenotype.txt \\
        --ped        cohort.ped \\
        --vcf-dir    /path/to/all/reannotated/vcfs \\
        --outdir     results/stripy/

Version : 1.0.0
"""

import argparse
import json
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path


# ── Constants ─────────────────────────────────────────────────────────────────

STRIPY_API     = "https://api.stripy.org"
API_TIMEOUT    = 15
API_RETRY      = 3
API_RETRY_WAIT = 5

NO_CALL_TOKENS = {'.', './.', ''}

LOCUS_ID_MAP = {
    'RFC1:AAGGG': 'RFC1',
    'RFC1:ACAGG': 'RFC1',
}

OUTPUT_COLUMNS = [
    # Section A - positional + EH
    "SampleID",
    "Chr", "Start", "End", "LocationCoordinates",
    "Ref", "Alt",
    "GT", "LocusCov", "RepeatLength",
    "FilterStatus", "CallStatus", "CallQuality",
    # Section A - locus identity
    "Locus", "LocationRegion", "RepeatType", "Motif",
    # Section A - clinical ranges + disease
    "Normal Range (Min)", "Normal Range (Max)",
    "IntermediateRange (Min)", "IntermediateRange (Max)",
    "PathogenicCutoff",
    "Disease Name", "Disease OMIM",
    "Inheritance", "Onset",
    "Stripy_Gene", "HPO_Gene",
    "HPO Phenotype",
    "HPO_Filtered",
    # Section B - per-allele repeat counts
    "Allele1_Repeats", "Allele2_Repeats",
    # Section B - literature
    "Literature:DiseaseName",
    "Literature:Inheritance",
    "Allele1_Range", "Allele2_Range",
    # Section B - per-allele population stats
    "Allele1_Zscore", "Allele2_Zscore",
    "Allele1_Outlier", "Allele2_Outlier",
    "Allele1_Direction", "Allele2_Direction",
    # Section C - per-allele DeNovo
    "Allele1_DeNovo", "Allele1_DeNovo_Direction",
    "Allele2_DeNovo", "Allele2_DeNovo_Direction",
    "DeNovo_Status",
    "Caller",
]


# ── Helpers ───────────────────────────────────────────────────────────────────

def normalise_locus_id(locus_id: str) -> str:
    return LOCUS_ID_MAP.get(locus_id, locus_id)


def fmt(val) -> str:
    """Serialise value for TSV output - handles list, None, bool."""
    if val is None:
        return '.'
    if isinstance(val, list):
        return '/'.join('.' if v is None else str(v) for v in val)
    return str(val)


# ── Stripy /compare API ───────────────────────────────────────────────────────

def api_get(endpoint: str) -> dict:
    url = f"{STRIPY_API}{endpoint}"
    for attempt in range(1, API_RETRY + 1):
        try:
            with urllib.request.urlopen(url, timeout=API_TIMEOUT) as resp:
                return json.loads(resp.read())
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return {"_error": "404_not_found"}
            if attempt < API_RETRY:
                time.sleep(API_RETRY_WAIT)
                continue
            return {"_error": f"HTTP_{e.code}"}
        except Exception as exc:
            if attempt < API_RETRY:
                time.sleep(API_RETRY_WAIT)
                continue
            return {"_error": str(exc)}
    return {"_error": "max_retries_exceeded"}


def fetch_compare(locus_id: str, rep1: str, rep2) -> dict:
    api_id = normalise_locus_id(locus_id)
    if rep2 is not None and rep2 != rep1:
        return api_get(f"/compare/{api_id}/{rep1}/{rep2}")
    return api_get(f"/compare/{api_id}/{rep1}")


# ── VCF parsing ───────────────────────────────────────────────────────────────

def parse_info(info_str: str) -> dict:
    info = {}
    for field in info_str.split(';'):
        if '=' in field:
            k, v = field.split('=', 1)
            info[k] = v
        else:
            info[field] = True
    return info


def parse_format(fmt_keys: str, fmt_vals: str) -> dict:
    keys = fmt_keys.split(':')
    vals = fmt_vals.split(':')
    vals += ['.'] * (len(keys) - len(vals))
    return dict(zip(keys, vals))


def parse_repcn(repcn_str: str):
    """
    Returns (rep1, rep2) strings or (None, None) for no-call.
    Both rep1 and rep2 are strings representing integer repeat counts.
    """
    if repcn_str in NO_CALL_TOKENS:
        return None, None
    parts = repcn_str.split('/')
    for p in parts:
        if p == '.' or p == '' or not p.lstrip('-').isdigit():
            return None, None
    if len(parts) == 2:
        return parts[0], parts[1]
    return parts[0], None


def read_vcf(vcf_path: str):
    """Generator - yields parsed record dicts from a VCF file."""
    with open(vcf_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            cols = line.strip().split('\t')
            if len(cols) < 10:
                continue
            chrom, pos, _, ref, alt, _, filt, info_str, fmt_keys, fmt_vals = \
                cols[:10]
            yield {
                'chrom':  chrom,
                'pos':    pos,
                'ref':    ref,
                'alt':    alt,
                'filter': filt,
                'info':   parse_info(info_str),
                'fmt':    parse_format(fmt_keys, fmt_vals),
            }


def parse_multi_disease(info: dict) -> list:
    """
    Expand multi-disease INFO fields into list of per-disease dicts.
    Returns single empty entry if locus not in Stripy DB (no DISID).
    """
    if 'DISID' not in info:
        return [{"dis_id": None, "dis_name": None,
                 "disinher": None, "dis_range": None}]

    dis_ids    = info.get('DISID',    '').split(',')
    dis_names  = [n.replace('_', ' ')
                  for n in info.get('DISNAME', '').split(',')]
    dis_inher  = info.get('DISINHER', '').split(',')
    dis_ranges = info.get('DISRANGE', '').split(',')

    n = len(dis_ids)

    def pad(lst):
        return lst + [None] * (n - len(lst))

    return [
        {
            "dis_id":    dis_ids[i],
            "dis_name":  pad(dis_names)[i],
            "disinher":  pad(dis_inher)[i],
            "dis_range": pad(dis_ranges)[i],
        }
        for i in range(n)
    ]


# ── Locus reference ───────────────────────────────────────────────────────────

def load_locus_reference(ref_path: str) -> dict:
    p = Path(ref_path)
    if not p.exists():
        print(f"[ERROR] Locus reference not found: {ref_path}", file=sys.stderr)
        sys.exit(1)
    with open(p) as f:
        data = json.load(f)
    return {entry['Locus']: entry
            for entry in data
            if isinstance(entry, dict) and 'Locus' in entry}


def get_locus_meta(locus_id: str, locus_ref: dict) -> dict:
    return locus_ref.get(normalise_locus_id(locus_id), {})


# ── HPO lookup ────────────────────────────────────────────────────────────────

def load_hpo_file(hpo_path: str) -> dict:
    """
    Load genes_to_phenotype.txt.

    Indexed by disease_id ONLY (e.g. "OMIM:143100") - NOT by gene symbol.
    This handles cases where Stripy's gene differs from HPO's gene for the
    same disease (e.g. ABCD3 vs LRP12 for OMIM:164310).

    Structure:
        {
          "OMIM:143100": {"gene": "HTT",   "terms": ["HP:0001250|Seizure", ...]},
          "OMIM:164310": {"gene": "LRP12", "terms": ["HP:0002460|Distal...", ...]},
        }
    """
    p = Path(hpo_path)
    if not p.exists():
        print(f"[ERROR] HPO file not found: {hpo_path}", file=sys.stderr)
        sys.exit(1)
    lookup = {}
    with open(p) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#') or \
               line.startswith('ncbi_gene_id'):
                continue
            parts = line.split('\t')
            if len(parts) < 6:
                continue
            gene_symbol = parts[1].strip()
            hpo_id      = parts[2].strip()
            hpo_name    = parts[3].strip()
            disease_id  = parts[5].strip()   # e.g. "OMIM:143100"
            pair        = f"{hpo_id}|{hpo_name}"
            if disease_id not in lookup:
                lookup[disease_id] = {"gene": gene_symbol, "terms": []}
            if pair not in lookup[disease_id]["terms"]:
                lookup[disease_id]["terms"].append(pair)
    return lookup


def get_hpo_terms(disease_omim: str, hpo_lookup: dict) -> tuple:
    """
    Look up HPO terms by disease OMIM ID alone - gene symbol not used for matching.

    Returns (hpo_terms_string, hpo_gene_symbol):
        hpo_terms_string : '; '-separated "HP:XXXXXXX|Name" pairs, or '.'
        hpo_gene_symbol  : gene from HPO file that matched, or '.' if no match

    Tries OMIM: prefix first, then ORPHA:.
    """
    if not disease_omim or disease_omim == '.':
        return '.', '.'
    for prefix in ('OMIM', 'ORPHA'):
        entry = hpo_lookup.get(f"{prefix}:{disease_omim}")
        if entry:
            return '; '.join(entry["terms"]), entry["gene"]
    return '.', '.'


# ── HPO filter ────────────────────────────────────────────────────────────────

def load_hpo_filter_file(filter_path: str) -> set:
    """
    Load a project-specific HPO filter file.
    Returns a set of HP: term IDs to keep (e.g. {"HP:0001250", "HP:0001251"}).
    File format: HP:XXXXXXX<tab>Name<tab>Category  (comments # ignored)
    Returns empty set if path is None or file not found (no filtering applied).
    """
    if not filter_path:
        return set()
    p = Path(filter_path)
    if not p.exists():
        print(f"[WARN] HPO filter file not found: {filter_path} - HPO_Filtered will be '.'",
              file=sys.stderr)
        return set()
    terms = set()
    with open(p) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split('	')
            hp_id = parts[0].strip()
            if hp_id.startswith('HP:'):
                terms.add(hp_id)
    print(f"[INFO] HPO filter: loaded {len(terms)} terms from {p.name}",
          file=sys.stderr)
    return terms


def get_hpo_filtered(hpo_terms_str: str, filter_terms: set) -> str:
    """
    Filter HPO terms string to only those in filter_terms set.
    hpo_terms_str: '; '-separated "HP:XXXXXXX|Name" pairs
    Returns filtered '; '-separated string, or '.' if no matches or no filter.
    """
    if not filter_terms or not hpo_terms_str or hpo_terms_str == '.':
        return '.'
    matches = []
    for term in hpo_terms_str.split('; '):
        term = term.strip()
        if not term:
            continue
        hp_id = term.split('|')[0].strip()
        if hp_id in filter_terms:
            matches.append(term)
    return '; '.join(matches) if matches else '.'


# ── Z-score direction ─────────────────────────────────────────────────────────

def zscore_direction(zscore) -> str:
    if zscore is None:
        return '.'
    try:
        if isinstance(zscore, list):
            return '/'.join(
                '.' if z is None else ('HIGH' if float(z) > 0 else 'LOW')
                for z in zscore
            )
        return 'HIGH' if float(zscore) > 0 else 'LOW'
    except (TypeError, ValueError):
        return '.'



def zscore_direction_single(zscore) -> str:
    """Return HIGH/LOW/. for a single scalar Z-score value."""
    if zscore is None:
        return '.'
    try:
        return 'HIGH' if float(zscore) > 0 else 'LOW'
    except (TypeError, ValueError):
        return '.'


# ── PED file parsing ──────────────────────────────────────────────────────────

def load_ped(ped_path: str) -> dict:
    """
    Parse PED file into dict keyed by sample_id.
    Returns:
        {
          "LH0014A": {
              "family_id": "FAM001",
              "father_id": "LH0014B",
              "mother_id": "LH0014C",
              "sex": "1"
          }, ...
        }
    Father/mother = "0" means not available.
    """
    p = Path(ped_path)
    if not p.exists():
        print(f"[ERROR] PED file not found: {ped_path}", file=sys.stderr)
        sys.exit(1)
    ped = {}
    with open(p) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            cols = line.split()
            if len(cols) < 6:
                continue
            fam_id, sample_id, father_id, mother_id, sex, _ = cols[:6]
            ped[sample_id] = {
                "family_id": fam_id,
                "father_id": father_id,
                "mother_id": mother_id,
                "sex":       sex,
            }
    return ped


def get_family_type(sample_id: str, ped: dict) -> str:
    """Return TRIO / DUO / SINGLETON based on parent availability."""
    entry = ped.get(sample_id, {})
    has_father = entry.get('father_id', '0') not in ('0', '', None)
    has_mother = entry.get('mother_id', '0') not in ('0', '', None)
    if has_father and has_mother:
        return 'TRIO'
    if has_father or has_mother:
        return 'DUO'
    return 'SINGLETON'


# ── Parent REPCN lookup ───────────────────────────────────────────────────────

def load_parent_repcn(vcf_dir: str, parent_id: str) -> dict:
    """
    Load REPCN values for a parent sample from their Reannotated VCF.
    Returns dict: {locus_id: (rep1, rep2)} - (None, None) for no-calls.
    Returns {} if parent VCF not found.
    """
    if not parent_id or parent_id == '0':
        return {}

    vcf_dir_path = Path(vcf_dir)
    # Try common naming patterns
    candidates = [
        vcf_dir_path / f"{parent_id}.Reannotated.vcf",
        vcf_dir_path / f"{parent_id}_Reannotated.vcf",
        vcf_dir_path / f"{parent_id}.stripy.vcf",
    ]
    vcf_path = None
    for c in candidates:
        if c.exists():
            vcf_path = c
            break

    if vcf_path is None:
        return {}

    repcn_map = {}
    for record in read_vcf(str(vcf_path)):
        locus_id = record['info'].get('REPID') or record['info'].get('VARID')
        if not locus_id:
            continue
        repcn    = record['fmt'].get('REPCN', '.')
        rep1, rep2 = parse_repcn(repcn)
        repcn_map[locus_id] = (rep1, rep2, repcn)
    return repcn_map


# ── DeNovo logic ──────────────────────────────────────────────────────────────

def assess_allele_denovo(
    allele_val:   int,
    direction:    str,
    nr_min:       int,
    nr_max:       int,
    father_repcn: dict,
    mother_repcn: dict,
    locus_id:     str,
    family_type:  str,
    ped_entry:    dict,
) -> str:
    """
    Assess DeNovo status for a SINGLE allele value.

    Returns: Yes / No / Possible / .
    DeNovo_Status (shared parent no-call flag) handled separately.

    Logic:
      - Check whether this allele value appears in parent(s)
      - If allele is outside normal range in proband but within range
        in both parents → Yes (de novo)
    """
    if family_type == 'SINGLETON':
        return '.'

    father_id = ped_entry.get('father_id', '0')
    mother_id = ped_entry.get('mother_id', '0')

    def allele_in_parent(parent_repcn_map, parent_id):
        """
        Check if this specific allele value exists in parent's REPCN.
        Returns:
          'same_direction' → parent has allele outside range in same direction
          'within'         → parent allele within normal range
          'no_call'        → parent REPCN = ./.
          'unavailable'    → parent VCF not found
        """
        if not parent_id or parent_id == '0':
            return 'unavailable'
        if not parent_repcn_map:
            return 'unavailable'
        entry = parent_repcn_map.get(locus_id)
        if entry is None:
            return 'unavailable'
        p_rep1, p_rep2, _ = entry
        if p_rep1 is None:
            return 'no_call'
        try:
            p_alleles = [int(p_rep1)]
            if p_rep2 is not None:
                p_alleles.append(int(p_rep2))
            for pa in p_alleles:
                if direction == 'HIGH' and pa > nr_max:
                    return 'same_direction'
                if direction == 'LOW' and pa < nr_min:
                    return 'same_direction'
            return 'within'
        except ValueError:
            return 'no_call'

    father_status = allele_in_parent(father_repcn, father_id)
    mother_status = allele_in_parent(mother_repcn, mother_id)

    # ── TRIO ──────────────────────────────────────────────────────────────
    if family_type == 'TRIO':
        if father_status == 'no_call' or mother_status == 'no_call':
            return '.'   # Undetermined - handled by shared DeNovo_Status
        if father_status == 'same_direction' or            mother_status == 'same_direction':
            return 'No'
        if father_status == 'within' and mother_status == 'within':
            return 'Yes'
        return '.'

    # ── DUO ───────────────────────────────────────────────────────────────
    if family_type == 'DUO':
        available = father_status if father_id != '0' else mother_status
        if available == 'no_call':
            return '.'
        if available == 'same_direction':
            return 'No'
        if available == 'within':
            return 'Possible'
        return '.'

    return '.'


def assess_denovo_status(
    proband_rep1: str,
    proband_rep2: str,
    nr_min:       int,
    nr_max:       int,
    father_repcn: dict,
    mother_repcn: dict,
    locus_id:     str,
) -> str:
    """
    Assess shared DeNovo_Status - returns Undetermined_ParentNoCall
    if any relevant parent has a no-call at this locus.
    Only checked when proband has at least one outlier allele.
    """
    father_id = '0'
    mother_id = '0'

    def is_no_call(parent_repcn_map, parent_id):
        if not parent_id or parent_id == '0':
            return False
        if not parent_repcn_map:
            return False
        entry = parent_repcn_map.get(locus_id)
        if entry is None:
            return False
        p_rep1, _, _ = entry
        return p_rep1 is None

    if is_no_call(father_repcn, father_id) or        is_no_call(mother_repcn, mother_id):
        return 'Undetermined_ParentNoCall'
    return '.'


def assess_denovo(
    proband_rep1: str,
    proband_rep2: str,
    normal_min,
    normal_max,
    father_repcn: dict,
    mother_repcn: dict,
    locus_id: str,
    family_type: str,
    ped_entry: dict,
) -> tuple:
    """
    Per-allele DeNovo assessment.

    Returns:
        (allele1_denovo, allele1_denovo_dir,
         allele2_denovo, allele2_denovo_dir,
         denovo_status)

    Each allele assessed independently:
      IF allele > NormalRange.Max → direction = HIGH → assess
      IF allele < NormalRange.Min → direction = LOW  → assess
      IF within range             → direction = '.'  → skip
    """
    dot5 = ('.', '.', '.', '.', '.')

    if proband_rep1 is None:
        return dot5
    if normal_min in ('.', None) or normal_max in ('.', None):
        return dot5

    try:
        nr_min = int(normal_min)
        nr_max = int(normal_max)
    except (ValueError, TypeError):
        return dot5

    def allele_direction(allele_str):
        if allele_str is None:
            return None, None
        try:
            a = int(allele_str)
            if a > nr_max:
                return a, 'HIGH'
            if a < nr_min:
                return a, 'LOW'
            return a, None   # within range - not a candidate
        except (ValueError, TypeError):
            return None, None

    a1_val, a1_dir = allele_direction(proband_rep1)
    a2_val, a2_dir = allele_direction(proband_rep2)

    # Neither allele is outside normal range
    if a1_dir is None and a2_dir is None:
        return dot5

    # Assess each allele independently
    a1_denovo = '.'
    a2_denovo = '.'

    if a1_dir is not None:
        a1_denovo = assess_allele_denovo(
            allele_val   = a1_val,
            direction    = a1_dir,
            nr_min       = nr_min,
            nr_max       = nr_max,
            father_repcn = father_repcn,
            mother_repcn = mother_repcn,
            locus_id     = locus_id,
            family_type  = family_type,
            ped_entry    = ped_entry,
        )

    if a2_dir is not None:
        a2_denovo = assess_allele_denovo(
            allele_val   = a2_val,
            direction    = a2_dir,
            nr_min       = nr_min,
            nr_max       = nr_max,
            father_repcn = father_repcn,
            mother_repcn = mother_repcn,
            locus_id     = locus_id,
            family_type  = family_type,
            ped_entry    = ped_entry,
        )

    # Shared DeNovo_Status - parent no-call flag
    # Only relevant if at least one allele is a candidate
    father_id = ped_entry.get('father_id', '0')
    mother_id = ped_entry.get('mother_id', '0')

    def is_no_call(parent_map, pid):
        if not pid or pid == '0' or not parent_map:
            return False
        entry = parent_map.get(locus_id)
        if entry is None:
            return False
        return entry[0] is None

    denovo_status = '.'
    if family_type in ('TRIO', 'DUO'):
        if is_no_call(father_repcn, father_id) or            is_no_call(mother_repcn, mother_id):
            denovo_status = 'Undetermined_ParentNoCall'

    return (
        a1_denovo,
        a1_dir if a1_dir else '.',
        a2_denovo,
        a2_dir if a2_dir else '.',
        denovo_status,
    )


# ── Core row builder ──────────────────────────────────────────────────────────

def build_rows(
    sample_id:    str,
    record:       dict,
    locus_ref:    dict,
    hpo_lookup:   dict,
    filter_terms: set,
    ped:          dict,
    father_repcn: dict,
    mother_repcn: dict,
) -> list:
    """Build output rows for one VCF record - one row per linked disease."""

    chrom = record['chrom']
    pos   = record['pos']
    ref   = record['ref']
    alt   = record['alt']
    filt  = record['filter']
    info  = record['info']
    frmt  = record['fmt']

    gt    = frmt.get('GT',    '.')
    repcn = frmt.get('REPCN', '.')
    lc    = frmt.get('LC',    '.')
    rl    = info.get('RL',    '.')    # RepeatLength in bp
    ru    = info.get('RU',    '.')    # Motif

    locus_id   = info.get('REPID') or info.get('VARID', '.')
    end        = info.get('END',  '.')
    loc_coords = f"{chrom}:{pos}-{end}"

    rep1, rep2  = parse_repcn(repcn)
    call_status = "CALLED" if rep1 is not None else "NO_CALL"

    meta    = get_locus_meta(locus_id, locus_ref)
    meta_ok = bool(meta) and '_error' not in meta

    loc_region = meta.get('LocationRegion', '.') if meta_ok else '.'
    rep_type   = meta.get('RepeatType',     '.') if meta_ok else '.'
    gene       = meta.get('Gene',           '.') if meta_ok else '.'

    disease_entries = parse_multi_disease(info)
    ped_entry       = ped.get(sample_id, {})
    family_type     = get_family_type(sample_id, ped)

    rows = []

    for dis in disease_entries:
        dis_id   = dis['dis_id']
        dis_name = dis['dis_name']
        disinher = dis['disinher']

        normal_min = normal_max = '.'
        inter_min  = inter_max  = '.'
        path_cutoff = disease_omim = onset = '.'

        if meta_ok and dis_id and 'Diseases' in meta:
            d = meta['Diseases'].get(dis_id, {})

            # Extract PathogenicCutoff first - needed for IR.Max derivation below
            path_cutoff  = d.get('PathogenicCutoff', '.')
            disease_omim = d.get('DiseaseOMIM',      '.')
            onset        = d.get('Onset',            '.')

            # NormalRange - Max=0 with Min>0 is a Stripy DB quirk → keep as '.'
            # Cannot safely derive NormalRange.Max from PathogenicCutoff alone
            nr = d.get('NormalRange') or {}
            nr_min = nr.get('Min', '.')
            nr_max = nr.get('Max', '.')
            if nr_min not in ('.', None) and nr_max not in ('.', None):
                if int(nr_max) < int(nr_min):
                    nr_max = '.'   # cannot derive safely - set to missing
            normal_min = nr_min
            normal_max = nr_max

            # IntermediateRange - Max < Min is a Stripy DB quirk
            # Fix:
            #   IF PathogenicCutoff available → IR.Max = PathogenicCutoff - 1
            #   IF PathogenicCutoff NOT available → keep IR.Max = 0 (flag for analyst)
            ir = d.get('IntermediateRange')
            if ir:
                ir_min = ir.get('Min', '.')
                ir_max = ir.get('Max', '.')
                if ir_min not in ('.', None) and ir_max not in ('.', None):
                    if int(ir_max) < int(ir_min):
                        if path_cutoff not in ('.', None):
                            ir_max = int(path_cutoff) - 1   # derive from PC
                        else:
                            ir_max = 0   # cannot derive - keep as 0 (analyst flag)
                inter_min = ir_min
                inter_max = ir_max

        hpo, hpo_gene = get_hpo_terms(str(disease_omim), hpo_lookup)

        # /compare API - called loci only
        # Per-allele columns - explicit, filterable, no packed '/' values
        allele1_repeats = allele2_repeats = '.'
        allele1_zscore  = allele2_zscore  = '.'
        allele1_outlier = allele2_outlier = '.'
        allele1_dir     = allele2_dir     = '.'
        allele1_range   = allele2_range   = '.'
        lit_dis_name    = lit_inher       = '.'

        # Per-allele repeat counts from REPCN
        if rep1 is not None:
            allele1_repeats = rep1
            allele2_repeats = rep2 if rep2 is not None else '.'

        # Per-allele direction - REPCN vs NormalRange (not Z-score)
        # HIGH: above normal max  LOW: below normal min  '.': within or unavailable
        def _rdir(a_str, nm, nx):
            if a_str in (None, '.') or nm in ('.', None) or nx in ('.', None):
                return '.'
            try:
                a, lo, hi = int(a_str), int(nm), int(nx)
                return 'HIGH' if a > hi else ('LOW' if a < lo else '.')
            except (ValueError, TypeError):
                return '.'
        allele1_dir = _rdir(rep1, normal_min, normal_max)
        allele2_dir = _rdir(rep2, normal_min, normal_max)

        # CallQuality - flag ambiguous calls
        # Nested/Replaced loci with GT!=0/0 but REPCN=0/0 are unreliable:
        # EH detected a variant but could not count the pathogenic motif
        call_quality = 'OK'
        if call_status == 'CALLED':
            is_zero_repcn = (rep1 == '0' and (rep2 is None or rep2 == '0'))
            is_nonref_gt  = gt not in ('0/0', '0', './.', '.')
            is_nested     = rep_type in ('Nested', 'Replaced', 'Imperfect GCN')
            if is_zero_repcn and is_nonref_gt and is_nested:
                call_quality = 'AmbiguousNestedCall'
        elif call_status == 'NO_CALL':
            call_quality = 'NoCall'

        if call_status == "CALLED":
            compare = fetch_compare(locus_id, rep1, rep2)
            if '_error' not in compare:
                pop   = compare.get('Population', {})
                raw_z = pop.get('Zscore')
                raw_o = pop.get('IsOutlier')

                # Per-allele Z-score
                if isinstance(raw_z, list):
                    allele1_zscore = fmt(raw_z[0]) if len(raw_z) > 0 else '.'
                    allele2_zscore = fmt(raw_z[1]) if len(raw_z) > 1 else '.'
                elif raw_z is not None:
                    allele1_zscore = fmt(raw_z)
                    allele2_zscore = fmt(raw_z) if rep2 is not None else '.'

                # Per-allele outlier
                if isinstance(raw_o, list):
                    allele1_outlier = fmt(raw_o[0]) if len(raw_o) > 0 else '.'
                    allele2_outlier = fmt(raw_o[1]) if len(raw_o) > 1 else '.'
                elif raw_o is not None:
                    allele1_outlier = fmt(raw_o)
                    allele2_outlier = fmt(raw_o) if rep2 is not None else '.'

                # Literature - per-allele range
                lit = compare.get('Literature', {})
                lit_entry = {}
                if dis_id and dis_id in lit:
                    lit_entry = lit[dis_id]
                elif lit:
                    lit_entry = next(iter(lit.values()))

                if lit_entry:
                    lit_dis_name = lit_entry.get('DiseaseName', '.')
                    lit_inher    = lit_entry.get('Inheritance',  '.')
                    raw_range    = lit_entry.get('Range',        '.')
                    if isinstance(raw_range, list):
                        allele1_range = str(raw_range[0]) if len(raw_range) > 0 else '.'
                        allele2_range = str(raw_range[1]) if len(raw_range) > 1 else '.'
                    else:
                        allele1_range = str(raw_range)
                        allele2_range = str(raw_range) if rep2 is not None else '.'
            else:
                print(f"    [WARN] /compare error {locus_id}: "
                      f"{compare['_error']}", file=sys.stderr)

        # DeNovo assessment - per allele independently
        a1_denovo, a1_denovo_dir, \
        a2_denovo, a2_denovo_dir, \
        denovo_status = assess_denovo(
            proband_rep1  = rep1,
            proband_rep2  = rep2,
            normal_min    = normal_min,
            normal_max    = normal_max,
            father_repcn  = father_repcn,
            mother_repcn  = mother_repcn,
            locus_id      = locus_id,
            family_type   = family_type,
            ped_entry     = ped_entry,
        )

        # Rule 1: If DeNovo='.', direction is meaningless → clear it
        if a1_denovo == '.':
            a1_denovo_dir = '.'
        if a2_denovo == '.':
            a2_denovo_dir = '.'

        # Rule 2: Homozygous (rep1==rep2) - Allele2 is not independent
        # Suppress ALL Allele2 DeNovo fields to avoid duplicating Allele1
        if rep1 is not None and rep2 is not None and str(rep1) == str(rep2):
            a2_denovo     = '.'
            a2_denovo_dir = '.'

        rows.append({
            "SampleID":                sample_id,
            "Chr":                     chrom,
            "Start":                   pos,
            "End":                     end,
            "LocationCoordinates":     loc_coords,
            "Ref":                     ref,
            "Alt":                     alt,
            "GT":                      gt,
            "LocusCov":                lc,
            "RepeatLength":            rl,
            "FilterStatus":            filt,
            "CallStatus":              call_status,
            "CallQuality":             call_quality,
            "Locus":                   locus_id,
            "LocationRegion":          loc_region,
            "RepeatType":              rep_type,
            "Motif":                   ru,
            "Normal Range (Min)":      fmt(normal_min),
            "Normal Range (Max)":      fmt(normal_max),
            "IntermediateRange (Min)": fmt(inter_min),
            "IntermediateRange (Max)": fmt(inter_max),
            "PathogenicCutoff":        fmt(path_cutoff),
            "Disease Name":            dis_name   or '.',
            "Disease OMIM":            fmt(disease_omim),
            "Inheritance":             disinher   or '.',
            "Onset":                   onset,
            "Stripy_Gene":             gene,
            "HPO_Gene":                hpo_gene,
            "HPO Phenotype":           hpo,
            "HPO_Filtered":            get_hpo_filtered(hpo, filter_terms),
            "Allele1_Repeats":         allele1_repeats,
            "Allele2_Repeats":         allele2_repeats,
            "Literature:DiseaseName":  lit_dis_name,
            "Literature:Inheritance":  lit_inher,
            "Allele1_Range":           allele1_range,
            "Allele2_Range":           allele2_range,
            "Allele1_Zscore":          allele1_zscore,
            "Allele2_Zscore":          allele2_zscore,
            "Allele1_Outlier":         allele1_outlier,
            "Allele2_Outlier":         allele2_outlier,
            "Allele1_Direction":       allele1_dir,
            "Allele2_Direction":       allele2_dir,
            "Allele1_DeNovo":          a1_denovo,
            "Allele1_DeNovo_Direction": a1_denovo_dir,
            "Allele2_DeNovo":          a2_denovo,
            "Allele2_DeNovo_Direction": a2_denovo_dir,
            "DeNovo_Status":           denovo_status,
            "Caller":                  "ExpansionHunter",
        })

    return rows


# ── TSV writer ────────────────────────────────────────────────────────────────

def write_tsv(rows: list, out_path: Path):
    with open(out_path, 'w') as f:
        f.write('\t'.join(OUTPUT_COLUMNS) + '\n')
        for row in rows:
            f.write('\t'.join(
                str(row.get(col, '.')) for col in OUTPUT_COLUMNS
            ) + '\n')


# ── Summary ───────────────────────────────────────────────────────────────────

def print_summary(sample_id: str, rows: list, family_type: str):
    total    = len(rows)
    called   = sum(1 for r in rows if r['CallStatus'] == 'CALLED')
    no_call  = total - called
    # Genuine expansion outliers - either allele HIGH
    high_out = [r for r in rows
                if ('True' in str(r.get('Allele1_Outlier','')) or
                    'True' in str(r.get('Allele2_Outlier','')))
                and ('HIGH' in str(r.get('Allele1_Direction','')) or
                     'HIGH' in str(r.get('Allele2_Direction','')))
                and r['CallStatus'] == 'CALLED']
    denovo   = [r for r in rows if r.get('Allele1_DeNovo') in ('Yes','Possible') or r.get('Allele2_DeNovo') in ('Yes','Possible')]

    print(f"\n{'='*65}")
    print(f"  SAMPLE              : {sample_id}")
    print(f"  Family type         : {family_type}")
    print(f"  Total rows          : {total}")
    print(f"  Called loci rows    : {called}")
    print(f"  No-call rows        : {no_call}")
    print(f"  Expansion outliers  : {len(high_out)}")
    print(f"  DeNovo candidates   : {len(denovo)}")
    if high_out:
        print()
        for r in high_out:
            print(f"    [OUTLIER] {r['Locus']:<20} "
                  f"A1={r['Allele1_Repeats']} A2={r['Allele2_Repeats']} "
                  f"Z1={r['Allele1_Zscore']} Z2={r['Allele2_Zscore']}")
    if denovo:
        print()
        for r in denovo:
            a1 = r.get('Allele1_DeNovo','.')
            a2 = r.get('Allele2_DeNovo','.')
            label = f'A1:{a1}' if a1 not in ('.','') else ''
            label += f' A2:{a2}' if a2 not in ('.','') else ''
            print(f"    [DENOVO] {r['Locus']:<20} "
                  f"REPCN={r['GT']:<10} {label}")
    print(f"{'='*65}\n")


# ── CLI ───────────────────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(
        description="Build per-sample STRipy post-processing TSV report"
    )
    p.add_argument('--vcf',        required=True,
                   help="STRipy-annotated VCF (*.Reannotated.vcf)")
    p.add_argument('--sample-id',  required=True,
                   help="Sample ID e.g. LH0014A")
    p.add_argument('--locus-ref',  required=True,
                   help="Stripy locus reference JSON")
    p.add_argument('--hpo-file',   required=True,
                   help="HPO genes_to_phenotype.txt")
    p.add_argument('--hpo-filter-file', default=None,
                   help="Project HPO filter terms file (optional). "
                        "E.g. hpo_epilepsy_terms.txt. "
                        "If omitted HPO_Filtered column will be '.'")
    p.add_argument('--ped',        required=True,
                   help="PED file for family relationships")
    p.add_argument('--vcf-dir',    required=True,
                   help="Directory containing all Reannotated VCFs "
                        "(for parent REPCN lookup)")
    p.add_argument('--outdir',     default='.',
                   help="Output directory")
    return p.parse_args()


def main():
    args   = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    print(f"[INFO] Sample     : {args.sample_id}")
    print(f"[INFO] VCF        : {args.vcf}")
    print(f"[INFO] Locus ref  : {args.locus_ref}")
    print(f"[INFO] HPO file   : {args.hpo_file}")
    print(f"[INFO] HPO filter : "
          f"{getattr(args,'hpo_filter_file',None) or 'disabled'}")
    print(f"[INFO] PED        : {args.ped}")
    print(f"[INFO] VCF dir    : {args.vcf_dir}")
    print(f"[INFO] Output dir : {outdir}")

    locus_ref    = load_locus_reference(args.locus_ref)
    hpo_lookup   = load_hpo_file(args.hpo_file)
    filter_terms = load_hpo_filter_file(
        getattr(args, 'hpo_filter_file', None))
    ped          = load_ped(args.ped)

    print(f"[INFO] Locus ref  : {len(locus_ref)} loci")
    print(f"[INFO] HPO lookup : {len(hpo_lookup)} gene-disease pairs")
    print(f"[INFO] HPO filter : {len(filter_terms)} terms "
          f"({'active' if filter_terms else 'disabled - HPO_Filtered will be .'})")
    print(f"[INFO] PED        : {len(ped)} samples")

    # Load parent VCFs for DeNovo
    ped_entry   = ped.get(args.sample_id, {})
    family_type = get_family_type(args.sample_id, ped)
    father_id   = ped_entry.get('father_id', '0')
    mother_id   = ped_entry.get('mother_id', '0')

    print(f"[INFO] Family type: {family_type}  "
          f"(Father: {father_id}  Mother: {mother_id})")

    father_repcn = load_parent_repcn(args.vcf_dir, father_id)
    mother_repcn = load_parent_repcn(args.vcf_dir, mother_id)

    if father_repcn:
        print(f"[INFO] Father VCF : loaded {len(father_repcn)} loci")
    if mother_repcn:
        print(f"[INFO] Mother VCF : loaded {len(mother_repcn)} loci")
    print()

    all_rows = []
    for i, record in enumerate(read_vcf(args.vcf), 1):
        locus_id = (record['info'].get('REPID') or
                    record['info'].get('VARID', f'LOCUS_{i}'))
        repcn    = record['fmt'].get('REPCN', '.')
        rep1, _  = parse_repcn(repcn)
        status   = "CALLED" if rep1 is not None else "NO_CALL"

        print(f"  [{i:>3}] {locus_id:<25} REPCN={repcn:<12} {status}",
              end='', flush=True)

        rows = build_rows(
            sample_id    = args.sample_id,
            record       = record,
            locus_ref    = locus_ref,
            hpo_lookup   = hpo_lookup,
            filter_terms = filter_terms,
            ped          = ped,
            father_repcn = father_repcn,
            mother_repcn = mother_repcn,
        )
        all_rows.extend(rows)
        print(f" → {len(rows)} row(s)")

    out_tsv = outdir / f"{args.sample_id}.stripy_report.tsv"
    write_tsv(all_rows, out_tsv)
    print(f"\n[INFO] Written : {out_tsv}")
    print_summary(args.sample_id, all_rows, family_type)


if __name__ == '__main__':
    main()
