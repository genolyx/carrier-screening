import os
import re
import csv
import gzip
import shutil
import tempfile
import sqlite3
import subprocess
from pathlib import Path
from typing import List, Dict, Any, Set, Optional, Tuple
from collections import defaultdict, deque

import pysam
from intervaltree import IntervalTree

from fastapi import FastAPI, Request, UploadFile, File, Form
from fastapi.responses import HTMLResponse, JSONResponse
from fastapi.templating import Jinja2Templates


# ============================================
# DEBUG CONFIG + LOGGER
# ============================================

DEBUG = os.environ.get("PORTAL_DEBUG", "0").strip().lower() in ("1", "true", "yes", "y")
SHOW_DEBUG_IN_RESULTS = os.environ.get("SHOW_DEBUG_IN_RESULTS", "0").strip().lower() in ("1", "true", "yes", "y")
DEBUG_MAX_LINES = int(os.environ.get("DEBUG_MAX_LINES", "200"))

_DEBUG_BUF = deque(maxlen=DEBUG_MAX_LINES)

def dlog(*args):
    """Debug logger: prints + stores messages for optional on-page display."""
    if not DEBUG:
        return
    msg = "[DEBUG] " + " ".join(str(a) for a in args)
    _DEBUG_BUF.append(msg)
    print(msg, flush=True)

def get_debug_dump() -> str:
    if not DEBUG:
        return ""
    return "\n".join(_DEBUG_BUF)


# ============================================
# CONFIGURATION
# ============================================

HPO_GENE_FILE = "genes_to_phenotype.txt"
GENE_COORDS_FILE = "genes_hg38.bed"
MANE_REFSEQ_GFF_FILE = "MANE.GRCh38.v1.3.ensembl_genomic.gff"

# snpEff
SNPEFF_JAR = "snpEff/snpEff.jar"
SNPEFF_DB = "GRCh38.86"

# Local curated variant DB (optional stub)
DB_PATH = "variants.db"

# ClinVar (override on Linux with env vars!)
DEFAULT_CLINVAR_DIR = "clinvar"
CLINVAR_DIR = os.environ.get("CLINVAR_DIR", DEFAULT_CLINVAR_DIR)
CLINVAR_VCF_GZ = os.environ.get("CLINVAR_VCF_GZ", "")
if not CLINVAR_VCF_GZ:
    CLINVAR_VCF_GZ = str(Path(CLINVAR_DIR) / "clinvar.vcf.gz")
if CLINVAR_VCF_GZ and os.path.isdir(CLINVAR_VCF_GZ):
    CLINVAR_VCF_GZ = str(Path(CLINVAR_VCF_GZ) / "clinvar.vcf.gz")

# ClinGen dosage TSV (override on Linux with env vars!)
DEFAULT_CLINGEN_DIR = "clingen"
CLINGEN_DIR = os.environ.get("CLINGEN_DIR", DEFAULT_CLINGEN_DIR)
CLINGEN_GENE_TSV = os.environ.get(
    "CLINGEN_GENE_TSV",
    str(Path(CLINGEN_DIR) / "ClinGen_gene_curation_list_GRCh38.tsv")
)

app = FastAPI()
templates = Jinja2Templates(directory="templates")


# ============================================
# SQLITE DB (OPTIONAL CURATED VARIANTS)
# ============================================

def get_db_connection():
    conn = sqlite3.connect(DB_PATH)
    conn.row_factory = sqlite3.Row
    return conn

def init_db():
    conn = get_db_connection()
    cur = conn.cursor()
    cur.execute(
        """
        CREATE TABLE IF NOT EXISTS curated_variants (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            gene TEXT NOT NULL,
            hgvsc TEXT,
            hgvsp TEXT,
            classification TEXT NOT NULL,
            source TEXT,
            notes TEXT
        )
        """
    )
    conn.commit()
    conn.close()

def lookup_curated_variant(gene: str, hgvsc: str, hgvsp: str):
    if not gene:
        return None
    conn = get_db_connection()
    cur = conn.cursor()
    cur.execute(
        """
        SELECT classification, source, notes
        FROM curated_variants
        WHERE gene = ?
          AND (
              (hgvsc IS NOT NULL AND hgvsc <> '' AND hgvsc = ?)
              OR
              (hgvsp IS NOT NULL AND hgvsp <> '' AND hgvsp = ?)
          )
        LIMIT 1
        """,
        (gene, hgvsc or "", hgvsp or "")
    )
    row = cur.fetchone()
    conn.close()
    if row:
        return {"classification": row["classification"], "source": row["source"], "notes": row["notes"]}
    return None

init_db()


# ============================================
# HPO GENE MAP
# ============================================

def load_hpo_gene_map(path: str = HPO_GENE_FILE):
    if not os.path.exists(path):
        raise RuntimeError(f"HPO gene file '{path}' not found.")

    gene_map: Dict[str, Set[str]] = {}
    name_map: Dict[str, str] = {}

    with open(path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue

            gene = parts[1].strip()
            hpo_id = parts[2].strip()
            hpo_name = parts[3].strip()

            gene_map.setdefault(hpo_id, set()).add(gene)
            name_map.setdefault(hpo_id, hpo_name)

    return gene_map, name_map

HPO_GENE_MAP, HPO_NAME_MAP = load_hpo_gene_map()
ALL_GENE_SYMBOLS: List[str] = sorted({g for genes in HPO_GENE_MAP.values() for g in genes})

def normalize_hpo_terms(raw: str) -> List[str]:
    if not raw:
        return []
    tokens = [t.strip() for t in raw.replace("\n", ",").split(",") if t.strip()]
    ids: List[str] = []
    pat = re.compile(r"(HP:\d{7})")
    for t in tokens:
        m = pat.search(t)
        if m:
            ids.append(m.group(1))
        elif t.startswith("HP:"):
            ids.append(t.split()[0])
    return ids


# ============================================
# GENE INTERVALS (BED)
# ============================================

def load_gene_intervals(path: str = GENE_COORDS_FILE):
    if not os.path.exists(path):
        raise RuntimeError(f"Gene BED '{path}' not found.")

    trees: Dict[str, IntervalTree] = defaultdict(IntervalTree)
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            chrom, start, end, gene = parts[0], int(parts[1]), int(parts[2]), parts[3].strip()
            trees[chrom].addi(start, end + 1, gene)
    return trees

GENE_INTERVALS = load_gene_intervals()

def positional_gene_lookup(chrom: str, pos: int) -> Optional[str]:
    tree = GENE_INTERVALS.get(chrom)
    if not tree:
        return None
    hits = tree[pos]
    if hits:
        return next(iter(hits)).data
    return None


# ============================================
# MANE (ENST + NM) MAP
# ============================================

def load_mane_refseq_map(path: str = MANE_REFSEQ_GFF_FILE):
    mane: Dict[str, Dict[str, str]] = {}
    if not os.path.exists(path):
        print(f"WARNING: MANE file '{path}' missing.")
        return mane

    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            if parts[2] != "transcript":
                continue

            attrs: Dict[str, str] = {}
            for kv in parts[8].split(";"):
                if "=" in kv:
                    k, v = kv.split("=", 1)
                    attrs[k.strip()] = v.strip()

            gene = attrs.get("gene_name")
            enst = attrs.get("transcript_id")

            nm = ""
            dbx = attrs.get("Dbxref", "")
            if dbx:
                for item in dbx.split(","):
                    item = item.strip()
                    if item.startswith("RefSeq:NM_"):
                        nm = item.replace("RefSeq:", "")
                        break

            if not gene:
                continue
            entry = mane.setdefault(gene, {})
            if enst and "enst" not in entry:
                entry["enst"] = enst
            if nm and "nm" not in entry:
                entry["nm"] = nm

    return mane

MANE_MAP = load_mane_refseq_map()


# ============================================
# CLINGEN DOSAGE TSV (ROBUST PARSER)
# ============================================

def load_clingen_gene_curation(path: str) -> Dict[str, Dict[str, Any]]:
    if not path or not os.path.exists(path):
        print(f"WARNING: ClinGen TSV not found: {path}")
        return {}

    def norm(x: str) -> str:
        return (x or "").strip()

    def norm_gene(x: str) -> str:
        return norm(x).upper()

    def to_int(x):
        try:
            return int(str(x).strip())
        except Exception:
            return None

    out: Dict[str, Dict[str, Any]] = {}

    with open(path, "r", encoding="utf-8", errors="replace") as f:
        header_cols = None

        for line in f:
            if not line.strip():
                continue
            if line.startswith("#Gene Symbol\t"):
                header_cols = line.lstrip("#").rstrip("\n").split("\t")
                break

        if not header_cols:
            print("WARNING: ClinGen header not found (expected '#Gene Symbol\\t...').")
            return {}

        reader = csv.DictReader(f, delimiter="\t", fieldnames=header_cols)

        for row in reader:
            gene = norm_gene(row.get("Gene Symbol") or row.get("geneSymbol") or row.get("Gene") or "")
            if not gene:
                continue

            hi_score = to_int(row.get("Haploinsufficiency Score") or row.get("HI Score") or row.get("haploScore"))
            ts_score = to_int(row.get("Triplosensitivity Score") or row.get("TS Score") or row.get("triploScore"))

            candidate = {
                "hgnc_id": norm(row.get("HGNC ID") or row.get("HGNCID") or ""),
                "ensembl_gene_id": norm(row.get("Ensembl Gene ID") or ""),
                "omim_gene_id": norm(row.get("OMIM Gene ID") or ""),

                "hi_score": hi_score,
                "ts_score": ts_score,
                "hi_desc": norm(row.get("Haploinsufficiency Description") or row.get("haploDescription") or ""),
                "ts_desc": norm(row.get("Triplosensitivity Description") or row.get("triploDescription") or ""),
                "last_eval": norm(row.get("Date Last Evaluated") or row.get("dateLastEvaluated") or ""),

                "hi_disease_id": norm(row.get("Haploinsufficiency Disease ID") or ""),
                "ts_disease_id": norm(row.get("Triplosensitivity Disease ID") or ""),

                "url": norm(row.get("url") or row.get("URL") or ""),
            }

            existing = out.get(gene)
            if existing is None:
                out[gene] = candidate
            else:
                old_rank = (existing.get("hi_score") is not None, existing.get("ts_score") is not None)
                new_rank = (candidate.get("hi_score") is not None, candidate.get("ts_score") is not None)
                if new_rank > old_rank:
                    out[gene] = candidate

    print("ClinGen loaded genes:", len(out))
    return out

CLINGEN_MAP = load_clingen_gene_curation(CLINGEN_GENE_TSV)


# ============================================
# VCF FORMAT CLEANER (FIXES '.' IN GQ/SB ETC)
# ============================================

TROUBLE_FORMATS = {"GQ", "SB"}

def clean_vcf_remove_formats(input_vcf: str, output_vcf: str):
    opener_in = gzip.open if input_vcf.endswith(".gz") else open
    with opener_in(input_vcf, "rt", encoding="utf-8", errors="replace") as fin, \
            open(output_vcf, "w", encoding="utf-8") as fout:
        for line in fin:
            line = line.rstrip("\n")

            if line.startswith("##FORMAT=<ID="):
                if any(f"ID={fid}" in line for fid in TROUBLE_FORMATS):
                    continue
                fout.write(line + "\n")
                continue

            if line.startswith("#"):
                fout.write(line + "\n")
                continue

            parts = line.split("\t")
            if len(parts) <= 8:
                fout.write(line + "\n")
                continue

            fmt_fields = parts[8].split(":")
            trouble_indices = [i for i, name in enumerate(fmt_fields) if name in TROUBLE_FORMATS]
            if not trouble_indices:
                fout.write(line + "\n")
                continue

            keep_indices = [i for i in range(len(fmt_fields)) if i not in trouble_indices]
            parts[8] = ":".join(fmt_fields[i] for i in keep_indices) if keep_indices else "."

            for i in range(9, len(parts)):
                sample = parts[i]
                if sample in (".", ""):
                    continue
                sf = sample.split(":")
                if len(sf) != len(fmt_fields):
                    continue
                parts[i] = ":".join(sf[j] for j in keep_indices) if keep_indices else "."

            fout.write("\t".join(parts) + "\n")


# ============================================
# snpEff ANNOTATION
# ============================================

def run_snpeff(input_vcf: str, output_vcf: str):
    if not os.path.exists(SNPEFF_JAR):
        raise RuntimeError(f"snpEff jar not found at '{SNPEFF_JAR}'")
    cmd = ["java", "-Xmx4g", "-jar", SNPEFF_JAR, SNPEFF_DB, "-hgvs", input_vcf]
    dlog("Running snpEff:", " ".join(cmd))
    with open(output_vcf, "w", encoding="utf-8") as out:
        res = subprocess.run(cmd, stdout=out, stderr=subprocess.PIPE, text=True)
    if res.returncode != 0:
        raise RuntimeError(f"snpEff error (exit {res.returncode}):\n{res.stderr}")


# ============================================
# HPO AUTOCOMPLETE API
# ============================================

@app.get("/api/hpo_search")
async def hpo_search(q: str = ""):
    q = q.lower().strip()
    if not q:
        return []
    out = []
    for hpo_id, name in HPO_NAME_MAP.items():
        if q in hpo_id.lower() or q in name.lower():
            out.append({"id": hpo_id, "name": name})
            if len(out) >= 20:
                break
    return out


# ============================================
# gnomAD LOCAL LOOKUP (THIS IS WHAT WE USE FOR AF)
# ============================================

DEFAULT_GNOMAD_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "gnomad"))
GNOMAD_DIR = os.environ.get("GNOMAD_DIR", DEFAULT_GNOMAD_DIR)

# matches: gnomad.genomes.v3.1.2.sites.chr1.vcf.bgz
# matches: gnomad.genomes.v3.1.2.sites.chr1.vcf.bgz OR gnomad.genomes.v3.1.2.sites.vcf.bgz
GNOMAD_GENOMES_GLOB = os.environ.get("GNOMAD_GENOMES_GLOB", "gnomad.genomes.v*.sites*.bgz")
GNOMAD_EXOMES_GLOB  = os.environ.get("GNOMAD_EXOMES_GLOB",  "gnomad.exomes.v*.sites*.bgz")

_GNOMAD_VCF_CACHE: Dict[str, pysam.VariantFile] = {}
_GNOMAD_HAS_CHR: Dict[str, bool] = {}
_GNOMAD_PICK_CACHE: Dict[Tuple[str, str], Optional[str]] = {}

# gnomAD sites VCF AF keys to try (varies by file/version)
_GNOMAD_AF_KEYS = ["AF", "AF_total", "AF_joint", "AF_popmax", "AF_POPMAX"]

def _first_numeric(val) -> Optional[float]:
    if val is None:
        return None
    if isinstance(val, (list, tuple)):
        val = val[0] if val else None
    if val is None:
        return None
    try:
        return float(val)
    except Exception:
        try:
            return float(str(val))
        except Exception:
            return None

def _detect_vcf_has_chr(vcf_path: str) -> bool:
    try:
        vf = pysam.VariantFile(vcf_path)
        contigs = list(vf.header.contigs)
        vf.close()
        has_chr = any(str(c).startswith("chr") for c in contigs)
        dlog("gnomAD contigs chr-prefix?", Path(vcf_path).name, "=>", has_chr)
        return has_chr
    except Exception as e:
        dlog("Could not inspect gnomAD contigs; assuming chr:", Path(vcf_path).name, "err:", e)
        return True

def _norm_chrom_for_vcf(chrom: str, vcf_has_chr: bool) -> str:
    chrom = str(chrom)
    if vcf_has_chr:
        return chrom if chrom.startswith("chr") else "chr" + chrom
    return chrom[3:] if chrom.startswith("chr") else chrom

def _open_gnomad_vcf(vcf_path: str) -> pysam.VariantFile:
    vf = _GNOMAD_VCF_CACHE.get(vcf_path)
    if vf is None:
        dlog("Opening gnomAD VCF:", vcf_path)
        vf = pysam.VariantFile(vcf_path)
        _GNOMAD_VCF_CACHE[vcf_path] = vf
        _GNOMAD_HAS_CHR[vcf_path] = _detect_vcf_has_chr(vcf_path)
    return vf

def _pick_gnomad_file(dataset: str, chrom: str) -> Optional[str]:
    """
    Pick the correct per-chromosome gnomAD file.
    Your files are like: gnomad.genomes.v3.1.2.sites.chr1.vcf.bgz
    """
    key = (dataset, str(chrom))
    if key in _GNOMAD_PICK_CACHE:
        return _GNOMAD_PICK_CACHE[key]

    if not GNOMAD_DIR or not os.path.isdir(GNOMAD_DIR):
        dlog("GNOMAD_DIR missing/not a dir:", GNOMAD_DIR)
        _GNOMAD_PICK_CACHE[key] = None
        return None

    c = str(chrom)
    c_nochr = c[3:] if c.startswith("chr") else c
    token_chr = f"chr{c_nochr}".lower()     # chr1, chr2 ...
    token_no = c_nochr.lower()              # 1, 2 ...

    glob_pat = GNOMAD_EXOMES_GLOB if dataset == "exomes" else GNOMAD_GENOMES_GLOB
    candidates = list(Path(GNOMAD_DIR).glob(glob_pat))
    dlog("gnomAD candidates", dataset, "glob:", glob_pat, "count:", len(candidates))

    best = None
    dlog("DEBUG: Looking for chrom:", chrom, "token_chr:", token_chr, "token_no:", token_no)
    
    # prefer explicit chr match
    for p in candidates:
        name = p.name.lower()
        if token_chr in name:
            dlog("DEBUG: MATCH found (token_chr):", name)
            best = str(p)
            break
        else:
             # Very verbose, but helpful for debugging once
             dlog("DEBUG: check failed:", name, "did not contain", token_chr)

    # fallback: match bare number (less safe but okay)
    if best is None:
        dlog("DEBUG: No token_chr match. Trying bare token_no:", token_no)
        for p in candidates:
            name = p.name.lower()
            if f".{token_no}." in name or f"_{token_no}_" in name:
                dlog("DEBUG: MATCH found (bare token):", name)
                best = str(p)
                break
    # fallback: if only one file, use it
    if best is None and len(candidates) == 1:
        dlog("DEBUG: Only 1 candidate found, effectively assuming it covers all:", candidates[0])
        best = str(candidates[0])

    if best:
        dlog("Picked gnomAD file:", dataset, chrom, "=>", best)
    else:
        dlog("FAILED to pick gnomAD file for", dataset, chrom, "Candidates were:", [p.name for p in candidates])

    _GNOMAD_PICK_CACHE[key] = best
    return best

def _extract_alt_specific_info(rec: pysam.VariantRecord, alt: str, key: str) -> Optional[float]:
    if key not in rec.info:
        return None
    val = rec.info[key]
    if isinstance(val, (list, tuple)):
        if not rec.alts:
            return _first_numeric(val)
        try:
            idx = list(rec.alts).index(alt)
        except ValueError:
            idx = 0
        if idx < len(val):
            return _first_numeric(val[idx])
        return _first_numeric(val[0])
    return _first_numeric(val)

def gnomad_lookup_af(chrom: str, pos: int, ref: str, alt: str) -> Dict[str, Any]:
    """
    Returns:
      {
        "af": max(exomes_af, genomes_af),
        "exomes_af": ...,
        "genomes_af": ...,
        "source": "gnomAD-local" if any found else ""
      }
    """
    out = {"af": None, "exomes_af": None, "genomes_af": None, "source": ""}

    if not GNOMAD_DIR or not os.path.isdir(GNOMAD_DIR):
        dlog("GNOMAD_DIR not set/invalid:", GNOMAD_DIR)
        return out

    def lookup_one(dataset: str) -> Optional[float]:
        vcf_path = _pick_gnomad_file(dataset, chrom)
        if not vcf_path or not os.path.exists(vcf_path):
            dlog("No gnomAD VCF for", dataset, chrom, "->", vcf_path)
            return None

        has_index = (os.path.exists(vcf_path + ".tbi") or os.path.exists(vcf_path + ".csi"))
        if not has_index:
            dlog("gnomAD file missing index (.tbi/.csi):", vcf_path)
            return None

        try:
            vf = _open_gnomad_vcf(vcf_path)
            vcf_has_chr = _GNOMAD_HAS_CHR.get(vcf_path, True)
            qchrom = _norm_chrom_for_vcf(chrom, vcf_has_chr)

            dlog("gnomAD fetch", dataset, "query:", qchrom, pos - 1, pos, "target:", chrom, pos, ref, alt)

            any_region = False
            for r in vf.fetch(qchrom, pos - 1, pos):
                any_region = True

                if r.pos != pos:
                    continue
                if r.ref != ref:
                    continue
                if not r.alts or alt not in r.alts:
                    continue

                # MATCH! now find AF key
                for k in _GNOMAD_AF_KEYS:
                    af = _extract_alt_specific_info(r, alt, k)
                    if af is not None:
                        dlog("gnomAD HIT", dataset, chrom, pos, ref, alt, "key:", k, "AF:", af)
                        return af

                dlog("gnomAD MATCH but no AF key among", _GNOMAD_AF_KEYS,
                     "available keys sample:", list(r.info.keys())[:50])

            if not any_region:
                dlog("gnomAD region returned 0 records:", qchrom, pos)
        except Exception as e:
            dlog("gnomAD lookup error", dataset, chrom, pos, ref, alt, "err:", e)
            return None

        return None

    # Try both (even if you only have genomes; exomes will just be None)
    ex_af = lookup_one("exomes")
    gn_af = lookup_one("genomes")

    out["exomes_af"] = ex_af
    out["genomes_af"] = gn_af

    vals = [x for x in (ex_af, gn_af) if isinstance(x, (int, float))]
    if vals:
        out["af"] = max(vals)
        out["source"] = "gnomAD-local"

    return out


@app.get("/api/gnomad_lookup")
async def api_gnomad_lookup(chrom: str, pos: int, ref: str, alt: str):
    """
    Quick test endpoint:
      /api/gnomad_lookup?chrom=chr1&pos=123&ref=A&alt=G
    """
    res = gnomad_lookup_af(chrom, pos, ref, alt)
    return JSONResponse({"query": {"chrom": chrom, "pos": pos, "ref": ref, "alt": alt}, "result": res, "debug": get_debug_dump()})


# ============================================
# SAMPLE METRICS (DP / AD / VAF / GT)
# ============================================

def _safe_int(x) -> Optional[int]:
    try:
        if x is None:
            return None
        if isinstance(x, bool):
            return None
        return int(x)
    except Exception:
        try:
            return int(float(str(x)))
        except Exception:
            return None

def _safe_float(x) -> Optional[float]:
    try:
        if x is None:
            return None
        if isinstance(x, bool):
            return None
        return float(x)
    except Exception:
        try:
            return float(str(x))
        except Exception:
            return None

def _alts_as_str_list(rec) -> List[str]:
    return [str(a) for a in (rec.alts or [])]

def _pick_first_sample(rec):
    if not rec.samples:
        return None, None
    name = next(iter(rec.samples.keys()))
    return name, rec.samples[name]

def _zygosity_from_gt(gt: Tuple[Optional[int], ...]) -> Tuple[str, str]:
    if not gt:
        return "", ""
    gt_str = "/".join("." if a is None else str(a) for a in gt)
    if len(gt) < 2:
        return "Other", gt_str
    a0, a1 = gt[0], gt[1]
    if a0 is None or a1 is None:
        return "Other", gt_str
    if (a0 == 0 and a1 == 1) or (a0 == 1 and a1 == 0):
        return "Het", gt_str
    if a0 == 1 and a1 == 1:
        return "Hom", gt_str
    if a0 == 0 and a1 == 0:
        return "Ref", gt_str
    return "Other", gt_str

def get_sample_metrics(rec, alt_allele: str) -> Dict[str, Any]:
    dp = None
    ref_depth = None
    alt_depth = None
    vaf = None
    gt_str = ""
    zygosity = ""
    source = ""

    alts = _alts_as_str_list(rec)
    alt_index = 0
    if alts:
        try:
            alt_index = alts.index(str(alt_allele))
        except ValueError:
            alt_index = 0

    sample_name, s = _pick_first_sample(rec)

    if s is not None:
        if "DP" in s and s["DP"] is not None:
            dp = _safe_int(s["DP"])
            if dp is not None:
                source = source or "FORMAT:DP"

        gt = s.get("GT")
        if gt:
            zygosity, gt_str = _zygosity_from_gt(gt)

        if "AD" in s and s["AD"] is not None:
            try:
                ad = list(s["AD"])
                if len(ad) >= 2:
                    rd = _safe_int(ad[0])
                    ai = 1 + alt_index
                    ad_alt = _safe_int(ad[ai]) if ai < len(ad) else None
                    if rd is not None:
                        ref_depth = rd
                    if ad_alt is not None:
                        alt_depth = ad_alt
                    if ref_depth is not None and alt_depth is not None:
                        denom = ref_depth + alt_depth
                        if denom > 0:
                            vaf = alt_depth / denom
                            source = source or "FORMAT:AD"
            except Exception:
                pass

    if dp is None and "DP" in rec.info and rec.info["DP"] is not None:
        val = rec.info["DP"]
        if isinstance(val, (list, tuple)):
            val = val[0] if val else None
        dp = _safe_int(val)
        if dp is not None:
            source = source or "INFO:DP"

    if vaf is None and ref_depth is not None and alt_depth is not None:
        denom = ref_depth + alt_depth
        if denom > 0:
            vaf = alt_depth / denom

    if dp is None and ref_depth is not None and alt_depth is not None:
        dp = ref_depth + alt_depth

    if zygosity == "Ref":
        zygosity = "Other"

    return {
        "sample_name": sample_name or "",
        "dp": dp,
        "ref_depth": ref_depth,
        "alt_depth": alt_depth,
        "vaf": vaf,
        "gt": gt_str,
        "zygosity": zygosity,
        "metrics_source": source,
    }
# ============================================
# dbSNP URL helper
# ============================================

def make_dbsnp_url(rsid: str) -> str:
    if rsid is None:
        return ""

    s = str(rsid).strip()
    if not s or s in (".", "0"):
        return ""

    # If it contains something like rs3007429 anywhere, use that number
    m = re.search(r"\brs(\d+)\b", s, flags=re.IGNORECASE)
    if m:
        return f"https://www.ncbi.nlm.nih.gov/snp/?term={m.group(1)}"

    # Otherwise if it contains a bare number anywhere, use the first one
    m2 = re.search(r"\b(\d+)\b", s)
    if m2:
        return f"https://www.ncbi.nlm.nih.gov/snp/?term={m2.group(1)}"

    return ""

# ============================================
# CLINVAR LOOKUP (LOCAL clinvar.vcf.gz + .tbi)
# ============================================

from collections import Counter

_CLINVAR_HAS_CHR: Optional[bool] = None

def _detect_has_chr(vcf_path: str) -> Optional[bool]:
    try:
        vf = pysam.VariantFile(vcf_path)
        contigs = list(vf.header.contigs)
        vf.close()
        if not contigs:
            return None
        return any(str(c).startswith("chr") for c in contigs)
    except Exception:
        return None

def _normalize_chrom(chrom: str, target_has_chr: bool) -> str:
    chrom = str(chrom)
    if target_has_chr:
        return chrom if chrom.startswith("chr") else "chr" + chrom
    return chrom[3:] if chrom.startswith("chr") else chrom

def _first_val(info, key: str) -> str:
    """Return first value for an INFO key (works for list/tuple scalars)."""
    if key not in info:
        return ""
    v = info[key]
    if v is None:
        return ""
    if isinstance(v, (list, tuple)):
        return str(v[0]) if v else ""
    return str(v)

def _all_vals_joined(info, key: str, sep: str = "|") -> str:
    """Return all values joined (useful for multi-valued ClinVar fields)."""
    if key not in info:
        return ""
    v = info[key]
    if v is None:
        return ""
    if isinstance(v, (list, tuple)):
        vals = [str(x) for x in v if x not in (None, "", ".")]
        return sep.join(vals)
    s = str(v).strip()
    return "" if s in ("", ".") else s

def _clinvar_sig_primary(clnsig_raw: str) -> str:
    if not clnsig_raw:
        return ""
    s = clnsig_raw.replace("_", " ").strip().lower()
    if "conflict" in s:
        return "Conflicting"
    if "pathogenic" in s and "likely" in s:
        return "Likely pathogenic"
    if "pathogenic" in s:
        return "Pathogenic"
    if "uncertain" in s or "vus" in s:
        return "VUS"
    if "benign" in s and "likely" in s:
        return "Likely benign"
    if "benign" in s:
        return "Benign"
    return clnsig_raw.strip()

def _clinvar_is_conflicting(clnsig_raw: str) -> bool:
    if not clnsig_raw:
        return False
    s = clnsig_raw.lower()
    return ("conflicting" in s) or ("conflict" in s) or (("benign" in s) and ("pathogenic" in s))

def _clinvar_review_stars(revstat_raw: str) -> int:
    if not revstat_raw:
        return 0
    s = revstat_raw.replace("_", " ").lower()
    if "practice guideline" in s:
        return 4
    if "reviewed by expert panel" in s:
        return 3
    if "multiple submitters" in s and "no conflicts" in s:
        return 2
    if "single submitter" in s:
        return 1
    return 0

def _normalize_sig_label(label: str) -> str:
    """Normalize ClinVar labels into your display buckets."""
    l = (label or "").strip().replace("_", " ")
    low = l.lower()
    if "likely benign" in low:
        return "Likely benign"
    if low == "benign":
        return "Benign"
    if "likely pathogenic" in low:
        return "Likely pathogenic"
    if low == "pathogenic":
        return "Pathogenic"
    if "uncertain" in low or "vus" in low:
        return "VUS"
    if "conflict" in low:
        return "Conflicting"
    return l.strip()

def _parse_clnsigconf_counts(clnsigconf_raw: str) -> Counter:
    """
    Parse CLNSIGCONF like:
      'Pathogenic (1),Uncertain_significance (1)'
    Returns Counter({ 'Pathogenic': 1, 'VUS': 1 })
    """
    c = Counter()
    if not clnsigconf_raw:
        return c

    # Split on comma or pipe; ClinVar annotations vary by source/tools.
    parts = re.split(r"[|,]", clnsigconf_raw)
    for p in parts:
        p = p.strip()
        if not p:
            continue

        # Match "Label (N)" or just "Label"
        m = re.match(r"^(.*?)(?:\s*\(\s*(\d+)\s*\))?$", p)
        if not m:
            continue

        label = _normalize_sig_label(m.group(1))
        n = int(m.group(2)) if m.group(2) else 1
        if label:
            c[label] += n

    # Don’t count the literal word "Conflicting" as one of the sides
    if "Conflicting" in c:
        del c["Conflicting"]

    return c

def _format_conflict_counts(counter: Counter) -> str:
    """
    Format as:
      'Benign (3) / Likely benign (1) vs Pathogenic (2)'
    """
    if not counter:
        return ""

    benignish = ["Benign", "Likely benign"]
    pathish = ["Pathogenic", "Likely pathogenic"]

    def fmt_group(keys):
        items = [f"{k} ({counter[k]})" for k in keys if counter.get(k, 0) > 0]
        return " / ".join(items)

    left = fmt_group(benignish)
    right = fmt_group(pathish)

    if left and right:
        return f"{left} vs {right}"

    # fallback: show everything we have (stable order-ish)
    return " / ".join(f"{k} ({v})" for k, v in counter.most_common())

def clinvar_lookup_variant(chrom: str, pos: int, ref: str, alt: str) -> Optional[Dict[str, Any]]:
    global _CLINVAR_HAS_CHR

    if not CLINVAR_VCF_GZ or not os.path.exists(CLINVAR_VCF_GZ):
        return None

    if _CLINVAR_HAS_CHR is None:
        _CLINVAR_HAS_CHR = _detect_has_chr(CLINVAR_VCF_GZ)
        if _CLINVAR_HAS_CHR is None:
            _CLINVAR_HAS_CHR = True

    cv_chrom = _normalize_chrom(chrom, target_has_chr=_CLINVAR_HAS_CHR)

    try:
        cv = pysam.VariantFile(CLINVAR_VCF_GZ)

        for rec in cv.fetch(cv_chrom, pos - 1, pos):
            if rec.pos != pos:
                continue
            if rec.ref != ref:
                continue
            if not rec.alts or alt not in rec.alts:
                continue

            info = rec.info

            # Pull ALL values where possible
            clnsig = _all_vals_joined(info, "CLNSIG", sep="|")
            clnsigconf = _all_vals_joined(info, "CLNSIGCONF", sep="|")  # <- THIS is the key for counts
            revstat = _first_val(info, "CLNREVSTAT")
            clndn = _first_val(info, "CLNDN")

            # Variation ID = VCF ID (often VCV########)
            variation_id = rec.id or ""

            # rs numeric ID (ClinVar VCFs often store numeric in INFO/RS)
            rs_raw = _first_val(info, "RS")
            rsid = ""
            if rs_raw and rs_raw not in (".", "0"):
                m = re.search(r"(\d+)", str(rs_raw))
                if m:
                    rsid = f"rs{m.group(1)}"

            # Build conflict detail (counts preferred)
            conflict_detail = ""
            if _clinvar_is_conflicting(clnsig):
                counts = _parse_clnsigconf_counts(clnsigconf)
                conflict_detail = _format_conflict_counts(counts)

            out = {
                "clnsig": clnsig,
                "clnsig_primary": _clinvar_sig_primary(clnsig),
                "conflicting": _clinvar_is_conflicting(clnsig),

                # what your template prints under “Details …”
                "conflict_summary": conflict_detail,

                "revstat": revstat,
                "stars": _clinvar_review_stars(revstat),
                "clndn": clndn,
                "clnvc": _first_val(info, "CLNVC"),
                "allele_id": _first_val(info, "ALLELEID"),
                "variation_id": variation_id,
                "rsid": rsid,
            }

            cv.close()
            return out

        cv.close()
        return None

    except Exception as e:
        print(f"ClinVar lookup error {chrom}:{pos} {ref}>{alt}: {e}")
        return None
# ============================================
# ANNOTATION LAYOUT (ANN/CSQ)
# ============================================

def get_annotation_layout(vcf):
    if "ANN" in vcf.header.info:
        desc = vcf.header.info["ANN"].description or ""
        fields_part = desc.split("Format:")[-1].strip().strip('"')
        fields = [f.strip() for f in fields_part.split("|")] if fields_part else []
        uf = [f.upper() for f in fields]

        def idx(names: Set[str]) -> Optional[int]:
            for n in names:
                if n in uf:
                    return uf.index(n)
            return None

        gene_idx = idx({"GENE_NAME", "GENE"})
        tx_idx = idx({"FEATURE_ID", "TRANSCRIPT"})
        c_idx = idx({"HGVSC", "HGVS.C"})
        p_idx = idx({"HGVSP", "HGVS.P"})
        eff_idx = idx({"ANNOTATION"})

        if not fields or gene_idx is None or tx_idx is None or c_idx is None or p_idx is None or eff_idx is None:
            gene_idx = 3
            eff_idx = 1
            tx_idx = 6
            c_idx = 9
            p_idx = 10

        return "ANN", gene_idx, tx_idx, c_idx, p_idx, eff_idx

    if "CSQ" in vcf.header.info:
        desc = vcf.header.info["CSQ"].description or ""
        fields_part = desc.split("Format:")[-1].strip().strip('"')
        fields = [f.strip() for f in fields_part.split("|")] if fields_part else []
        uf = [f.upper() for f in fields]

        def idx(names: Set[str]) -> Optional[int]:
            for n in names:
                if n in uf:
                    return uf.index(n)
            return None

        gene_idx = idx({"SYMBOL", "GENE"})
        tx_idx = idx({"FEATURE", "TRANSCRIPT"})
        c_idx = idx({"HGVSC", "HGVS_C"})
        p_idx = idx({"HGVSP", "HGVS_P"})
        eff_idx = idx({"CONSEQUENCE"})
        return "CSQ", gene_idx, tx_idx, c_idx, p_idx, eff_idx

    return None, None, None, None, None, None

def extract_variant_info(rec, tag, gi, ti, ci, pi, ei):
    gene = transcript = hgvsc = hgvsp = effect = ""

    if tag and tag in rec.info:
        raw = rec.info[tag]
        if isinstance(raw, (list, tuple)):
            raw = raw[0]
        parts = str(raw).split("|")

        if gi is not None and gi < len(parts):
            g = parts[gi].strip()
            if g and g != ".":
                gene = g

        if ti is not None and ti < len(parts):
            t = parts[ti].strip()
            if t and t != ".":
                transcript = t

        if ci is not None and ci < len(parts):
            cv = parts[ci].strip()
            if cv and cv != ".":
                if ":" in cv:
                    tx, cval = cv.split(":", 1)
                    if not transcript:
                        transcript = tx.strip()
                    hgvsc = cval.strip()
                else:
                    hgvsc = cv

        if pi is not None and pi < len(parts):
            pv = parts[pi].strip()
            if pv and pv != ".":
                if ":" in pv:
                    hgvsp = pv.split(":", 1)[1].strip()
                else:
                    hgvsp = pv

        if ei is not None and ei < len(parts):
            ev = parts[ei].strip()
            if ev and ev != ".":
                effect = ev

    if not gene:
        gene = positional_gene_lookup(rec.chrom, rec.pos) or ""

    return gene, transcript, hgvsc, hgvsp, effect


# ============================================
# BASIC FILTERS
# ============================================

PROTEIN_ALTERING = {
    "missense_variant", "stop_gained", "stop_lost",
    "frameshift_variant",
    "splice_acceptor_variant", "splice_donor_variant", "splice_region_variant",
    "start_lost",
    "inframe_insertion", "inframe_deletion",
    "protein_altering_variant",
}

def is_protein_altering(effect: str) -> bool:
    if not effect:
        return True
    parts = re.split(r"[&,]", effect)
    return any(p.strip() in PROTEIN_ALTERING for p in parts)

def is_nonref(rec) -> bool:
    if not rec.samples:
        return True
    for s in rec.samples.values():
        gt = s.get("GT")
        if gt and any(a != 0 for a in gt if a is not None):
            return True
    return False


# ============================================
# ACMG (STUB)
# ============================================

def dummy_acmg(rec, gene: str, hgvsc: str, hgvsp: str) -> Dict[str, Any]:
    hit = lookup_curated_variant(gene, hgvsc, hgvsp)
    if hit:
        return {
            "classification": hit["classification"] or "VUS",
            "evidence": [
                f"Curated classification (source: {hit.get('source') or 'internal'}).",
                hit.get("notes") or "",
            ],
        }
    return {"classification": "VUS", "evidence": ["Placeholder ACMG logic."]}


# ============================================
# ROUTES
# ============================================

@app.get("/", response_class=HTMLResponse)
async def index(request: Request):
    return templates.TemplateResponse(
        "index.html",
        {"request": request, "gene_symbols": ALL_GENE_SYMBOLS},
    )


@app.post("/analyze", response_class=HTMLResponse)
async def analyze(
    request: Request,
    vcf_file: UploadFile = File(...),
    hpo_terms: str = Form(""),
    annotate: str = Form("no"),
    acmg_filter: List[str] = Form([]),
    gene_filter: str = Form(""),
    max_af: str = Form(""),
    clinvar_filter: List[str] = Form([]),
    exclude_clinvar_conflicts: str = Form("no"),
):
    tmp = tempfile.mkdtemp()
    raw_vcf = os.path.join(tmp, vcf_file.filename)
    cleaned_vcf = os.path.join(tmp, "cleaned_" + vcf_file.filename.replace(".gz", ""))
    ann_vcf = os.path.join(tmp, "annotated_" + vcf_file.filename.replace(".gz", ""))

    try:
        with open(raw_vcf, "wb") as f:
            shutil.copyfileobj(vcf_file.file, f)

        annotation_warning = None

        hpo_display_terms = [t.strip() for t in hpo_terms.replace("\n", ",").split(",") if t.strip()]
        hpo_ids = normalize_hpo_terms(hpo_terms) if hpo_terms.strip() else []

        # HPO genes
        hpo_genes: Set[str] = set()
        for h in hpo_ids:
            hpo_genes |= HPO_GENE_MAP.get(h, set())

        # Manual gene filter
        gene_filter_set: Set[str] = set()
        if gene_filter.strip():
            gene_filter_set = {
                g.strip().upper()
                for g in gene_filter.replace("\r", "\n").replace("\n", ",").split(",")
                if g.strip()
            }
        gene_filter_display: List[str] = sorted(gene_filter_set)

        # Max AF (will apply to gnomAD AF)
        max_af_value: Optional[float] = None
        if max_af.strip():
            try:
                max_af_value = float(max_af)
            except ValueError:
                max_af_value = None

        # Clean FORMAT fields to avoid pysam parse errors
        base_vcf = raw_vcf
        try:
            clean_vcf_remove_formats(raw_vcf, cleaned_vcf)
            base_vcf = cleaned_vcf
        except Exception as e:
            annotation_warning = (annotation_warning or "") + f" FORMAT cleanup failed ({e}). "

        # Annotate with snpEff if requested
        vcf_to_read = base_vcf
        if annotate == "yes":
            try:
                run_snpeff(base_vcf, ann_vcf)
                vcf_to_read = ann_vcf
            except Exception as e:
                annotation_warning = (annotation_warning or "") + f" snpEff failed ({e}). Using unannotated VCF. "

        dlog("Reading VCF:", vcf_to_read)
        dlog("GNOMAD_DIR:", GNOMAD_DIR)
        dlog("GNOMAD genomes glob:", GNOMAD_GENOMES_GLOB)
        dlog("GNOMAD exomes  glob:", GNOMAD_EXOMES_GLOB)

        vcf = pysam.VariantFile(vcf_to_read)
        tag, gi, ti, ci, pi, ei = get_annotation_layout(vcf)
        if tag is None:
            annotation_warning = (annotation_warning or "") + " No ANN/CSQ detected; c./p. may be missing. "

        selected_acmg = set(acmg_filter) if acmg_filter else None

        filtered: List[Dict[str, Any]] = []
        parsed_gene_count = 0
        parse_error = None

        try:
            for rec in vcf:
                gene, transcript, hgvsc, hgvsp, effect = extract_variant_info(rec, tag, gi, ti, ci, pi, ei)
                if gene:
                    parsed_gene_count += 1

                if not is_nonref(rec):
                    continue
                if not is_protein_altering(effect):
                    continue

                if hpo_genes and gene not in hpo_genes:
                    continue
                if gene_filter_set and gene.upper() not in gene_filter_set:
                    continue

                alts = list(rec.alts or [])
                if not alts:
                    continue

                mane = MANE_MAP.get(gene, {})
                clingen = CLINGEN_MAP.get((gene or "").strip().upper(), {})

                hgnc_id = clingen.get("hgnc_id", "")
                clingen_gene_url = ""
                if hgnc_id:
                    clingen_gene_url = f"https://search.clinicalgenome.org/kb/genes/{hgnc_id.replace(':', '%3A')}"
                elif clingen.get("url"):
                    clingen_gene_url = clingen["url"]

                for alt in alts:
                    alt_str = str(alt)

                    # sample metrics
                    sample_metrics = get_sample_metrics(rec, alt_str) or {}
                    sample_metrics.setdefault("dp", None)
                    sample_metrics.setdefault("ref_depth", None)
                    sample_metrics.setdefault("alt_depth", None)
                    sample_metrics.setdefault("vaf", None)
                    sample_metrics.setdefault("gt", "")
                    sample_metrics.setdefault("zygosity", "")

                    # ✅ ALWAYS use gnomAD local AF (skip VCF AF entirely)
                    gnomad = gnomad_lookup_af(rec.chrom, rec.pos, rec.ref, alt_str)
                    gnomad_af = (gnomad or {}).get("af")
                    
                    # Only log if we actually found something to reduce noise
                    if gnomad_af is not None:
                         dlog("Variant", rec.chrom, rec.pos, rec.ref, alt_str,
                              "-> gnomAD:", gnomad)

                    # Apply max_af filter using gnomAD AF
                    if max_af_value is not None and gnomad_af is not None and gnomad_af > max_af_value:
                        continue

                    clinvar = clinvar_lookup_variant(rec.chrom, rec.pos, rec.ref, alt_str)
                    has_clinvar = bool(clinvar and (clinvar.get("clnsig") or clinvar.get("allele_id") or clinvar.get("rsid")))
                    cv_primary = (clinvar or {}).get("clnsig_primary", "")
                    cv_conflict = bool((clinvar or {}).get("conflicting", False))
                    
                    # dbSNP rsID: prefer ClinVar rsID if present, else uploaded VCF ID
                    dbsnp_rsid = (clinvar or {}).get("rsid") or (rec.id or "")
                    dbsnp_url = make_dbsnp_url(dbsnp_rsid)

                    selected_cv = set(clinvar_filter or [])
                    apply_cv_filter = bool(selected_cv) and ("any" not in selected_cv)

                    if apply_cv_filter:
                        matches = False
                        if "has" in selected_cv and has_clinvar:
                            matches = True
                        if "plp" in selected_cv and cv_primary in ("Pathogenic", "Likely pathogenic"):
                            matches = True
                        if "vus" in selected_cv and cv_primary == "VUS":
                            matches = True
                        if "blb" in selected_cv and cv_primary in ("Benign", "Likely benign"):
                            matches = True
                        if "conflicting" in selected_cv and cv_conflict:
                            matches = True
                        if "not_conflicting" in selected_cv and has_clinvar and not cv_conflict:
                            matches = True
                        if not matches:
                            continue

                    if exclude_clinvar_conflicts == "yes" and cv_conflict:
                        continue

                    # -------------------------------------------------------------
                    # ✅ LIVE CLASSIFICATION (Lite Mode)
                    # Use the shared ACMG Classifier to check ClinVar/Freq/PVS1 locally.
                    # -------------------------------------------------------------
                    
                    # Lazy Import to avoid circular deps or path issues at top level
                    if 'ACMGClassifier' not in globals():
                        import sys
                        sys.path.append('../genetic_reporter')
                        try:
                            from acmg_classifier import ACMGClassifier
                            # Initialize without API key for Lite Mode (safe)
                            global_classifier = ACMGClassifier(api_key="sk-dummy") 
                        except ImportError:
                            global_classifier = None
                    
                    acmg_class = "VUS" # Default
                    acmg_evidence = []
                    
                    if 'global_classifier' in locals() and global_classifier:
                         # Construct Context for Classifier
                         # We need to pass: user_af, user_clinvar, user_effect, chrom, pos, ref, alt
                         ctx = {
                             'af': gnomad_af,
                             # Use Primary Signal (e.g. "Pathogenic"), fallback to legacy 'clnsig' if primary is missing
                             'clinvar': cv_primary or (clinvar or {}).get("clnsig", ""),
                             'effect': effect,
                             'chrom': rec.chrom,
                             'pos': rec.pos,
                             'ref': rec.ref,
                             'alt': alt_str
                         }
                         try:
                             # RUN LITE MODE (No AI, just Rules)
                             # Disable local lookup to prevent PySAM vs CyVCF2 conflicts in the Portal process
                             result = global_classifier.classify(gene, hgvsc, context=ctx, lite_mode=True, do_local_lookup=False)
                             acmg_class = result.acmg_classification
                             acmg_evidence = result.acmg_evidence
                         except Exception as e:
                             print(f"Classifier Error for {gene}: {e}")
                             acmg_class = "VUS"

                    # If user selected specific classes to view, filter now
                    if selected_acmg and acmg_class not in selected_acmg:
                        continue

                    filtered.append({
                        "chrom": rec.chrom,
                        "pos": rec.pos,
                        "ref": rec.ref,
                        "alt": alt_str,

                        "gene": gene,
                        "transcript": transcript,
                        "canonical_enst": mane.get("enst", ""),
                        "clinical_nm": mane.get("nm", ""),
                        "hgvsc": hgvsc,
                        "hgvsp": hgvsp,
                        "effect": effect,

                        "dp": sample_metrics["dp"],
                        "ref_depth": sample_metrics["ref_depth"],
                        "alt_depth": sample_metrics["alt_depth"],
                        "vaf": sample_metrics["vaf"],
                        "gt": sample_metrics["gt"],
                        "zygosity": sample_metrics["zygosity"],

                        # ✅ This is what your results column should display now:
                        "gnomad_af": gnomad_af,
                        "gnomad_exomes_af": (gnomad or {}).get("exomes_af"),
                        "gnomad_genomes_af": (gnomad or {}).get("genomes_af"),
                        "gnomad_source": (gnomad or {}).get("source", ""),

                        "clinvar_sig": (clinvar or {}).get("clnsig", ""),
                        "clinvar_sig_primary": (clinvar or {}).get("clnsig_primary", ""),
                        "clinvar_conflicting": bool((clinvar or {}).get("conflicting", False)),
                        "clinvar_conflict_detail": (clinvar or {}).get("conflict_summary", ""),
                        "clinvar_revstat": (clinvar or {}).get("revstat", ""),
                        "clinvar_stars": int((clinvar or {}).get("stars", 0) or 0),
                        "clinvar_dn": (clinvar or {}).get("clndn", ""),
                        "dbsnp_rsid": dbsnp_rsid,
                        "dbsnp_url": dbsnp_url,

                        # optional but useful if you want to display it somewhere:
                        "clinvar_variation_id": (clinvar or {}).get("variation_id", ""),

                        "clingen_hi_score": clingen.get("hi_score"),
                        "clingen_ts_score": clingen.get("ts_score"),
                        "clingen_hi_desc": clingen.get("hi_desc", ""),
                        "clingen_ts_desc": clingen.get("ts_desc", ""),
                        "clingen_last_eval": clingen.get("last_eval", ""),
                        "clingen_hi_disease_id": clingen.get("hi_disease_id", ""),
                        "clingen_ts_disease_id": clingen.get("ts_disease_id", ""),
                        "clingen_hgnc_id": hgnc_id,
                        "clingen_url": clingen_gene_url,

                        "acmg_class": acmg_class,
                        "acmg_evidence": acmg_evidence,
                    })

        except OSError as e:
            parse_error = f"Error parsing VCF: {e}"

        vcf.close()

        if parse_error:
            annotation_warning = (annotation_warning or "") + " " + parse_error

        # Notes
        clinvar_note = ""
        if not os.path.exists(CLINVAR_VCF_GZ):
            clinvar_note = f"ClinVar VCF not found at: {CLINVAR_VCF_GZ} (set CLINVAR_DIR or CLINVAR_VCF_GZ)."

        clingen_note = ""
        if not os.path.exists(CLINGEN_GENE_TSV):
            clingen_note = f"ClinGen TSV not found at: {CLINGEN_GENE_TSV} (set CLINGEN_DIR or CLINGEN_GENE_TSV)."

        # Optional: show debug buffer inside results page
        debug_text = get_debug_dump() if (DEBUG and SHOW_DEBUG_IN_RESULTS) else ""

        return templates.TemplateResponse(
            "results.html",
            {
                "request": request,
                "variants": filtered,
                "hpo_terms_display": hpo_display_terms,
                "hpo_terms": hpo_ids,
                "hpo_genes": sorted(hpo_genes),
                "gene_filter_display": gene_filter_display,

                "annotation_warning": annotation_warning,
                "clinvar_note": clinvar_note,
                "clingen_note": clingen_note,
                "parsed_gene_count": parsed_gene_count,
                "matched_variant_count": len(filtered),

                # Debug context (optional)
                "debug_enabled": DEBUG,
                "debug_text": debug_text,
            },
        )

    finally:
        shutil.rmtree(tmp, ignore_errors=True)
