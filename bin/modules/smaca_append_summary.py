#!/usr/bin/env python3
"""
Append ##DARK_GENE_PIPELINE_SUMMARY + one line after raw SMAca TSV (called from Nextflow SMACA_RUN).
Kept as a separate file so Groovy does not parse Python f-strings inside process script \"\"\" blocks.
"""
from __future__ import annotations

import ast
import os
import re
import sys

import pysam


def pick_chr5(bam: pysam.AlignmentFile) -> str | None:
    for c in ("chr5", "5"):
        if c in bam.references:
            return c
    return None


def complement_dna(b: str) -> str:
    m = {"A": "T", "T": "A", "C": "G", "G": "C"}
    return m.get(b.upper(), "N")


def read_base_on_reference_forward(aln: pysam.AlignedSegment, query_pos: int) -> str:
    """Base aligned to the forward reference strand at query_pos (IUPAC)."""
    b = aln.query_sequence[query_pos].upper()
    if aln.is_reverse:
        return complement_dna(b)
    return b


def base_in_hgvs_coding_space(ref_forward_base: str, minus_strand_gene: bool) -> str:
    """Map forward-reference base to HGVS coding-strand letter for the gene."""
    if not minus_strand_gene:
        return ref_forward_base
    return complement_dna(ref_forward_base)


def pileup_acgt_at_1bp(
    bam_path: str,
    pos_1based: int,
    *,
    hgvs_coding_minus_strand: bool = False,
    contig: str | None = None,
):
    # Count A/C/G/T/N at one reference position (1-based).
    # contig: if set (e.g. unified SMN mini-ref), use it; else chr5/5 on hg38 BAM.
    # When hgvs_coding_minus_strand=True (SMN c.840), bases are tallied in transcript
    # (HGVS coding) sense for the minus-strand SMN genes.
    bam = pysam.AlignmentFile(bam_path, "rb")
    chrom = contig if contig else pick_chr5(bam)
    if not chrom:
        bam.close()
        return None
    start0 = pos_1based - 1
    end0 = pos_1based
    counts = {"A": 0, "C": 0, "G": 0, "T": 0, "N": 0}
    try:
        for col in bam.pileup(
            contig=chrom,
            start=start0,
            stop=end0,
            truncate=True,
            min_mapping_quality=0,
            min_base_quality=0,
            max_depth=1_000_000,
        ):
            for pr in col.pileups:
                if pr.is_del or pr.is_refskip:
                    continue
                if pr.query_position is None:
                    continue
                aln = pr.alignment
                qp = pr.query_position
                rf = read_base_on_reference_forward(aln, qp)
                b = base_in_hgvs_coding_space(rf, hgvs_coding_minus_strand)
                if b in counts:
                    counts[b] += 1
                else:
                    counts["N"] += 1
    except Exception:
        bam.close()
        return None
    bam.close()
    return counts


def fmt_norm_ct_fields(ct: dict | None, key_prefix: str) -> str:
    # c.840 discriminator-style C vs T: merge G→C and A→T (complement pairs on minus strand).
    if not ct:
        return (
            f"|{key_prefix}_norm_C=NA|{key_prefix}_norm_T=NA|"
            f"{key_prefix}_norm_C_frac_CplusT=NA"
        )
    c_n = ct["C"] + ct["G"]
    t_n = ct["T"] + ct["A"]
    s = c_n + t_n
    frac = f"{c_n / s:.3f}" if s > 0 else "?"
    return (
        f"|{key_prefix}_norm_C={c_n}|{key_prefix}_norm_T={t_n}|"
        f"{key_prefix}_norm_C_frac_CplusT={frac}"
    )


def fmt_c840_locus_both(tag: str, ct_hgvs: dict | None, ct_fwd: dict | None) -> str:
    # One hg38 coordinate: SMN1 vs SMN2 paralog site (c.840 discriminator at that alignment column).
    pref = f"c840_{tag}"
    parts = []
    if not ct_hgvs:
        parts.append(
            f"|{pref}_hgvs_A=NA|{pref}_hgvs_C=NA|{pref}_hgvs_G=NA|{pref}_hgvs_T=NA|"
            f"{pref}_hgvs_depth=NA|{pref}_hgvs_C_frac_CplusT_coding=NA"
        )
        parts.append(fmt_norm_ct_fields(None, f"{pref}_hgvs"))
    else:
        a, c, g, t = ct_hgvs["A"], ct_hgvs["C"], ct_hgvs["G"], ct_hgvs["T"]
        dep = sum(ct_hgvs.values())
        ct_sum = c + t
        c_frac = f"{c / ct_sum:.3f}" if ct_sum > 0 else "?"
        parts.append(
            f"|{pref}_hgvs_A={a}|{pref}_hgvs_C={c}|{pref}_hgvs_G={g}|{pref}_hgvs_T={t}|"
            f"{pref}_hgvs_depth={dep}|{pref}_hgvs_C_frac_CplusT_coding={c_frac}"
        )
        parts.append(fmt_norm_ct_fields(ct_hgvs, f"{pref}_hgvs"))
    if not ct_fwd:
        parts.append(
            f"|{pref}_fwd_A=NA|{pref}_fwd_C=NA|{pref}_fwd_G=NA|{pref}_fwd_T=NA|"
            f"{pref}_fwd_depth=NA|{pref}_fwd_C_frac_CplusT=NA"
        )
        parts.append(fmt_norm_ct_fields(None, f"{pref}_fwd"))
    else:
        a, c, g, t = ct_fwd["A"], ct_fwd["C"], ct_fwd["G"], ct_fwd["T"]
        dep = sum(ct_fwd.values())
        ct_sum = c + t
        c_frac = f"{c / ct_sum:.3f}" if ct_sum > 0 else "?"
        parts.append(
            f"|{pref}_fwd_A={a}|{pref}_fwd_C={c}|{pref}_fwd_G={g}|{pref}_fwd_T={t}|"
            f"{pref}_fwd_depth={dep}|{pref}_fwd_C_frac_CplusT={c_frac}"
        )
        parts.append(fmt_norm_ct_fields(ct_fwd, f"{pref}_fwd"))
    return "".join(parts)


def fmt_c840_merged_both(ct_hgvs: dict | None, ct_fwd: dict | None) -> str:
    # Single c.840 column (unified BAM only): genome forward + HGVS coding ACGT.
    parts = []
    if not ct_hgvs:
        parts.append(
            "|c840_merged_hgvs_A=NA|c840_merged_hgvs_C=NA|c840_merged_hgvs_G=NA|c840_merged_hgvs_T=NA|"
            "c840_merged_hgvs_depth=NA|c840_merged_hgvs_C_frac_CplusT_coding=NA"
        )
        parts.append(fmt_norm_ct_fields(None, "c840_merged_hgvs"))
    else:
        a, c, g, t = ct_hgvs["A"], ct_hgvs["C"], ct_hgvs["G"], ct_hgvs["T"]
        dep = sum(ct_hgvs.values())
        ct_sum = c + t
        c_frac = f"{c / ct_sum:.3f}" if ct_sum > 0 else "?"
        parts.append(
            f"|c840_merged_hgvs_A={a}|c840_merged_hgvs_C={c}|c840_merged_hgvs_G={g}|c840_merged_hgvs_T={t}|"
            f"c840_merged_hgvs_depth={dep}|c840_merged_hgvs_C_frac_CplusT_coding={c_frac}"
        )
        parts.append(fmt_norm_ct_fields(ct_hgvs, "c840_merged_hgvs"))
    if not ct_fwd:
        parts.append(
            "|c840_merged_fwd_A=NA|c840_merged_fwd_C=NA|c840_merged_fwd_G=NA|c840_merged_fwd_T=NA|"
            "c840_merged_fwd_depth=NA|c840_merged_fwd_C_frac_CplusT=NA"
        )
        parts.append(fmt_norm_ct_fields(None, "c840_merged_fwd"))
    else:
        a, c, g, t = ct_fwd["A"], ct_fwd["C"], ct_fwd["G"], ct_fwd["T"]
        dep = sum(ct_fwd.values())
        ct_sum = c + t
        c_frac = f"{c / ct_sum:.3f}" if ct_sum > 0 else "?"
        parts.append(
            f"|c840_merged_fwd_A={a}|c840_merged_fwd_C={c}|c840_merged_fwd_G={g}|c840_merged_fwd_T={t}|"
            f"c840_merged_fwd_depth={dep}|c840_merged_fwd_C_frac_CplusT={c_frac}"
        )
        parts.append(fmt_norm_ct_fields(ct_fwd, "c840_merged_fwd"))
    return "".join(parts)


def unified_bam_contig(bam_path: str) -> str | None:
    bam = pysam.AlignmentFile(bam_path, "rb")
    refs = list(bam.references)
    bam.close()
    if len(refs) == 1:
        return refs[0]
    return None


def silent_carrier_dup_g27134_allele_counts(raw: str):
    # g.27134T>G DUP_MARK cell (silent carrier / duplication context), not c.840 exon 7.
    if not raw or raw == "-":
        return None
    m = re.search(r"\[.*\]", raw, re.DOTALL)
    if not m:
        return None
    try:
        arr = ast.literal_eval(m.group(0))
    except Exception:
        return None
    if not isinstance(arr, list) or len(arr) != 4:
        return None
    if not all(isinstance(r, list) for r in arr):
        return None
    A = int(sum(arr[0]))
    C = int(sum(arr[1]))
    G = int(sum(arr[2]))
    T = int(sum(arr[3]))
    tg = T + G
    ct = C + T
    t_frac_tg = f"{T / tg:.3f}" if tg > 0 else "?"
    c_frac_ct = f"{C / ct:.3f}" if ct > 0 else "?"
    return {
        "A": A,
        "C": C,
        "G": G,
        "T": T,
        "T_frac_TG": t_frac_tg,
        "C_frac_CT": c_frac_ct,
    }


def silent_carrier_summary_label(raw_val: str, acgt: dict | None) -> str:
    # Silent carrier (g.27134T>G dup context): use SMAca G read counts when parsed;
    # otherwise fall back to legacy "G" in the text before the [[...]] array.
    if raw_val == "-":
        return "False"
    raw_head = raw_val.split("[", 1)[0]
    if acgt is not None:
        g = acgt["G"]
        call = "True" if g > 0 else "False"
        return f"{call} (g27134_G={g}; Raw: {raw_head})"
    # Unparsed cell: avoid "G" in "T>G" false positives before [[...]]
    if raw_head.strip().startswith("G"):
        return f"True (g27134_G=NA; Raw: {raw_head})"
    return f"False (g27134_G=NA; Raw: {raw_head})"


def depth_tot(ct) -> float | None:
    if not ct:
        return None
    return float(sum(ct.values()))


def fmt_hgvs_c_t(ct: dict | None) -> str:
    """HGVS-coding strand: literal C and T read counts at c.840 pileup."""
    if not ct:
        return "C=NA T=NA depth=NA"
    dep = int(sum(ct.values()))
    return f"C={ct['C']} T={ct['T']} depth={dep}"


def fmt_hgvs_norm_disc(ct: dict | None) -> str:
    """HGVS pileup: discriminator-style C vs T after G→C and A→T merge."""
    if not ct:
        return "norm_C=NA norm_T=NA C/(C+T)_norm=NA"
    c_n = ct["C"] + ct["G"]
    t_n = ct["T"] + ct["A"]
    s = c_n + t_n
    frac = f"{c_n / s:.3f}" if s > 0 else "?"
    return f"norm_C={c_n} norm_T={t_n} C/(C+T)_norm={frac}"


def line_c840_ct_readable(
    unified_used: bool,
    ct_u: dict | None,
    ct_s1: dict | None,
    ct_s2: dict | None,
) -> str:
    """One grep-friendly line with C and T counts (HGVS coding)."""
    if unified_used and ct_u is not None:
        return f"##C840_C_T_PILEUP mode=unified_SMN1_ref HGVS_coding {fmt_hgvs_c_t(ct_u)}"
    return (
        "##C840_C_T_PILEUP mode=hg38_two_positions "
        f"SMN1_HGVS_coding {fmt_hgvs_c_t(ct_s1)} | SMN2_HGVS_coding {fmt_hgvs_c_t(ct_s2)}"
    )


def line_c840_norm_readable(
    unified_used: bool,
    ct_u: dict | None,
    ct_s1: dict | None,
    ct_s2: dict | None,
) -> str:
    """G→C, A→T normalized C:T discriminator line (HGVS coding)."""
    if unified_used and ct_u is not None:
        return (
            "##C840_C_T_NORM (G→C, A→T) mode=unified_SMN1_ref HGVS_coding "
            f"{fmt_hgvs_norm_disc(ct_u)}"
        )
    return (
        "##C840_C_T_NORM (G→C, A→T) mode=hg38_two_positions "
        f"SMN1 {fmt_hgvs_norm_disc(ct_s1)} | SMN2 {fmt_hgvs_norm_disc(ct_s2)}"
    )


def main() -> None:
    if len(sys.argv) < 6:
        print(
            "Usage: smaca_append_summary.py <sample_id> <bam_path> "
            "<smn_cnv_est_total_copies> <smn_c840_sm1_1bp> <smn_c840_sm2_1bp> "
            "[<unified_smn_bam|-> [<unified_c840_1bp_in_unified_ref>]]",
            file=sys.stderr,
        )
        sys.exit(1)
    sample_id, bam_path, tc_est_s, pos_s1_s, pos_s2_s = sys.argv[1:6]
    tc_est = int(tc_est_s)
    pos_s1 = int(pos_s1_s)
    pos_s2 = int(pos_s2_s)
    unified_bam: str | None = None
    unified_pos: int | None = None
    if len(sys.argv) >= 8 and sys.argv[7] not in ("", "-"):
        unified_bam = sys.argv[7]
        unified_pos = int(sys.argv[8]) if len(sys.argv) > 8 else None

    path = f"{sample_id}_smaca.txt"
    with open(path) as f:
        lines = f.read().splitlines()
    if any(line.strip().startswith("##DARK_GENE_PIPELINE_SUMMARY") for line in lines):
        sys.exit(0)

    header = []
    data = []
    for line in lines:
        line = line.strip()
        if line.startswith("# id") or line.startswith("id"):
            header = [h.strip() for h in line.replace("#", "").strip().split(",")]
        elif not line.startswith("#") and len(line) > 5:
            if "|" in line:
                data = [d.strip() for d in line.split("|")]
            elif "," in line:
                data = [d.strip() for d in line.split(",")]

    nl = chr(10)
    if not (header and data):
        with open(path, "a") as f:
            f.write(nl + "##DARK_GENE_PIPELINE_SUMMARY" + nl)
            f.write("Parsing Failed" + nl)
        sys.exit(0)

    min_len = min(len(header), len(data))
    row = {header[i]: data[i] for i in range(min_len)}
    s1 = row.get("avg_cov_SMN1", "?")
    s2 = row.get("avg_cov_SMN2", "?")
    s1_val = None
    s2_val = None
    try:
        s1_val = float(s1)
        s1 = f"{s1_val:.2f}"
    except Exception:
        pass
    try:
        s2_val = float(s2)
        s2 = f"{s2_val:.2f}"
    except Exception:
        pass

    variant_col = "g.27134T>G"
    raw_val = row.get(variant_col, "-")
    acgt = silent_carrier_dup_g27134_allele_counts(raw_val) if raw_val != "-" else None
    silent_carrier = silent_carrier_summary_label(raw_val, acgt)

    c_ratio = "?"
    try:
        if isinstance(s1_val, float) and isinstance(s2_val, float):
            total_cov = s1_val + s2_val
            if total_cov > 0:
                c_ratio = f"{s1_val / total_cov:.3f}"
            else:
                c_ratio = "0.00"
    except Exception:
        pass

    cn_est1 = "?"
    cn_est2 = "?"
    try:
        if isinstance(s1_val, float) and isinstance(s2_val, float):
            tot = s1_val + s2_val
            if tot > 0:
                sm1_frac = s1_val / tot
                cn_est1 = int(round(tc_est * sm1_frac))
                cn_est2 = tc_est - cn_est1
    except Exception:
        pass

    pi_b = row.get("Pi_b", "?")
    cov_s1b = row.get("cov_SMN1_b", "?")
    cov_s2b = row.get("cov_SMN2_b", "?")
    exon7 = f"|exon7_c840_Pi_b={pi_b}|cov_SMN1_b={cov_s1b}|cov_SMN2_b={cov_s2b}"
    silent_allelic = ""
    if acgt:
        silent_allelic = (
            f"|silent_g27134_A={acgt['A']}|silent_g27134_C={acgt['C']}|silent_g27134_G={acgt['G']}|silent_g27134_T={acgt['T']}"
            f"|silent_g27134_T_frac_TplusG={acgt['T_frac_TG']}|silent_g27134_C_frac_CplusT={acgt['C_frac_CT']}"
        )
    else:
        silent_allelic = "|silent_g27134_alleles=NA"

    ct_s1 = None
    ct_s2 = None
    ct_u = None
    ct_u_fwd = None
    unified_used = False
    if (
        unified_bam
        and os.path.isfile(unified_bam)
        and os.path.getsize(unified_bam) > 800
        and unified_pos is not None
    ):
        uc = unified_bam_contig(unified_bam)
        if uc:
            ct_u = pileup_acgt_at_1bp(
                unified_bam,
                unified_pos,
                hgvs_coding_minus_strand=True,
                contig=uc,
            )
            ct_u_fwd = pileup_acgt_at_1bp(
                unified_bam,
                unified_pos,
                hgvs_coding_minus_strand=False,
                contig=uc,
            )
            if ct_u is not None:
                unified_used = True

    if unified_used and ct_u is not None:
        c840_ct = fmt_c840_merged_both(ct_u, ct_u_fwd)
        pair_frac = "NA"
        delta_smaca_minus_c840 = "NA"
        delta_cov_minus_merged_c_frac = "?"
        try:
            if ct_u_fwd and (ct_u_fwd["C"] + ct_u_fwd["T"]) > 0:
                mcf = ct_u_fwd["C"] / (ct_u_fwd["C"] + ct_u_fwd["T"])
                cr2 = float(c_ratio) if c_ratio != "?" else None
                if cr2 is not None:
                    delta_cov_minus_merged_c_frac = f"{cr2 - mcf:.3f}"
        except Exception:
            pass
        c840_ct += (
            f"|c840_pileup_mode=unified_SMN1_ref|"
            f"c840_SMN1_depth_frac_of_pair=NA|SMAca_cov_frac_minus_c840_depth_frac=NA|"
            f"SMAca_cov_frac_minus_c840_merged_C_frac_CplusT={delta_cov_minus_merged_c_frac}"
        )
    else:
        ct_s1 = pileup_acgt_at_1bp(bam_path, pos_s1, hgvs_coding_minus_strand=True)
        ct_s2 = pileup_acgt_at_1bp(bam_path, pos_s2, hgvs_coding_minus_strand=True)
        ct_s1f = pileup_acgt_at_1bp(bam_path, pos_s1, hgvs_coding_minus_strand=False)
        ct_s2f = pileup_acgt_at_1bp(bam_path, pos_s2, hgvs_coding_minus_strand=False)
        c840_ct = fmt_c840_locus_both("SMN1", ct_s1, ct_s1f) + fmt_c840_locus_both(
            "SMN2", ct_s2, ct_s2f
        )

        d1 = depth_tot(ct_s1)
        d2 = depth_tot(ct_s2)
        pair_frac = "?"
        if d1 is not None and d2 is not None and (d1 + d2) > 0:
            pair_frac = f"{d1 / (d1 + d2):.3f}"

        delta_smaca_minus_c840 = "?"
        try:
            cr = float(c_ratio) if c_ratio != "?" else None
            pr = float(pair_frac) if pair_frac != "?" else None
            if cr is not None and pr is not None:
                delta_smaca_minus_c840 = f"{cr - pr:.3f}"
        except Exception:
            pass

        # Compare SMAca SMN1 fraction to HGVS C/(C+T) at the SMN1 c.840 column only (not summed with SMN2).
        delta_cov_minus_smn1_hgvs_c_frac = "?"
        try:
            if ct_s1 and (ct_s1["C"] + ct_s1["T"]) > 0:
                mcf = ct_s1["C"] / (ct_s1["C"] + ct_s1["T"])
                cr2 = float(c_ratio) if c_ratio != "?" else None
                if cr2 is not None:
                    delta_cov_minus_smn1_hgvs_c_frac = f"{cr2 - mcf:.3f}"
        except Exception:
            pass

        c840_ct += (
            f"|c840_pileup_mode=hg38_two_positions|"
            f"c840_SMN1_depth_frac_of_pair={pair_frac}|SMAca_cov_frac_minus_c840_depth_frac={delta_smaca_minus_c840}|"
            f"SMAca_cov_frac_minus_c840_SMN1_hgvs_C_frac_CplusT={delta_cov_minus_smn1_hgvs_c_frac}"
        )

    summary = (
        f"SMN1_cov_frac={c_ratio}|C_Ratio={c_ratio}"
        f"|Cov(1,2)={s1},{s2}|SMN1_CN_est={cn_est1}|SMN2_CN_est={cn_est2}|CNV_est_total_assumed={tc_est}"
        f"|SilentCarrier={silent_carrier}{exon7}{silent_allelic}{c840_ct}"
    )

    ct_line = line_c840_ct_readable(
        unified_used and ct_u is not None,
        ct_u,
        ct_s1,
        ct_s2,
    )
    norm_line = line_c840_norm_readable(
        unified_used and ct_u is not None,
        ct_u,
        ct_s1,
        ct_s2,
    )

    with open(path, "a") as f:
        f.write(nl + "##DARK_GENE_PIPELINE_SUMMARY" + nl)
        f.write(summary + nl)
        f.write(ct_line + nl)
        f.write(norm_line + nl)
    print(ct_line)
    print(norm_line)


if __name__ == "__main__":
    main()
