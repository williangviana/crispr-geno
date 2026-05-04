"""Read-level haplotype phasing for long-read CRISPR amplicon data.

For each primary read in a BAM, derive a per-read haplotype tuple of
(allele_at_g1, allele_at_g2, ...). Cluster reads by tuple to recover
the underlying chromosomes, then make a zygosity call per sample
(WT / homozygous / heterozygous / biallelic / mosaic / lowN).

This is a separate pass from the per-guide caller in analysis.py. It is
the substrate for the consolidated zygosity column in results.xlsx.
"""
from __future__ import annotations

import collections
import re
from dataclasses import dataclass, field
from typing import Optional

import pysam

from .analysis import AlleleCall, Guide, read_signature


# ---- Types --------------------------------------------------------------

# Slot label vocabulary used as the canonical clustering key.
#   "WT"     no edit at this guide
#   "+K"     K-bp insertion (e.g. "+1", "+5")
#   "-K"     K-bp deletion  (e.g. "-1", "-3")
#   "SV"     this guide is inside a detected structural variant
#   "NC"     read does not span this guide's window — cannot phase
GuideAllele = str
HaplotypeTuple = tuple[GuideAllele, ...]


@dataclass(frozen=True)
class SlotDetail:
    """Per-slot detail kept alongside the slot label, for rendering."""
    label: str
    size: int = 0
    ins_seq: str = ""
    ref_pos: int = -1


@dataclass
class Haplotype:
    """One inferred chromosome."""
    rank: int                       # 1 = top, 2 = second
    tuple_: HaplotypeTuple
    details: list[SlotDetail]
    support_reads: int
    support_frac: float             # 0..1 of phased reads
    description: str                # e.g. "SV:-5411bp" or "gRNA2:+1bp"


@dataclass
class PhasingResult:
    sample: str
    n_phased: int                   # reads contributing to a non-NC tuple
    n_partial: int                  # reads dropped because of NC slots
    haplotypes: list[Haplotype]     # 0, 1, or 2 entries
    zygosity_call: str              # see _zygosity_call()
    is_mosaic: bool                 # >=3 dominant clusters
    notes: list[str]
    all_tuples: list[tuple[HaplotypeTuple, int]]  # post-merge, sorted desc


# ---- Per-read classification --------------------------------------------

def classify_read_at_guide(read, cut_pos: int, window: int,
                           ref_seq: str | None = None) -> SlotDetail:
    """Classify one read's allele at one guide's cut site.

    Returns 'NC' if the read does not span the cut window (cannot phase
    this slot). Otherwise reuses read_signature() to find the largest
    indel event near the cut. When ref_seq is provided, indels are
    left-normalized so equivalent placements (e.g. +T at pos N vs N+1
    inside a TTTT homopolymer) cluster together.
    """
    if read.reference_start > cut_pos - window or read.reference_end < cut_pos + window:
        return SlotDetail("NC")
    events, ins_seqs = read_signature(read, cut_pos, window, ref_seq)
    if not events:
        return SlotDetail("WT")
    # If multiple events near the same cut, keep the largest by size.
    largest = max(events, key=lambda e: e[1])
    t, sz, p = largest
    sign = "+" if t == "I" else "-"
    label = f"{sign}{sz}"
    ins_seq = ins_seqs.get(largest, "") if t == "I" else ""
    return SlotDetail(label, sz, ins_seq, p)


def _guide_in_sv_span(guide_idx: int, sv_pair: tuple[int, int],
                      guides: list[Guide]) -> bool:
    """A guide is inside the SV span if its cut position lies between the
    pair's two cut positions (inclusive of both endpoints — the endpoint
    guides are themselves at the deletion breakpoints)."""
    i, j = sv_pair
    lo = min(guides[i].cut_pos, guides[j].cut_pos)
    hi = max(guides[i].cut_pos, guides[j].cut_pos)
    return lo <= guides[guide_idx].cut_pos <= hi


def _build_read_to_sv(sv_results: dict[tuple[int, int], dict]) -> dict[str, tuple[int, int]]:
    """Map each SV-supporting query_name to its (i,j) pair. If a read supports
    more than one SV pair, keep the largest deletion."""
    read_to_sv: dict[str, tuple[int, int]] = {}
    for pair, info in sv_results.items():
        for qn in info.get("supporting_read_names", ()):
            existing = read_to_sv.get(qn)
            if existing is None or info["median_size"] > sv_results[existing]["median_size"]:
                read_to_sv[qn] = pair
    return read_to_sv


def build_read_haplotypes(
    bam_path: str,
    ref_name: str,
    guides: list[Guide],
    sv_results: dict[tuple[int, int], dict],
    window: int,
    min_read_len: int,
    imputation: list[Optional[SlotDetail]] | None = None,
    ref_seq: str | None = None,
    min_mapq: int = 0,
) -> tuple[collections.Counter, dict[HaplotypeTuple, list[SlotDetail]], int, int]:
    """Single pass over primary reads, producing:
      - a Counter of HaplotypeTuple
      - a dict of HaplotypeTuple -> representative per-slot details (first seen)
      - the number of reads dropped due to NC slots (partial coverage)
      - the number of reads rescued from partial-coverage by imputation

    `imputation` is an optional per-guide pre-imputed SlotDetail. When a
    read has NC at slot i and imputation[i] is provided, the NC is filled
    with imputation[i] instead of dropping the read. This rescues reads
    that span only some of the guides — common for long amplicons or
    reads from a chromosome with a large lesion that prevents spanning.

    Imputation is only performed at slots where the per-guide caller
    is unambiguously confident (computed by `per_guide_imputation()`),
    so we never fabricate evidence at uncertain sites.
    """
    read_to_sv = _build_read_to_sv(sv_results)
    use_impute = imputation is not None

    bam = pysam.AlignmentFile(bam_path, "rb")
    counter: collections.Counter = collections.Counter()
    details_store: dict[HaplotypeTuple, list[SlotDetail]] = {}
    n_partial = 0
    n_imputed = 0
    try:
        for read in bam.fetch(ref_name):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            if read.mapping_quality < min_mapq:
                continue
            if read.query_length is not None and read.query_length < min_read_len:
                continue
            sv_pair = read_to_sv.get(read.query_name)
            details: list[SlotDetail] = []
            had_nc = False
            for i, g in enumerate(guides):
                if sv_pair is not None and _guide_in_sv_span(i, sv_pair, guides):
                    details.append(SlotDetail("SV"))
                else:
                    d = classify_read_at_guide(read, g.cut_pos, window, ref_seq)
                    if d.label == "NC":
                        if use_impute and imputation[i] is not None:
                            details.append(imputation[i])
                            had_nc = True
                        else:
                            details.append(d)
                    else:
                        details.append(d)
            if any(d.label == "NC" for d in details):
                n_partial += 1
                continue
            if had_nc:
                n_imputed += 1
            tup = tuple(d.label for d in details)
            counter[tup] += 1
            details_store.setdefault(tup, details)
    finally:
        bam.close()
    return counter, details_store, n_partial, n_imputed


def per_guide_imputation(
    plant_calls: list[AlleleCall],
    guides: list[Guide],
    ref_seq: str,
) -> list[Optional[SlotDetail]]:
    """Per-guide SlotDetail to use when imputing NC slots in reads.

    Returns one entry per guide. Conservative: only impute at guides
    where the per-guide caller is **high confidence** AND the call is
    a single allele (not a het, biallelic, or lowN). For high-confidence
    WT calls we impute "WT"; for high-confidence single-edit calls we
    impute the corresponding phased slot label.

    Anything ambiguous (het, biallelic, mosaic-like, low coverage) gets
    `None` — those slots stay NC and reads NC at them stay dropped.
    """
    out: list[Optional[SlotDetail]] = []
    for i, ac in enumerate(plant_calls):
        if ac.confidence != "high":
            out.append(None)
            continue
        call = ac.call.strip()
        if call == "WT":
            out.append(SlotDetail("WT"))
            continue
        if "/" in call or call in ("lowN", ""):
            out.append(None)
            continue
        # Single-allele edit call (e.g. "+T", "-CCC") — convert to a
        # phased slot label and synthesize a SlotDetail. ref_pos is set
        # to -1 (unknown) so the renderer falls back to size-only
        # notation ('+Tbp', '-3bp') if no observed read populates the
        # tuple first; we don't fabricate a position that might point
        # at the wrong reference bases.
        phased = _per_guide_label_to_phased(call)
        if phased is None:
            out.append(None)
            continue
        size = int(phased[1:])
        ins_seq = call[1:] if call.startswith("+") else ""
        out.append(SlotDetail(phased, size, ins_seq, -1))
    return out


# ---- Noise merging ------------------------------------------------------

def _hamming(a: HaplotypeTuple, b: HaplotypeTuple) -> int:
    return sum(1 for x, y in zip(a, b) if x != y)


def merge_noise_tuples(
    counter: collections.Counter,
    noise_frac: float = 0.05,
    noise_min_reads: int = 2,
) -> collections.Counter:
    """Merge minor tuples into the nearest dominant tuple by Hamming distance.
    Minor = support < max(noise_min_reads, noise_frac * total).
    Ties broken by support of the candidate dominant tuple.
    """
    if not counter:
        return counter
    total = sum(counter.values())
    threshold = max(noise_min_reads, int(noise_frac * total))
    dominants = [t for t, n in counter.items() if n >= threshold]
    minors = [t for t, n in counter.items() if n < threshold]
    if not dominants:
        return counter  # nothing to merge into
    merged: collections.Counter = collections.Counter(
        {t: counter[t] for t in dominants}
    )
    for m in minors:
        # nearest dominant by Hamming, tie-break on highest support
        best = min(dominants,
                   key=lambda d: (_hamming(m, d), -counter[d]))
        merged[best] += counter[m]
    return merged


# ---- Haplotype picking --------------------------------------------------

def pick_haplotypes(
    merged: collections.Counter,
    het_threshold_frac: float = 0.20,
    mosaic_threshold_frac: float = 0.10,
) -> tuple[list[HaplotypeTuple], bool]:
    """Pick the top 1-2 haplotypes from the merged Counter.

    Returns (top_tuples, is_mosaic). top_tuples has length 0, 1, or 2.
    is_mosaic is True if either 3+ tuples cross mosaic_threshold_frac, or
    one chromosome reads as WT and the non-WT mass is fragmented across
    3+ distinct lesions with no single edited tuple clearing the het
    threshold (mosaic on one chromosome).
    """
    if not merged:
        return [], False
    total = sum(merged.values())
    ranked = merged.most_common()

    top = [(t, n) for t, n in ranked if n / total >= het_threshold_frac]
    is_mosaic = sum(1 for _, n in ranked if n / total >= mosaic_threshold_frac) >= 3

    # Mosaic-on-one-chromosome: WT reads cleanly on one chromosome but the
    # other chromosome's reads are fragmented across many distinct lesions
    # (no single edited haplotype clears the het threshold). The 3-above-
    # threshold check above misses this — typically only WT and the tallest
    # edited tuple cross it. Detect it from the non-WT mass directly.
    if not is_mosaic:
        non_wt = [(t, n) for t, n in ranked if not all(s == "WT" for s in t)]
        if non_wt:
            non_wt_total = sum(n for _, n in non_wt)
            non_wt_top = max(n for _, n in non_wt)
            n_distinct = sum(1 for _, n in non_wt
                             if n / total >= mosaic_threshold_frac / 2)
            if (non_wt_total / total >= het_threshold_frac
                    and non_wt_top / total < het_threshold_frac
                    and n_distinct >= 3):
                is_mosaic = True

    if not top:
        return [], is_mosaic
    if len(top) == 1:
        return [top[0][0]], is_mosaic
    # Two or more clear haplotypes; take the top two.
    return [top[0][0], top[1][0]], is_mosaic


# ---- Zygosity call ------------------------------------------------------

def _zygosity_call(
    tuples: list[HaplotypeTuple], is_homozygous: bool, is_mosaic: bool,
) -> str:
    """Synthesize the sample-level zygosity call from inferred haplotypes.

    Vocabulary (DNA-level, no protein-function inference):
      WT           — both chromosomes wildtype at every guide
      homozygous   — both chromosomes carry the SAME edit
      biallelic    — both chromosomes are edited but with DIFFERENT lesions
      heterozygous — one chromosome is edited, the other is wildtype
      mosaic       — 3+ haplotypes detected (Cas9 still active or chimera)
      lowN         — fewer than --min-spanning reads phased
    """
    if is_mosaic:
        return "mosaic"
    if len(tuples) == 0:
        return "lowN"
    if len(tuples) == 1:
        return "WT" if all(s == "WT" for s in tuples[0]) else "homozygous"
    a_wt = all(s == "WT" for s in tuples[0])
    b_wt = all(s == "WT" for s in tuples[1])
    if a_wt and b_wt:
        return "WT"
    if a_wt or b_wt:
        return "heterozygous"
    if is_homozygous:
        return "homozygous"
    return "biallelic"


# ---- Allele description rendering ---------------------------------------

def _slot_description(detail: SlotDetail, guide_name: str, ref_seq: str) -> str:
    """Render a single guide's edit on one chromosome.

    'WT' returns ''. SV slots aren't rendered here — collapsed at the
    haplotype level by render_haplotype_description.
    """
    if detail.label in ("WT", "NC", "SV"):
        return ""
    if detail.label.startswith("+"):
        if detail.ins_seq:
            return f"{guide_name}:+{detail.ins_seq}"
        return f"{guide_name}:{detail.label}bp"
    if detail.label.startswith("-"):
        if 0 <= detail.ref_pos and detail.size > 0:
            seq = ref_seq[detail.ref_pos:detail.ref_pos + detail.size]
            if seq:
                return f"{guide_name}:-{seq}"
        return f"{guide_name}:{detail.label}bp"
    return f"{guide_name}:{detail.label}"


def render_haplotype_description(
    tuple_: HaplotypeTuple,
    details: list[SlotDetail],
    guides: list[Guide],
    ref_seq: str,
    sv_results: dict[tuple[int, int], dict],
) -> str:
    """Compose a chromosome's edit set into a single human-readable string.

    Examples: 'WT', 'SV:-5411bp', 'gRNA2:+T', 'gRNA1:+T, gRNA3:-GATC'.
    """
    sv_indices = {i for i, s in enumerate(tuple_) if s == "SV"}
    if sv_indices:
        # Find the SV pair that explains these indices. Prefer the pair that
        # spans exactly the SV indices; otherwise take the largest deletion.
        chosen: Optional[tuple[int, int]] = None
        for pair, info in sv_results.items():
            if all(_guide_in_sv_span(i, pair, guides) for i in sv_indices):
                if chosen is None or info["median_size"] > sv_results[chosen]["median_size"]:
                    chosen = pair
        if chosen is None:
            # Should not happen if we labelled SV from sv_results, but degrade gracefully.
            return f"SV (guides {','.join(guides[i].name for i in sorted(sv_indices))})"
        size = sv_results[chosen]["median_size"]
        sv_str = f"SV:-{size}bp"
        # Render any non-SV edits on this chromosome alongside the SV label.
        extras = []
        for i, d in enumerate(details):
            if d.label in ("WT", "NC", "SV"):
                continue
            extras.append(_slot_description(d, guides[i].name, ref_seq))
        if extras:
            return sv_str + ", " + ", ".join(extras)
        return sv_str

    parts = []
    for i, d in enumerate(details):
        s = _slot_description(d, guides[i].name, ref_seq)
        if s:
            parts.append(s)
    return ", ".join(parts) if parts else "WT"


def _slot_short(detail: SlotDetail, ref_seq: str) -> str:
    """Compact per-guide cell for the per-guide A1/A2 columns ('SV', '+T', '-GATC', 'WT')."""
    if detail.label == "SV":
        return "SV"
    if detail.label == "WT":
        return "WT"
    if detail.label == "NC":
        return "?"
    if detail.label.startswith("+"):
        if detail.ins_seq:
            return f"+{detail.ins_seq}"
        return detail.label + "bp"
    if detail.label.startswith("-"):
        if 0 <= detail.ref_pos and detail.size > 0:
            seq = ref_seq[detail.ref_pos:detail.ref_pos + detail.size]
            if seq:
                return f"-{seq}"
        return detail.label + "bp"
    return detail.label


# ---- Orchestrator -------------------------------------------------------

def phase_sample(
    sample_id: str,
    bam_path: str,
    ref_name: str,
    ref_seq: str,
    guides: list[Guide],
    sv_results: dict[tuple[int, int], dict],
    window: int,
    min_spanning: int,
    min_read_len: int,
    noise_frac: float = 0.05,
    noise_min_reads: int = 2,
    het_threshold_frac: float = 0.20,
    mosaic_threshold_frac: float = 0.10,
    imputation: list[Optional[SlotDetail]] | None = None,
    min_mapq: int = 0,
) -> PhasingResult:
    """Run the full phasing pass for one sample.

    `imputation` (optional, one entry per guide) lets reads that don't
    span every guide still be phased, by filling in NC slots at sites
    where the per-guide caller is unambiguously confident. See
    `per_guide_imputation()` for how to derive it.
    """
    counter, details_store, n_partial, n_imputed = build_read_haplotypes(
        bam_path, ref_name, guides, sv_results, window, min_read_len,
        imputation=imputation, ref_seq=ref_seq, min_mapq=min_mapq,
    )
    n_phased = sum(counter.values())
    notes: list[str] = []
    if n_partial:
        notes.append(f"partial-coverage:{n_partial}")
    if n_imputed:
        notes.append(f"imputed:{n_imputed}")

    if n_phased < min_spanning:
        return PhasingResult(
            sample=sample_id, n_phased=n_phased, n_partial=n_partial,
            haplotypes=[], zygosity_call="lowN", is_mosaic=False,
            notes=notes, all_tuples=[],
        )

    merged = merge_noise_tuples(counter, noise_frac, noise_min_reads)
    n_merged_in = sum(counter.values()) - sum(
        v for k, v in counter.items() if k in merged
    )
    if n_merged_in:
        notes.append(f"noise-merged:{n_merged_in}")

    top_tuples, is_mosaic = pick_haplotypes(
        merged, het_threshold_frac, mosaic_threshold_frac,
    )
    # If the noise-merge step absorbed a large fraction of phased reads,
    # the pre-merge counter held many distinct minor tuples (each below
    # 5%) — i.e. real haplotype diversity, not basecalling drift. The
    # post-merge tuple-counting can't see that anymore, so flag mosaic
    # directly off the merged-ratio.
    if n_phased and n_merged_in / n_phased > 0.30:
        is_mosaic = True
    if is_mosaic:
        notes.append("mosaic")

    haplotypes: list[Haplotype] = []
    total = sum(merged.values())
    for rank, tup in enumerate(top_tuples, start=1):
        details = details_store[tup]
        support = merged[tup]
        haplotypes.append(Haplotype(
            rank=rank,
            tuple_=tup,
            details=details,
            support_reads=support,
            support_frac=support / total if total else 0.0,
            description=render_haplotype_description(
                tup, details, guides, ref_seq, sv_results,
            ),
        ))

    is_homozygous = (
        len(haplotypes) == 2 and haplotypes[0].tuple_ == haplotypes[1].tuple_
    )
    zygosity_call = _zygosity_call(
        [h.tuple_ for h in haplotypes], is_homozygous, is_mosaic,
    )

    all_tuples = sorted(merged.items(), key=lambda kv: -kv[1])

    return PhasingResult(
        sample=sample_id, n_phased=n_phased, n_partial=n_partial,
        haplotypes=haplotypes, zygosity_call=zygosity_call,
        is_mosaic=is_mosaic, notes=notes, all_tuples=all_tuples,
    )


# ---- Per-guide A1/A2 cells (for the primary results.xlsx sheet) --------

def per_guide_cells(result: PhasingResult, guides: list[Guide], ref_seq: str) -> list[str]:
    """Return one cell per guide describing what each chromosome carries.

    Homozygous (both chromosomes identical) collapses to a single label
    ('+T', 'WT', '-CC') matching the legacy notation. Heterozygous keeps
    both halves separated by ' / ' ('SV / +T', '+T / WT').
    """
    cells = []
    for gi, _g in enumerate(guides):
        if not result.haplotypes:
            cells.append("")
            continue
        h1 = result.haplotypes[0]
        h2 = result.haplotypes[1] if len(result.haplotypes) >= 2 else h1
        a = _slot_short(h1.details[gi], ref_seq)
        b = _slot_short(h2.details[gi], ref_seq)
        cells.append(a if a == b else f"{a} / {b}")
    return cells


# ---- Cross-check: per-guide caller vs phased haplotypes ----------------

_INDEL_RE = re.compile(r"^([+-])([ACGT]+)$", re.IGNORECASE)


def _per_guide_label_to_phased(label: str) -> str | None:
    """Map a per-guide caller's allele label ('+T', '-CCC') to a phased
    slot label ('+1', '-3'). Returns None for SV / WT / unparseable."""
    m = _INDEL_RE.match(label)
    if m:
        sign, seq = m.groups()
        return f"{sign}{len(seq)}"
    return None


def residual_non_wt_phased(
    result: "PhasingResult",
    slot_signals: list[float] | None = None,
    slot_noise_floors: list[float] | None = None,
) -> tuple[float, int]:
    """(percent, absolute reads) of phased reads carrying a haplotype
    that introduces a *new* edit not present in any called chromosome.

    A "new edit" means a slot where the residual haplotype is non-WT
    but every called chromosome is WT. Variants of an existing called
    edit (e.g. a +1 insertion at a homopolymer site getting basecalled
    as WT in some reads, producing a residual that drops the +1) do
    NOT count — they reflect per-base error around a fixed allele,
    not new editing activity.

    Optional cross-sample noise floors: if slot_signals (this sample's
    per-guide top non-WT frac) and slot_noise_floors (per-guide batch
    background, derived from confidently-WT samples) are both provided,
    a "new edit" at slot i only counts when slot_signals[i] is above
    slot_noise_floors[i]. This suppresses false positives at sites with
    consistent batch noise — e.g. a homopolymer-prone guide showing 5%
    del1 in every WT sample shouldn't mark a residual as Cas9 activity.

    Examples (without floors):
      called=(WT,WT,WT,WT), residual=(-1,WT,-24,WT)
        → slots 0 and 2 carry edits the call lacks → counts (real Cas9+)

      called=(WT,-2,WT,+1), residual=(WT,-2,WT,WT)
        → only differs at slot 3 by *removing* +1 → no new edit → noise

      called=(WT,-2,WT,+1), residual=(WT,-2,WT,+2)
        → slot 3 is non-WT in both call and residual (variant of +1)
        → no new edit → noise

    Returns the absolute read count alongside the percent so callers
    can require both a fractional bar AND a minimum absolute support
    — at low n_phased a 5% bar can be tripped by 2-3 noise reads.
    """
    if not result.n_phased:
        return 0.0, 0
    called_tuples = [h.tuple_ for h in result.haplotypes]
    called_set = set(called_tuples)
    use_floors = slot_signals is not None and slot_noise_floors is not None

    def has_new_edit(tup: HaplotypeTuple) -> bool:
        for i, allele in enumerate(tup):
            if allele == "WT":
                continue
            # New iff every called chromosome is WT at this slot. With
            # no called chromosomes (lowN edge case) this is vacuously
            # true, so any non-WT residual still counts.
            if not all(c[i] == "WT" for c in called_tuples):
                continue
            # Cross-sample noise gate: even a "new edit" slot is
            # discarded if its per-guide signal in this sample is
            # within batch noise.
            if use_floors and i < len(slot_signals) and i < len(slot_noise_floors):
                if slot_signals[i] <= slot_noise_floors[i]:
                    continue
            return True
        return False

    residual = sum(
        n for tup, n in result.all_tuples
        if tup not in called_set and has_new_edit(tup)
    )
    return residual / result.n_phased * 100.0, residual


def detect_phasing_mismatch(
    plant_calls: list[AlleleCall],
    result: PhasingResult,
    guides: list[Guide],
    threshold_pct: float = 15.0,
) -> list[str]:
    """Find guides where the per-guide caller saw a non-WT, non-SV allele
    above threshold_pct that no phased haplotype represents.

    The phasing pass needs every read to span every guide cut window. If
    one chromosome carries a large SV, reads spanning the SV breakpoint
    trivially cover all guides — but the *other* chromosome's reads may
    need to span the full amplicon and often get dropped as partial-
    coverage. The per-guide caller doesn't have that constraint, so it
    can see edits the phased view missed. This helper surfaces that gap
    so the call doesn't silently collapse to a single-haplotype answer.

    Returns one note per affected guide, e.g.
    'phasing-mismatch:gRNA2:+T@41%'.
    """
    if not result.haplotypes:
        return []
    n = len(result.haplotypes[0].tuple_)
    notes: list[str] = []
    for i in range(min(n, len(plant_calls), len(guides))):
        ac = plant_calls[i]
        if ac.confidence == "lowN":
            continue
        if ac.top_frac < threshold_pct:
            continue
        slots = {h.tuple_[i] for h in result.haplotypes}
        # Try to attribute the mismatch to a specific allele from the call
        found = False
        for part in ac.call.split("/"):
            part = part.strip()
            if part in ("WT", "", "lowN") or part.startswith("SV"):
                continue
            phased_eq = _per_guide_label_to_phased(part)
            if phased_eq is None or phased_eq in slots:
                continue
            # Same-lesion guard: a guide inside an SV span shows 'SV' in
            # the phased view, but the per-guide caller (which doesn't
            # know about SV) sees the same deletion as a long indel.
            # Skip when phased has SV here and the per-guide indel is
            # large enough to be the SV itself.
            if "SV" in slots and phased_eq.startswith("-"):
                try:
                    if int(phased_eq[1:]) >= 50:
                        continue
                except ValueError:
                    pass
            notes.append(
                f"phasing-mismatch:{guides[i].name}:{part}@{ac.top_frac:.0f}%"
            )
            found = True
            break
        if found:
            continue
        # Generic fallback only when the call truly has no specific
        # non-WT non-SV allele to attribute the signal to. If the call
        # named an allele but the loop above skipped it (e.g. it matched
        # a phased slot, or it was the SV showing through as a long
        # indel), the signal is already explained — don't emit a
        # spurious generic 'edit' note on top.
        has_named_allele = any(
            p.strip() not in ("WT", "", "lowN") and not p.strip().startswith("SV")
            for p in ac.call.split("/")
        )
        if has_named_allele:
            continue
        if all(s in ("WT", "SV") for s in slots):
            notes.append(
                f"phasing-mismatch:{guides[i].name}:edit@{ac.top_frac:.0f}%"
            )
    return notes
