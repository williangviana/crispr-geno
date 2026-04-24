"""Core analysis: align reads with minimap2, count indels near cut sites,
call genotypes and extract dominant allele sequences."""
from __future__ import annotations

import collections
import os
import shutil
import subprocess
from dataclasses import dataclass
from typing import Optional

import pysam


@dataclass
class Guide:
    name: str
    sequence: str
    cut_pos: int = -1
    strand: str = "+"


@dataclass
class GuideCounts:
    spanning: int
    ins1: int
    del1: int
    ins_ge2: int
    del_ge2: int

    def pct(self, attr: str) -> float:
        n = getattr(self, attr)
        return 100.0 * n / self.spanning if self.spanning else 0.0


@dataclass
class AlleleCall:
    """Dominant non-WT allele at a cut site, as a readable string."""
    call: str            # 'WT', '+T', '-GATC', '+T/WT', etc., 'lowN'
    spanning: int
    top_frac: float      # % of reads supporting the dominant non-WT signature
    second_frac: float   # % of reads supporting the second-most-common non-WT signature
    wt_frac: float       # % of reads with no indel near the cut
    confidence: str      # 'high' | 'med' | 'low' | 'ambiguous' | 'lowN'


def rc(s: str) -> str:
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def find_cut(amp: str, guide: str) -> tuple[Optional[int], str]:
    """Locate a guide's Cas9 cut site in the amplicon. Returns (cut_pos, strand)."""
    if guide in amp:
        return amp.find(guide) + 17, "+"
    g_rc = rc(guide)
    if g_rc in amp:
        return amp.find(g_rc) + 3, "-"
    return None, ""


def resolve_guides(amp: str, guides: list[Guide]) -> list[Guide]:
    """Return a new list with cut_pos and strand filled. Skip ones not found."""
    resolved = []
    for g in guides:
        cut, strand = find_cut(amp, g.sequence.upper())
        if cut is None:
            print(f"  WARN: guide {g.name} ({g.sequence}) not found in amplicon, skipping")
            continue
        resolved.append(Guide(g.name, g.sequence.upper(), cut, strand))
    return resolved


def check_tools() -> None:
    for tool in ("minimap2", "samtools"):
        if shutil.which(tool) is None:
            raise RuntimeError(
                f"{tool} not found on PATH. Activate an environment with minimap2 and samtools installed."
            )


def align(ref_fa: str, fastq: str, bam_out: str, threads: int = 4) -> None:
    """Align a fastq to the amplicon with minimap2 map-ont preset, sort, index."""
    sam = bam_out + ".sam"
    with open(sam, "w") as s, open(os.devnull, "w") as null:
        subprocess.run(
            ["minimap2", "-ax", "map-ont", "-t", str(threads), ref_fa, fastq],
            stdout=s, stderr=null, check=True,
        )
    subprocess.run(
        ["samtools", "sort", "-@", str(threads), "-o", bam_out, sam],
        check=True, stderr=subprocess.DEVNULL,
    )
    subprocess.run(["samtools", "index", bam_out], check=True)
    os.remove(sam)


def count_indels(bam_path: str, ref_name: str, cut_pos: int,
                 window: int, min_read_len: int) -> GuideCounts:
    """Count reads spanning the cut and their indels, broken out by size class."""
    bam = pysam.AlignmentFile(bam_path, "rb")
    span_start, span_end = cut_pos - window, cut_pos + window
    spanning = ins1 = del1 = ins_ge2 = del_ge2 = 0
    for read in bam.fetch(ref_name, max(0, span_start), span_end):
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        if read.query_length is not None and read.query_length < min_read_len:
            continue
        if read.reference_start > span_start or read.reference_end < span_end:
            continue
        spanning += 1
        pos = read.reference_start
        r_ins1 = r_del1 = r_ins_ge2 = r_del_ge2 = False
        for op, length in read.cigartuples:
            if op in (0, 7, 8):
                pos += length
            elif op == 2:  # D
                if pos < span_end and pos + length > span_start:
                    if length == 1:
                        r_del1 = True
                    else:
                        r_del_ge2 = True
                pos += length
            elif op == 1:  # I
                if span_start <= pos <= span_end:
                    if length == 1:
                        r_ins1 = True
                    else:
                        r_ins_ge2 = True
            elif op == 3:
                pos += length
        if r_ins1: ins1 += 1
        if r_del1: del1 += 1
        if r_ins_ge2: ins_ge2 += 1
        if r_del_ge2: del_ge2 += 1
    bam.close()
    return GuideCounts(spanning, ins1, del1, ins_ge2, del_ge2)


def _read_signature(read, cut_pos: int, window: int):
    """Tuple of (('I'|'D', size, ref_pos), ...) near the cut and inserted seqs dict."""
    pos = read.reference_start
    qpos = 0
    events = []
    ins_seqs = {}
    for op, length in read.cigartuples:
        if op in (0, 7, 8):
            pos += length
            qpos += length
        elif op == 2:  # D
            if pos < cut_pos + window and pos + length > cut_pos - window:
                # keep all >=2bp; 1bp only if right at cut (suppresses Nanopore noise)
                if length >= 2 or abs(pos - cut_pos) <= 1:
                    events.append(("D", length, pos))
            pos += length
        elif op == 1:  # I
            if cut_pos - window <= pos <= cut_pos + window:
                if length >= 2 or abs(pos - cut_pos) <= 2:
                    ev = ("I", length, pos)
                    events.append(ev)
                    if read.query_sequence:
                        ins_seqs[ev] = read.query_sequence[qpos:qpos + length]
            qpos += length
        elif op == 4:
            qpos += length
        elif op == 3:
            pos += length
    return tuple(events), ins_seqs


def call_allele(bam_path: str, ref_name: str, ref_seq: str, cut_pos: int,
                window: int, min_spanning: int, homo_threshold: float,
                het_threshold: float, min_read_len: int) -> AlleleCall:
    """Extract the dominant allele sequence at a cut site."""
    bam = pysam.AlignmentFile(bam_path, "rb")
    win_s, win_e = cut_pos - window, cut_pos + window
    spanning = 0
    sig_counts: collections.Counter = collections.Counter()
    ins_seq_store: dict = {}
    for read in bam.fetch(ref_name, max(0, win_s), win_e):
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        if read.query_length is not None and read.query_length < min_read_len:
            continue
        if read.reference_start > win_s or read.reference_end < win_e:
            continue
        spanning += 1
        events, seqs = _read_signature(read, cut_pos, window)
        sig_counts[events] += 1
        for ev, s in seqs.items():
            ins_seq_store.setdefault(ev, s)
    bam.close()

    if spanning < min_spanning:
        return AlleleCall("lowN", spanning, 0.0, 0.0, 0.0, "lowN")

    wt_n = sig_counts.get((), 0)
    wt_frac = 100.0 * wt_n / spanning
    non_wt = [(sig, n) for sig, n in sig_counts.items() if sig]
    non_wt.sort(key=lambda x: -x[1])
    top_sig, top_n = (non_wt[0] if non_wt else (None, 0))
    top_frac = 100.0 * top_n / spanning if spanning else 0.0
    second_n = non_wt[1][1] if len(non_wt) >= 2 else 0
    second_frac = 100.0 * second_n / spanning if spanning else 0.0

    # Aggregate non-WT load across ALL distinct signatures (not just top two).
    # Used to catch cases where many small deletions add up to a big fraction
    # that no single signature captures.
    non_wt_total_frac = sum(100.0 * n / spanning for _, n in non_wt) if spanning else 0.0

    def confidence_of(call_frac: float, call_is_het_vs_wt: bool = False) -> str:
        # Ambiguous rule 1: a second single signature is a significant competitor.
        if second_frac >= 15.0:
            return "ambiguous"
        # Ambiguous rule 2: het/WT call but WT reads are near-absent (the "WT half"
        # is really a pile of different mutant alleles).
        if call_is_het_vs_wt and wt_frac < 10.0:
            return "ambiguous"
        # Ambiguous rule 3: homozygous call but lots of other non-WT noise.
        # If the call is the dominant allele and the non-WT total exceeds it
        # by >=20 points, the rest is probably meaningful structure.
        if call_frac and (non_wt_total_frac - call_frac) >= 20.0:
            return "ambiguous"
        if spanning < 50:
            return "low"
        # Homo calls barely crossing the threshold
        if call_frac and call_frac < 80.0:
            return "low"
        if spanning >= 100 and call_frac >= 90.0:
            return "high"
        return "med"

    # WT case: no non-WT signal, or top non-WT under het threshold
    if not non_wt or top_frac < het_threshold:
        # For WT, confidence comes from how clean the WT signal is.
        if second_frac >= 15.0:
            conf = "ambiguous"
        elif wt_frac < 80.0:
            # Called WT but the wildtype reads are a minority — aggregate editing
            # load is spread across many alleles (mosaic / still-active Cas9).
            conf = "ambiguous"
        elif spanning < 50:
            conf = "low"
        elif spanning >= 100 and wt_frac >= 90.0:
            conf = "high"
        else:
            conf = "med"
        return AlleleCall("WT", spanning, top_frac, second_frac, wt_frac, conf)

    def event_str(ev):
        t, sz, p = ev
        if t == "I":
            seq = ins_seq_store.get(ev, "N" * sz)
            return f"+{seq}"
        else:
            seq = ref_seq[p:p + sz]
            return f"-{seq}"

    def sig_str(sig):
        return ";".join(event_str(e) for e in sig)

    top_str = sig_str(top_sig)

    # Homozygous mutant
    if top_frac >= homo_threshold:
        return AlleleCall(top_str, spanning, top_frac, second_frac, wt_frac,
                          confidence_of(top_frac))

    # Biallelic: two different mutant alleles each above het threshold, WT minor
    if len(non_wt) >= 2:
        second_sig, _ = non_wt[1]
        if second_frac >= het_threshold and wt_frac < het_threshold:
            return AlleleCall(f"{top_str}/{sig_str(second_sig)}",
                              spanning, top_frac, second_frac, wt_frac,
                              confidence_of(top_frac))

    # Heterozygous vs WT
    return AlleleCall(f"{top_str}/WT", spanning, top_frac, second_frac, wt_frac,
                      confidence_of(top_frac, call_is_het_vs_wt=True))


def find_large_deletions(
    bam_path: str,
    ref_name: str,
    guides: list[Guide],
    min_size: int = 50,
    breakpoint_window: int = 20,
    min_read_len: int = 200,
    homo_threshold: float = 70.0,
    het_threshold: float = 20.0,
    min_supporting: int = 5,
) -> dict[tuple[int, int], dict]:
    """Find large deletions that span two cut sites (structural variants).

    Scans all primary alignments for deletion CIGAR ops ≥ min_size whose
    start AND end fall within breakpoint_window bp of two different guide
    cut sites. Returns a dict keyed by (guide_i_index, guide_j_index) where
    i < j, mapping to:
      {
        'supporting_reads': int,      # reads carrying this SV
        'spanning_reads':  int,       # reads spanning both cut sites
                                       # (so we can compute a fraction)
        'median_size':     int,
        'call':            str,       # 'SV:-<size>bp' | 'SV:-<size>bp/WT' | 'WT' | 'lowN'
        'top_frac':        float,
      }
    Only pairs with ≥ min_supporting reads are returned.
    """
    bam = pysam.AlignmentFile(bam_path, "rb")

    # supporting reads per pair: key = (i, j), value = list of deletion sizes
    sv_sizes: dict[tuple[int, int], list[int]] = {}
    # spanning reads per pair: reads whose alignment covers both cut sites
    # (needed to compute the correct denominator)
    span_reads: dict[tuple[int, int], int] = {}
    # pre-compute all pairs
    cuts = [(i, g.cut_pos) for i, g in enumerate(guides) if g.cut_pos >= 0]
    pairs = [
        (i, j, ci, cj) for (i, ci) in cuts for (j, cj) in cuts
        if i < j
    ]

    for read in bam.fetch(ref_name):
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        if read.query_length is not None and read.query_length < min_read_len:
            continue

        # which cut-site pairs does this read *span*?
        r_start = read.reference_start
        r_end = read.reference_end
        for i, j, ci, cj in pairs:
            lo, hi = (ci, cj) if ci < cj else (cj, ci)
            # read must cover from upstream breakpoint-window to downstream breakpoint-window
            if r_start <= lo - breakpoint_window and r_end >= hi + breakpoint_window:
                span_reads[(i, j)] = span_reads.get((i, j), 0) + 1

        # scan CIGAR for large deletions
        pos = read.reference_start
        for op, length in read.cigartuples:
            if op in (0, 7, 8):
                pos += length
            elif op == 2:  # D
                if length >= min_size:
                    del_start = pos
                    del_end = pos + length
                    for i, j, ci, cj in pairs:
                        lo_c, hi_c = (ci, cj) if ci < cj else (cj, ci)
                        if (abs(del_start - lo_c) <= breakpoint_window and
                                abs(del_end - hi_c) <= breakpoint_window):
                            sv_sizes.setdefault((i, j), []).append(length)
                pos += length
            elif op == 1:
                pass  # I doesn't consume ref
            elif op == 3:
                pos += length

    bam.close()

    results: dict[tuple[int, int], dict] = {}
    for pair, sizes in sv_sizes.items():
        n = len(sizes)
        if n < min_supporting:
            continue
        span = span_reads.get(pair, 0)
        median = sorted(sizes)[len(sizes) // 2]
        if span < 20:
            call = "lowN"
            frac = 0.0
        else:
            frac = 100.0 * n / span
            if frac >= homo_threshold:
                call = f"SV:-{median}bp"
            elif frac >= het_threshold:
                call = f"SV:-{median}bp/WT"
            else:
                call = "WT"
        results[pair] = {
            "supporting_reads": n,
            "spanning_reads": span,
            "median_size": median,
            "call": call,
            "top_frac": frac,
        }
    return results


def classify_counts(c: GuideCounts, homo_threshold: float, het_threshold: float,
                    min_spanning: int) -> str:
    """Call a coarse genotype from aggregate indel percentages."""
    if c.spanning < min_spanning:
        return "lowN"
    ins1 = c.pct("ins1")
    delge2 = c.pct("del_ge2")
    ins_homo = ins1 >= homo_threshold
    ins_het = ins1 >= het_threshold
    del_homo = delge2 >= homo_threshold
    del_het = delge2 >= het_threshold
    if ins_homo and del_homo: return "T+del"
    if ins_homo: return "T"
    if del_homo: return "del"
    if ins_het and del_het: return "T/+ & del/+"
    if ins_het: return "T/+"
    if del_het: return "del/+"
    return "WT"
