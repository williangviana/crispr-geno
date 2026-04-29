"""Unit tests for the pure-Python phasing logic.

The BAM-traversal part of phasing is exercised end-to-end against the
crd2 sample in the validation pass; these tests cover the pieces that
can be reasoned about without real read data.
"""
from __future__ import annotations

import collections

from crispr_geno.analysis import AlleleCall, Guide
from crispr_geno.phasing import (
    SlotDetail,
    _slot_description,
    _zygosity_call,
    detect_phasing_mismatch,
    merge_noise_tuples,
    per_guide_cells,
    PhasingResult,
    Haplotype,
    pick_haplotypes,
    render_haplotype_description,
)


# ---- merge_noise_tuples -------------------------------------------------

def test_merge_noise_tuples_collapses_minor_into_nearest_dominant():
    counter = collections.Counter({
        ("WT", "+1", "WT", "WT"): 100,        # dominant haplotype A
        ("SV", "SV", "SV", "SV"): 90,         # dominant haplotype B
        ("WT", "+2", "WT", "WT"): 1,          # minor — closest to A
        ("SV", "SV", "WT", "SV"): 1,          # minor — closest to B
    })
    merged = merge_noise_tuples(counter, noise_frac=0.05, noise_min_reads=2)
    assert merged[("WT", "+1", "WT", "WT")] == 101
    assert merged[("SV", "SV", "SV", "SV")] == 91
    # Minor tuples should not appear as separate entries.
    assert ("WT", "+2", "WT", "WT") not in merged
    assert ("SV", "SV", "WT", "SV") not in merged


def test_merge_noise_tuples_no_dominants_returns_input_unchanged():
    counter = collections.Counter({
        ("WT", "WT"): 1,
        ("SV", "SV"): 1,
    })
    # No tuple meets noise_min_reads=2, so no merging happens.
    merged = merge_noise_tuples(counter, noise_frac=0.05, noise_min_reads=2)
    assert merged == counter


# ---- pick_haplotypes ----------------------------------------------------

def test_pick_haplotypes_returns_two_when_both_above_het_threshold():
    merged = collections.Counter({
        ("SV", "SV"): 70,
        ("WT", "+1"): 30,
    })
    top, mosaic = pick_haplotypes(merged, het_threshold_frac=0.20,
                                  mosaic_threshold_frac=0.15)
    assert ("SV", "SV") in top
    assert ("WT", "+1") in top
    assert mosaic is False


def test_pick_haplotypes_flags_mosaic_when_three_clusters_pass_threshold():
    merged = collections.Counter({
        ("WT",): 35,
        ("+1",): 35,
        ("-1",): 30,
    })
    top, mosaic = pick_haplotypes(merged, het_threshold_frac=0.20,
                                  mosaic_threshold_frac=0.15)
    assert mosaic is True
    # Top-2 still returned even when mosaic.
    assert len(top) == 2


def test_pick_haplotypes_returns_one_when_homozygous():
    merged = collections.Counter({
        ("+1", "+1"): 100,
        ("WT", "WT"): 5,  # below het threshold
    })
    top, mosaic = pick_haplotypes(merged, het_threshold_frac=0.20,
                                  mosaic_threshold_frac=0.15)
    assert top == [("+1", "+1")]
    assert mosaic is False


# ---- _zygosity_call -----------------------------------------------------

def test_zygosity_homozygous_when_two_identical_edited_haplotypes():
    tup = ("WT", "+1", "WT", "WT")
    assert _zygosity_call([tup, tup], is_homozygous=True, is_mosaic=False) == "homozygous"


def test_zygosity_homozygous_when_inframe_indel_on_both_chromosomes():
    # The motivating case: -6/WT/-3 on both chromosomes is homozygous,
    # regardless of frameshift vs in-frame.
    tup = ("-6", "WT", "-3")
    assert _zygosity_call([tup, tup], is_homozygous=True, is_mosaic=False) == "homozygous"


def test_zygosity_biallelic_when_two_different_edited_haplotypes():
    a = ("WT", "+1", "WT", "WT")
    b = ("SV", "SV", "SV", "SV")
    assert _zygosity_call([a, b], is_homozygous=False, is_mosaic=False) == "biallelic"


def test_zygosity_heterozygous_when_one_haplotype_is_all_wt():
    edited = ("WT", "+1", "WT", "WT")
    wt = ("WT", "WT", "WT", "WT")
    assert _zygosity_call([edited, wt], is_homozygous=False, is_mosaic=False) == "heterozygous"
    assert _zygosity_call([wt, edited], is_homozygous=False, is_mosaic=False) == "heterozygous"


def test_zygosity_wt_when_both_haplotypes_are_all_wt():
    wt = ("WT", "WT", "WT", "WT")
    assert _zygosity_call([wt, wt], is_homozygous=True, is_mosaic=False) == "WT"


def test_zygosity_single_haplotype_treated_as_homozygous():
    edited = ("WT", "+1", "WT", "WT")
    wt = ("WT", "WT", "WT", "WT")
    assert _zygosity_call([edited], is_homozygous=False, is_mosaic=False) == "homozygous"
    assert _zygosity_call([wt], is_homozygous=False, is_mosaic=False) == "WT"


def test_zygosity_lowN_when_no_haplotypes():
    assert _zygosity_call([], is_homozygous=False, is_mosaic=False) == "lowN"


def test_zygosity_mosaic_overrides_everything():
    tup = ("WT", "+1", "WT", "WT")
    assert _zygosity_call([tup, tup], is_homozygous=True, is_mosaic=True) == "mosaic"


# ---- _slot_description / render_haplotype_description -------------------

def test_slot_description_renders_insertion_with_seq():
    d = SlotDetail(label="+1", size=1, ins_seq="T", ref_pos=1639)
    assert _slot_description(d, "gRNA2", "ACGT") == "gRNA2:+T"


def test_slot_description_renders_deletion_from_ref():
    d = SlotDetail(label="-4", size=4, ins_seq="", ref_pos=2)
    # ref[2:6] = 'CDEF'
    assert _slot_description(d, "gRNA1", "ABCDEFGH") == "gRNA1:-CDEF"


def test_slot_description_returns_empty_for_wt_and_sv_and_nc():
    assert _slot_description(SlotDetail("WT"), "gRNA1", "ACGT") == ""
    assert _slot_description(SlotDetail("NC"), "gRNA1", "ACGT") == ""
    assert _slot_description(SlotDetail("SV"), "gRNA1", "ACGT") == ""


def _g(name: str, cut: int) -> Guide:
    return Guide(name=name, sequence="A" * 20, cut_pos=cut, strand="+")


def test_render_haplotype_collapses_sv_span():
    guides = [_g("gRNA1", 100), _g("gRNA2", 200), _g("gRNA3", 300), _g("gRNA4", 400)]
    sv_results = {(0, 3): {"median_size": 5411,
                           "supporting_read_names": set()}}
    tup = ("SV", "SV", "SV", "SV")
    details = [SlotDetail("SV")] * 4
    out = render_haplotype_description(tup, details, guides, "X" * 500, sv_results)
    assert out == "SV:-5411bp"


def test_render_haplotype_renders_per_guide_edits():
    guides = [_g("gRNA1", 100), _g("gRNA2", 200)]
    sv_results = {}
    tup = ("WT", "+1")
    details = [SlotDetail("WT"), SlotDetail("+1", 1, "T", 200)]
    out = render_haplotype_description(tup, details, guides, "Z" * 500, sv_results)
    assert out == "gRNA2:+T"


def test_render_haplotype_handles_pure_wt():
    guides = [_g("gRNA1", 100), _g("gRNA2", 200)]
    tup = ("WT", "WT")
    details = [SlotDetail("WT"), SlotDetail("WT")]
    out = render_haplotype_description(tup, details, guides, "Z" * 500, {})
    assert out == "WT"


# ---- per_guide_cells ----------------------------------------------------

def test_per_guide_cells_collapses_homozygous_to_single_label():
    guides = [_g("gRNA1", 100), _g("gRNA2", 200)]
    plus_t = SlotDetail("+1", 1, "T", 200)
    h1 = Haplotype(
        rank=1, tuple_=("WT", "+1"),
        details=[SlotDetail("WT"), plus_t],
        support_reads=80, support_frac=0.5,
        description="gRNA2:+T",
    )
    h2 = Haplotype(
        rank=2, tuple_=("WT", "+1"),
        details=[SlotDetail("WT"), plus_t],
        support_reads=80, support_frac=0.5,
        description="gRNA2:+T",
    )
    result = PhasingResult(
        sample="hom", n_phased=160, n_partial=0,
        haplotypes=[h1, h2], zygosity_call="homozygous",
        is_mosaic=False, notes=[], all_tuples=[],
    )
    cells = per_guide_cells(result, guides, "X" * 500)
    # both halves match → collapse; non-matching halves stay split.
    assert cells == ["WT", "+T"]


def test_per_guide_cells_renders_a1_slash_a2():
    guides = [_g("gRNA1", 100), _g("gRNA2", 200), _g("gRNA3", 300), _g("gRNA4", 400)]
    h1 = Haplotype(
        rank=1, tuple_=("SV", "SV", "SV", "SV"),
        details=[SlotDetail("SV")] * 4,
        support_reads=411, support_frac=0.7,
        description="SV:-5411bp",
    )
    h2 = Haplotype(
        rank=2, tuple_=("WT", "+1", "WT", "WT"),
        details=[SlotDetail("WT"), SlotDetail("+1", 1, "T", 200),
                 SlotDetail("WT"), SlotDetail("WT")],
        support_reads=174, support_frac=0.3,
        description="gRNA2:+T",
    )
    result = PhasingResult(
        sample="crd2", n_phased=585, n_partial=0,
        haplotypes=[h1, h2], zygosity_call="biallelic",
        is_mosaic=False, notes=[], all_tuples=[],
    )
    cells = per_guide_cells(result, guides, "X" * 500)
    assert cells == ["SV / WT", "SV / +T", "SV / WT", "SV / WT"]


# ---- detect_phasing_mismatch -------------------------------------------

def _make_sv_phasing_result(n_guides: int = 4) -> PhasingResult:
    """Build a phasing result with a single dominant SV/SV/.../SV haplotype.
    Mirrors the WV301-1-style case where partial-coverage stranded the
    second chromosome's reads and the phased view collapsed to one SV
    haplotype."""
    tup = tuple(["SV"] * n_guides)
    h = Haplotype(
        rank=1, tuple_=tup, details=[SlotDetail("SV")] * n_guides,
        support_reads=150, support_frac=1.0, description="SV:-5412bp",
    )
    return PhasingResult(
        sample="WV301-1", n_phased=150, n_partial=2800,
        haplotypes=[h], zygosity_call="homozygous", is_mosaic=False,
        notes=["partial-coverage:2800"], all_tuples=[],
    )


def test_phasing_mismatch_flags_per_guide_edit_missing_from_phased():
    # Real-world WV301-1 case: per-guide caller flagged 'ambiguous'
    # because of competing +T and SV-deletion at gRNA2; phased view only
    # sees SV. Mismatch detection should fire on the +T.
    guides = [_g("gRNA1", 100), _g("gRNA2", 200), _g("gRNA3", 300), _g("gRNA4", 400)]
    result = _make_sv_phasing_result(n_guides=4)
    plant_calls = [
        AlleleCall("WT",     500, 2.0,  1.0, 95.0, "high"),
        AlleleCall("+T/WT",  424, 41.0, 5.0, 9.0,  "ambiguous"),
        AlleleCall("WT",     500, 3.0,  1.0, 90.0, "high"),
        AlleleCall("WT",     500, 1.0,  0.0, 95.0, "high"),
    ]
    notes = detect_phasing_mismatch(plant_calls, result, guides)
    assert notes == ["phasing-mismatch:gRNA2:+T@41%"]


def test_phasing_mismatch_silent_when_phased_already_has_the_edit():
    guides = [_g("gRNA1", 100), _g("gRNA2", 200)]
    h = Haplotype(
        rank=1, tuple_=("WT", "+1"), details=[SlotDetail("WT"), SlotDetail("+1", 1, "T", 200)],
        support_reads=200, support_frac=1.0, description="gRNA2:+T",
    )
    result = PhasingResult(
        sample="x", n_phased=200, n_partial=0,
        haplotypes=[h], zygosity_call="homozygous", is_mosaic=False,
        notes=[], all_tuples=[],
    )
    plant_calls = [
        AlleleCall("WT",   200, 2.0,  1.0, 95.0, "high"),
        AlleleCall("+T",   200, 95.0, 1.0, 4.0,  "high"),  # already in phased
    ]
    assert detect_phasing_mismatch(plant_calls, result, guides) == []


def test_phasing_mismatch_skips_low_fraction_signals():
    guides = [_g("gRNA1", 100), _g("gRNA2", 200)]
    result = _make_sv_phasing_result(n_guides=2)
    plant_calls = [
        AlleleCall("WT",    500, 2.0,  1.0, 95.0, "high"),
        AlleleCall("+T/WT", 500, 8.0,  2.0, 88.0, "low"),  # below 15% threshold
    ]
    assert detect_phasing_mismatch(plant_calls, result, guides) == []


def test_phasing_mismatch_skips_lowN_per_guide():
    guides = [_g("gRNA1", 100), _g("gRNA2", 200)]
    result = _make_sv_phasing_result(n_guides=2)
    plant_calls = [
        AlleleCall("WT",   500, 1.0,  0.0, 98.0, "high"),
        AlleleCall("lowN", 5,   0.0,  0.0, 0.0,  "lowN"),
    ]
    assert detect_phasing_mismatch(plant_calls, result, guides) == []


def test_phasing_mismatch_flags_unnamed_edit_when_call_is_WT_but_top_frac_high():
    # WV301-3-style case: per-guide call is 'WT' (signal didn't reach
    # the caller's own threshold for naming an allele), but top_frac is
    # 17% — above the cross-check threshold. Phased view is SV-only at
    # this slot, so we should still flag it (with a generic 'edit' note).
    guides = [_g("gRNA1", 100), _g("gRNA2", 200)]
    result = _make_sv_phasing_result(n_guides=2)
    plant_calls = [
        AlleleCall("WT", 200, 5.0,  3.0,  90.0, "high"),
        AlleleCall("WT", 190, 17.0, 10.0, 2.6,  "ambiguous"),
    ]
    notes = detect_phasing_mismatch(plant_calls, result, guides)
    assert notes == ["phasing-mismatch:gRNA2:edit@17%"]


def test_phasing_mismatch_silent_when_phased_already_has_a_real_edit_at_slot():
    # If the phased view already has a non-WT non-SV slot at this guide,
    # an unnamed per-guide signal is presumed to be that same edit and
    # we don't re-flag it.
    guides = [_g("gRNA1", 100), _g("gRNA2", 200)]
    h = Haplotype(
        rank=1, tuple_=("WT", "+1"),
        details=[SlotDetail("WT"), SlotDetail("+1", 1, "T", 200)],
        support_reads=200, support_frac=1.0, description="gRNA2:+T",
    )
    result = PhasingResult(
        sample="x", n_phased=200, n_partial=0,
        haplotypes=[h], zygosity_call="homozygous", is_mosaic=False,
        notes=[], all_tuples=[],
    )
    plant_calls = [
        AlleleCall("WT", 200, 2.0,  1.0,  95.0, "high"),
        AlleleCall("WT", 200, 17.0, 5.0,  70.0, "ambiguous"),
    ]
    assert detect_phasing_mismatch(plant_calls, result, guides) == []


def test_phasing_mismatch_returns_empty_when_no_haplotypes():
    guides = [_g("gRNA1", 100)]
    result = PhasingResult(
        sample="x", n_phased=0, n_partial=0,
        haplotypes=[], zygosity_call="lowN", is_mosaic=False,
        notes=[], all_tuples=[],
    )
    plant_calls = [AlleleCall("+T", 200, 90.0, 0.0, 5.0, "high")]
    assert detect_phasing_mismatch(plant_calls, result, guides) == []
