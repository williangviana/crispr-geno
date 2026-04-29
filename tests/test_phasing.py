"""Unit tests for the pure-Python phasing logic.

The BAM-traversal part of phasing is exercised end-to-end against the
crd2 sample in the validation pass; these tests cover the pieces that
can be reasoned about without real read data.
"""
from __future__ import annotations

import collections

from crispr_geno.analysis import Guide
from crispr_geno.phasing import (
    SlotDetail,
    _slot_description,
    _zygosity_call,
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
