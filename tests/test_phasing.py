"""Unit tests for the pure-Python phasing logic.

The BAM-traversal part of phasing is exercised end-to-end against the
crd2 sample in the validation pass; these tests cover the pieces that
can be reasoned about without real read data.
"""
from __future__ import annotations

import collections

from crispr_geno.analysis import (
    AlleleCall,
    Guide,
    _left_normalize_deletion,
    _left_normalize_insertion,
    detect_per_guide_mosaic,
    per_guide_noise_floors,
)
from crispr_geno.phasing import (
    SlotDetail,
    _slot_description,
    _zygosity_call,
    detect_phasing_mismatch,
    merge_noise_tuples,
    per_guide_cells,
    per_guide_imputation,
    PhasingResult,
    Haplotype,
    pick_haplotypes,
    render_haplotype_description,
    residual_non_wt_phased,
)


# ---- indel left-normalization -------------------------------------------

def test_left_normalize_deletion_in_homopolymer():
    # Reference: "AAATTTTCGT". Deleting any single T at positions 3, 4,
    # 5, or 6 produces the same biological outcome — should all
    # normalize to position 3 (the leftmost T).
    ref = "AAATTTTCGT"
    for pos in (3, 4, 5, 6):
        assert _left_normalize_deletion(pos, 1, ref) == 3


def test_left_normalize_deletion_unique_position_unchanged():
    # Non-repetitive context: nothing to shift.
    ref = "AAACGTACGT"
    assert _left_normalize_deletion(4, 2, ref) == 4   # del 'GT'


def test_left_normalize_insertion_in_homopolymer():
    # Reference: "AAATTTTCGT". Inserting +T at positions 3-7 inside the
    # TTTT run all normalize to position 3 with seq "T".
    ref = "AAATTTTCGT"
    for pos in (3, 4, 5, 6, 7):
        norm_pos, norm_seq = _left_normalize_insertion(pos, "T", ref)
        assert norm_pos == 3
        assert norm_seq == "T"


def test_left_normalize_insertion_rotates_multibase_seq():
    # Inserting "TC" at the end of a TTTT homopolymer — rotates as it
    # shifts left. After full rotation the result is unchanged at the
    # start of TTTT, since inserting "TC" before the run vs. after is
    # only equivalent for cyclic-rotation-friendly seqs.
    ref = "AAATTTTCGT"
    norm_pos, norm_seq = _left_normalize_insertion(7, "TC", ref)
    # 'TC' inserted at pos 7: ref[6]='T'==seq[-1]='C'? No → no shift.
    assert (norm_pos, norm_seq) == (7, "TC")


def test_left_normalize_insertion_no_seq_returns_unchanged():
    # Empty inserted sequence → no rotation possible.
    ref = "AAATTTTCGT"
    assert _left_normalize_insertion(5, "", ref) == (5, "")


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


def test_pick_haplotypes_flags_mosaic_when_wt_dominant_but_non_wt_fragmented():
    # WV151-2-4-7 pattern: WT at ~52%, plus 6 distinct edited haplotypes
    # each at 5-12%. No single edited tuple clears het threshold, but the
    # non-WT mass is ~48% spread across many lesions — mosaic on the
    # second chromosome, not a clean het.
    merged = collections.Counter({
        ("WT", "WT", "WT", "WT"): 339,
        ("-1", "WT", "-32", "-12"): 81,
        ("-1", "WT", "-16", "-12"): 54,
        ("-6", "WT", "-10", "-15"): 52,
        ("-7", "WT", "WT", "-10"): 47,
        ("-1", "WT", "-26", "-49"): 45,
        ("-45", "WT", "WT", "WT"): 38,
    })
    _, mosaic = pick_haplotypes(merged, het_threshold_frac=0.20,
                                mosaic_threshold_frac=0.10)
    assert mosaic is True


def test_pick_haplotypes_clean_het_not_flagged_as_mosaic():
    # Clean het: WT and one edited haplotype both around 50%. The
    # mosaic-on-one-chromosome check must not fire here.
    merged = collections.Counter({
        ("WT", "WT"): 50,
        ("+1", "WT"): 50,
    })
    _, mosaic = pick_haplotypes(merged, het_threshold_frac=0.20,
                                mosaic_threshold_frac=0.10)
    assert mosaic is False


def test_pick_haplotypes_het_with_minor_noise_not_flagged_as_mosaic():
    # WT + dominant edit + a couple of minor edited tuples below het
    # threshold. The dominant edit clears het threshold, so this is a
    # clean het — not mosaic.
    merged = collections.Counter({
        ("WT", "WT"): 50,
        ("+1", "WT"): 40,
        ("-1", "WT"): 5,
        ("-2", "WT"): 5,
    })
    _, mosaic = pick_haplotypes(merged, het_threshold_frac=0.20,
                                mosaic_threshold_frac=0.10)
    assert mosaic is False


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


# ---- detect_per_guide_mosaic --------------------------------------------

def _ac(top_frac: float, wt_frac: float, conf: str = "med") -> AlleleCall:
    return AlleleCall(call="x", spanning=1000, top_frac=top_frac,
                      second_frac=0.0, wt_frac=wt_frac, confidence=conf)


def test_per_guide_mosaic_fires_on_two_fragmented_guides():
    # WV151-2-4-7 pattern: gRNA3 and gRNA4 both have lots of editing
    # but no dominant edited allele.
    plant_calls = [
        _ac(top_frac=30.35, wt_frac=42.22),  # gRNA1 — has dominant edit
        _ac(top_frac=0.34, wt_frac=97.90),   # gRNA2 — clean WT
        _ac(top_frac=7.87, wt_frac=67.22),   # gRNA3 — fragmented
        _ac(top_frac=16.61, wt_frac=54.54),  # gRNA4 — fragmented
    ]
    assert detect_per_guide_mosaic(plant_calls, het_threshold=20.0) is True


def test_per_guide_mosaic_does_not_fire_on_single_fragmented_guide():
    # WV161-2-6-8 pattern: only gRNA4 is fragmented in the per-guide
    # view; gRNA1 and gRNA3 both have dominant edited alleles.
    plant_calls = [
        _ac(top_frac=32.50, wt_frac=51.26),  # gRNA1 — dominant edit
        _ac(top_frac=0.18, wt_frac=98.66),   # gRNA2 — clean WT
        _ac(top_frac=35.84, wt_frac=51.73),  # gRNA3 — dominant edit
        _ac(top_frac=12.37, wt_frac=66.10),  # gRNA4 — fragmented
    ]
    assert detect_per_guide_mosaic(plant_calls, het_threshold=20.0) is False


def test_per_guide_mosaic_ignores_lowN_guides():
    # A lowN guide with weird-looking numbers must not contribute to
    # the fragmented-guide count.
    plant_calls = [
        _ac(top_frac=5.0, wt_frac=50.0, conf="lowN"),
        _ac(top_frac=5.0, wt_frac=50.0, conf="lowN"),
    ]
    assert detect_per_guide_mosaic(plant_calls, het_threshold=20.0) is False


def test_per_guide_mosaic_clean_het_at_every_guide_not_mosaic():
    # All guides have a clean ~50/50 het signal — dominant edited
    # allele clears the threshold everywhere. Not mosaic.
    plant_calls = [
        _ac(top_frac=48.0, wt_frac=50.0),
        _ac(top_frac=48.0, wt_frac=50.0),
        _ac(top_frac=48.0, wt_frac=50.0),
    ]
    assert detect_per_guide_mosaic(plant_calls, het_threshold=20.0) is False


# ---- residual_non_wt_phased ---------------------------------------------

def _phasing(
    n_phased: int,
    haplotypes: list[tuple],
    all_tuples: list[tuple],
) -> PhasingResult:
    haps = [Haplotype(rank=i + 1, tuple_=t, details=[],
                      support_reads=0, support_frac=0.0, description="")
            for i, t in enumerate(haplotypes)]
    return PhasingResult(
        sample="x", n_phased=n_phased, n_partial=0, haplotypes=haps,
        zygosity_call="x", is_mosaic=False, notes=[], all_tuples=all_tuples,
    )


def test_residual_non_wt_zero_for_clean_het():
    # WV155-2-9 pattern: WT chromosome + edited chromosome, both called.
    # No residual signal — Cas9 fired once and is done.
    p = _phasing(
        n_phased=706,
        haplotypes=[("WT", "WT", "WT", "WT"), ("+1", "WT", "WT", "WT")],
        all_tuples=[(("WT", "WT", "WT", "WT"), 442),
                    (("+1", "WT", "WT", "WT"), 264)],
    )
    frac, reads = residual_non_wt_phased(p)
    assert frac == 0.0
    assert reads == 0


def test_residual_non_wt_flags_small_subpopulation_in_wt_dominant():
    # WV161-2-6-2 pattern: WT alone is the called chromosome, the
    # 7% edited tuple is residual (Cas9 still active in some lineages).
    # 78 reads of residual at n_phased=1115 — clears both 5% and 10-read floors.
    p = _phasing(
        n_phased=1115,
        haplotypes=[("WT", "WT", "WT", "WT")],
        all_tuples=[(("WT", "WT", "WT", "WT"), 1037),
                    (("-1", "WT", "-24", "WT"), 78)],
    )
    frac, reads = residual_non_wt_phased(p)
    assert abs(frac - 78 / 1115 * 100.0) < 1e-9
    assert reads == 78


def test_residual_non_wt_low_n_phased_returns_few_reads():
    # 24-297-8 pattern: only 28 phased reads (95%+ went to partial-coverage),
    # 3 of them carry an extra noise indel at gRNA4. Frac is 10.7% but only
    # 3 absolute reads — caller (cli) requires both >=5% AND >=10 reads, so
    # this should NOT trip Cas9+ at the call site.
    p = _phasing(
        n_phased=28,
        haplotypes=[("-25", "WT", "WT", "WT")],
        all_tuples=[(("-25", "WT", "WT", "WT"), 25),
                    (("-25", "WT", "WT", "-2"), 3)],
    )
    frac, reads = residual_non_wt_phased(p)
    assert abs(frac - 3 / 28 * 100.0) < 1e-9
    assert reads == 3


def test_residual_non_wt_zero_when_no_phased_reads():
    p = _phasing(n_phased=0, haplotypes=[], all_tuples=[])
    frac, reads = residual_non_wt_phased(p)
    assert frac == 0.0
    assert reads == 0


def test_residual_non_wt_ignores_called_edited_tuples():
    # A biallelic call where a minor tuple introduces an edit at a slot
    # the called chromosomes don't touch (both calls are WT at slot 1).
    # The minor tuple is a real new edit, not a noise variant.
    p = _phasing(
        n_phased=1000,
        haplotypes=[("+1", "WT"), ("-1", "WT")],
        all_tuples=[(("+1", "WT"), 480),
                    (("-1", "WT"), 470),
                    (("+1", "+1"), 50)],   # +1 at slot 1 is new
    )
    frac, reads = residual_non_wt_phased(p)
    assert abs(frac - 5.0) < 1e-9
    assert reads == 50


def test_residual_non_wt_ignores_basecaller_variants_of_called_edit():
    # crd4-plant4 pattern: called chromosome is WT/-2/WT/+1 (homozygous).
    # The two minor tuples drop or swap the +1 at gRNA4 — both are
    # nanopore homopolymer noise around the called +T, not new edits.
    # Neither residual introduces a non-WT slot where the call is WT.
    p = _phasing(
        n_phased=101,
        haplotypes=[("WT", "-2", "WT", "+1")],
        all_tuples=[(("WT", "-2", "WT", "+1"), 91),
                    (("WT", "-2", "WT", "WT"), 5),
                    (("WT", "-2", "WT", "+2"), 5)],
    )
    frac, reads = residual_non_wt_phased(p)
    assert frac == 0.0
    assert reads == 0


def test_residual_non_wt_ignores_dropped_edit_from_called():
    # crd4-plant3 pattern: called WT/-2/WT/+1, residual WT/-2/WT/WT.
    # Residual is the called chromosome with the +1 at gRNA4 dropped —
    # a basecaller miscall at the homopolymer, not a real edit event.
    p = _phasing(
        n_phased=137,
        haplotypes=[("WT", "-2", "WT", "+1")],
        all_tuples=[(("WT", "-2", "WT", "+1"), 126),
                    (("WT", "-2", "WT", "WT"), 11)],
    )
    frac, reads = residual_non_wt_phased(p)
    assert frac == 0.0
    assert reads == 0


def test_residual_non_wt_noise_floor_suppresses_within_floor_signal():
    # Residual claims a "new edit" at slot 0 (call is WT there) — but
    # this sample's per-guide signal at slot 0 is 4%, below the batch
    # noise floor of 8%. Discard as batch noise.
    p = _phasing(
        n_phased=1000,
        haplotypes=[("WT", "WT")],
        all_tuples=[(("WT", "WT"), 920),
                    (("-1", "WT"), 80)],
    )
    frac, reads = residual_non_wt_phased(
        p,
        slot_signals=[4.0, 0.0],          # 4% at gRNA1 — below floor
        slot_noise_floors=[8.0, 0.0],     # batch noise at gRNA1 is 8%
    )
    assert frac == 0.0
    assert reads == 0


def test_residual_non_wt_noise_floor_passes_above_floor_signal():
    # Same shape as previous, but per-guide signal (15%) clears the
    # batch noise floor (8%) — the residual is real Cas9 activity.
    p = _phasing(
        n_phased=1000,
        haplotypes=[("WT", "WT")],
        all_tuples=[(("WT", "WT"), 920),
                    (("-1", "WT"), 80)],
    )
    frac, reads = residual_non_wt_phased(
        p,
        slot_signals=[15.0, 0.0],
        slot_noise_floors=[8.0, 0.0],
    )
    assert reads == 80
    assert abs(frac - 8.0) < 1e-9


def test_residual_non_wt_noise_floor_optional():
    # When floors aren't provided, behavior is unchanged from the
    # no-floor version — the residual at slot 0 still counts.
    p = _phasing(
        n_phased=1000,
        haplotypes=[("WT", "WT")],
        all_tuples=[(("WT", "WT"), 920),
                    (("-1", "WT"), 80)],
    )
    frac, reads = residual_non_wt_phased(p)
    assert reads == 80


# ---- per_guide_noise_floors ---------------------------------------------

def test_noise_floors_returns_zero_when_too_few_clean_wt_samples():
    # Only 2 high-confidence WT samples at gRNA1 — below the
    # min_clean_samples=3 default → no floor learned.
    plant_calls_per_sample = [
        [AlleleCall("WT", 500, 1.0, 0.5, 99.0, "high")],
        [AlleleCall("WT", 500, 2.0, 1.0, 98.0, "high")],
        [AlleleCall("+T", 500, 95.0, 1.0, 4.0, "high")],
    ]
    floors = per_guide_noise_floors(plant_calls_per_sample, n_guides=1)
    assert floors == [0.0]


def test_noise_floors_learns_from_clean_wt_samples():
    # Five high-confidence WT samples at gRNA1 with top_frac
    # 1, 2, 3, 4, 6 → 90th percentile ≈ 6 → floor = 6 * 1.5 = 9.0.
    plant_calls_per_sample = [
        [AlleleCall("WT", 500, top, 0.5, 99.0 - top, "high")]
        for top in (1.0, 2.0, 3.0, 4.0, 6.0)
    ]
    floors = per_guide_noise_floors(plant_calls_per_sample, n_guides=1)
    assert abs(floors[0] - 6.0 * 1.5) < 1e-9


def test_noise_floors_ignores_non_wt_calls_and_low_confidence():
    # Six samples at gRNA1: only the 3 high-confidence WT samples
    # contribute. The +T calls and the med/low WT calls are excluded.
    plant_calls_per_sample = [
        [AlleleCall("WT", 500, 1.0, 0.0, 99.0, "high")],
        [AlleleCall("WT", 500, 2.0, 0.0, 98.0, "high")],
        [AlleleCall("WT", 500, 3.0, 0.0, 97.0, "high")],
        [AlleleCall("+T", 500, 90.0, 0.0, 9.0, "high")],
        [AlleleCall("WT", 30, 5.0, 0.0, 95.0, "low")],
        [AlleleCall("WT", 500, 8.0, 0.0, 92.0, "med")],
    ]
    floors = per_guide_noise_floors(plant_calls_per_sample, n_guides=1)
    # 90th percentile of [1, 2, 3] is 3 → floor = 3 * 1.5 = 4.5.
    assert abs(floors[0] - 3.0 * 1.5) < 1e-9


def test_noise_floors_independent_per_guide():
    # gRNA4 has noisy WT (homopolymer site); gRNA1 is quiet.
    plant_calls_per_sample = [
        [AlleleCall("WT", 500, 1.0, 0.0, 99.0, "high"),
         AlleleCall("WT", 500, 5.0, 0.0, 95.0, "high")],
        [AlleleCall("WT", 500, 2.0, 0.0, 98.0, "high"),
         AlleleCall("WT", 500, 6.0, 0.0, 94.0, "high")],
        [AlleleCall("WT", 500, 3.0, 0.0, 97.0, "high"),
         AlleleCall("WT", 500, 7.0, 0.0, 93.0, "high")],
    ]
    floors = per_guide_noise_floors(plant_calls_per_sample, n_guides=2)
    assert abs(floors[0] - 3.0 * 1.5) < 1e-9
    assert abs(floors[1] - 7.0 * 1.5) < 1e-9


# ---- per_guide_imputation -----------------------------------------------

def test_per_guide_imputation_wt_high_confidence():
    guides = [_g("gRNA1", 100)]
    plant_calls = [AlleleCall("WT", 500, 0.5, 0.0, 99.5, "high")]
    out = per_guide_imputation(plant_calls, guides, "ACGT" * 100)
    assert out[0] is not None
    assert out[0].label == "WT"


def test_per_guide_imputation_homozygous_insertion_high_confidence():
    guides = [_g("gRNA2", 200)]
    plant_calls = [AlleleCall("+T", 500, 96.0, 0.5, 1.0, "high")]
    out = per_guide_imputation(plant_calls, guides, "ACGT" * 100)
    assert out[0] is not None
    assert out[0].label == "+1"
    assert out[0].ins_seq == "T"


def test_per_guide_imputation_homozygous_deletion_high_confidence():
    guides = [_g("gRNA1", 100)]
    plant_calls = [AlleleCall("-CCC", 500, 95.0, 0.5, 2.0, "high")]
    out = per_guide_imputation(plant_calls, guides, "ACGT" * 100)
    assert out[0] is not None
    assert out[0].label == "-3"


def test_per_guide_imputation_het_call_returns_none():
    # Het calls are ambiguous — we can't impute one of two chromosomes.
    guides = [_g("gRNA1", 100)]
    plant_calls = [AlleleCall("+T/WT", 500, 36.0, 0.5, 60.0, "low")]
    out = per_guide_imputation(plant_calls, guides, "ACGT" * 100)
    assert out[0] is None


def test_per_guide_imputation_low_confidence_returns_none():
    # Even a single-allele call gets dropped if confidence isn't high.
    guides = [_g("gRNA1", 100)]
    plant_calls = [AlleleCall("WT", 500, 5.0, 0.5, 95.0, "med")]
    out = per_guide_imputation(plant_calls, guides, "ACGT" * 100)
    assert out[0] is None


def test_per_guide_imputation_lowN_returns_none():
    guides = [_g("gRNA1", 100)]
    plant_calls = [AlleleCall("lowN", 5, 0.0, 0.0, 0.0, "lowN")]
    out = per_guide_imputation(plant_calls, guides, "ACGT" * 100)
    assert out[0] is None


def test_per_guide_imputation_independent_per_guide():
    # Mixed: gRNA1 het (no impute), gRNA2 high-conf homozygous +T (impute).
    guides = [_g("gRNA1", 100), _g("gRNA2", 200)]
    plant_calls = [
        AlleleCall("+T/WT", 500, 36.0, 0.5, 60.0, "low"),
        AlleleCall("+T", 500, 96.0, 0.5, 1.0, "high"),
    ]
    out = per_guide_imputation(plant_calls, guides, "ACGT" * 100)
    assert out[0] is None
    assert out[1] is not None
    assert out[1].label == "+1"


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


def test_phasing_mismatch_silent_when_per_guide_large_deletion_is_the_SV():
    # ME034v-pWGV242-4-5-style: SV is a 102 bp deletion spanning gRNA1
    # and gRNA2. Phased view shows 'SV' at those slots; the per-guide
    # caller (which doesn't know about SVs) reports the same 102 bp
    # deletion as a long indel ('-AATATGCTGAT...AT'). Same lesion —
    # don't flag it as a mismatch.
    guides = [_g("gRNA1", 100), _g("gRNA2", 200), _g("gRNA3", 300)]
    h = Haplotype(
        rank=1, tuple_=("SV", "SV", "-1"),
        details=[SlotDetail("SV"), SlotDetail("SV"), SlotDetail("-1", 1, "", 300)],
        support_reads=49, support_frac=1.0, description="SV:-102bp",
    )
    result = PhasingResult(
        sample="x", n_phased=49, n_partial=329,
        haplotypes=[h], zygosity_call="homozygous", is_mosaic=False,
        notes=[], all_tuples=[],
    )
    long_del = "-" + "A" * 102  # 102 bp deletion = the SV size
    plant_calls = [
        AlleleCall(long_del, 99, 99.0, 1.0, 0.0, "med"),
        AlleleCall(long_del, 99, 99.0, 1.0, 0.0, "med"),
        AlleleCall("-G",     228, 96.9, 0.9, 0.9, "high"),
    ]
    assert detect_phasing_mismatch(plant_calls, result, guides) == []


def test_phasing_mismatch_still_flags_small_deletion_when_phased_has_SV():
    # Same SV setup, but the per-guide caller sees a small deletion
    # (10 bp, below the SV threshold) at a guide that the phased view
    # marks as SV — that's a different lesion from a real second
    # chromosome, still a mismatch.
    guides = [_g("gRNA1", 100), _g("gRNA2", 200)]
    h = Haplotype(
        rank=1, tuple_=("SV", "SV"),
        details=[SlotDetail("SV"), SlotDetail("SV")],
        support_reads=200, support_frac=1.0, description="SV:-5kbp",
    )
    result = PhasingResult(
        sample="x", n_phased=200, n_partial=0,
        haplotypes=[h], zygosity_call="homozygous", is_mosaic=False,
        notes=[], all_tuples=[],
    )
    plant_calls = [
        AlleleCall("WT",          200, 1.0,  0.5, 98.0, "high"),
        AlleleCall("-AAAAAAAAAA", 200, 25.0, 5.0, 70.0, "ambiguous"),  # 10 bp
    ]
    notes = detect_phasing_mismatch(plant_calls, result, guides)
    assert notes == ["phasing-mismatch:gRNA2:-AAAAAAAAAA@25%"]


def test_phasing_mismatch_returns_empty_when_no_haplotypes():
    guides = [_g("gRNA1", 100)]
    result = PhasingResult(
        sample="x", n_phased=0, n_partial=0,
        haplotypes=[], zygosity_call="lowN", is_mosaic=False,
        notes=[], all_tuples=[],
    )
    plant_calls = [AlleleCall("+T", 200, 90.0, 0.0, 5.0, "high")]
    assert detect_phasing_mismatch(plant_calls, result, guides) == []
