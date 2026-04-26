"""Write color-coded XLSX versions of the genotype and sequence tables."""
from __future__ import annotations

from openpyxl import Workbook
from openpyxl.styles import Alignment, Border, Font, PatternFill, Side


_HF = Font(name="Arial", size=10, bold=True, color="FFFFFF")
_HFILL = PatternFill("solid", start_color="305496")
_CEN = Alignment(horizontal="center", vertical="center")
_LFT = Alignment(horizontal="left", vertical="center")
_THIN = Side(border_style="thin", color="000000")
_BD = Border(left=_THIN, right=_THIN, top=_THIN, bottom=_THIN)
_AMB = Side(border_style="medium", color="CC0000")
_BD_AMBIG = Border(left=_AMB, right=_AMB, top=_AMB, bottom=_AMB)
_ARIAL = Font(name="Arial", size=10)
_MONO = Font(name="Arial", size=10)

_GENO_FILLS = {
    "WT":            PatternFill("solid", start_color="E7E6E6"),
    "T":             PatternFill("solid", start_color="C6E0B4"),
    "T/+":           PatternFill("solid", start_color="E2EFDA"),
    "del":           PatternFill("solid", start_color="F8CBAD"),
    "del/+":         PatternFill("solid", start_color="FCE4D6"),
    "T+del":         PatternFill("solid", start_color="FFD966"),
    "T/+ & del/+":   PatternFill("solid", start_color="FFE699"),
    "lowN":          PatternFill("solid", start_color="D9D9D9"),
}


def _style_header(ws) -> None:
    for c in ws[1]:
        c.font = _HF; c.fill = _HFILL; c.alignment = _CEN; c.border = _BD


def genotype_xlsx(header: list[str], rows: list[list[str]], xlsx: str) -> None:
    wb = Workbook(); ws = wb.active; ws.title = "Genotypes"
    ws.append(header)
    for row in rows:
        ws.append(row)
    _style_header(ws)
    for row in ws.iter_rows(min_row=2):
        for i, c in enumerate(row):
            c.font = _ARIAL; c.alignment = _CEN; c.border = _BD
            if i > 0 and c.value in _GENO_FILLS:
                c.fill = _GENO_FILLS[c.value]
    ws.column_dimensions["A"].width = 16
    for col_idx in range(2, ws.max_column + 1):
        ws.column_dimensions[ws.cell(1, col_idx).column_letter].width = 18
    ws.freeze_panes = "A2"
    wb.save(xlsx)


def _pick_seq_fill(value: str) -> PatternFill | None:
    if value == "WT":
        return PatternFill("solid", start_color="E7E6E6")
    if value == "lowN":
        return PatternFill("solid", start_color="F4CCCC")
    # Plant-level status values
    if value == "stable":
        return PatternFill("solid", start_color="D9E1F2")  # pale blue
    if value == "Cas9+":
        return PatternFill("solid", start_color="F4CCCC")  # light red — needs attention
    # Structural variant (large deletion between cuts)
    if value.startswith("SV:") and value.endswith("/WT"):
        return PatternFill("solid", start_color="D5A6BD")  # light purple/pink (het)
    if value.startswith("SV:"):
        return PatternFill("solid", start_color="B084CB")  # purple (homo)
    # het vs WT: /WT suffix
    if value.endswith("/WT"):
        if value.startswith("+"):
            return PatternFill("solid", start_color="E2EFDA")  # light green
        if value.startswith("-"):
            return PatternFill("solid", start_color="FCE4D6")  # light orange
    # biallelic: two different mutant alleles separated by '/' (no WT suffix)
    if "/" in value and not value.endswith("/WT"):
        return PatternFill("solid", start_color="FFE699")  # yellow
    # homozygous
    if value.startswith("+"):
        return PatternFill("solid", start_color="C6E0B4")  # green
    if value.startswith("-"):
        return PatternFill("solid", start_color="F8CBAD")  # orange
    return PatternFill("solid", start_color="FFF2CC")


def sequence_xlsx(
    header: list[str],
    rows: list[list[str]],
    xlsx: str,
    conf_rows: list[list[str]] | None = None,
) -> None:
    """Write color-coded sequence table with an inline legend.

    If conf_rows is provided (same shape as rows without the sample-id column),
    cells flagged 'ambiguous' get a red border.
    """
    wb = Workbook(); ws = wb.active; ws.title = "Sequences"
    ws.append(header)
    for row in rows:
        ws.append(row)
    _style_header(ws)

    for ridx, row in enumerate(ws.iter_rows(min_row=2), start=0):
        row[0].font = _ARIAL; row[0].alignment = _CEN; row[0].border = _BD
        conf = conf_rows[ridx] if conf_rows and ridx < len(conf_rows) else None
        for ci, c in enumerate(row[1:]):
            c.font = _MONO; c.alignment = _LFT
            c.border = _BD_AMBIG if (conf and ci < len(conf) and conf[ci] == "ambiguous") else _BD
            fill = _pick_seq_fill(c.value or "")
            if fill is not None:
                c.fill = fill

    ws.column_dimensions["A"].width = 16
    for col_idx in range(2, ws.max_column + 1):
        ws.column_dimensions[ws.cell(1, col_idx).column_letter].width = 36
    ws.freeze_panes = "A2"

    # --- Legend: inline, same sheet, one column of buffer to the right of the table ---
    last_col = ws.max_column          # last column of the data table
    buffer_col = last_col + 1         # empty column separator
    color_col = last_col + 2          # color swatch column
    label_col = last_col + 3          # allele / notation label
    meaning_col = last_col + 4        # explanation

    # Legend header (row 1, aligned with the data table's header row)
    for col, title in ((color_col, "Color"),
                       (label_col, "Allele / Notation"),
                       (meaning_col, "Meaning")):
        cell = ws.cell(1, col, value=title)
        cell.font = _HF; cell.fill = _HFILL; cell.alignment = _CEN; cell.border = _BD

    # Generic notation entries (explanations of what the symbols mean)
    generic_entries = [
        ("stable",   "Status: all sites resolve cleanly — Cas9 likely inactive or absent",
            PatternFill("solid", start_color="D9E1F2")),
        ("Cas9+",    "Status: at least one site shows mosaic / residual editing — Cas9 still active",
            PatternFill("solid", start_color="F4CCCC")),
        ("WT",       "wildtype — no edit detected at this site",
            PatternFill("solid", start_color="E7E6E6")),
        ("lowN",     "fewer than 10 reads spanning the cut — cannot call",
            PatternFill("solid", start_color="F4CCCC")),
        ("+XYZ",     "homozygous insertion of XYZ at cut site",
            PatternFill("solid", start_color="C6E0B4")),
        ("+XYZ/WT",  "heterozygous: one chromosome has +XYZ, the other is wildtype",
            PatternFill("solid", start_color="E2EFDA")),
        ("-XYZ",     "homozygous deletion of XYZ from reference",
            PatternFill("solid", start_color="F8CBAD")),
        ("-XYZ/WT",  "heterozygous: one chromosome has -XYZ, the other is wildtype",
            PatternFill("solid", start_color="FCE4D6")),
        ("A/B",      "biallelic: one chromosome has allele A, the other has a different mutant allele B",
            PatternFill("solid", start_color="FFE699")),
        ("SV:-Nbp",  "homozygous large deletion between two cut sites (N bp removed)",
            PatternFill("solid", start_color="B084CB")),
        ("SV:-Nbp/WT", "heterozygous large deletion between two cut sites",
            PatternFill("solid", start_color="D5A6BD")),
        ("ambiguous", "competing allele present at ≥15% — check editing_rates.csv (shown as a red border on the cell)",
            None),
    ]
    row_idx = 2
    for label, meaning, fill in generic_entries:
        if fill is not None:
            ws.cell(row_idx, color_col).fill = fill
            ws.cell(row_idx, color_col).border = _BD
        else:
            # "red border" entry: show the border itself as the swatch
            ws.cell(row_idx, color_col).border = _BD_AMBIG
        lbl = ws.cell(row_idx, label_col, value=label)
        lbl.font = _ARIAL; lbl.alignment = _LFT
        mng = ws.cell(row_idx, meaning_col, value=meaning)
        mng.font = _ARIAL; mng.alignment = _LFT
        for col in (label_col, meaning_col):
            ws.cell(row_idx, col).border = _BD
        row_idx += 1

    # Sizing for the legend columns
    ws.column_dimensions[ws.cell(1, buffer_col).column_letter].width = 3
    ws.column_dimensions[ws.cell(1, color_col).column_letter].width = 8
    ws.column_dimensions[ws.cell(1, label_col).column_letter].width = 22
    ws.column_dimensions[ws.cell(1, meaning_col).column_letter].width = 75

    wb.save(xlsx)
