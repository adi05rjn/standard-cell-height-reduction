# Prints the measured spacings of every built model next to the value drawn in the papers.
# Plain Python, no SEMulator3D needed:  python summarize_all.py
import json, os

here = os.path.dirname(os.path.abspath(__file__))
out = os.path.join(here, "out")

def load(name):
    with open(os.path.join(out, name)) as f:
        return json.load(f)

def runs(prof, names):
    """Merge neighbouring segments whose material is in names -> list of (start, end)."""
    r = []
    for n, a, b in prof:
        if n in names:
            if r and abs(r[-1][1] - a) < 1e-6:
                r[-1] = (r[-1][0], b)
            else:
                r.append((a, b))
    return r

def gaps(r):
    return [round(r[i + 1][0] - r[i][1], 2) for i in range(len(r) - 1)]

def seg_runs(spans):
    r = []
    for a, b in spans:
        if r and abs(r[-1][1] - a) < 1e-6:
            r[-1] = (r[-1][0], b)
        else:
            r.append((a, b))
    return r

rows = []
def row(model, what, paper, measured):
    rows.append((model, what, paper, measured))

METAL = ("W", "Ti")
for t, tag, pn_tip, tap in ((6.5, "T6p5", "31.5 (41.5 - 2x5)", "16.5 (43/2 - 5)"),
                            (5.5, "T5p5", "10.5", "16.5"),
                            (4.5, "T4p5", "30", "17.25")):
    j = load("stage2_%s.json" % tag)
    c0 = gaps(seg_runs(j["contact_y_spans"]))
    c1 = gaps(seg_runs(j["contact_y_spans_col1"]))
    row("Stage 2 %sT" % t, "contact gaps along column 1 (no power tap): N-P, border, ...", pn_tip + " / 33 or 34.5", c1)
    row("Stage 2 %sT" % t, "contact gaps along column 0 (middle row taps power)", "border with tap: " + tap, c0)
    row("Stage 2 %sT" % t, "epi gaps along column 0", "-", [g for g in j["epi_gaps"]])

j = load("stage3_T4p5.json")
p = j["profiles"]
lk0 = [round(b - a, 2) for a, b in runs(p["y_col0_wing_level_row0"], ("SiO2_SiCOH",))]
lk1 = [round(b - a, 2) for a, b in runs(p["y_col1_wing_level_row0"], ("SiO2_SiCOH",))]
row("Stage 3 4.5T", "low-k widths along column 0 (N wing): edge, N tip to P contact, edge", "10.5 tip-to-tip", lk0)
row("Stage 3 4.5T", "low-k widths along column 1 (P wing): edge, N contact to P tip, edge", "10.5 tip-to-tip", lk1)
row("Stage 3 4.5T", "N wing metal span (row 0)", "-", runs(p["y_col0_wing_level_row0"], METAL)[0])
row("Stage 3 4.5T", "P wing metal span (row 0)", "-", runs(p["y_col1_wing_level_row0"], METAL)[-1])
row("Stage 3 4.5T", "gate metal to contact along x (low-k spacer)", "-", [round(b - a, 2) for a, b in runs(p["x_through_N_fin_wing_level"], ("SiO2_SiCOH",))])
row("Stage 3 4.5T", "bridge metal along x", "-", runs(p["b_x_bridge_level"], ("Co",)))
row("Stage 3 4.5T", "via on N wing, along y", "-", runs(p["b_y_col0_via_level_row0"], ("Co",)))
row("Stage 3 4.5T", "via on P wing, along y", "-", runs(p["b_y_col1_via_level_row0"], ("Co",)))
row("Stage 3 4.5T", "stack at N-wing via (bottom to top)", "-", [s[0] for s in p["b_z_col0_via"]])
row("Stage 3 4.5T", "stack under N wing", "-", [s[0] for s in p["z_col0_under_N_wing"]])

j = load("stage4_T5p0_BPR.json")
p = j["profiles"]
row("Stage 4 5T BPR", "contact gaps along signal column: centre, border, ...", "18.75 centre / 39.75 border", gaps(runs(p["y_contact_level_signal_col"], METAL)))
row("Stage 4 5T BPR", "contact gaps along power column", "18.75 centre, merged at border", gaps(runs(p["y_contact_level_power_col"], METAL)))
yr = p["y_rail_level_gatecol"]
fins = runs(yr, ("Si_Xtal",)); rail = runs(yr, ("Co",)); lin = runs(yr, ("SiO2_HTO",))
row("Stage 4 5T BPR", "fin positions at rail level", "PN 37, cell-to-cell 58", gaps(fins))
row("Stage 4 5T BPR", "fin edge to rail trench (first full rail)", "13.25 (> 10)", [round(lin[1][0] - fins[1][1], 2), round(fins[2][0] - lin[2][1], 2)])
row("Stage 4 5T BPR", "rail trench width incl. liner", "31.5 (3 CD)", round(lin[2][1] - lin[1][0], 2))
row("Stage 4 5T BPR", "rail metal height", "-", [j["z_rail"][0], j["z_rail"][1]])
row("Stage 4 5T BPR", "stack at power via (bottom to top)", "-", [s[0] for s in p["z_power_col_at_border"]])
row("Stage 4 5T BPR", "stack over rail in signal column", "-", [s[0] for s in p["z_signal_col_at_border"]])

for m, w, pa, me in rows:
    print("%-15s | %-70s | paper: %-28s | model: %s" % (m, w, pa, me))
