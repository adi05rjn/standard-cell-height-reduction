# Design-level checks that involve no process modelling at all.
#   A. Arithmetic: every paper dimension recomputed by formula from the pitches, and compared with
#      what cellheight/cells.py produces.
#   B. Layout: the GDS files in manager/ are read back with a separate reader and the mask gaps are
#      measured in 2D, to 0.01 nm.
# Plain Python, no SEMulator3D license:   python check_design.py
import os, sys, struct, json

here = os.path.dirname(os.path.abspath(__file__))
root = os.path.abspath(os.path.join(here, ".."))
sys.path.insert(0, os.path.join(root, "cellheight"))
import cells

results, fails = [], 0
def check(group, what, got, want, tol=0.005):
    global fails
    ok = abs(got - want) <= tol
    fails += (not ok)
    results.append(dict(group=group, what=what, got=round(got, 3), want=want, ok=ok))
    print("%-4s %-11s %-62s got %8.2f  want %8.2f" % ("ok" if ok else "FAIL", group, what, got, want))

# ---- A. arithmetic -----------------------------------------------------------------
MP, FW, FP = 21.0, 5.0, 21.0            # metal pitch, fin width, fin pitch (papers)
COFE = 5.0                              # contact extension over fin (papers)
# Paper values: (tracks, fins per device, N-P fin gap, border fin gap, N-P contact tip gap, border gap beside a power tap)
PAPER = [
    (6.5, 2, 41.5, 43.0, None, None),   # 2018 fig. 2
    (5.5, 2, 20.5, 43.0, 10.5, 16.5),   # 2018 fig. 4 (fin gap labelled 21.5 in the figure, see RESULTS.md section 3)
    (4.5, 1, 40.0, 44.5, 30.0, 17.25),  # 2018 figs. 6, 7
    (5.0, 1, 37.0, 58.0, None, None),   # 2019 fig. 6
]
for t, nf, pn, out, tip, tap in PAPER:
    H = t * MP
    span = (nf - 1) * FP + FW
    g = cells.cell(t)
    check("arithmetic", "%.1fT cell height = tracks x 21" % t, g["H"], H)
    check("arithmetic", "%.1fT heights close: 2 x device + N-P gap + border gap" % t, 2 * span + pn + out, H)
    check("arithmetic", "%.1fT border fin gap from cells.py" % t, g["out"], out)
    n, p = cells.fins_in_cell(t)
    check("arithmetic", "%.1fT N-P fin gap from cells.py" % t, p[0][0] - n[-1][1], pn)
    if tip is not None:
        check("arithmetic", "%.1fT N-P contact gap = fin gap - 2 x 5" % t, pn - 2 * COFE, tip)
    if tap is not None:
        check("arithmetic", "%.1fT border gap beside power tap = border gap / 2 - 5" % t, out / 2 - COFE, tap)
# 5T buried rail (2019 fig. 6)
RAIL_W = 1.5 * MP
CUME = MP / 2 - FW / 2 + MP / 4
check("arithmetic", "5T fins two pitches apart: 2 x 21 - 5", 2 * FP - FW, 37.0)
check("arithmetic", "5T border: three pitches minus fin: 3 x 21 - 5", 3 * FP - FW, 58.0)
check("arithmetic", "5T rail width = 3 CD = 1.5 x 21", RAIL_W, 31.5)
check("arithmetic", "5T fin to rail = border gap / 2 - rail / 2", 58.0 / 2 - RAIL_W / 2, 13.25)
check("arithmetic", "5T contact reach under next track (CuMe)", CUME, 13.25)
check("arithmetic", "5T contact gap at centre = 37 - 5 - 13.25", 37.0 - COFE - CUME, 18.75)
check("arithmetic", "5T contact gap at border = 58 - 5 - 13.25", 58.0 - COFE - CUME, 39.75)
check("arithmetic", "4.5T handshake gap = 2 x 5.25 low-k", 2 * 5.25, 10.5)

# ---- B. layout read-back ---------------------------------------------------------------
def read_gds(path):
    """Minimal GDSII reader: returns {layer: [(x0, y0, x1, y1) in nm]} and the database unit in nm."""
    data = open(path, "rb").read()
    pos, layer, dbu_nm, out = 0, None, None, {}
    def real8(b):
        sign = -1 if b[0] & 0x80 else 1
        return sign * (int.from_bytes(b[1:], "big") / float(1 << 56)) * 16.0 ** ((b[0] & 0x7F) - 64)
    while pos < len(data):
        length, rtype, dtype = struct.unpack(">HBB", data[pos:pos + 4])
        body = data[pos + 4:pos + length]
        if rtype == 0x03:   dbu_nm = real8(body[8:16]) * 1e9
        elif rtype == 0x0D: layer = struct.unpack(">h", body)[0]
        elif rtype == 0x10:
            xy = struct.unpack(">%di" % (len(body) // 4), body)
            xs, ys = xy[0::2], xy[1::2]
            out.setdefault(layer, []).append((min(xs) * dbu_nm, min(ys) * dbu_nm, max(xs) * dbu_nm, max(ys) * dbu_nm))
        pos += length
    return out, dbu_nm

def layers(name):
    folder = os.path.join(root, "manager", name)
    shapes, dbu = read_gds(os.path.join(folder, name + ".gds"))
    names = {}
    for line in open(os.path.join(folder, name + "-layermap.txt")):
        parts = line.split()
        if len(parts) >= 2: names[parts[0]] = int(parts[1])
    return dict((n, sorted(shapes.get(l, []))) for n, l in names.items()), dbu

def ygaps(rects, x):
    """Gaps along y between shapes that cover the given x."""
    ys = sorted((r[1], r[3]) for r in rects if r[0] <= x <= r[2])
    merged = []
    for a, b in ys:
        if merged and a <= merged[-1][1] + 1e-9: merged[-1][1] = max(merged[-1][1], b)
        else: merged.append([a, b])
    return [round(merged[i + 1][0] - merged[i][1], 4) for i in range(len(merged) - 1)]

def xcols(rects):
    return sorted(set((r[0], r[2]) for r in rects))

for name, t, pn, out, tip, border_plain, tap, cpp in (
        ("cell_6p5T", 6.5, 41.5, 43.0, 31.5, 33.0, 16.5, 42.0),
        ("cell_5p5T", 5.5, 20.5, 43.0, 10.5, 33.0, 16.5, 42.0),
        ("cell_4p5T", 4.5, 40.0, 44.5, 30.0, 34.5, 17.25, 42.0)):
    L, dbu = layers(name)
    check("layout", "%s database unit (nm)" % name, dbu, 0.01, 1e-9)
    fins = ygaps(L["FinAll"], 1.0)
    check("layout", "%s smallest N-P fin gap on the Fin layer" % name, min(g for g in fins if abs(g - pn) < 1.0), pn)
    check("layout", "%s border fin gap on the Fin layer" % name, min(g for g in fins if abs(g - out) < 1.0), out)
    gates = xcols(L["Gate"])
    check("layout", "%s gate length" % name, gates[0][1] - gates[0][0], 16.0)
    check("layout", "%s gate pitch" % name, gates[1][0] - gates[0][0], cpp)
    cut = L["GateCut"][0]
    check("layout", "%s gate cut width" % name, cut[3] - cut[1], 18.0)
    cols = xcols(L["Contact"])
    xa, xb = (cols[0][0] + cols[0][1]) / 2, (cols[-1][0] + cols[-1][1]) / 2
    ga, gb = ygaps(L["Contact"], xa), ygaps(L["Contact"], xb)
    check("layout", "%s N-P contact gap, plain column" % name, min(gb), tip)
    check("layout", "%s border contact gap, plain column" % name, max(gb), border_plain)
    check("layout", "%s border gap beside power tap" % name, min(g for g in ga if g > tip + 1 or tip > 20) if tip < 20 else min(ga), tap)

L, dbu = layers("bpr_5T")
fins = ygaps(L["FinAll"], 1.0)
check("layout", "bpr_5T N-P fin gap", min(fins), 37.0)
check("layout", "bpr_5T border fin gap", max(fins), 58.0)
gates = xcols(L["Gate"])
check("layout", "bpr_5T gate pitch", gates[1][0] - gates[0][0], 45.0)
rail = [r for r in L["BPR"] if r[1] > 1 and r[3] < 300][0]
check("layout", "bpr_5T rail width", rail[3] - rail[1], 31.5)
finr = [r for r in L["FinAll"] if r[3] <= rail[1]][-1]
check("layout", "bpr_5T fin to rail", rail[1] - finr[3], 13.25)
cols = xcols(L["Contact"])
gs = ygaps(L["Contact"], (cols[-1][0] + cols[-1][1]) / 2)
check("layout", "bpr_5T contact gap at cell centre, signal column", min(gs), 18.75)
check("layout", "bpr_5T contact gap at cell border, signal column", max(gs), 39.75)
gp = ygaps(L["Contact"], (cols[0][0] + cols[0][1]) / 2)
check("layout", "bpr_5T power contacts merged at borders: only centre gaps remain", max(gp), 18.75)
via = [r for r in L["VBPR"] if r[1] > 1 and r[3] < 300][0]
check("layout", "bpr_5T rail via inset from rail edge", via[1] - rail[1], 5.0)

L, dbu = layers("handshake_4p5T")
n_deep = sorted(set((r[1], r[3]) for r in L["DeepN"]))[0]
p_deep = sorted(set((r[1], r[3]) for r in L["DeepP"]))[0]
n_wing = sorted(set((r[1], r[3]) for r in L["WingN"]))[0]
p_wing = sorted(set((r[1], r[3]) for r in L["WingP"]))[0]
check("layout", "handshake N wing starts where the N cavity ends", n_wing[0] - n_deep[1], 0.0)
check("layout", "handshake N wing ends where the P cavity starts", p_deep[0] - n_wing[1], 0.0)
check("layout", "handshake P wing spans from N cavity to P cavity", (p_wing[0] - n_deep[1]) + (p_deep[0] - p_wing[1]), 0.0)
check("layout", "handshake cavity reaches 10.25 nm past the fin (5 + 5.25)", 22.25 - n_deep[0], 10.25)
check("layout", "handshake metal tip gap implied: 2 x 5.25 low-k", 2 * 5.25, 10.5)
br = sorted(L["Bridge"])[0]
check("layout", "handshake bridge width (half the metal pitch)", br[3] - br[1], 10.5)
check("layout", "handshake bridge centred on the cell", (br[1] + br[3]) / 2, 94.5 / 2)

# ---- nanosheet cell (2018 fig. 8) ---------------------------------------------------------
gn = cells.cell("4.5NS")
check("arithmetic", "nanosheet heights close: 2 x 11 + 40 + 32.5", 2 * 11.0 + 40.0 + 32.5, 94.5)
check("arithmetic", "nanosheet border sheet gap from cells.py", gn["out"], 32.5)
check("arithmetic", "nanosheet N-P contact gap = 40 - 2 x 5", 40.0 - 2 * COFE, 30.0)
check("arithmetic", "nanosheet plain border contact gap = 32.5 - 2 x 5", 32.5 - 2 * COFE, 22.5)
L, dbu = layers("nanosheet_4p5T")
sh = ygaps(L["Sheet"], 1.0)
one = L["Sheet"][0]
check("layout", "nanosheet sheet width", one[3] - one[1], 11.0)
check("layout", "nanosheet N-P sheet gap", max(sh), 40.0)
check("layout", "nanosheet border sheet gap", min(sh), 32.5)
cols = xcols(L["Contact"])
ga = ygaps(L["Contact"], (cols[0][0] + cols[0][1]) / 2)
gb = ygaps(L["Contact"], (cols[-1][0] + cols[-1][1]) / 2)
check("layout", "nanosheet N-P contact gap", max(gb), 30.0)
check("layout", "nanosheet border contact gap, plain column", min(gb), 22.5)
check("layout", "nanosheet border gap beside power tap (paper 10.5)", min(ga), 10.5)

print("\n%d checks, %d failed" % (len(results), fails))
with open(os.path.join(here, "check_design.json"), "w") as f:
    json.dump(dict(checks=len(results), failed=fails, results=results), f, indent=1)
sys.exit(1 if fails else 0)
