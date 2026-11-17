#! python
# Stage 3: 4.5T mid-track handshake (Sherazi 2018, sec. 4.1).
#   - S/D contacts are made BEFORE the replacement metal gate, N first, then P.
#   - Each contact is patterned twice: one deep cavity over the fins, one shallow "wing"
#     toward the cell centre (the staircase), so the other device's fin stays buried.
#   - The nitride gate spacer is swapped for a low-k spacer in each module. The same film lines
#     the cavity tip, so where an N cavity meets a P cavity there is a double low-k spacer.
# Usage (inside SEMulator3D shell):  python build_stage3.py 4.5 [stopAfter]
_SEMulator3D_API_VERSION = '13.1'
import sys, os, json
from voxelModeler import *
from CModeler import mpiInitialize, mpiFinalize
from CDocument import ECoordinates
mpiInitialize(len(sys.argv), sys.argv)

here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, here)
import cells

tracks = float(sys.argv[1]) if len(sys.argv) > 1 else 4.5
stop_after = sys.argv[2] if len(sys.argv) > 2 else "end"
tag = ("T%.1f" % tracks).replace(".", "p")
outdir = os.path.join(here, "out")
os.makedirs(outdir, exist_ok=True)
zam = os.path.join(outdir, "stage3_%s.zam" % tag)

# ---- geometry (nm) -------------------------------------------------------
g = cells.cell(tracks)
nfins, pfins = cells.fins_stack(tracks, rows=3)
devs = cells.devices_stack(tracks, rows=3)
H = g["H"]; H3 = 3 * H
CPP = cells.CPP
X0, X1 = 0.0, 3 * CPP
LG, SP = 16.0, 6.0                      # gate length, nitride spacer thickness (assumptions)
gate_cx = [CPP / 2 + CPP * k for k in range(3)]
gate_x = [(c - LG / 2, c + LG / 2) for c in gate_cx]
sd_x = [(gate_x[k][1] + SP, gate_x[k + 1][0] - SP) for k in range(2)]
CUT = 18.0                              # gate cut width
CEXT = 5.0                              # contact metal extension past the fin (paper: M0A extension 5 nm)
LK = 5.25                               # low-k replacement spacer: two of them give the paper's 10.5 nm
                                        # handshake tip-to-tip (2018 paper, figs. 9 and 19)
DEX = CEXT + LK                         # cavity edge past the fin; metal ends CEXT past the fin once lined
gx = [(gate_x[k][1], gate_x[k + 1][0]) for k in range(2)]            # gate edge to gate edge
strips = [[(gate_x[k][1], sd_x[k][0]), (sd_x[k][1], gate_x[k + 1][0])] for k in range(2)]  # spacer footprints
WINGCOL = {"N": 0, "P": 1}             # staggered handshake: N wing in S/D column 0, P wing in column 1
for d in devs:
    q = [e for e in devs if e["row"] == d["row"] and e["type"] != d["type"]][0]
    d["deep"] = (d["y0"] - DEX, d["y1"] + DEX)
    if q["y0"] > d["y0"]:
        d["wing"] = (d["y1"] + DEX, q["y0"] - DEX)
    else:
        d["wing"] = (q["y1"] + DEX, d["y0"] - DEX)
Z_SI, HM, TRENCH, REVEAL = 150.0, 25.0, 100.0, 45.0
POLY, CAP, EPI = 60.0, 20.0, 8.0

res = 1.0
doc = VoxelDocument(zam, res, NANOMETERS, OVERWRITE)
doc.GetProperties().SetModelArea(ECoordinates)
doc.GetProperties().SetBoundaryCoordinates((X0, 0.0, X1, H3))
modeler = VoxelModeler(doc)
modeler.SetDefaultMaskUnitsPerMicron(1000)

def mask(name, rects):
    m = Mask(name, 1000)
    for (xa, ya, xb, yb) in rects:
        m.AddRectangle(xa, ya, xb, yb, NANOMETERS)
    modeler.AddActiveMask(m)
    return m

mask("FinAll", [(X0, a, X1, b) for a, b in sorted(nfins + pfins)])
mask("Gate", [(a, 0.0, b, H3) for a, b in gate_x])
mask("GateCut", [(X0, k * H - CUT / 2, X1, k * H + CUT / 2) for k in (1, 2)])
for T in ("N", "P"):
    dd = [d for d in devs if d["type"] == T]
    wc = WINGCOL[T]
    mask("Deep" + T, [(gx[k][0], d["deep"][0], gx[k][1], d["deep"][1]) for d in dd for k in range(2)])
    mask("Wing" + T, [(gx[wc][0], d["wing"][0], gx[wc][1], d["wing"][1]) for d in dd])
    mask("Open" + T, [(sd_x[k][0], d["deep"][0], sd_x[k][1], d["deep"][1]) for d in dd for k in range(2)])
    mask("SpDeep" + T, [(a_, d["deep"][0], b_, d["deep"][1]) for d in dd for k in range(2) for (a_, b_) in strips[k]])
    mask("SpWing" + T, [(a_, d["wing"][0], b_, d["wing"][1]) for d in dd for (a_, b_) in strips[wc]])
M0_T = 20.0                             # bridge metal thickness (assumption)
TRK = [(r * H + H / 2 - cells.MP / 4, r * H + H / 2 + cells.MP / 4) for r in range(3)]   # centre track, half-pitch wide
mask("Bridge", [(sd_x[0][0], a, sd_x[1][1], b) for a, b in TRK])
mask("BridgeVia", [(sd_x[k][0], a, sd_x[k][1], b) for a, b in TRK for k in range(2)])
modeler.CacheBuildBoundary()

# ---- materials -----------------------------------------------------------
matDoc = CMaterialData.MaterialDocument()
matDoc.openWithPasswords(os.path.join(here, "common", "16nmFINFET.vmpd"), "", "", False)
for mat in matDoc.getAllMaterials():
    doc.GetMaterials().Add(mat, DUP_ADD_ALLOW)
for dop in matDoc.getAllDopants():
    doc.GetDopants().Add(dop, DUP_ADD_ALLOW)

SI, SIN, STI = "/Silicon/Si_Xtal", "/Nitrides/Si3N4_LPCVD", "/Oxides/SiO2_HDP"
APOLY, ILD = "/Silicon/Si_Amorph", "/Oxides/SiO2_TEOS"
NEPI, PEPI = "/Silicon/SiC", "/Silicon/SiGe"
HIK, WF, W = "/Oxides/HfOx", "/Metals/TiN", "/Metals/W"
LINER = "/Metals/Ti"
LOWK, CCAP = "/Oxides/SiO2_SiCOH", "/Nitrides/Si3N4_PECVD"   # low-k spacer, contact cap
BRIDGE = [m.GetName() for m in doc.GetMaterials() if m.GetName().split("/")[-1] == "Co"][0]
print("BRIDGE material path:", BRIDGE)
seq = doc.GetModel()
def out(n): return OutputParams(seq, n, -1)

def mat_obj(path):
    nm = path.split("/")[-1]
    for m in doc.GetMaterials():
        if m.GetName().split("/")[-1] == nm:
            return m
    raise RuntimeError("material not found: " + path)

wafer = None
def report(label):
    info = wafer.GetMaterialBounds()
    parts = []
    for m in info.MaterialsFound():
        parts.append("%s=%.0f" % (m.GetName().split("/")[-1], info.VolumeForMaterial(m)))
    print("VOL[%s] %s" % (label, " ".join(sorted(parts))))

def done(label):
    if stop_after == label:
        print("STOP after", label)
        doc.GetProperties().SetModelBuildComplete(True)
        doc.Save()
        mpiFinalize()
        sys.exit(0)

def etch(mats, depth, mname, label, polarity=0):
    p = StepProperties("Etch"); p.type = "straight"; p.materials = mats; p.depth = depth; p.waferSide = "top"
    if hasattr(p, "selectivity"): p.selectivity = False
    if mname:
        p.advanced.setEnabled(True); p.mask.name = mname; p.mask.polarity = polarity
    wafer.Etch(p, out(label))

def deposit(kind, mat, thk, label):
    p = StepProperties("Deposit"); p.type = kind; p.material = mat; p.thickness = thk; p.waferSide = "top"
    wafer.Deposit(p, out(label))

def polish(stops, label):
    p = StepProperties("Polish"); p.type = "stopOnMaterial"; p.waferSide = "top"
    p.type.stopOnMaterial.materials = stops; p.type.stopOnMaterial.overpolish = 0
    wafer.Polish(p, out(label))

# ---- fins + STI ------------------------------------------------------------
p = StepProperties("WaferSetup"); p.wafer = "Wafer1"; p.material = SI; p.thickness = Z_SI
modeler.WaferSetup(p, out("Wafer setup"))
wafer = modeler.GetActiveWafers()["Wafer1"]
deposit("conformal", SIN, HM, "Hard mask nitride")
etch((SIN, SI), HM + TRENCH, "FinAll", "Fin etch", polarity=1)
deposit("conformal", STI, 90, "STI fill")
polish((SIN,), "STI CMP")
p = StepProperties("Etch"); p.type = "basic"; p.materials = (STI,); p.depth = HM + REVEAL; p.waferSide = "top"
if hasattr(p, "selectivity"): p.selectivity = False
wafer.Etch(p, out("STI recess"))
p = StepProperties("RemoveMaterials"); p.materials = (SIN,)
wafer.RemoveMaterials(p, out("Strip hard mask"))
report("fins"); done("fins")
full = wafer.GetMaterialBounds()
z_fin_top = full.BoundsForMaterial(mat_obj(SI)).GetBounds()[5]
z_sti_top = full.BoundsForMaterial(mat_obj(STI)).GetBounds()[5]

# ---- dummy gate + gate cut --------------------------------------------------
deposit("planarizing", APOLY, POLY, "Dummy gate poly")
deposit("straight", SIN, CAP, "Gate cap")
z_gate_top = wafer.GetMaterialBounds().BoundsForMaterial(mat_obj(APOLY)).GetBounds()[5]
etch((SIN, APOLY), CAP + (z_gate_top - z_sti_top) + 5, "Gate", "Gate etch", polarity=1)
etch((SIN, APOLY), CAP + (z_gate_top - z_sti_top) + 5, "GateCut", "Gate cut")
report("gate"); done("gate")

# ---- nitride spacer kept over the fins + ILD (paper fig. 10) -----------------
deposit("conformal", SIN, SP, "Nitride spacer (kept over fins)")
deposit("planarizing", ILD, 40, "ILD0")
polish((SIN,), "ILD CMP to nitride")
report("ild")
z_top = wafer.GetMaterialBounds().BoundsForMaterial(mat_obj(ILD)).GetBounds()[5]
Z_SHELF = z_fin_top + 15.0              # staircase floor: above the epi, keeps the other fin buried
Z_METAL = z_gate_top - 25.0             # contact metal top after recess (below the poly-open CMP level)
SH = z_top - Z_SHELF
RC = z_top - Z_METAL
print("LEVELS z_sti_top=%g z_fin_top=%g z_gate_top=%g z_top=%g shelf=%g metal_top=%g SH=%g RC=%g" % (
    z_sti_top, z_fin_top, z_gate_top, z_top, Z_SHELF, Z_METAL, SH, RC))
done("ild")

# ---- contact module, run once for N and once for P (paper figs. 11-16) -------
def contact_module(T, epi_mat):
    etch((ILD,), 300, "Deep" + T, T + " contact deep etch")
    etch((ILD,), SH, "Wing" + T, T + " contact shallow etch (staircase)")
    etch((SIN,), 300, "Open" + T, T + " fin open (nitride)")
    report(T + "-open"); done(T + "open")
    p = StepProperties("Growth"); p.type = "selectiveEpitaxy"; p.material = epi_mat; p.thickness = EPI; p.waferSide = "top"
    p.type.selectiveEpitaxy.targetMaterials = (SI,)
    wafer.Growth(p, out(T + " epi"))
    etch((SIN,), 300, "SpDeep" + T, T + " nitride spacer strip (deep)")
    etch((SIN,), SH, "SpWing" + T, T + " nitride spacer strip (wing)")
    deposit("conformal", LOWK, LK, T + " low-k spacer deposit")
    etch((LOWK,), LK + 2, None, T + " low-k spacer etch-back")
    report(T + "-lowk"); done(T + "lowk")
    deposit("conformal", LINER, 2, T + " contact liner")
    deposit("conformal", W, 10, T + " contact W")
    polish((ILD, SIN), T + " contact CMP")
    etch((W, LINER), RC, None, T + " contact recess")
    deposit("conformal", CCAP, 10, T + " contact cap")
    polish((ILD, SIN), T + " contact cap CMP")
    report(T + "-contact"); done(T + "contact")

contact_module("N", NEPI)
contact_module("P", PEPI)

# ---- replacement metal gate, after the contacts ---------------------------------
polish((APOLY,), "Poly-open CMP")
etch((APOLY,), 300, None, "Dummy gate removal")
report("poly-removed"); done("polyremoved")
deposit("conformal", HIK, 2, "High-k")
deposit("conformal", WF, 3, "Work function")
deposit("conformal", W, 8, "Gate W")
polish((ILD,), "RMG CMP")
etch((W, WF, HIK), 25, None, "Gate recess")
deposit("conformal", SIN, 14, "SAC cap")
polish((ILD,), "SAC cap CMP")
report("rmg"); done("rmg")

# ---- mid-track handshake bridge (paper fig. 17) --------------------------------
# Dual damascene on the centre routing track: trench in a new dielectric, vias through the
# contact caps (self-aligned: only the cap material is etched), cobalt fill, polish.
deposit("straight", ILD, M0_T, "Bridge dielectric")
etch((ILD,), M0_T, "Bridge", "Bridge trench etch")
etch((CCAP,), 60, "BridgeVia", "Via etch through contact caps")
deposit("conformal", BRIDGE, 8, "Bridge metal (Co)")
polish((ILD,), "Bridge CMP")
report("bridge")
modeler.SaveFinalStep()

# ---- measurements ---------------------------------------------------------------
USED = (SI, SIN, STI, APOLY, ILD, NEPI, PEPI, HIK, WF, W, LINER, LOWK, CCAP, BRIDGE)
objs_all = [(u.split("/")[-1], mat_obj(u)) for u in USED]

def dominant(box):
    info = wafer.GetMaterialBounds(box)
    best, bv = "air", 0.0
    for nm, m in objs_all:
        v = info.VolumeForMaterial(m)
        if v > bv:
            best, bv = nm, v
    return best

def profile(axis, x, y, z, lo, hi):
    out_, cur, start = [], None, lo
    v = lo
    while v < hi:
        if axis == "x":   box = [v, y, z - 0.25, v + res, y + res, z + 0.25]
        elif axis == "y": box = [x, v, z - 0.25, x + res, v + res, z + 0.25]
        else:             box = [x, y, v, x + res, y + res, v + res]
        d = dominant(box)
        if d != cur:
            if cur is not None:
                out_.append((cur, start, v))
            cur, start = d, v
        v += res
    out_.append((cur, start, hi))
    return out_

cx0 = (sd_x[0][0] + sd_x[0][1]) / 2.0 - 0.5
cx1 = (sd_x[1][0] + sd_x[1][1]) / 2.0 - 0.5
zw = (Z_SHELF + Z_METAL) / 2.0          # mid-height of the wing metal
n0 = [d for d in devs if d["type"] == "N" and d["row"] == 0][0]
p0 = [d for d in devs if d["type"] == "P" and d["row"] == 0][0]
yn = int((n0["y0"] + n0["y1"]) / 2.0)
yp = int((p0["y0"] + p0["y1"]) / 2.0)
ymid = int((n0["wing"][0] + n0["wing"][1]) / 2.0)
res_out = dict(tracks=tracks, cell_height=H, z_sti_top=z_sti_top, z_fin_top=z_fin_top, z_gate_top=z_gate_top,
               z_shelf=Z_SHELF, z_metal_top=Z_METAL, z_wing_probe=zw, gate_x=gate_x, sd_x=sd_x,
               devices=[dict(type=d["type"], row=d["row"], fin=(d["y0"], d["y1"]), deep=d["deep"], wing=d["wing"]) for d in devs])
prof = {}
prof["y_col0_wing_level_row0"] = profile("y", cx0, 0, zw, 0.0, H)
prof["y_col1_wing_level_row0"] = profile("y", cx1, 0, zw, 0.0, H)
prof["y_col0_fin_level_row0"] = profile("y", cx0, 0, z_fin_top - 10.0, 0.0, H)
prof["y_col1_fin_level_row0"] = profile("y", cx1, 0, z_fin_top - 10.0, 0.0, H)
prof["y_col0_wing_level_all"] = profile("y", cx0, 0, zw, 0.0, H3)
prof["y_col1_wing_level_all"] = profile("y", cx1, 0, zw, 0.0, H3)
prof["z_col0_N_fin"] = profile("z", cx0, yn, 0, z_sti_top - 5, z_gate_top + 1)
prof["z_col0_under_N_wing"] = profile("z", cx0, ymid, 0, z_sti_top - 5, z_gate_top + 1)
prof["z_col1_P_fin"] = profile("z", cx1, yp, 0, z_sti_top - 5, z_gate_top + 1)
prof["z_col1_under_P_wing"] = profile("z", cx1, ymid, 0, z_sti_top - 5, z_gate_top + 1)
prof["x_through_N_fin_wing_level"] = profile("x", 0, yn, zw, X0, X1)
prof["x_through_wing_mid"] = profile("x", 0, ymid, zw, X0, X1)
prof["y_gate_wing_level_all"] = profile("y", int(gate_cx[1]), 0, zw, 0.0, H3)
yt = int(H / 2)                          # centre of the handshake track, row 0
zv = (Z_METAL + z_gate_top) / 2.0       # via level, inside the contact cap
zb = z_gate_top + M0_T / 2.0            # bridge level
prof["b_z_col0_via"] = profile("z", cx0, yt - 2, 0, Z_SHELF - 5, z_gate_top + M0_T + 1)
prof["b_z_col1_via"] = profile("z", cx1, yt + 2, 0, Z_SHELF - 5, z_gate_top + M0_T + 1)
prof["b_z_over_gate"] = profile("z", int(gate_cx[1]), yt, 0, Z_SHELF - 5, z_gate_top + M0_T + 1)
prof["b_x_bridge_level"] = profile("x", 0, yt, zb, X0, X1)
prof["b_x_via_level"] = profile("x", 0, yt, zv, X0, X1)
prof["b_y_col0_bridge_level_row0"] = profile("y", cx0, 0, zb, 0.0, H)
prof["b_y_col0_via_level_row0"] = profile("y", cx0, 0, zv, 0.0, H)
prof["b_y_col1_via_level_row0"] = profile("y", cx1, 0, zv, 0.0, H)
res_out["profiles"] = prof
with open(os.path.join(outdir, "stage3_%s.json" % tag), "w") as f:
    json.dump(res_out, f, indent=1)
for k in sorted(prof):
    print("PROF", k, prof[k])

doc.GetProperties().SetModelBuildComplete(True)
doc.Save()
mpiFinalize()
