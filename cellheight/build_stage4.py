#! python
# Stage 4: 5T cell with buried power rail (BPR), Sherazi et al., SPIE 10962 (2019), section 2.
#   - CPP 45 nm, metal pitch 21 nm, cell height 105 nm, one fin per device.
#   - Power rail buried in the STI / silicon under the cell border, 3 CD wide (31.5 nm).
#   - Source contacts (M0A) extend to the border and drop a via onto the rail; contacts from the
#     two abutting cells merge there (allowed by the paper because both go to the same rail).
# Usage (inside SEMulator3D shell):  python build_stage4.py [stopAfter]
_SEMulator3D_API_VERSION = '13.1'
import sys, os, json
from voxelModeler import *
from CModeler import mpiInitialize, mpiFinalize
from CDocument import ECoordinates
mpiInitialize(len(sys.argv), sys.argv)

here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, here)
import cells

tracks = 5.0
stop_after = sys.argv[1] if len(sys.argv) > 1 else "end"
outdir = os.path.join(here, "out")
os.makedirs(outdir, exist_ok=True)
zam = os.path.join(outdir, "stage4_T5p0_BPR.zam")

# ---- geometry (nm) -------------------------------------------------------
g = cells.cell(tracks)
nfins, pfins = cells.fins_stack(tracks, rows=3)
devs = cells.devices_stack(tracks, rows=3)
H = g["H"]; H3 = 3 * H
CPP = 45.0                              # 2019 paper relaxes CPP from 42 to 45 nm
X0, X1 = 0.0, 3 * CPP
LG, SP = 16.0, 6.0                      # gate length, spacer thickness (assumptions)
gate_cx = [CPP / 2 + CPP * k for k in range(3)]
gate_x = [(c - LG / 2, c + LG / 2) for c in gate_cx]
sd_x = [(gate_x[k][1] + SP, gate_x[k + 1][0] - SP) for k in range(2)]
CUT = 18.0                              # gate cut width (EUV, paper)
CEXT = 5.0                              # contact extension over fin (CoFe)
CUME = cells.MP / 2 - cells.FIN_W / 2 + cells.MP / 4   # 13.25 nm: reach to the far edge of the next MINT track
BPR_W = 1.5 * cells.MP                  # 3 CD rail = 31.5 nm (trench width, liner included)
BPR_LINER = 3.0                         # dielectric liner around the rail (assumption)
BPR_BELOW = 20.0                        # rail depth below the fin-trench bottom (assumption)
BPR_CAP = 15.0                          # dielectric left over the rail under the STI surface (assumption)
VIA_IN = 5.0                            # via edge inset from the rail trench edge (assumption)
PWRCOL = 0                              # S/D column whose contacts are tied to the rails
borders = [k * H for k in range(4)]     # y = 0, H, 2H, 3H : rails alternate VSS / VDD
Z_SI, HM, TRENCH, REVEAL = 150.0, 25.0, 100.0, 45.0
POLY, CAP, EPI = 60.0, 20.0, 8.0

for d in devs:                          # which border each device's source contact runs to
    c = (d["y0"] + d["y1"]) / 2.0
    d["border"] = min(borders, key=lambda b: abs(b - c))
    # Signal contacts (2019 paper, figs. 6 and 16): CoFe = 5 nm over the fin on one side, CuMe = 13.25 nm
    # on the other so the contact sits fully under the next MINT track. This gives the paper's
    # 18.75 nm tip-to-tip at the cell centre and 39.75 nm at the border.
    d["con"] = (d["y0"] - CEXT, d["y1"] + CUME)
    d["pwr"] = (min(d["con"][0], d["border"]), max(d["con"][1], d["border"]))

res = 1.0
doc = VoxelDocument(zam, res, NANOMETERS, OVERWRITE)
doc.GetProperties().SetModelArea(ECoordinates)
doc.GetProperties().SetBoundaryCoordinates((X0, 0.0, X1, H3))
modeler = VoxelModeler(doc)
modeler.SetDefaultMaskUnitsPerMicron(1000)

def clipy(a, b):
    return max(0.0, a), min(H3, b)

def mask(name, rects):
    m = Mask(name, 1000)
    for (xa, ya, xb, yb) in rects:
        m.AddRectangle(xa, ya, xb, yb, NANOMETERS)
    modeler.AddActiveMask(m)
    return m

mask("FinAll", [(X0, a, X1, b) for a, b in sorted(nfins + pfins)])
mask("BPR", [(X0,) + (clipy(b - BPR_W / 2, b + BPR_W / 2)[0],) + (X1,) + (clipy(b - BPR_W / 2, b + BPR_W / 2)[1],) for b in borders])
mask("Gate", [(a, 0.0, b, H3) for a, b in gate_x])
mask("GateCut", [(X0, k * H - CUT / 2, X1, k * H + CUT / 2) for k in (1, 2)])
mask("SD", [(a, 0.0, b, H3) for a, b in sd_x])
mask("NBlock", [(X0, d["y0"] - 10, X1, d["y1"] + 10) for d in devs if d["type"] == "N"])
mask("PBlock", [(X0, d["y0"] - 10, X1, d["y1"] + 10) for d in devs if d["type"] == "P"])
con_rects = []
for d in devs:
    for k, (a, b) in enumerate(sd_x):
        ya, yb = d["pwr"] if k == PWRCOL else d["con"]
        con_rects.append((a - 2, ya, b + 2, yb))
mask("Contact", con_rects)
a, b = sd_x[PWRCOL]
mask("VBPR", [(a, clipy(bd - BPR_W / 2 + VIA_IN, bd + BPR_W / 2 - VIA_IN)[0], b, clipy(bd - BPR_W / 2 + VIA_IN, bd + BPR_W / 2 - VIA_IN)[1]) for bd in borders])
modeler.CacheBuildBoundary()

# ---- materials -----------------------------------------------------------
matDoc = CMaterialData.MaterialDocument()
matDoc.openWithPasswords(os.path.join(here, "common", "16nmFINFET.vmpd"), "", "", False)
for mat in matDoc.getAllMaterials():
    doc.GetMaterials().Add(mat, DUP_ADD_ALLOW)
for dop in matDoc.getAllDopants():
    doc.GetDopants().Add(dop, DUP_ADD_ALLOW)

def path(leaf):
    for m in doc.GetMaterials():
        if m.GetName().split("/")[-1] == leaf:
            return m.GetName()
    raise RuntimeError("material not found: " + leaf)

SI, SIN, STI = path("Si_Xtal"), path("Si3N4_LPCVD"), path("SiO2_HDP")
APOLY, ILD, BLOCK = path("Si_Amorph"), path("SiO2_TEOS"), path("SiO2_LTO")
NEPI, PEPI = path("SiC"), path("SiGe")
HIK, WF, W = path("HfOx"), path("TiN"), path("W")
LINER = path("Ti")
RLINER, RAIL = path("SiO2_HTO"), path("Co")     # rail liner oxide; rail metal kept distinct from contact W
seq = doc.GetModel()
def out(n): return OutputParams(seq, n, -1)

def mat_obj(p_):
    for m in doc.GetMaterials():
        if m.GetName() == p_:
            return m
    raise RuntimeError("material not found: " + p_)

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

def etch(mats, depth, mname, label, polarity=0, kind="straight"):
    p = StepProperties("Etch"); p.type = kind; p.materials = mats; p.depth = depth; p.waferSide = "top"
    if hasattr(p, "selectivity"): p.selectivity = False
    if mname:
        p.advanced.setEnabled(True); p.mask.name = mname; p.mask.polarity = polarity
    wafer.Etch(p, out(label))

def deposit(kind, mat, thk, label, mname=None):
    p = StepProperties("Deposit"); p.type = kind; p.material = mat; p.thickness = thk; p.waferSide = "top"
    if mname:
        p.advanced.setEnabled(True); p.mask.name = mname; p.mask.polarity = 0
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

# ---- buried power rail, before fin reveal --------------------------------------
etch((STI, SI), HM + TRENCH + BPR_BELOW, "BPR", "BPR trench etch")
deposit("conformal", RLINER, BPR_LINER, "BPR liner oxide")
deposit("conformal", RAIL, 16, "BPR metal fill")
polish((SIN,), "BPR CMP")
etch((RAIL,), HM + REVEAL + BPR_CAP, None, "BPR metal recess")
deposit("conformal", STI, 16, "BPR cap oxide")
polish((SIN,), "BPR cap CMP")
report("bpr"); done("bpr")

etch((STI, RLINER), HM + REVEAL, None, "STI recess (fin reveal)", kind="basic")
p = StepProperties("RemoveMaterials"); p.materials = (SIN,)
wafer.RemoveMaterials(p, out("Strip hard mask"))
report("fins"); done("fins")
full = wafer.GetMaterialBounds()
z_fin_top = full.BoundsForMaterial(mat_obj(SI)).GetBounds()[5]
z_sti_top = full.BoundsForMaterial(mat_obj(STI)).GetBounds()[5]
rb = full.BoundsForMaterial(mat_obj(RAIL)).GetBounds()
z_rail_bot, z_rail_top = rb[2], rb[5]
print("LEVELS z_sti_top=%g z_fin_top=%g rail z=%g..%g" % (z_sti_top, z_fin_top, z_rail_bot, z_rail_top))

# ---- dummy gate + gate cut --------------------------------------------------
deposit("planarizing", APOLY, POLY, "Dummy gate poly")
deposit("straight", SIN, CAP, "Gate cap")
z_gate_top = wafer.GetMaterialBounds().BoundsForMaterial(mat_obj(APOLY)).GetBounds()[5]
etch((SIN, APOLY), CAP + (z_gate_top - z_sti_top) + 5, "Gate", "Gate etch", polarity=1)
etch((SIN, APOLY), CAP + (z_gate_top - z_sti_top) + 5, "GateCut", "Gate cut")
report("gate"); done("gate")

# ---- spacers + epi -------------------------------------------------------------
deposit("conformal", SIN, SP, "Spacer deposit")
etch((SIN,), 200, "SD", "Spacer open in S/D")
def epi(block_mask, epi_mat, label):
    deposit("conformal", BLOCK, 8, label + " block", mname=block_mask)
    p = StepProperties("Growth"); p.type = "selectiveEpitaxy"; p.material = epi_mat; p.thickness = EPI; p.waferSide = "top"
    p.type.selectiveEpitaxy.targetMaterials = (SI,)
    wafer.Growth(p, out(label + " epi"))
    p = StepProperties("RemoveMaterials"); p.materials = (BLOCK,)
    wafer.RemoveMaterials(p, out(label + " unblock"))
epi("PBlock", NEPI, "N")
epi("NBlock", PEPI, "P")
report("epi"); done("epi")

# ---- ILD + replacement metal gate ----------------------------------------------
deposit("planarizing", ILD, 40, "ILD0")
polish((APOLY,), "ILD CMP to poly")
etch((APOLY,), 300, None, "Dummy gate removal")
deposit("conformal", HIK, 2, "High-k")
deposit("conformal", WF, 3, "Work function")
deposit("conformal", W, 8, "Gate W")
polish((ILD,), "RMG CMP")
etch((W, WF, HIK), 25, None, "Gate recess")
deposit("conformal", SIN, 14, "SAC cap")
polish((ILD,), "SAC cap CMP")
report("rmg"); done("rmg")

# ---- contacts, with the power contacts dropping a via onto the rail ----------------
etch((ILD,), 250, "Contact", "Contact etch (self-aligned)")
etch((STI,), BPR_CAP + 25, "VBPR", "Via to buried rail (through rail cap)")
report("contact-open"); done("contactopen")
deposit("conformal", LINER, 2, "Contact liner")
deposit("conformal", W, 9, "Contact W")
polish((ILD, SIN), "Contact CMP")
report("contact")
modeler.SaveFinalStep()

# ---- measurements ---------------------------------------------------------------
USED = (SI, SIN, STI, ILD, NEPI, PEPI, HIK, WF, W, LINER, RLINER, RAIL)
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

cxp = int((sd_x[PWRCOL][0] + sd_x[PWRCOL][1]) / 2.0)        # power column
cxs = int((sd_x[1 - PWRCOL][0] + sd_x[1 - PWRCOL][1]) / 2.0)  # signal column
zr = (z_rail_bot + z_rail_top) / 2.0
zc = z_fin_top + 20.0
yb = int(H)                              # border between row 0 and row 1 (P-P, VDD rail)
zlo, zhi = z_rail_bot - 8, z_gate_top + 1
res_out = dict(tracks=tracks, cell_height=H, cpp=CPP, z_sti_top=z_sti_top, z_fin_top=z_fin_top, z_gate_top=z_gate_top,
               z_rail=(z_rail_bot, z_rail_top), gate_x=gate_x, sd_x=sd_x, borders=borders,
               devices=[dict(type=d["type"], row=d["row"], fin=(d["y0"], d["y1"]), border=d["border"]) for d in devs])
prof = {}
prof["y_rail_level_gatecol"] = profile("y", int(gate_cx[1]), 0, zr, 0.0, H3)
prof["y_contact_level_power_col"] = profile("y", cxp, 0, zc, 0.0, H3)
prof["y_contact_level_signal_col"] = profile("y", cxs, 0, zc, 0.0, H3)
prof["y_epi_level_power_col"] = profile("y", cxp, 0, z_fin_top + 3.0, 0.0, H3)
prof["z_power_col_at_border"] = profile("z", cxp, yb, 0, zlo, zhi)
prof["z_signal_col_at_border"] = profile("z", cxs, yb, 0, zlo, zhi)
prof["z_gate_over_rail"] = profile("z", int(gate_cx[1]), yb - 12, 0, zlo, zhi)
prof["z_power_col_between_fin_and_rail"] = profile("z", cxp, yb - 22, 0, zlo, zhi)
prof["x_rail_level_at_border"] = profile("x", 0, yb, zr, X0, X1)
prof["x_contact_level_through_P_fin"] = profile("x", 0, int((devs[1]["y0"] + devs[1]["y1"]) / 2.0), zc, X0, X1)
res_out["profiles"] = prof
with open(os.path.join(outdir, "stage4_T5p0_BPR.json"), "w") as f:
    json.dump(res_out, f, indent=1)
for k in sorted(prof):
    print("PROF", k, prof[k])

doc.GetProperties().SetModelBuildComplete(True)
doc.Save()
mpiFinalize()
