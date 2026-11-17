#! python
# Stage 5: 4.5T cell with stacked nanosheets (gate-all-around), Sherazi et al., SPIE 10588 (2018), fig. 8.
#   - One stack of three silicon sheets per device, 11 nm wide, in place of the single fin.
#   - N-to-P separation 40 nm, which leaves 32.5 nm between the sheets of neighbouring cells.
#   - Sacrificial SiGe between the sheets, inner spacers at the source/drain, sheets released after
#     the dummy gate is removed, gate stack wrapped all the way around each sheet.
# Usage (inside the SEMulator3D shell):  python build_stage5.py [stopAfter]
_SEMulator3D_API_VERSION = '13.1'
import sys, os, json
from voxelModeler import *
from CModeler import mpiInitialize, mpiFinalize
from CDocument import ECoordinates
mpiInitialize(len(sys.argv), sys.argv)

here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, here)
import cells

KEY = "4.5NS"
stop_after = sys.argv[1] if len(sys.argv) > 1 else "end"
outdir = os.path.join(here, "out")
os.makedirs(outdir, exist_ok=True)
zam = os.path.join(outdir, "stage5_T4p5_nanosheet.zam")

# ---- geometry (nm) -------------------------------------------------------
g = cells.cell(KEY)
nsheets, psheets = cells.fins_stack(KEY, rows=3)
devs = cells.devices_stack(KEY, rows=3)
H = g["H"]; H3 = 3 * H
CPP = cells.CPP
X0, X1 = 0.0, 3 * CPP
LG, SP, CUT, CEXT = 16.0, 6.0, 18.0, 5.0        # gate length, spacer, gate cut, contact extension past the sheet
T_SHEET, T_SAC, NSHEET = 6.0, 10.0, 3           # sheet thickness, sacrificial thickness, sheets per stack (assumptions)
STACK = NSHEET * (T_SHEET + T_SAC)              # 48 nm
INDENT = 7.0                                    # sacrificial indent: 6 nm under the spacer plus 1 nm
ISP_T = 5.0                                     # inner spacer film; the etch-back leaves a plug about 5 nm thick
TAP_GAP = 10.5                                  # paper: tip-to-tip at the cell boundary beside a power tap
TAP = g["out"] - CEXT - TAP_GAP                 # 17 nm: power-tap contact reach past the sheet
gate_cx = [CPP / 2 + CPP * k for k in range(3)]
gate_x = [(c - LG / 2, c + LG / 2) for c in gate_cx]
sd_x = [(gate_x[k][1] + SP, gate_x[k + 1][0] - SP) for k in range(2)]
sd_all = [(X0, gate_x[0][0] - SP)] + sd_x + [(gate_x[2][1] + SP, X1)]   # with the half columns at the model edges
HM, BELOW = 25.0, 52.0                          # hard mask, trench depth below the stack
POLY, CAP = 60.0, 20.0
EPI = 10.0                                      # 10 nm, not 8: growth from the sheet ends must merge across the 10 nm gaps between sheets

def con_y(d, k):
    ya, yb = d["y0"] - CEXT, d["y1"] + CEXT
    if k == 0 and d["row"] == 1:                # power tap: middle row, column 0
        c = (d["y0"] + d["y1"]) / 2.0
        if abs(c - H) < abs(c - 2 * H): ya = d["y0"] - TAP if c > H else ya; yb = d["y1"] + TAP if c < H else yb
        else:                           ya = d["y0"] - TAP if c > 2 * H else ya; yb = d["y1"] + TAP if c < 2 * H else yb
    return ya, yb

res = 1.0
UPM = 100000                                    # mask units per micron: 0.01 nm, so quarter-nanometre edges are exact
doc = VoxelDocument(zam, res, NANOMETERS, OVERWRITE)
doc.GetProperties().SetModelArea(ECoordinates)
doc.GetProperties().SetBoundaryCoordinates((X0, 0.0, X1, H3))
modeler = VoxelModeler(doc)
modeler.SetDefaultMaskUnitsPerMicron(UPM)

def mask(name, rects):
    m = Mask(name, UPM)
    for (xa, ya, xb, yb) in rects:
        m.AddRectangle(xa, ya, xb, yb, NANOMETERS)
    modeler.AddActiveMask(m)

mask("Sheet", [(X0, a, X1, b) for a, b in sorted(nsheets + psheets)])
mask("Gate", [(a, 0.0, b, H3) for a, b in gate_x])
mask("GateCut", [(X0, k * H - CUT / 2, X1, k * H + CUT / 2) for k in (1, 2)])
mask("SD", [(a, 0.0, b, H3) for a, b in sd_all])
mask("NBlock", [(X0, d["y0"] - 12, X1, d["y1"] + 12) for d in devs if d["type"] == "N"])
mask("PBlock", [(X0, d["y0"] - 12, X1, d["y1"] + 12) for d in devs if d["type"] == "P"])
mask("Contact", [(a - 2, con_y(d, k)[0], b + 2, con_y(d, k)[1]) for d in devs for k, (a, b) in enumerate(sd_x)])
modeler.CacheBuildBoundary()

# ---- materials -----------------------------------------------------------
matDoc = CMaterialData.MaterialDocument()
matDoc.openWithPasswords(os.path.join(here, "common", "cellheight_nanosheet.vmpd"), "", "", False)
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
HIK, WF, W, LINER = path("HfOx"), path("TiN"), path("W"), path("Ti")
SAC, ISP = path("SiGe_Sacrificial"), path("Si3N4_MLD")          # sacrificial layers, inner spacer
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
    parts = ["%s=%.0f" % (m.GetName().split("/")[-1], info.VolumeForMaterial(m)) for m in info.MaterialsFound()]
    print("VOL[%s] %s" % (label, " ".join(sorted(parts))))

def done(label):
    if stop_after == label:
        print("STOP after", label)
        doc.GetProperties().SetModelBuildComplete(True)
        doc.Save(); mpiFinalize(); sys.exit(0)

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

def remove(mats, label):
    p = StepProperties("RemoveMaterials"); p.materials = mats
    wafer.RemoveMaterials(p, out(label))

USED = (SI, SIN, STI, APOLY, ILD, NEPI, PEPI, HIK, WF, W, LINER, SAC, ISP)
def dominant(box):
    info = wafer.GetMaterialBounds(box)
    best, bv = "air", 0.0
    for u in USED:
        v = info.VolumeForMaterial(mat_obj(u))
        if v > bv:
            best, bv = u.split("/")[-1], v
    return best
def profile(axis, x, y, z, lo, hi):
    out_, cur, start, v = [], None, lo, lo
    while v < hi:
        if axis == "x":   box = [v, y, z - 0.25, v + res, y + res, z + 0.25]
        elif axis == "y": box = [x, v, z - 0.25, x + res, v + res, z + 0.25]
        else:             box = [x, y, v, x + res, y + res, v + res]
        d = dominant(box)
        if d != cur:
            if cur is not None: out_.append((cur, start, v))
            cur, start = d, v
        v += res
    out_.append((cur, start, hi))
    return out_

# ---- sheet stack, patterning and STI ---------------------------------------------
p = StepProperties("WaferSetup"); p.wafer = "Wafer1"; p.material = SI; p.thickness = 150.0
modeler.WaferSetup(p, out("Wafer setup"))
wafer = modeler.GetActiveWafers()["Wafer1"]
for i in range(NSHEET):
    deposit("straight", SAC, T_SAC, "Sacrificial SiGe %d" % (i + 1))
    deposit("straight", SI, T_SHEET, "Silicon sheet %d" % (i + 1))
z_top_si = wafer.GetMaterialBounds().BoundsForMaterial(mat_obj(SI)).GetBounds()[5]     # top of the stack
z_stack0 = z_top_si - STACK                                                           # bottom of the stack
deposit("conformal", SIN, HM, "Hard mask nitride")
etch((SIN, SI, SAC), HM + STACK + BELOW, "Sheet", "Sheet stack etch", polarity=1)
deposit("conformal", STI, 90, "STI fill")
polish((SIN,), "STI CMP")
etch((STI,), HM + STACK, None, "STI recess (stack reveal)", kind="basic")
remove((SIN,), "Strip hard mask")
report("stack"); done("stack")
z_sti_top = wafer.GetMaterialBounds().BoundsForMaterial(mat_obj(STI)).GetBounds()[5]
print("LEVELS stack bottom=%g stack top=%g sti top=%g" % (z_stack0, z_top_si, z_sti_top))

# ---- dummy gate ------------------------------------------------------------------
deposit("planarizing", APOLY, POLY, "Dummy gate poly")
deposit("straight", SIN, CAP, "Gate cap")
z_gate_top = wafer.GetMaterialBounds().BoundsForMaterial(mat_obj(APOLY)).GetBounds()[5]
etch((SIN, APOLY), CAP + (z_gate_top - z_sti_top) + 5, "Gate", "Gate etch", polarity=1)
etch((SIN, APOLY), CAP + (z_gate_top - z_sti_top) + 5, "GateCut", "Gate cut")
report("gate"); done("gate")

# ---- spacer, source/drain recess, inner spacers -----------------------------------
deposit("conformal", SIN, SP, "Spacer deposit")
etch((SIN,), 200, "SD", "Spacer open in S/D")
etch((SI, SAC), STACK, "SD", "S/D recess (cut the stack)")
etch((SAC,), INDENT, None, "Sacrificial indent", kind="basic")
n0 = [d for d in devs if d["type"] == "N" and d["row"] == 0][0]
p0 = [d for d in devs if d["type"] == "P" and d["row"] == 0][0]
yn, yp = int((n0["y0"] + n0["y1"]) / 2.0), int((p0["y0"] + p0["y1"]) / 2.0)
z_sac2 = z_stack0 + (T_SAC + T_SHEET) + T_SAC / 2.0            # middle of the second sacrificial layer
z_sh2 = z_stack0 + (T_SAC + T_SHEET) + T_SAC + T_SHEET / 2.0   # middle of the second sheet
print("CHK after indent, x at sacrificial level:", profile("x", 0, yn, z_sac2, X0, X1))
print("CHK after indent, x at sheet level:      ", profile("x", 0, yn, z_sh2, X0, X1))
done("indent")
deposit("conformal", ISP, ISP_T, "Inner spacer deposit")
etch((ISP,), ISP_T, None, "Inner spacer etch-back", kind="basic")
print("CHK after inner spacer, x at sacrificial level:", profile("x", 0, yn, z_sac2, X0, X1))
print("CHK after inner spacer, x at sheet level:      ", profile("x", 0, yn, z_sh2, X0, X1))
report("innerspacer"); done("innerspacer")

# ---- source/drain epi ---------------------------------------------------------------
def epi(block_mask, epi_mat, label):
    deposit("conformal", BLOCK, 8, label + " block", mname=block_mask)
    p = StepProperties("Growth"); p.type = "selectiveEpitaxy"; p.material = epi_mat; p.thickness = EPI; p.waferSide = "top"
    p.type.selectiveEpitaxy.targetMaterials = (SI,)
    wafer.Growth(p, out(label + " epi"))
    remove((BLOCK,), label + " unblock")
epi("PBlock", NEPI, "N")
epi("NBlock", PEPI, "P")
cx0, cx1 = (sd_x[0][0] + sd_x[0][1]) / 2.0 - 0.5, (sd_x[1][0] + sd_x[1][1]) / 2.0 - 0.5
print("CHK epi, z through N cavity:", profile("z", cx0, yn, 0, z_stack0 - 8, z_top_si + 14))
print("CHK epi, x at sheet level:  ", profile("x", 0, yn, z_sh2, X0, X1))
report("epi"); done("epi")

# ---- dummy gate removal, sheet release, gate all around --------------------------------
deposit("planarizing", ILD, 40, "ILD0")
polish((APOLY,), "ILD CMP to poly")
etch((APOLY,), 300, None, "Dummy gate removal")
remove((SAC,), "Sheet release (remove sacrificial SiGe)")
print("CHK released, z under gate:", profile("z", int(gate_cx[1]), yn, 0, z_stack0 - 4, z_top_si + 6))
report("release"); done("release")
deposit("conformal", HIK, 2, "High-k")
deposit("conformal", WF, 3, "Work function")
deposit("conformal", W, 8, "Gate W")
polish((ILD,), "RMG CMP")
etch((W, WF, HIK), 25, None, "Gate recess")
deposit("conformal", SIN, 14, "SAC cap")
polish((ILD,), "SAC cap CMP")
report("rmg"); done("rmg")

# ---- contacts ---------------------------------------------------------------------------
etch((ILD,), 250, "Contact", "Contact etch (self-aligned)")
deposit("conformal", LINER, 2, "Contact liner")
deposit("conformal", W, 8, "Contact W")
polish((ILD, SIN), "Contact CMP")
report("contact")
modeler.SaveFinalStep()

# ---- measurements ---------------------------------------------------------------
def runs(prof, names):
    r = []
    for n, a, b in prof:
        if n in names:
            if r and abs(r[-1][1] - a) < 1e-6: r[-1] = (r[-1][0], b)
            else: r.append((a, b))
    return r
def gaps(r): return [round(r[i + 1][0] - r[i][1], 1) for i in range(len(r) - 1)]

gx = int(gate_cx[1])
zc = z_top_si + 20.0
res_out = dict(cell=KEY, cell_height=H, sheet_width=g["w"], sheets=NSHEET, t_sheet=T_SHEET, t_sacrificial=T_SAC,
               z_stack=(z_stack0, z_top_si), z_sti_top=z_sti_top, z_gate_top=z_gate_top, gate_x=gate_x, sd_x=sd_x,
               devices=[dict(type=d["type"], row=d["row"], sheet=(d["y0"], d["y1"])) for d in devs])
prof = {}
prof["z_under_gate_through_N_sheets"] = profile("z", gx, yn, 0, z_stack0 - 4, z_top_si + 14)
prof["y_across_N_sheet_under_gate"] = profile("y", gx, 0, z_sh2, n0["y0"] - 12, n0["y1"] + 12)
prof["y_between_sheets_under_gate"] = profile("y", gx, 0, z_sac2, n0["y0"] - 12, n0["y1"] + 12)
prof["x_at_sheet_level_N"] = profile("x", 0, yn, z_sh2, X0, X1)
prof["x_between_sheets_N"] = profile("x", 0, yn, z_sac2, X0, X1)
prof["x_at_sheet_level_P"] = profile("x", 0, yp, z_sh2, X0, X1)
prof["z_through_N_source_drain"] = profile("z", cx0, yn, 0, z_stack0 - 8, z_gate_top + 1)
prof["y_sheets_at_sheet_level"] = profile("y", gx, 0, z_sh2, 0.0, H3)
prof["y_contact_level_tap_column"] = profile("y", cx0, 0, zc, 0.0, H3)
prof["y_contact_level_plain_column"] = profile("y", cx1, 0, zc, 0.0, H3)
prof["y_epi_level_plain_column"] = profile("y", cx1, 0, z_sh2, 0.0, H3)
res_out["profiles"] = prof
res_out["sheet_gaps"] = gaps(runs(prof["y_sheets_at_sheet_level"], ("Si_Xtal",)))
res_out["contact_gaps_tap_column"] = gaps(runs(prof["y_contact_level_tap_column"], ("W", "Ti")))
res_out["contact_gaps_plain_column"] = gaps(runs(prof["y_contact_level_plain_column"], ("W", "Ti")))
res_out["epi_gaps"] = gaps(runs(prof["y_epi_level_plain_column"], ("SiC", "SiGe")))
with open(os.path.join(outdir, "stage5_T4p5_nanosheet.json"), "w") as f:
    json.dump(res_out, f, indent=1)
for k in sorted(prof):
    if not k.startswith("y_contact") and not k.startswith("y_epi") and not k.startswith("y_sheets"):
        print("PROF", k, prof[k])
for k in ("sheet_gaps", "contact_gaps_tap_column", "contact_gaps_plain_column", "epi_gaps"):
    print("RESULT", k, res_out[k])

doc.GetProperties().SetModelBuildComplete(True)
doc.Save()
mpiFinalize()
