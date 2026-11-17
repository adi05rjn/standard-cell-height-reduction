#! python
# Stage 2: fin + STI + dummy gate + gate cut + spacers + N/P epi + ILD + RMG + SAC cap + S/D contacts.
# Usage (inside SEMulator3D shell):  python build_stage2.py 6.5 [stopAfter]
_SEMulator3D_API_VERSION = '13.1'
import sys, os, json
from voxelModeler import *
from CModeler import mpiInitialize, mpiFinalize
from CDocument import ECoordinates
mpiInitialize(len(sys.argv), sys.argv)

here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, here)
import cells

tracks = float(sys.argv[1])
stop_after = sys.argv[2] if len(sys.argv) > 2 else "end"
tag = ("T%.1f" % tracks).replace(".", "p")
outdir = os.path.join(here, "out")
os.makedirs(outdir, exist_ok=True)
zam = os.path.join(outdir, "stage2_%s.zam" % tag)

# ---- geometry (nm) -------------------------------------------------------
g = cells.cell(tracks)
nfins, pfins = cells.fins_stack(tracks, rows=3)
devs = cells.devices_stack(tracks, rows=3)
H = g["H"]; H3 = 3 * H
CPP = cells.CPP
X0, X1 = 0.0, 3 * CPP
LG, SP = 16.0, 6.0                      # gate length, spacer thickness (assumptions)
gate_cx = [CPP / 2 + CPP * k for k in range(3)]
gate_x = [(c - LG / 2, c + LG / 2) for c in gate_cx]
sd_x = [(gate_x[k][1] + SP, gate_x[k + 1][0] - SP) for k in range(2)]
CUT = 18.0                              # gate cut width (EUV, paper 2019)
CEXT = 5.0                              # contact extension over fin (assumption; gives 10.5 nm N-P tip-to-tip at 5.5T)
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
mask("SD", [(a, 0.0, b, H3) for a, b in sd_x])
mask("NBlock", [(X0, d["y0"] - 10, X1, d["y1"] + 10) for d in devs if d["type"] == "N"])
mask("PBlock", [(X0, d["y0"] - 10, X1, d["y1"] + 10) for d in devs if d["type"] == "P"])
# Power tap (2018 paper, figs. 4 and 7): in S/D column 0 the middle row's contacts run out to the
# cell boundary, where the power rail sits. The neighbouring cells' contacts stay short, which
# gives the paper's border tip-to-tip (16.5 nm at 5.5T, 17.25 nm at 4.5T). Column 1 has no tap.
def con_y(d, k):
    ya, yb = d["y0"] - CEXT, d["y1"] + CEXT
    if k == 0 and d["row"] == 1:
        bd = min((H, 2 * H), key=lambda b_: abs(b_ - (d["y0"] + d["y1"]) / 2.0))
        ya, yb = min(ya, bd), max(yb, bd)
    return ya, yb
mask("Contact", [(a - 2, con_y(d, k)[0], b + 2, con_y(d, k)[1]) for d in devs for k, (a, b) in enumerate(sd_x)])
modeler.CacheBuildBoundary()

# ---- materials -----------------------------------------------------------
procDoc = CProcessData.ProcessDocument()
procDoc.openWithPasswords(os.path.join(here, "common", "16nmFINFET.vproc"), "", "", False)
matDoc = CMaterialData.MaterialDocument()
matDoc.openWithPasswords(os.path.join(here, "common", "16nmFINFET.vmpd"), "", "", False)
for mat in matDoc.getActiveMaterialsInProcess(procDoc):
    doc.GetMaterials().Add(mat, DUP_ADD_ALLOW)
for dop in matDoc.getActiveDopantsInProcess(procDoc):
    doc.GetDopants().Add(dop, DUP_ADD_ALLOW)

SI, SIN, STI = "/Silicon/Si_Xtal", "/Nitrides/Si3N4_LPCVD", "/Oxides/SiO2_HDP"
APOLY, ILD, BLOCK = "/Silicon/Si_Amorph", "/Oxides/SiO2_TEOS", "/Oxides/SiO2_LTO"
NEPI, PEPI = "/Silicon/SiC", "/Silicon/SiGe"
HIK, WF, W = "/Oxides/HfOx", "/Metals/TiN", "/Metals/W"
LINER = "/Metals/Ti"
seq = doc.GetModel()
def out(n): return OutputParams(seq, n, -1)

def mat_obj(path):
    for m in doc.GetMaterials():
        if m.GetName().endswith(path.split("/")[-1]):
            return m
    raise RuntimeError("material not found: " + path)

wafer = None
def report(label):
    info = wafer.GetMaterialBounds()
    parts = []
    for m in info.MaterialsFound():
        parts.append("%s=%.0f" % (m.GetName().split("/")[-1], info.VolumeForMaterial(m)))
    print("VOL[%s] %s" % (label, " ".join(sorted(parts))))

def intervals(axis, fixed_a, fixed_b, mats, lo, hi, step=res):
    """Scan along 'x' or 'y' at fixed other coords (z mid slab) -> intervals where any of mats present."""
    objs = [mat_obj(m) for m in mats]
    zc = fixed_b
    res_iv, inside, start = [], False, None
    v = lo
    while v < hi:
        if axis == "x":
            box = [v, fixed_a, zc - 0.25, v + step, fixed_a + step, zc + 0.25]
        else:
            box = [fixed_a, v, zc - 0.25, fixed_a + step, v + step, zc + 0.25]
        info = wafer.GetMaterialBounds(box)
        here_ = any(info.VolumeForMaterial(o) > 0 for o in objs)
        if here_ and not inside:
            inside, start = True, v
        elif not here_ and inside:
            inside = False
            res_iv.append((start, v))
        v += step
    if inside:
        res_iv.append((start, hi))
    return res_iv

def done(label):
    if stop_after == label:
        print("STOP after", label)
        doc.GetProperties().SetModelBuildComplete(True)
        doc.Save()
        mpiFinalize()
        sys.exit(0)

# ---- process: fins + STI (stage 1) ---------------------------------------
p = StepProperties("WaferSetup"); p.wafer = "Wafer1"; p.material = SI; p.thickness = Z_SI
modeler.WaferSetup(p, out("Wafer setup"))
wafer = modeler.GetActiveWafers()["Wafer1"]

p = StepProperties("Deposit"); p.type = "conformal"; p.material = SIN; p.thickness = HM; p.waferSide = "top"
wafer.Deposit(p, out("Hard mask nitride"))
p = StepProperties("Etch"); p.type = "straight"; p.materials = (SIN, SI); p.depth = HM + TRENCH; p.waferSide = "top"
if hasattr(p, "selectivity"): p.selectivity = False
p.advanced.setEnabled(True); p.mask.name = "FinAll"; p.mask.polarity = 1
wafer.Etch(p, out("Fin etch"))
p = StepProperties("Deposit"); p.type = "conformal"; p.material = STI; p.thickness = 90; p.waferSide = "top"
wafer.Deposit(p, out("STI fill"))
p = StepProperties("Polish"); p.type = "stopOnMaterial"; p.waferSide = "top"
p.type.stopOnMaterial.materials = (SIN,); p.type.stopOnMaterial.overpolish = 0
wafer.Polish(p, out("STI CMP"))
p = StepProperties("Etch"); p.type = "basic"; p.materials = (STI,); p.depth = HM + REVEAL; p.waferSide = "top"
if hasattr(p, "selectivity"): p.selectivity = False
wafer.Etch(p, out("STI recess"))
p = StepProperties("RemoveMaterials"); p.materials = (SIN,)
wafer.RemoveMaterials(p, out("Strip hard mask"))
report("fins"); done("fins")

full = wafer.GetMaterialBounds()
z_fin_top = full.BoundsForMaterial(mat_obj(SI)).GetBounds()[5]
z_sti_top = full.BoundsForMaterial(mat_obj(STI)).GetBounds()[5]
z_fin_mid = (z_fin_top + z_sti_top) / 2.0

# ---- dummy gate ------------------------------------------------------------
p = StepProperties("Deposit"); p.type = "planarizing"; p.material = APOLY; p.thickness = POLY; p.waferSide = "top"
wafer.Deposit(p, out("Dummy gate poly"))
p = StepProperties("Deposit"); p.type = "straight"; p.material = SIN; p.thickness = CAP; p.waferSide = "top"
wafer.Deposit(p, out("Gate cap"))
report("poly+cap")
z_gate_top = wafer.GetMaterialBounds().BoundsForMaterial(mat_obj(APOLY)).GetBounds()[5]
p = StepProperties("Etch"); p.type = "straight"; p.materials = (SIN, APOLY); p.depth = CAP + (z_gate_top - z_sti_top) + 5; p.waferSide = "top"
if hasattr(p, "selectivity"): p.selectivity = False
p.advanced.setEnabled(True); p.mask.name = "Gate"; p.mask.polarity = 1
wafer.Etch(p, out("Gate etch"))
p = StepProperties("Etch"); p.type = "straight"; p.materials = (SIN, APOLY); p.depth = CAP + (z_gate_top - z_sti_top) + 5; p.waferSide = "top"
if hasattr(p, "selectivity"): p.selectivity = False
p.advanced.setEnabled(True); p.mask.name = "GateCut"; p.mask.polarity = 0
wafer.Etch(p, out("Gate cut"))
report("gate")
_zhi = z_gate_top - 10
print("CHK gate poly x @mid-N:", intervals("x", (devs[0]["y0"]+devs[0]["y1"])/2.0, _zhi, (APOLY,), X0, X1))
print("CHK gate poly y @x=%g:" % gate_cx[1], intervals("y", gate_cx[1], _zhi, (APOLY,), 0.0, H3))
done("gate")

# ---- spacers ------------------------------------------------------------
p = StepProperties("Deposit"); p.type = "conformal"; p.material = SIN; p.thickness = SP; p.waferSide = "top"
wafer.Deposit(p, out("Spacer deposit"))
p = StepProperties("Etch"); p.type = "straight"; p.materials = (SIN,); p.depth = 200; p.waferSide = "top"
if hasattr(p, "selectivity"): p.selectivity = False
p.advanced.setEnabled(True); p.mask.name = "SD"; p.mask.polarity = 0
wafer.Etch(p, out("Spacer open in S/D"))
report("spacer")
_ym = (devs[0]["y0"]+devs[0]["y1"])/2.0
print("CHK spacer+poly x @z_hi:", intervals("x", _ym, z_gate_top - 10, (SIN, APOLY), X0, X1))
print("CHK nitride on fin level x @z_fin_mid:", intervals("x", devs[0]["y0"]+2.0, z_fin_mid, (SIN,), X0, X1))
done("spacer")

# ---- N / P epitaxy ------------------------------------------------------
def epi(block_mask, epi_mat, label):
    p = StepProperties("Deposit"); p.type = "conformal"; p.material = BLOCK; p.thickness = 8; p.waferSide = "top"
    p.advanced.setEnabled(True); p.mask.name = block_mask; p.mask.polarity = 0
    wafer.Deposit(p, out(label + " block"))
    p = StepProperties("Growth"); p.type = "selectiveEpitaxy"; p.material = epi_mat; p.thickness = EPI; p.waferSide = "top"
    p.type.selectiveEpitaxy.targetMaterials = (SI,)
    wafer.Growth(p, out(label + " epi"))
    p = StepProperties("RemoveMaterials"); p.materials = (BLOCK,)
    wafer.RemoveMaterials(p, out(label + " unblock"))

epi("PBlock", NEPI, "N")
epi("NBlock", PEPI, "P")
report("epi")
for _k, _d in enumerate(devs[:4]):
    _em = NEPI if _d["type"] == "N" else PEPI
    print("CHK epi", _d["type"], "y-span", intervals("y", (sd_x[0][0]+sd_x[0][1])/2.0, z_fin_mid + 8, (_em,), _d["y0"]-15, _d["y1"]+15),
          "fins", (round(_d["y0"],1), round(_d["y1"],1)))
done("epi")

# ---- ILD + replacement metal gate ----------------------------------------------
p = StepProperties("Deposit"); p.type = "planarizing"; p.material = ILD; p.thickness = 40; p.waferSide = "top"
wafer.Deposit(p, out("ILD0"))
p = StepProperties("Polish"); p.type = "stopOnMaterial"; p.waferSide = "top"
p.type.stopOnMaterial.materials = (APOLY,); p.type.stopOnMaterial.overpolish = 0
wafer.Polish(p, out("ILD CMP to poly"))
p = StepProperties("Etch"); p.type = "straight"; p.materials = (APOLY,); p.depth = 300; p.waferSide = "top"
if hasattr(p, "selectivity"): p.selectivity = False
wafer.Etch(p, out("Dummy gate removal"))
report("poly-removed"); done("polyremoved")

p = StepProperties("Deposit"); p.type = "conformal"; p.material = HIK; p.thickness = 2; p.waferSide = "top"
wafer.Deposit(p, out("High-k"))
p = StepProperties("Deposit"); p.type = "conformal"; p.material = WF; p.thickness = 3; p.waferSide = "top"
wafer.Deposit(p, out("Work function"))
p = StepProperties("Deposit"); p.type = "conformal"; p.material = W; p.thickness = 8; p.waferSide = "top"
wafer.Deposit(p, out("Gate W"))
p = StepProperties("Polish"); p.type = "stopOnMaterial"; p.waferSide = "top"
p.type.stopOnMaterial.materials = (ILD,); p.type.stopOnMaterial.overpolish = 0
wafer.Polish(p, out("RMG CMP"))
report("rmg"); done("rmg")

p = StepProperties("Etch"); p.type = "straight"; p.materials = (W, WF, HIK); p.depth = 25; p.waferSide = "top"
if hasattr(p, "selectivity"): p.selectivity = False
wafer.Etch(p, out("Gate recess"))
p = StepProperties("Deposit"); p.type = "conformal"; p.material = SIN; p.thickness = 14; p.waferSide = "top"
wafer.Deposit(p, out("SAC cap"))
p = StepProperties("Polish"); p.type = "stopOnMaterial"; p.waferSide = "top"
p.type.stopOnMaterial.materials = (ILD,); p.type.stopOnMaterial.overpolish = 0
wafer.Polish(p, out("SAC cap CMP"))
report("sac"); done("sac")

# ---- S/D contacts ----------------------------------------------------------------
p = StepProperties("Etch"); p.type = "straight"; p.materials = (ILD,); p.depth = 250; p.waferSide = "top"
if hasattr(p, "selectivity"): p.selectivity = False
p.advanced.setEnabled(True); p.mask.name = "Contact"; p.mask.polarity = 0
wafer.Etch(p, out("Contact etch (self-aligned)"))
p = StepProperties("Deposit"); p.type = "conformal"; p.material = LINER; p.thickness = 2; p.waferSide = "top"
wafer.Deposit(p, out("Contact liner"))
p = StepProperties("Deposit"); p.type = "conformal"; p.material = W; p.thickness = 8; p.waferSide = "top"
wafer.Deposit(p, out("Contact W"))
p = StepProperties("Polish"); p.type = "stopOnMaterial"; p.waferSide = "top"
p.type.stopOnMaterial.materials = (ILD, SIN); p.type.stopOnMaterial.overpolish = 0
wafer.Polish(p, out("Contact CMP"))
report("contact")
modeler.SaveFinalStep()

# ---- measurements ---------------------------------------------------------------
z_m = z_fin_top + 20.0                   # above the epi, inside gate metal and contact metal
objs_all = [(m.GetName().split("/")[-1], m) for m in doc.GetMaterials()]

def dominant(box):
    info = wafer.GetMaterialBounds(box)
    best, bv = "air", 0.0
    for nm, m in objs_all:
        v = info.VolumeForMaterial(m)
        if v > bv:
            best, bv = nm, v
    return best

def profile(axis, fixed_a, z, lo, hi):
    out_, cur, start = [], None, lo
    v = lo
    while v < hi:
        box = [v, fixed_a, z - 0.25, v + res, fixed_a + res, z + 0.25] if axis == "x" else \
              [fixed_a, v, z - 0.25, fixed_a + res, v + res, z + 0.25]
        d = dominant(box)
        if d != cur:
            if cur is not None:
                out_.append((cur, start, v))
            cur, start = d, v
        v += res
    out_.append((cur, start, hi))
    return out_

def spans(prof, names):
    return [(a, b) for (n, a, b) in prof if n in names]

metal = ("W", "TiN", "HfOx", "Ti")
res_out = dict(tracks=tracks, cell_height=H, z_fin_top=z_fin_top, z_probe=z_m,
               gate_x=gate_x, sd_x=sd_x)
cx = (sd_x[0][0] + sd_x[0][1]) / 2.0
d0 = devs[0]
res_out["x_profile_through_device0"] = profile("x", (d0["y0"] + d0["y1"]) / 2.0, z_m, X0, X1)
res_out["x_profile_at_cell_border"] = profile("x", H, z_m, X0, X1)
yp_gate = profile("y", gate_cx[1], z_m, 0.0, H3)
res_out["gate_metal_y_spans"] = spans(yp_gate, ("W", "TiN", "HfOx"))
yp_con = profile("y", cx, z_m, 0.0, H3)
con = spans(yp_con, ("W", "Ti"))
res_out["contact_y_spans"] = con
res_out["contact_tip_to_tip"] = [round(con[i + 1][0] - con[i][1], 1) for i in range(len(con) - 1)]
cx1 = (sd_x[1][0] + sd_x[1][1]) / 2.0
con1 = spans(profile("y", cx1, z_m, 0.0, H3), ("W", "Ti"))
res_out["contact_y_spans_col1"] = con1
res_out["contact_tip_to_tip_col1"] = [round(con1[i + 1][0] - con1[i][1], 1) for i in range(len(con1) - 1)]
yp_epi = profile("y", cx, z_fin_top + 3.0, 0.0, H3)
epi = spans(yp_epi, ("SiC", "SiGe"))
res_out["epi_y_spans"] = epi
res_out["epi_gaps"] = [round(epi[i + 1][0] - epi[i][1], 1) for i in range(len(epi) - 1)]
# gate-metal to contact lateral separation through device 0 (should equal spacer thickness)
xp = res_out["x_profile_through_device0"]
gm = spans(xp, ("W", "TiN", "HfOx"))
cm = spans(xp, ("Ti", "W"))
res_out["gate_metal_x_spans"] = gm
with open(os.path.join(outdir, "stage2_%s.json" % tag), "w") as f:
    json.dump(res_out, f, indent=1)
print("RESULT", json.dumps(res_out))

doc.GetProperties().SetModelBuildComplete(True)
doc.Save()
mpiFinalize()
