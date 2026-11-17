#! python
# Stage 1: fin + STI cross-section for a given standard-cell track height.
# Usage (inside SEMulator3D shell):  python build_stage1.py 6.5
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
tag = ("T%.1f" % tracks).replace(".", "p")
outdir = os.path.join(here, "out")
os.makedirs(outdir, exist_ok=True)
zam = os.path.join(outdir, "stage1_%s.zam" % tag)

# ---- geometry (nm) -------------------------------------------------------
g = cells.cell(tracks)
nfins, pfins = cells.fins_stack(tracks, rows=3)
H3 = 3 * g["H"]
X0, X1 = 0.0, 60.0                      # fin length shown
PAD = 0.0
Z_SI = 150.0                            # silicon thickness
HM = 25.0                               # nitride hard mask
TRENCH = 100.0                          # fin etch depth into Si
REVEAL = 45.0                           # fin height above STI

# ---- document ------------------------------------------------------------
res = 1.0
doc = VoxelDocument(zam, res, NANOMETERS, OVERWRITE)
doc.GetProperties().SetModelArea(ECoordinates)
doc.GetProperties().SetBoundaryCoordinates((X0, 0.0, X1, H3))
modeler = VoxelModeler(doc)
modeler.SetDefaultMaskUnitsPerMicron(1000)

def rect_mask(name, ranges):
    m = Mask(name, 1000)
    for (a, b) in ranges:
        m.AddRectangle(X0, a, X1, b, NANOMETERS)
    return m

finAll = rect_mask("FinAll", sorted(nfins + pfins))
modeler.AddActiveMask(finAll)
modeler.AddActiveMask(rect_mask("FinN", nfins))
modeler.AddActiveMask(rect_mask("FinP", pfins))
# trench mask = everything except fins
full = Mask("Full", 1000)
full.AddRectangle(X0, 0.0, X1, H3, NANOMETERS)
trench = full - Mask(finAll)
trench.SetName("FinTrench")
modeler.AddActiveMask(trench)
modeler.CacheBuildBoundary()

# ---- materials (borrowed from the STI example database) -----------------
procDoc = CProcessData.ProcessDocument()
procDoc.openWithPasswords(os.path.join(here, "common", "STI_basic_prebuilt.vproc"), "", "", False)
matDoc = CMaterialData.MaterialDocument()
matDoc.openWithPasswords(os.path.join(here, "common", "STI.vmpd"), "", "", False)
for mat in matDoc.getActiveMaterialsInProcess(procDoc):
    doc.GetMaterials().Add(mat, DUP_ADD_ALLOW)
for dop in matDoc.getActiveDopantsInProcess(procDoc):
    doc.GetDopants().Add(dop, DUP_ADD_ALLOW)

SI, SIN, OX = "/Silicon/Si_Xtal", "/Nitrides/Si3N4_LPCVD", "/Oxides/SiO2_HDP"
seq = doc.GetModel()
def out(name):
    return OutputParams(seq, name, -1)

# ---- process -------------------------------------------------------------
p = StepProperties("WaferSetup")
p.wafer = "Wafer1"; p.material = SI; p.thickness = Z_SI
modeler.WaferSetup(p, out("Wafer setup"))
wafer = modeler.GetActiveWafers()["Wafer1"]

p = StepProperties("Deposit")
p.type = "conformal"; p.material = SIN; p.thickness = HM; p.waferSide = "top"
wafer.Deposit(p, out("Hard mask nitride"))

p = StepProperties("Etch")
p.type = "straight"; p.materials = (SIN, SI); p.depth = HM + TRENCH; p.waferSide = "top"
if hasattr(p, "selectivity"): p.selectivity = False
p.advanced.setEnabled(True)
p.mask.name = "FinAll"; p.mask.polarity = 1
wafer.Etch(p, out("Fin etch"))

p = StepProperties("Deposit")
p.type = "conformal"; p.material = OX; p.thickness = 90; p.waferSide = "top"
wafer.Deposit(p, out("STI fill"))

p = StepProperties("Polish")
p.type = "stopOnMaterial"; p.waferSide = "top"
p.type.stopOnMaterial.materials = (SIN,); p.type.stopOnMaterial.overpolish = 0
wafer.Polish(p, out("STI CMP"))

p = StepProperties("Etch")
p.type = "basic"; p.materials = (OX,); p.depth = HM + REVEAL; p.waferSide = "top"
if hasattr(p, "selectivity"): p.selectivity = False
wafer.Etch(p, out("STI recess"))

p = StepProperties("RemoveMaterials")
p.materials = (SIN,)
wafer.RemoveMaterials(p, out("Strip hard mask"))

modeler.SaveFinalStep()

# ---- measurement: scan fin positions from the finished 3D model ----------
def find_mat(suffix):
    for m in doc.GetMaterials():
        if m.GetName().endswith(suffix):
            return m
    raise RuntimeError("material not found: " + suffix)

si_mat, ox_mat = find_mat("Si_Xtal"), find_mat("SiO2_HDP")

def has_si(y0, y1, z0, z1):
    info = wafer.GetMaterialBounds([X0, y0, z0, X1, y1, z1])
    return info.VolumeForMaterial(si_mat) > 0

# find top of silicon / STI for a z reference
full_info = wafer.GetMaterialBounds()
si_b = full_info.BoundsForMaterial(si_mat).GetBounds()
ox_b = full_info.BoundsForMaterial(ox_mat).GetBounds()
z_si_top, z_ox_top = si_b[5], ox_b[5]
zmid = z_ox_top + REVEAL / 2.0

fin_edges, in_fin, start = [], False, None
y = 0.0
while y < H3:
    present = has_si(y, y + res, zmid - 0.25, zmid + 0.25)
    if present and not in_fin:
        in_fin, start = True, y
    elif not present and in_fin:
        in_fin = False
        fin_edges.append((start, y))
    y += res
if in_fin:
    fin_edges.append((start, H3))

exp = sorted(nfins + pfins)
ok = len(fin_edges) == len(exp) and all(
    abs(a - ea) <= res and abs(b - eb) <= res for (a, b), (ea, eb) in zip(fin_edges, exp))
meas_w = [round(b - a, 1) for a, b in fin_edges]
meas_gaps = [round(fin_edges[i + 1][0] - fin_edges[i][1], 1) for i in range(len(fin_edges) - 1)]

result = dict(
    pass_within_one_voxel=ok, measured_fin_widths=meas_w, measured_gaps=meas_gaps,
    tracks=tracks, cell_height=g["H"], expected_gaps=dict(pn=g["pn"], out=g["out"]),
    expected_fins=[[round(a, 2), round(b, 2)] for a, b in sorted(nfins + pfins)],
    measured_fins=[[round(a, 2), round(b, 2)] for a, b in fin_edges],
    z_si_top=z_si_top, z_ox_top=z_ox_top, fin_height=z_si_top - z_ox_top,
)
with open(os.path.join(outdir, "stage1_%s.json" % tag), "w") as f:
    json.dump(result, f, indent=1)
print("RESULT", json.dumps(result))

doc.GetProperties().SetModelBuildComplete(True)
doc.Save()
mpiFinalize()
