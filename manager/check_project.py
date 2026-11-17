#! python
# Rebuilds a Manager project from a fresh copy through its own generated build script, measures the
# result, and compares it with the script-built model of the same cell.
# Usage (inside the SEMulator3D shell):  python check_project.py <variant>     (6.5 | 5.5 | 4.5 | handshake | bpr)
_SEMulator3D_API_VERSION = '13.1'
import sys, os, json, io, shutil, subprocess

here = os.path.dirname(os.path.abspath(__file__))
variant = sys.argv[1] if len(sys.argv) > 1 else "5.5"
del sys.argv[1:]                                   # the generated script parses argv itself
name = subprocess.check_output([sys.executable, os.path.join(here, "make_project.py"), variant, "--name"]).decode().strip()
folder = os.path.join(here, name)
outref = os.path.join(here, "..", "cellheight", "out")

# A built model file will not rebuild in place, so work on a fresh unbuilt copy in a scratch folder.
scratch = os.path.join(folder, "_check")
shutil.rmtree(scratch, ignore_errors=True)
subprocess.check_call([sys.executable, os.path.join(here, "make_project.py"), variant], env=dict(os.environ, CELL_OUTDIR=scratch))
shutil.copyfile(os.path.join(folder, name + ".vproc"), os.path.join(scratch, name + ".vproc"))
src = io.open(os.path.join(folder, "__semu_%s.py" % name), encoding="utf-8").read()
src = src[:src.index("# Save Model Build")]         # build only: no save, no viewer
old_zam = (folder + "/" + name + ".zam").replace("\\", "/")
new_zam = (scratch + "/" + name + ".zam").replace("\\", "/")
assert old_zam in src, "model file path not found in generated script"
ns = {"__name__": "__main__"}
exec(compile(src.replace(old_zam, new_zam), "__semu_%s.py" % name, "exec"), ns)
if not ns.get("success"):
    print("CHECK build failed"); sys.exit(1)

modeler, doc = ns["modeler"], ns["voxDoc"]
wafer = modeler.GetActiveWafers()["Wafer1"]
res = 1.0
objs = [(m.GetName().split("/")[-1], m) for m in doc.GetMaterials()]
full = wafer.GetMaterialBounds()
vol = dict((u, round(full.VolumeForMaterial(m))) for u, m in objs if full.VolumeForMaterial(m) > 0)

def dominant(box):
    info = wafer.GetMaterialBounds(box)
    best, bv = "air", 0.0
    for nm, m in objs:
        v = info.VolumeForMaterial(m)
        if v > bv:
            best, bv = nm, v
    return best

def profile(axis, x, y, z, lo, hi):
    out_, cur, start, v = [], None, lo, lo
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

def runs(prof, names):
    r = []
    for n, a, b in prof:
        if n in names:
            if r and abs(r[-1][1] - a) < 1e-6: r[-1] = (r[-1][0], b)
            else: r.append((a, b))
    return r
def gaps(r):   return [round(r[i + 1][0] - r[i][1], 1) for i in range(len(r) - 1)]
def widths(r): return [round(b - a, 1) for a, b in r]
def seg_runs(spans): return runs([("m", a, b) for a, b in spans], ("m",))
METAL = ("W", "Ti")

rows = []          # (label, manager value, script value)
if variant == "nanosheet":
    ref = json.load(open(os.path.join(outref, "stage5_T4p5_nanosheet.json"))); p = ref["profiles"]
    H, sd_x = ref["cell_height"], ref["sd_x"]
    H3, X1 = 3 * H, 126.0
    z0, z1 = ref["z_stack"]
    ts, tc = ref["t_sheet"], ref["t_sacrificial"]
    z_sac2, z_sh2, zc = z0 + (tc + ts) + tc / 2.0, z0 + (tc + ts) + tc + ts / 2.0, z1 + 20.0
    cx0, cx1 = (sd_x[0][0] + sd_x[0][1]) / 2.0 - 0.5, (sd_x[1][0] + sd_x[1][1]) / 2.0 - 0.5
    gx = int((ref["gate_x"][1][0] + ref["gate_x"][1][1]) / 2.0)
    n0 = [d for d in ref["devices"] if d["type"] == "N" and d["row"] == 0][0]
    yn = int((n0["sheet"][0] + n0["sheet"][1]) / 2.0)
    rows.append(("contact gaps, power-tap column", gaps(runs(profile("y", cx0, 0, zc, 0.0, H3), METAL)), ref["contact_gaps_tap_column"]))
    rows.append(("contact gaps, plain column", gaps(runs(profile("y", cx1, 0, zc, 0.0, H3), METAL)), ref["contact_gaps_plain_column"]))
    rows.append(("sheet gaps under the gate", gaps(runs(profile("y", gx, 0, z_sh2, 0.0, H3), ("Si_Xtal",))), ref["sheet_gaps"]))
    rows.append(("stack under the gate through the sheets", [list(s) for s in profile("z", gx, yn, 0, z0 - 4, z1 + 14)], [list(s) for s in p["z_under_gate_through_N_sheets"]]))
    rows.append(("profile along the sheet", [list(s) for s in profile("x", 0, yn, z_sh2, 0.0, X1)], [list(s) for s in p["x_at_sheet_level_N"]]))
    rows.append(("profile between the sheets", [list(s) for s in profile("x", 0, yn, z_sac2, 0.0, X1)], [list(s) for s in p["x_between_sheets_N"]]))
    rows.append(("stack through the source/drain", [list(s) for s in profile("z", cx0, yn, 0, z0 - 8, ref["z_gate_top"] + 1)], [list(s) for s in p["z_through_N_source_drain"]]))
elif variant not in ("handshake", "bpr"):
    ref = json.load(open(os.path.join(outref, "stage2_%s.json" % ("T%.1f" % float(variant)).replace(".", "p"))))
    H3, X1, sd_x, zf = 3 * ref["cell_height"], 126.0, ref["sd_x"], ref["z_fin_top"]
    cx0, cx1, z_m = (sd_x[0][0] + sd_x[0][1]) / 2.0, (sd_x[1][0] + sd_x[1][1]) / 2.0, zf + 20.0
    yd = (ref["contact_y_spans"][0][0] + ref["contact_y_spans"][2][1]) / 2.0
    rows.append(("contact gaps, power-tap column", gaps(runs(profile("y", cx0, 0, z_m, 0.0, H3), METAL)), gaps(seg_runs(ref["contact_y_spans"]))))
    rows.append(("contact gaps, plain column", gaps(runs(profile("y", cx1, 0, z_m, 0.0, H3), METAL)), gaps(seg_runs(ref["contact_y_spans_col1"]))))
    rows.append(("epi gaps", gaps(runs(profile("y", cx0, 0, zf + 3.0, 0.0, H3), ("SiC", "SiGe"))), ref["epi_gaps"]))
    rows.append(("profile along gate direction", [list(s) for s in profile("x", 0, yd, z_m, 0.0, X1)], [list(s) for s in ref["x_profile_through_device0"]]))
elif variant == "handshake":
    ref = json.load(open(os.path.join(outref, "stage3_T4p5.json"))); p = ref["profiles"]
    H, sd_x, zw = ref["cell_height"], ref["sd_x"], ref["z_wing_probe"]
    cx0, cx1 = (sd_x[0][0] + sd_x[0][1]) / 2.0 - 0.5, (sd_x[1][0] + sd_x[1][1]) / 2.0 - 0.5
    yt, zv, zb = int(H / 2), (ref["z_metal_top"] + ref["z_gate_top"]) / 2.0, ref["z_gate_top"] + 10.0
    n0 = [d for d in ref["devices"] if d["type"] == "N" and d["row"] == 0][0]
    yn = int((n0["fin"][0] + n0["fin"][1]) / 2.0)
    m0, m1 = profile("y", cx0, 0, zw, 0.0, H), profile("y", cx1, 0, zw, 0.0, H)
    rows.append(("low-k widths, N-wing column", widths(runs(m0, ("SiO2_SiCOH",))), widths(runs(p["y_col0_wing_level_row0"], ("SiO2_SiCOH",)))))
    rows.append(("low-k widths, P-wing column", widths(runs(m1, ("SiO2_SiCOH",))), widths(runs(p["y_col1_wing_level_row0"], ("SiO2_SiCOH",)))))
    rows.append(("contact metal spans, N-wing column", runs(m0, METAL), runs(p["y_col0_wing_level_row0"], METAL)))
    rows.append(("contact metal spans, P-wing column", runs(m1, METAL), runs(p["y_col1_wing_level_row0"], METAL)))
    rows.append(("profile along gate direction", [list(s) for s in profile("x", 0, yn, zw, 0.0, 126.0)], [list(s) for s in p["x_through_N_fin_wing_level"]]))
    rows.append(("bridge metal along gate direction", runs(profile("x", 0, yt, zb, 0.0, 126.0), ("Co",)), runs(p["b_x_bridge_level"], ("Co",))))
    rows.append(("via on N wing", runs(profile("y", cx0, 0, zv, 0.0, H), ("Co",)), runs(p["b_y_col0_via_level_row0"], ("Co",))))
    rows.append(("via on P wing", runs(profile("y", cx1, 0, zv, 0.0, H), ("Co",)), runs(p["b_y_col1_via_level_row0"], ("Co",))))
    rows.append(("stack at N-wing via", [list(s) for s in profile("z", cx0, yt - 2, 0, ref["z_shelf"] - 5, ref["z_gate_top"] + 21)], [list(s) for s in p["b_z_col0_via"]]))
else:
    ref = json.load(open(os.path.join(outref, "stage4_T5p0_BPR.json"))); p = ref["profiles"]
    H, sd_x, zf = ref["cell_height"], ref["sd_x"], ref["z_fin_top"]
    H3, X1 = 3 * H, 3 * ref["cpp"]
    cxp, cxs = int((sd_x[0][0] + sd_x[0][1]) / 2.0), int((sd_x[1][0] + sd_x[1][1]) / 2.0)
    gx = int((ref["gate_x"][1][0] + ref["gate_x"][1][1]) / 2.0)
    zr, zc, yb = (ref["z_rail"][0] + ref["z_rail"][1]) / 2.0, zf + 20.0, int(H)
    zlo, zhi = ref["z_rail"][0] - 8, ref["z_gate_top"] + 1
    rows.append(("contact gaps, signal column", gaps(runs(profile("y", cxs, 0, zc, 0.0, H3), METAL)), gaps(runs(p["y_contact_level_signal_col"], METAL))))
    rows.append(("contact gaps, power column", gaps(runs(profile("y", cxp, 0, zc, 0.0, H3), METAL)), gaps(runs(p["y_contact_level_power_col"], METAL))))
    rows.append(("profile at rail level along the cell", [list(s) for s in profile("y", gx, 0, zr, 0.0, H3)], [list(s) for s in p["y_rail_level_gatecol"]]))
    rows.append(("stack at power via", [list(s) for s in profile("z", cxp, yb, 0, zlo, zhi)], [list(s) for s in p["z_power_col_at_border"]]))
    rows.append(("stack over rail, signal column", [list(s) for s in profile("z", cxs, yb, 0, zlo, zhi)], [list(s) for s in p["z_signal_col_at_border"]]))
    rows.append(("stack under gate over rail", [list(s) for s in profile("z", gx, yb - 12, 0, zlo, zhi)], [list(s) for s in p["z_gate_over_rail"]]))

def norm(v): return json.loads(json.dumps(v))
def maxdiff(a, b):
    """Largest numeric difference between two same-shaped nested lists, or None if shapes/materials differ."""
    if isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b): return None
        ds = [maxdiff(x, y) for x, y in zip(a, b)]
        return None if any(d is None for d in ds) else (max(ds) if ds else 0.0)
    if isinstance(a, str) or isinstance(b, str):
        return 0.0 if a == b else None
    return abs(a - b)

worst, struct_ok = 0.0, True
for label, mine, theirs in rows:
    d = maxdiff(norm(mine), norm(theirs))
    if d is None: struct_ok = False
    else: worst = max(worst, d)
    tag = "IDENTICAL" if d == 0 else ("DIFFERENT STRUCTURE" if d is None else "within %g nm" % d)
    print("CHECK %-38s %s" % (label, tag))
    if d != 0 and len(json.dumps(norm(mine))) < 400:
        print("      manager:", mine)
        print("      script :", theirs)
print("CHECK volumes (manager):", json.dumps(vol, sort_keys=True))
print("CHECK overall %s: %s" % (name, ("same structure, largest difference %g nm" % worst) if struct_ok else "STRUCTURE DIFFERS"))
with open(os.path.join(folder, name + "_check.json"), "w") as f:
    json.dump(dict(variant=variant, same_structure=struct_ok, largest_difference_nm=worst, volumes=vol,
                   rows=[dict(label=l, manager=norm(a), script=norm(b)) for l, a, b in rows]), f, indent=1)
from CModeler import mpiFinalize
mpiFinalize()
shutil.rmtree(scratch, ignore_errors=True)
