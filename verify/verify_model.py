#! python
# Whole-model verification of one cell. Rebuilds the Manager project from a fresh copy, pulls the
# full 3D voxel array out of the model and analyses every voxel, instead of probing along lines.
#   1. Connectivity: label every connected piece of metal; classify it; flag shorts.
#   2. Nets: label metal + epi together; count pieces that touch both N and P epi.
#   3. Separation: smallest axis-aligned dielectric gap between every pair of different metal pieces.
# Usage (inside the SEMulator3D shell):  python verify_model.py <variant> [resolution_nm]
#   variant: 6.5 | 5.5 | 4.5 | handshake | bpr      resolution: 1 (default) or 0.5
_SEMulator3D_API_VERSION = '13.1'
import sys, os, json, io, re, shutil, subprocess, time

here = os.path.dirname(os.path.abspath(__file__))
mgr = os.path.abspath(os.path.join(here, "..", "manager"))
variant = sys.argv[1] if len(sys.argv) > 1 else "5.5"
RES = float(sys.argv[2]) if len(sys.argv) > 2 else 1.0
del sys.argv[1:]
name = subprocess.check_output([sys.executable, os.path.join(mgr, "make_project.py"), variant, "--name"]).decode().strip()
folder = os.path.join(mgr, name)
restag = ("%g" % RES).replace(".", "p")
scratch = os.path.join(here, "_build_%s_%s" % (name, restag))
shutil.rmtree(scratch, ignore_errors=True)
subprocess.check_call([sys.executable, os.path.join(mgr, "make_project.py"), variant], env=dict(os.environ, CELL_OUTDIR=scratch))
shutil.copyfile(os.path.join(folder, name + ".vproc"), os.path.join(scratch, name + ".vproc"))
zam_path = os.path.join(scratch, name + ".zam")
if RES != 1.0:
    z = io.open(zam_path, encoding="utf-8").read().replace('ModelResolution="1.0"', 'ModelResolution="%g"' % RES)
    io.open(zam_path, "w", encoding="utf-8").write(z)

src = io.open(os.path.join(folder, "__semu_%s.py" % name), encoding="utf-8").read()
src = src[:src.index("# Save Model Build")]
old_zam = (folder + "/" + name + ".zam").replace("\\", "/")
assert old_zam in src
src = src.replace(old_zam, zam_path.replace("\\", "/"))
assert "modelResolution = 1\n" in src
src = src.replace("modelResolution = 1\n", "modelResolution = %r\n" % RES)
ns = {"__name__": "__main__"}
t0 = time.time()
exec(compile(src, "__semu_%s.py" % name, "exec"), ns)
if not ns.get("success"):
    print("VERIFY build failed"); sys.exit(1)
build_s = time.time() - t0
modeler, doc = ns["modeler"], ns["voxDoc"]
wafer = modeler.GetActiveWafers()["Wafer1"]

import numpy as np
arr = wafer.GetVoxelArray()
print("VERIFY array", arr.shape, arr.dtype, "resolution", RES, "build %.0f s" % build_s)
_si = [m for m in doc.GetMaterials() if m.GetName().split("/")[-1] == "Si_Xtal"][0]
Z0 = float(wafer.GetMaterialBounds().BoundsForMaterial(_si).GetBounds()[2])   # the array's z index 0 sits at this height (substrate bottom)
frac, m1, m2 = arr[..., 0], arr[..., 1], arr[..., 2]
dom = np.where(frac >= 64, m1, m2)            # dominant material id in each voxel
del arr

# material id -> name, from the model file the build just saved
doc.Save()
ztxt = io.open(zam_path, encoding="utf-8").read()
ids = dict((int(i), n.split("/")[-1]) for i, n in re.findall(r'<material MaterialID="(\d+)"[^>]*? name="([^"]+)"', ztxt))
present = dict((ids.get(int(i), "id%d" % i), int(c)) for i, c in zip(*np.unique(dom, return_counts=True)))
print("VERIFY materials present:", json.dumps(present, sort_keys=True))
def idset(*names): return [i for i, n in ids.items() if n in names]
def mask_of(*names): return np.isin(dom, idset(*names))

METALS = ("W", "Ti", "TiN", "Co")
metal = mask_of(*METALS)

try:
    from scipy import ndimage
    def label(mask):
        return ndimage.label(mask)            # 6-connectivity
    engine = "scipy"
except Exception:
    engine = "numpy"
    def label(mask):
        """6-connected labelling by repeated neighbour-maximum propagation."""
        L = np.zeros(mask.shape, np.int32)
        L[mask] = np.arange(1, int(mask.sum()) + 1, dtype=np.int32)
        while True:
            N = L.copy()
            for ax in range(3):
                a = np.swapaxes(N, 0, ax); s = np.swapaxes(L, 0, ax)
                np.maximum(a[1:], s[:-1], out=a[1:]); np.maximum(a[:-1], s[1:], out=a[:-1])
            N[~mask] = 0
            if np.array_equal(N, L): break
            L = N
        u, inv = np.unique(L, return_inverse=True)
        L = inv.reshape(L.shape).astype(np.int32)
        return L, len(u) - 1
print("VERIFY labelling engine:", engine)

# ---- optional: cross-section images straight from the voxel array -------------------------
# set VERIFY_SLICES, e.g. "x=63;x=42;y=21.5", to write figures/slice_<name>_<axis><pos>.png
COLORS = {"Si_Xtal": (150, 22, 22), "SiO2_HDP": (20, 150, 150), "SiO2_TEOS": (60, 130, 170), "Si3N4_LPCVD": (20, 110, 20),
          "Si3N4_MLD": (140, 225, 130), "Si3N4_PECVD": (0, 90, 80), "SiC": (225, 120, 25), "SiGe": (120, 40, 150),
          "HfOx": (235, 205, 50), "TiN": (125, 125, 125), "W": (45, 55, 62), "Ti": (190, 190, 190), "Co": (140, 120, 210),
          "SiO2_SiCOH": (120, 210, 215), "SiO2_HTO": (90, 180, 200), "SiGe_Sacrificial": (200, 150, 220)}
def write_png(path, rgb):
    import zlib, struct
    h, w, _ = rgb.shape
    raw = b"".join(b"\x00" + rgb[r].tobytes() for r in range(h))
    def chunk(t, d): return struct.pack(">I", len(d)) + t + d + struct.pack(">I", zlib.crc32(t + d) & 0xffffffff)
    open(path, "wb").write(b"\x89PNG\r\n\x1a\n" + chunk(b"IHDR", struct.pack(">IIBBBBB", w, h, 8, 2, 0, 0, 0)) +
                           chunk(b"IDAT", zlib.compress(raw, 6)) + chunk(b"IEND", b""))
if os.environ.get("VERIFY_SLICES"):
    lut = np.full((256, 3), 255, np.uint8)
    for i, nm in ids.items():
        if i < 256: lut[i] = COLORS.get(nm, (255, 255, 255))
    figdir = os.path.abspath(os.path.join(here, "..", "figures"))
    for spec in os.environ["VERIFY_SLICES"].split(";"):
        ax, posv = spec.split("="); k = int(round(float(posv) / RES))
        sl = dom[k, :, :] if ax == "x" else dom[:, k, :]          # (along, z)
        img = lut[sl.T[::-1]]                                     # rows: z from the top; columns: y or x
        scale = max(1, int(round(4 / RES / (1 / RES)))) if RES >= 1 else 2
        scale = 4 if RES >= 1 else 2
        img = np.repeat(np.repeat(img, scale, axis=0), scale, axis=1)
        fn = os.path.join(figdir, "slice_%s_%s%s.png" % (name, ax, posv.replace(".", "p")))
        write_png(fn, np.ascontiguousarray(img))
        print("VERIFY slice written:", fn, img.shape)

# ---- 0. voids: empty voxels not connected to the open top surface ----------------------
air = ~np.isin(dom, list(ids.keys()))
La, na = label(air)
open_ids = set(int(v) for v in np.unique(La[:, :, -1])) - {0}
void_list = []
for k in range(1, na + 1):
    if k in open_ids: continue
    vx, vy, vz = np.where(La == k)
    void_list.append(dict(volume_nm3=float(len(vx) * RES ** 3),
                          at=(float(vx.mean() * RES), float(vy.mean() * RES), float(vz.mean() * RES + Z0))))
void_list.sort(key=lambda v: -v["volume_nm3"])
print("VERIFY empty voxels: %d, open to the top: %d region(s), enclosed voids: %d, void volume %.1f nm3%s" % (
    int(air.sum()), len(open_ids), len(void_list), sum(v["volume_nm3"] for v in void_list),
    (", largest %.1f nm3 at x,y,z = %s" % (void_list[0]["volume_nm3"], [round(c, 1) for c in void_list[0]["at"]])) if void_list else ""))
del La

# ---- 1. connectivity -------------------------------------------------------------
L, n = label(metal)
pieces = []
for k in range(1, n + 1):
    pm = (L == k)
    xs, ys, zs = np.where(pm)
    mats = dict((ids[int(i)], int(c)) for i, c in zip(*np.unique(dom[pm], return_counts=True)))
    has = lambda m_: mats.get(m_, 0) > 0
    kind = ("SHORT gate-contact" if has("TiN") and has("Ti") else "gate" if has("TiN") else
            "contact+" + ("cobalt" if has("Co") else "") if has("Ti") and has("Co") else
            "contact" if has("Ti") else "cobalt only" if has("Co") else "tungsten only")
    kind = kind.rstrip("+")
    pieces.append(dict(id=k, kind=kind, voxels=int(pm.sum()),
                       x=(float(xs.min() * RES), float((xs.max() + 1) * RES)),
                       y=(float(ys.min() * RES), float((ys.max() + 1) * RES)),
                       z=(float(zs.min() * RES + Z0), float((zs.max() + 1) * RES + Z0)), materials=mats))
kinds = {}
for p in pieces: kinds[p["kind"]] = kinds.get(p["kind"], 0) + 1
tiny = [p for p in pieces if p["voxels"] * RES ** 3 < 50]
print("VERIFY metal pieces: %d  %s" % (n, json.dumps(kinds, sort_keys=True)))
print("VERIFY pieces under 50 nm3 (slivers):", [(p["kind"], p["voxels"], p["x"], p["y"], p["z"]) for p in tiny])

# ---- 2. nets: metal + epi -----------------------------------------------------------
nmask, pmask = mask_of("SiC"), mask_of("SiGe")
Ln, nn = label(metal | nmask | pmask)
both = n_only = p_only = neither = 0
for k in range(1, nn + 1):
    sel = (Ln == k)
    hn, hp = bool((sel & nmask).any()), bool((sel & pmask).any())
    both += hn and hp; n_only += hn and not hp; p_only += hp and not hn; neither += (not hn and not hp)
Le, ne = label(nmask | pmask)
epi_mixed = sum(1 for k in range(1, ne + 1) if (nmask[Le == k]).any() and (pmask[Le == k]).any())
print("VERIFY nets (metal+epi): %d  touching N only %d, P only %d, both N and P %d, no epi %d" % (nn, n_only, p_only, both, neither))
print("VERIFY epi pieces: %d, with N and P epi touching directly: %d" % (ne, epi_mixed))

# ---- 3. separation between different metal pieces --------------------------------------
kind_of = dict((p["id"], p["kind"]) for p in pieces)
gaps = {}                                   # (kindA, kindB, axis) -> [min gap nm, count at min, example position]
for ax, axn in enumerate("xyz"):
    Lm = np.moveaxis(L, ax, 0)
    npos = Lm.shape[0]
    pos = np.arange(npos, dtype=np.int32).reshape((-1, 1, 1))
    idx = np.where(Lm > 0, pos, -1)
    last = np.maximum.accumulate(idx, axis=0)           # index of latest metal voxel at or before
    prev = np.empty_like(last); prev[0] = -1; prev[1:] = last[:-1]
    prev_lab = np.take_along_axis(Lm, np.clip(prev, 0, None), axis=0)
    sel = (Lm > 0) & (prev >= 0) & (prev_lab != Lm) & (prev_lab > 0)
    a, b, c = np.where(sel)
    g = (a - prev[a, b, c] - 1).astype(np.float64) * RES
    la, lb = Lm[a, b, c], prev_lab[a, b, c]
    for i in np.argsort(g, kind="stable"):
        key = tuple(sorted((kind_of[int(la[i])], kind_of[int(lb[i])]))) + (axn,)
        if key not in gaps:
            p3 = [0, 0, 0]; p3[ax] = float(a[i] * RES); o = [d for d in range(3) if d != ax]
            p3[o[0]] = float(b[i] * RES); p3[o[1]] = float(c[i] * RES)
            p3[2] += Z0
            gaps[key] = [float(g[i]), int((g[(la == la[i]) | True] == g[i]).sum()), p3]
    del idx, last, prev, prev_lab, sel
rows = sorted(gaps.items(), key=lambda kv: kv[1][0])
for (ka, kb, axn), (gmin, cnt, where) in rows:
    print("VERIFY min gap %-18s to %-18s along %s: %5.1f nm   (e.g. at x,y,z = %s)" % (ka, kb, axn, gmin, where))
overall = rows[0] if rows else None
print("VERIFY smallest gap in the whole model: %s" % (("%.1f nm, %s to %s along %s" % (overall[1][0], overall[0][0], overall[0][1], overall[0][2])) if overall else "none"))

out = dict(variant=variant, name=name, resolution_nm=RES, engine=engine, shape=list(dom.shape), materials_present=present,
           metal_pieces=n, piece_kinds=kinds, slivers=len(tiny), pieces=pieces,
           nets=dict(total=int(nn), n_only=int(n_only), p_only=int(p_only), both=int(both), no_epi=int(neither)),
           epi_pieces=int(ne), epi_n_p_touching=int(epi_mixed),
           empty_voxels=int(air.sum()), enclosed_voids=len(void_list), voids=void_list[:20],
           min_gaps=[dict(a=k[0], b=k[1], axis=k[2], gap_nm=v[0], example_xyz=v[2]) for k, v in rows])
with open(os.path.join(here, "verify_%s_res%s.json" % (name, restag)), "w") as f:
    json.dump(out, f, indent=1)
from CModeler import mpiFinalize
mpiFinalize()
shutil.rmtree(scratch, ignore_errors=True)
