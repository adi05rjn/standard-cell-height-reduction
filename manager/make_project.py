# Generates the input files of a SEMulator3D Manager project:
#   <name>.gds            layout with one layer per mask
#   <name>-layermap.txt   mask name -> GDS layer number
#   <name>.json           process flow, converted to <name>.vproc by the JSON Translator (see build_project.bat)
#   <name>.zam            model file the Manager opens
#   <name>.vmpd           material file
# Usage (plain Python, no license needed):   python make_project.py <variant>
#   6.5 | 5.5 | 4.5   the cell at that track height        (same flow as cellheight/build_stage2.py)
#   handshake         4.5T mid-track handshake with bridge  (same flow as cellheight/build_stage3.py)
#   bpr               5T cell with buried power rail        (same flow as cellheight/build_stage4.py)
# Add --name to print the project name only.
# The geometry comes from cellheight/cells.py, the same source the build scripts use.
import os, sys, json, struct, datetime, shutil

here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(here, "..", "cellheight"))
import cells

variant = sys.argv[1] if len(sys.argv) > 1 else "5.5"
if variant == "handshake":
    tracks, name, title = 4.5, "handshake_4p5T", "4.5-track mid-track handshake, Sherazi et al. SPIE 10588 (2018) sec. 4.1"
elif variant == "bpr":
    tracks, name, title = 5.0, "bpr_5T", "5-track cell with buried power rail, Sherazi et al. SPIE 10962 (2019)"
elif variant == "nanosheet":
    tracks, name, title = 4.5, "nanosheet_4p5T", "4.5-track stacked-nanosheet cell, Sherazi et al. SPIE 10588 (2018) fig. 8"
else:
    tracks = float(variant)
    name = "cell_" + ("%.1fT" % tracks).replace(".", "p")
    title = "%.1f-track cell, Sherazi et al. SPIE 10588 (2018)" % tracks
if "--name" in sys.argv:
    print(name); sys.exit(0)
outdir = os.environ.get("CELL_OUTDIR") or os.path.join(here, name)   # CELL_OUTDIR: used by check_project.py
os.makedirs(outdir, exist_ok=True)

# ---- common geometry (nm) --------------------------------------------------------
ckey = "4.5NS" if variant == "nanosheet" else tracks      # key into cells.CELLS
g = cells.cell(ckey)
nfins, pfins = cells.fins_stack(ckey, rows=3)
devs = cells.devices_stack(ckey, rows=3)
H = g["H"]; H3 = 3 * H
CPP = 45.0 if variant == "bpr" else cells.CPP
X0, X1 = 0.0, 3 * CPP
LG, SP, CUT, CEXT = 16.0, 6.0, 18.0, 5.0
HM, TRENCH, REVEAL = 25.0, 100.0, 45.0
POLY, CAP, EPI = 60.0, 20.0, 8.0
gate_cx = [CPP / 2 + CPP * k for k in range(3)]
gate_x = [(c - LG / 2, c + LG / 2) for c in gate_cx]
sd_x = [(gate_x[k][1] + SP, gate_x[k + 1][0] - SP) for k in range(2)]
# Heights the scripts read back from the model; fixed here because a process file cannot measure.
Z_STI, Z_FIN, Z_GATE = 30.0, 75.0, 135.0          # STI top, fin top, dummy-gate top
GATE_ETCH = CAP + (Z_GATE - Z_STI) + 5

def clipy(a, b): return max(0.0, a), min(H3, b)

# ---- materials -------------------------------------------------------------------
SI, SIN, STI = "/Silicon/Si_Xtal", "/Nitrides/Si3N4_LPCVD", "/Oxides/SiO2_HDP"
APOLY, ILD, BLOCK = "/Silicon/Si_Amorph", "/Oxides/SiO2_TEOS", "/Oxides/SiO2_LTO"
NEPI, PEPI = "/Silicon/SiC", "/Silicon/SiGe"
HIK, WF, W, LINER = "/Oxides/HfOx", "/Metals/TiN", "/Metals/W", "/Metals/Ti"
LOWK, CCAP, CO, RLINER = "/Oxides/SiO2_SiCOH", "/Nitrides/Si3N4_PECVD", "/Metals/Co", "/Oxides/SiO2_HTO"

# ---- step helpers (JSON Translator syntax) -----------------------------------------
WAFER = "Wafer1"
def _mask(s, mask, polarity):
    if mask:
        s.update({"advanced": True, "mask.name": mask, "mask.polarity": bool(polarity)})
    return s
def dep(nm, kind, mat, thk, mask=None):
    return _mask({"stepType": "Deposit", "stepName": nm, "wafer": WAFER, "type": kind, "material": mat,
                  "thickness": str(thk)}, mask, 0)
def etch(nm, mats, depth, mask=None, polarity=0, kind="straight"):
    return _mask({"stepType": "Etch", "stepName": nm, "wafer": WAFER, "type": kind, "materials": ",".join(mats),
                  "selectivity": False, "depth": str(depth)}, mask, polarity)
def cmp_(nm, stops):  return {"stepType": "Polish", "stepName": nm, "wafer": WAFER, "type": "stopOnMaterial",
                              "type.stopOnMaterial.materials": ",".join(stops), "type.stopOnMaterial.overpolish": "0"}
def strip(nm, mats):  return {"stepType": "RemoveMaterials", "stepName": nm, "wafer": WAFER, "materials": ",".join(mats)}
def grow(nm, mat, thk=None):
    return {"stepType": "Growth", "stepName": nm, "wafer": WAFER, "type": "selectiveEpitaxy", "material": mat,
            "thickness": str(EPI if thk is None else thk), "type.selectiveEpitaxy.targetMaterials": SI}
def seq(nm, steps):   return {"stepType": "Sequence", "stepName": nm, "children": steps}

wafer_setup = {"stepType": "WaferSetup", "stepName": "Wafer setup", "wafer": WAFER, "material": SI, "thickness": "150"}
fin_etch_fill = [
    dep("Hard mask nitride", "conformal", SIN, HM),
    etch("Fin etch", (SIN, SI), HM + TRENCH, "FinAll", 1),
    dep("STI fill", "conformal", STI, 90),
    cmp_("STI CMP", (SIN,)),
]
dummy_gate = seq("Dummy gate", [
    dep("Dummy gate poly", "planarizing", APOLY, POLY),
    dep("Gate cap", "straight", SIN, CAP),
    etch("Gate etch", (SIN, APOLY), GATE_ETCH, "Gate", 1),
    etch("Gate cut", (SIN, APOLY), GATE_ETCH, "GateCut", 0),
])
spacer_epi = [
    seq("Spacer", [
        dep("Spacer deposit", "conformal", SIN, SP),
        etch("Spacer open in S/D", (SIN,), 200, "SD", 0),
    ]),
    seq("Source/drain epi", [
        dep("N block", "conformal", BLOCK, 8, "PBlock"),
        grow("N epi", NEPI),
        strip("N unblock", (BLOCK,)),
        dep("P block", "conformal", BLOCK, 8, "NBlock"),
        grow("P epi", PEPI),
        strip("P unblock", (BLOCK,)),
    ]),
]
def rmg(first):
    return seq("Replacement metal gate", first + [
        etch("Dummy gate removal", (APOLY,), 300),
        dep("High-k", "conformal", HIK, 2),
        dep("Work function", "conformal", WF, 3),
        dep("Gate W", "conformal", W, 8),
        cmp_("RMG CMP", (ILD,)),
        etch("Gate recess", (W, WF, HIK), 25),
        dep("SAC cap", "conformal", SIN, 14),
        cmp_("SAC cap CMP", (ILD,)),
    ])
rmg_gate_first = rmg([dep("ILD0", "planarizing", ILD, 40), cmp_("ILD CMP to poly", (APOLY,))])
base_masks = [
    ("BB",      [(X0, 0.0, X1, H3)]),
    ("FinAll",  [(X0, a, X1, b) for a, b in sorted(nfins + pfins)]),
    ("Gate",    [(a, 0.0, b, H3) for a, b in gate_x]),
    ("GateCut", [(X0, k * H - CUT / 2, X1, k * H + CUT / 2) for k in (1, 2)]),
]
sd_masks = [
    ("SD",      [(a, 0.0, b, H3) for a, b in sd_x]),
    ("NBlock",  [(X0, d["y0"] - 10, X1, d["y1"] + 10) for d in devs if d["type"] == "N"]),
    ("PBlock",  [(X0, d["y0"] - 10, X1, d["y1"] + 10) for d in devs if d["type"] == "P"]),
]

# ---- variant: plain cell (Stage 2) -------------------------------------------------
if variant not in ("handshake", "bpr", "nanosheet"):
    def con_y(d, k):
        ya, yb = d["y0"] - CEXT, d["y1"] + CEXT
        if k == 0 and d["row"] == 1:          # power tap: middle row, column 0, runs to the cell boundary
            bd = min((H, 2 * H), key=lambda b_: abs(b_ - (d["y0"] + d["y1"]) / 2.0))
            ya, yb = min(ya, bd), max(yb, bd)
        return ya, yb
    masks = base_masks + sd_masks + [
        ("Contact", [(a - 2, con_y(d, k)[0], b + 2, con_y(d, k)[1]) for d in devs for k, (a, b) in enumerate(sd_x)]),
    ]
    flow = [
        wafer_setup,
        seq("Fins and STI", fin_etch_fill + [
            etch("STI recess", (STI,), HM + REVEAL, kind="basic"),
            strip("Strip hard mask", (SIN,)),
        ]),
        dummy_gate,
    ] + spacer_epi + [
        rmg_gate_first,
        seq("Contacts", [
            etch("Contact etch (self-aligned)", (ILD,), 250, "Contact", 0),
            dep("Contact liner", "conformal", LINER, 2),
            dep("Contact W", "conformal", W, 8),
            cmp_("Contact CMP", (ILD, SIN)),
        ]),
    ]

# ---- variant: 4.5T mid-track handshake (Stage 3) -----------------------------------
elif variant == "handshake":
    LK = 5.25                                # low-k replacement spacer; two give the paper's 10.5 nm
    DEX = CEXT + LK
    M0_T = 20.0
    Z_TOP = Z_GATE + CAP + SP                # surface after the ILD polish stops on the spacer nitride
    Z_SHELF, Z_METAL = Z_FIN + 15.0, Z_GATE - 25.0
    SH, RC = Z_TOP - Z_SHELF, Z_TOP - Z_METAL
    gx = [(gate_x[k][1], gate_x[k + 1][0]) for k in range(2)]
    strips = [[(gate_x[k][1], sd_x[k][0]), (sd_x[k][1], gate_x[k + 1][0])] for k in range(2)]
    WINGCOL = {"N": 0, "P": 1}
    for d in devs:
        q = [e for e in devs if e["row"] == d["row"] and e["type"] != d["type"]][0]
        d["deep"] = (d["y0"] - DEX, d["y1"] + DEX)
        d["wing"] = (d["y1"] + DEX, q["y0"] - DEX) if q["y0"] > d["y0"] else (q["y1"] + DEX, d["y0"] - DEX)
    masks = list(base_masks)
    for T in ("N", "P"):
        dd = [d for d in devs if d["type"] == T]
        wc = WINGCOL[T]
        masks += [
            ("Deep" + T,   [(gx[k][0], d["deep"][0], gx[k][1], d["deep"][1]) for d in dd for k in range(2)]),
            ("Wing" + T,   [(gx[wc][0], d["wing"][0], gx[wc][1], d["wing"][1]) for d in dd]),
            ("Open" + T,   [(sd_x[k][0], d["deep"][0], sd_x[k][1], d["deep"][1]) for d in dd for k in range(2)]),
            ("SpDeep" + T, [(a_, d["deep"][0], b_, d["deep"][1]) for d in dd for k in range(2) for (a_, b_) in strips[k]]),
            ("SpWing" + T, [(a_, d["wing"][0], b_, d["wing"][1]) for d in dd for (a_, b_) in strips[wc]]),
        ]
    TRK = [(r * H + H / 2 - cells.MP / 4, r * H + H / 2 + cells.MP / 4) for r in range(3)]
    masks += [
        ("Bridge",    [(sd_x[0][0], a, sd_x[1][1], b) for a, b in TRK]),
        ("BridgeVia", [(sd_x[k][0], a, sd_x[k][1], b) for a, b in TRK for k in range(2)]),
    ]
    def contact_module(T, epi_mat):
        return seq(T + " contact module", [
            etch(T + " contact deep etch", (ILD,), 300, "Deep" + T),
            etch(T + " contact shallow etch (staircase)", (ILD,), SH, "Wing" + T),
            etch(T + " fin open (nitride)", (SIN,), 300, "Open" + T),
            grow(T + " epi", epi_mat),
            etch(T + " nitride spacer strip (deep)", (SIN,), 300, "SpDeep" + T),
            etch(T + " nitride spacer strip (wing)", (SIN,), SH, "SpWing" + T),
            dep(T + " low-k spacer deposit", "conformal", LOWK, LK),
            etch(T + " low-k spacer etch-back", (LOWK,), LK + 2),
            dep(T + " contact liner", "conformal", LINER, 2),
            dep(T + " contact W", "conformal", W, 10),
            cmp_(T + " contact CMP", (ILD, SIN)),
            etch(T + " contact recess", (W, LINER), RC),
            dep(T + " contact cap", "conformal", CCAP, 10),
            cmp_(T + " contact cap CMP", (ILD, SIN)),
        ])
    flow = [
        wafer_setup,
        seq("Fins and STI", fin_etch_fill + [
            etch("STI recess", (STI,), HM + REVEAL, kind="basic"),
            strip("Strip hard mask", (SIN,)),
        ]),
        dummy_gate,
        seq("Spacer kept over fins and ILD", [
            dep("Nitride spacer (kept over fins)", "conformal", SIN, SP),
            dep("ILD0", "planarizing", ILD, 40),
            cmp_("ILD CMP to nitride", (SIN,)),
        ]),
        contact_module("N", NEPI),
        contact_module("P", PEPI),
        rmg([cmp_("Poly-open CMP", (APOLY,))]),
        seq("Mid-track handshake bridge", [
            dep("Bridge dielectric", "straight", ILD, M0_T),
            etch("Bridge trench etch", (ILD,), M0_T, "Bridge"),
            etch("Via etch through contact caps", (CCAP,), 60, "BridgeVia"),
            dep("Bridge metal (Co)", "conformal", CO, 8),
            cmp_("Bridge CMP", (ILD,)),
        ]),
    ]

# ---- variant: 4.5T stacked nanosheets (Stage 5) ------------------------------------
elif variant == "nanosheet":
    SAC, ISP = "/SiGe_Sacrificial", "/Nitrides/Si3N4_MLD"
    T_SHEET, T_SAC, NSHEET = 6.0, 10.0, 3
    STACK = NSHEET * (T_SHEET + T_SAC)
    INDENT, ISP_T, EPI_NS, BELOW = 7.0, 5.0, 10.0, 52.0
    TAP = g["out"] - CEXT - 10.5                 # power-tap reach that leaves the paper's 10.5 nm at the border
    NS_Z_STI, NS_Z_GATE = 75.0, 75.0 + STACK + POLY   # STI top = stack bottom; dummy-gate top
    NS_GATE_ETCH = CAP + (NS_Z_GATE - NS_Z_STI) + 5
    sd_all = [(X0, gate_x[0][0] - SP)] + sd_x + [(gate_x[2][1] + SP, X1)]
    def con_y(d, k):
        ya, yb = d["y0"] - CEXT, d["y1"] + CEXT
        if k == 0 and d["row"] == 1:
            c = (d["y0"] + d["y1"]) / 2.0
            bd = min((H, 2 * H), key=lambda b_: abs(b_ - c))
            if c > bd: ya = d["y0"] - TAP
            else:      yb = d["y1"] + TAP
        return ya, yb
    masks = [
        ("BB",      [(X0, 0.0, X1, H3)]),
        ("Sheet",   [(X0, a, X1, b) for a, b in sorted(nfins + pfins)]),
        base_masks[2], base_masks[3],
        ("SD",      [(a, 0.0, b, H3) for a, b in sd_all]),
        ("NBlock",  [(X0, d["y0"] - 12, X1, d["y1"] + 12) for d in devs if d["type"] == "N"]),
        ("PBlock",  [(X0, d["y0"] - 12, X1, d["y1"] + 12) for d in devs if d["type"] == "P"]),
        ("Contact", [(a - 2, con_y(d, k)[0], b + 2, con_y(d, k)[1]) for d in devs for k, (a, b) in enumerate(sd_x)]),
    ]
    stack_steps = []
    for i in range(NSHEET):
        stack_steps += [dep("Sacrificial SiGe %d" % (i + 1), "straight", SAC, T_SAC),
                        dep("Silicon sheet %d" % (i + 1), "straight", SI, T_SHEET)]
    flow = [
        wafer_setup,
        seq("Sheet stack", stack_steps),
        seq("Stack patterning and STI", [
            dep("Hard mask nitride", "conformal", SIN, HM),
            etch("Sheet stack etch", (SIN, SI, SAC), HM + STACK + BELOW, "Sheet", 1),
            dep("STI fill", "conformal", STI, 90),
            cmp_("STI CMP", (SIN,)),
            etch("STI recess (stack reveal)", (STI,), HM + STACK, kind="basic"),
            strip("Strip hard mask", (SIN,)),
        ]),
        seq("Dummy gate", [
            dep("Dummy gate poly", "planarizing", APOLY, POLY),
            dep("Gate cap", "straight", SIN, CAP),
            etch("Gate etch", (SIN, APOLY), NS_GATE_ETCH, "Gate", 1),
            etch("Gate cut", (SIN, APOLY), NS_GATE_ETCH, "GateCut", 0),
        ]),
        seq("Spacer, S/D recess and inner spacers", [
            dep("Spacer deposit", "conformal", SIN, SP),
            etch("Spacer open in S/D", (SIN,), 200, "SD", 0),
            etch("S/D recess (cut the stack)", (SI, SAC), STACK, "SD", 0),
            etch("Sacrificial indent", (SAC,), INDENT, kind="basic"),
            dep("Inner spacer deposit", "conformal", ISP, ISP_T),
            etch("Inner spacer etch-back", (ISP,), ISP_T, kind="basic"),
        ]),
        seq("Source/drain epi", [
            dep("N block", "conformal", BLOCK, 8, "PBlock"),
            grow("N epi", NEPI, EPI_NS),
            strip("N unblock", (BLOCK,)),
            dep("P block", "conformal", BLOCK, 8, "NBlock"),
            grow("P epi", PEPI, EPI_NS),
            strip("P unblock", (BLOCK,)),
        ]),
        seq("Sheet release and gate all around", [
            dep("ILD0", "planarizing", ILD, 40),
            cmp_("ILD CMP to poly", (APOLY,)),
            etch("Dummy gate removal", (APOLY,), 300),
            strip("Sheet release (remove sacrificial SiGe)", (SAC,)),
            dep("High-k", "conformal", HIK, 2),
            dep("Work function", "conformal", WF, 3),
            dep("Gate W", "conformal", W, 8),
            cmp_("RMG CMP", (ILD,)),
            etch("Gate recess", (W, WF, HIK), 25),
            dep("SAC cap", "conformal", SIN, 14),
            cmp_("SAC cap CMP", (ILD,)),
        ]),
        seq("Contacts", [
            etch("Contact etch (self-aligned)", (ILD,), 250, "Contact", 0),
            dep("Contact liner", "conformal", LINER, 2),
            dep("Contact W", "conformal", W, 8),
            cmp_("Contact CMP", (ILD, SIN)),
        ]),
    ]

# ---- variant: 5T buried power rail (Stage 4) ---------------------------------------
else:
    CUME = cells.MP / 2 - cells.FIN_W / 2 + cells.MP / 4      # 13.25 nm, contact reach under the next MINT track
    BPR_W, BPR_LINER, BPR_BELOW, BPR_CAP, VIA_IN, PWRCOL = 1.5 * cells.MP, 3.0, 20.0, 15.0, 5.0, 0
    borders = [k * H for k in range(4)]
    for d in devs:
        c = (d["y0"] + d["y1"]) / 2.0
        d["border"] = min(borders, key=lambda b: abs(b - c))
        d["con"] = (d["y0"] - CEXT, d["y1"] + CUME)
        d["pwr"] = (min(d["con"][0], d["border"]), max(d["con"][1], d["border"]))
    con_rects = []
    for d in devs:
        for k, (a, b) in enumerate(sd_x):
            ya, yb = d["pwr"] if k == PWRCOL else d["con"]
            con_rects.append((a - 2, ya, b + 2, yb))
    a, b = sd_x[PWRCOL]
    masks = [base_masks[0], base_masks[1],
        ("BPR", [(X0,) + (clipy(bd - BPR_W / 2, bd + BPR_W / 2)[0],) + (X1,) + (clipy(bd - BPR_W / 2, bd + BPR_W / 2)[1],) for bd in borders]),
        base_masks[2], base_masks[3]] + sd_masks + [
        ("Contact", con_rects),
        ("VBPR", [(a, clipy(bd - BPR_W / 2 + VIA_IN, bd + BPR_W / 2 - VIA_IN)[0], b,
                   clipy(bd - BPR_W / 2 + VIA_IN, bd + BPR_W / 2 - VIA_IN)[1]) for bd in borders]),
    ]
    flow = [
        wafer_setup,
        seq("Fins and STI fill", fin_etch_fill),
        seq("Buried power rail", [
            etch("BPR trench etch", (STI, SI), HM + TRENCH + BPR_BELOW, "BPR"),
            dep("BPR liner oxide", "conformal", RLINER, BPR_LINER),
            dep("BPR metal fill", "conformal", CO, 16),
            cmp_("BPR CMP", (SIN,)),
            etch("BPR metal recess", (CO,), HM + REVEAL + BPR_CAP),
            dep("BPR cap oxide", "conformal", STI, 16),
            cmp_("BPR cap CMP", (SIN,)),
        ]),
        seq("Fin reveal", [
            etch("STI recess (fin reveal)", (STI, RLINER), HM + REVEAL, kind="basic"),
            strip("Strip hard mask", (SIN,)),
        ]),
        dummy_gate,
    ] + spacer_epi + [
        rmg_gate_first,
        seq("Contacts and rail via", [
            etch("Contact etch (self-aligned)", (ILD,), 250, "Contact", 0),
            etch("Via to buried rail (through rail cap)", (STI,), BPR_CAP + 25, "VBPR"),
            dep("Contact liner", "conformal", LINER, 2),
            dep("Contact W", "conformal", W, 9),
            cmp_("Contact CMP", (ILD, SIN)),
        ]),
    ]

masks = [(m, i + 1, r) for i, (m, r) in enumerate(masks)]     # GDS layer numbers 1, 2, 3 ...

# ---- GDSII writer ----------------------------------------------------------------
DBU_NM = 0.01                              # database unit: 0.01 nm, so quarter-nanometre edges are exact

def real8(v):
    if v == 0:
        return b"\x00" * 8
    sign = 0x80 if v < 0 else 0
    v = abs(v); e = 64
    while v >= 1.0:
        v /= 16.0; e += 1
    while v < 1.0 / 16.0:
        v *= 16.0; e -= 1
    m = int(round(v * (1 << 56)))
    return bytes([sign | e]) + m.to_bytes(7, "big")

def rec(rtype, dtype, data=b""):
    return struct.pack(">HBB", 4 + len(data), rtype, dtype) + data

def i16(*v): return struct.pack(">%dh" % len(v), *v)
def i32(*v): return struct.pack(">%di" % len(v), *v)
def txt(s):
    b = s.encode("ascii")
    return b + (b"\x00" if len(b) % 2 else b"")

TOP = name.upper()
now = datetime.datetime.now()
stamp = i16(*([now.year, now.month, now.day, now.hour, now.minute, now.second] * 2))
gds = rec(0x00, 2, i16(600)) + rec(0x01, 2, stamp) + rec(0x02, 6, txt("CELLHEIGHT"))
gds += rec(0x03, 5, real8(DBU_NM * 1e-3) + real8(DBU_NM * 1e-9))      # db unit in microns, db unit in metres
gds += rec(0x05, 2, stamp) + rec(0x06, 6, txt(TOP))
nrect = 0
for mname, layer, rects in masks:
    for (xa, ya, xb, yb) in rects:
        q = [int(round(c / DBU_NM)) for c in (xa, ya, xb, yb)]
        pts = [q[0], q[1], q[2], q[1], q[2], q[3], q[0], q[3], q[0], q[1]]
        gds += rec(0x08, 0) + rec(0x0D, 2, i16(layer)) + rec(0x0E, 2, i16(0)) + rec(0x10, 3, i32(*pts)) + rec(0x11, 0)
        nrect += 1
gds += rec(0x07, 0) + rec(0x04, 0)
with open(os.path.join(outdir, name + ".gds"), "wb") as f:
    f.write(gds)

with open(os.path.join(outdir, name + "-layermap.txt"), "w") as f:
    for mname, layer, rects in masks:
        f.write("%-27s %-5d 0\n" % (mname, layer))

with open(os.path.join(outdir, name + ".json"), "w") as f:      # input of the JSON Translator
    json.dump({"process": flow}, f, indent=1)

# ---- model file --------------------------------------------------------------------
zam = """
<SceneGraph ColorMode="Material">

  <DocumentProperties>
    <ModelProperties IsComplete="false" ModelResolution="1.0" ModelUnits="nm"/>
    <Version App="SEMulator3D" AppVersion="13.1.0.0" Number="2.2.12.0"/>
    <LayoutFile Filename="%(n)s.gds" TopCell="%(top)s"/>
    <ProcessFile Filename="%(n)s.vproc"/>
    <Zones Number="1"/>
    <Output Directory="" Enabled="true" ExtractFragments="false"/>
    <ManagerProperties CropOption="false" LayerMapFileName="%(n)s-layermap.txt" RunType="Start From Blank" StartViewer="false"/>
    <ModelArea BoundingMask="BB" Coordinates="0.0 0.0 0.0 0.0" PadBy="0.0" Type="Single Mask"/>
    <UserProperties Author="" AutoRev="true" Comments="%(c)s" CreatedOn="" LastModified="" Revision="1" Watermark=""/>
    <SecurityProperties FileAccessRestriction=""/>
  </DocumentProperties>

  <material blue="50" green="50" name="default" red="50" visibility="true"/>

  <Dopants/>

  <Views/>

  <Lights/>

  <Analyses Comments="" Name="" Transparency="0.0" Visibility="true"/>

</SceneGraph>
""" % dict(n=name, top=TOP, c=title)
with open(os.path.join(outdir, name + ".zam"), "w") as f:
    f.write(zam)

vmpd = "cellheight_nanosheet.vmpd" if variant == "nanosheet" else "16nmFINFET.vmpd"   # nanosheet adds SiGe_Sacrificial
shutil.copyfile(os.path.join(here, "..", "cellheight", "common", vmpd), os.path.join(outdir, name + ".vmpd"))

nsteps = sum(len(s.get("children", [1])) for s in flow)
print("wrote %s: %d masks, %d rectangles, %d steps, top cell %s" % (outdir, len(masks), nrect, nsteps, TOP))
