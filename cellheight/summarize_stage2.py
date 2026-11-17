import json, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cells

def merge(spans, tol=0.01):
    out = []
    for a, b in sorted(spans):
        if out and a <= out[-1][1] + tol:
            out[-1][1] = max(out[-1][1], b)
        else:
            out.append([a, b])
    return out

def gaps(iv):
    return [round(iv[i + 1][0] - iv[i][1], 1) for i in range(len(iv) - 1)]

rows = []
for t in (6.5, 5.5, 4.5):
    tag = ("T%.1f" % t).replace(".", "p")
    r = json.load(open(os.path.join("out", "stage2_%s.json" % tag)))
    g = cells.cell(t)
    devs = cells.devices_stack(t, 3)
    gate = merge(r["gate_metal_y_spans"])
    con = merge(r["contact_y_spans"])
    epi = merge(r["epi_y_spans"], 1.5)
    xp = r["x_profile_through_device0"]
    # lateral dielectric between gate metal and contact metal (nm)
    names = [(n, a, b) for n, a, b in xp]
    gm = [(a, b) for n, a, b in names if n in ("W", "TiN", "HfOx")]
    cm = [(a, b) for n, a, b in names if n in ("Ti",)]
    gm_m = merge([s for s in gm]); 
    # distance from first gate edge to first contact edge on its right
    sep = None
    for (ga, gb) in gm_m:
        right = [c for c in cm if c[0] >= gb]
        if right:
            sep = round(min(c[0] for c in right) - gb, 1); break
    exp_con = [(round(d["y0"] - 5, 1), round(d["y1"] + 5, 1)) for d in devs]
    ok_con = len(con) == len(exp_con) and all(abs(a - ea) <= 1.5 and abs(b - eb) <= 1.5 for (a, b), (ea, eb) in zip(con, exp_con))
    cuts = gaps(gate)
    rows.append(dict(
        tracks=t, H=g["H"], gate_cuts=cuts, contacts=len(con), contact_ok=ok_con,
        con_tip2tip_NP=[x for x in gaps(con)][::2][:1] + [x for x in gaps(con)][1::2][:1],
        con_gaps=gaps(con), epi_gaps=gaps(epi), gate_to_contact=sep,
        x_border=r["x_profile_at_cell_border"]))
for r in rows:
    print(json.dumps(r))
