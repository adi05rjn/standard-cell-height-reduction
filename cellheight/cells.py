# Cell-height geometry from Sherazi et al., SPIE 10588 (2018). All lengths in nm.
# MP = 21 nm, CPP = 42 nm. Cell height = tracks * MP.
# ASSUMPTIONS (not stated explicitly in the paper text): fin width 5 nm, fin pitch 21 nm,
# 5.5T PN separation 20.5 nm (derived so that out-gap equals the 6.5T value of 43 nm).
MP = 21.0
CPP = 42.0
FIN_W = 5.0
FIN_P = 21.0

CELLS = {
    # tracks: (fins per device, PN separation edge-to-edge nm)
    6.5: dict(nfin=2, pn=41.5),
    5.5: dict(nfin=2, pn=20.5),
    4.5: dict(nfin=1, pn=40.0),
    # 5T buried-power-rail cell, Sherazi et al., SPIE 10962 (2019), fin pitch 21 nm option:
    # PN separation = one dummy fin (2 fin pitches centre to centre -> 37 nm edge to edge),
    # cell-to-cell = two dummy fins (3 pitches). Fins sit midway between MINT tracks. CPP is 45 nm here.
    5.0: dict(nfin=1, pn=37.0),
    # 4.5T stacked-nanosheet cell, Sherazi et al. 2018, fig. 8: one 11 nm wide sheet stack per device,
    # N-P separation 40 nm, which leaves 32.5 nm between sheets of neighbouring cells.
    "4.5NS": dict(tracks=4.5, nfin=1, pn=40.0, w=11.0),
}

def cell(tracks):
    c = CELLS[tracks]
    w = c.get("w", FIN_W)                 # width of one fin, or of one nanosheet stack
    H = c.get("tracks", tracks) * MP
    span = (c["nfin"] - 1) * FIN_P + w
    out = H - 2 * span - c["pn"]
    return dict(tracks=c.get("tracks", tracks), H=H, nfin=c["nfin"], pn=c["pn"], span=span, out=out, w=w)

def fins_in_cell(tracks):
    """Return (nfins, pfins) as lists of (y0, y1) for an un-mirrored cell spanning y=[0,H]. N at bottom."""
    g = cell(tracks)
    y = g["out"] / 2.0
    n = [(y + i * FIN_P, y + i * FIN_P + g["w"]) for i in range(g["nfin"])]
    y2 = y + g["span"] + g["pn"]
    p = [(y2 + i * FIN_P, y2 + i * FIN_P + g["w"]) for i in range(g["nfin"])]
    return n, p

def fins_stack(tracks, rows=3):
    """Rows alternate normal / mirrored so abutting cells share the power-rail boundary."""
    g = cell(tracks)
    H = g["H"]
    nf, pf = [], []
    n0, p0 = fins_in_cell(tracks)
    for r in range(rows):
        base = r * H
        if r % 2 == 0:
            nf += [(base + a, base + b) for a, b in n0]
            pf += [(base + a, base + b) for a, b in p0]
        else:
            nf += [(base + H - b, base + H - a) for a, b in n0]
            pf += [(base + H - b, base + H - a) for a, b in p0]
    return sorted(nf), sorted(pf)

if __name__ == "__main__":
    for t in CELLS:
        g = cell(t)
        print(t, g)

def devices_stack(tracks, rows=3):
    """Per-row device groups: list of dict(type='N'|'P', y0, y1) covering the fins of one device."""
    g = cell(tracks)
    H = g["H"]
    n0, p0 = fins_in_cell(tracks)
    devs = []
    for r in range(rows):
        base = r * H
        for typ, fl in (("N", n0), ("P", p0)):
            if r % 2 == 0:
                a, b = base + fl[0][0], base + fl[-1][1]
            else:
                a, b = base + H - fl[-1][1], base + H - fl[0][0]
            devs.append(dict(type=typ, y0=a, y1=b, row=r))
    return sorted(devs, key=lambda d: d["y0"])
