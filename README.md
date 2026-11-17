# Standard Cell-Height Reduction

SEMulator3D 13.1 models of the standard-cell track-height options in two imec papers:

- Sherazi et al., SPIE 10588 (2018): 6.5T, 5.5T and 4.5T cells, and the 4.5T mid-track handshake.
- Sherazi et al., SPIE 10962 (2019): the 5T cell with a buried power rail.

## Folders

| Folder | Contents |
|---|---|
| `cellheight/` | The study: geometry (`cells.py`), the five build scripts, run helpers, and measured results in `out/*.json`. |
| `manager/` | All six cells as Manager projects (layout, layer map, process JSON) and the scripts that generate and check them. |
| `verify/` | Verification. `run_verify.bat <variant> [0.5]` checks connectivity, shorts, voids and minimum separation over the full voxel array and writes a `verify_*.json` per run. `check_design.py` runs 81 arithmetic and layout checks with no build. |
| `figures/` | Cross-section images of the built models. The `slice_*.png` files are rendered directly from the voxel data. |

Built models (`.zam`, Viewer data), generated process files (`.vproc`) and SEMulator3D example files are
not in the repository; the next section says how to regenerate or supply them.

## Setup

1. Install SEMulator3D 13.1 (default path `C:\Coventor\SEMulator3D13`) and make a license available:
   set the environment variable `COVLD_LICENSE_FILE` to `port@host` of your license server.
2. Create `cellheight/common/` and copy these vendor files into it from the SEMulator3D install
   (`Examples\Semiconductor`):
   - from `16nmFinFET`: `16nmFINFET.vmpd` and `16nmFINFET.vproc`
   - from `STI`: `STI.vmpd` and `STI_basic_prebuilt.vproc`
3. Run `python make_nanosheet_materials.py` inside the SEMulator3D shell once. It writes
   `cellheight/common/cellheight_nanosheet.vmpd`, which Stage 5 needs.

The two source papers are not included. Their full citations are above.

## Stages

| Stage | Script | Run | Model written to `cellheight/out/` |
|---|---|---|---|
| 1. Fins and isolation | `build_stage1.py` | `run_stage1.bat 6.5` (or 5.5, 4.5) | `stage1_T6p5.zam` etc. |
| 2. Full device: gate, spacers, epi, metal gate, contacts | `build_stage2.py` | `run_stage2.bat 6.5` (or 5.5, 4.5) | `stage2_T6p5.zam` etc. |
| 3. 4.5T mid-track handshake with double low-k spacer and centre-track bridge | `build_stage3.py` | `run_stage3.bat 4.5` | `stage3_T4p5.zam` |
| 4. 5T cell with buried power rail | `build_stage4.py` | `run_stage4.bat` | `stage4_T5p0_BPR.zam` |
| 5. 4.5T cell with stacked nanosheets (gate all around) | `build_stage5.py` | `run_stage5.bat` | `stage5_T4p5_nanosheet.zam` |

Run the `.bat` files from the `cellheight` folder. Each build takes about a minute and prints its
measurements, which are also saved as a `.json` next to the model. `summarize_stage2.py` prints the
Stage 2 comparison across the three track heights.

Stages 2 to 4 accept a stop label to end early at a checkpoint, for example `run_stage3.bat 4.5 rmg`
builds the handshake model without the bridge.

## Manager project (`manager/`)

Every cell also exists as a normal SEMulator3D Manager project, with a layout file and a process
file, so it can be opened, edited in the Process Editor and rebuilt the standard SEMulator3D way.

| Project folder in `manager/` | Cell | Steps | Masks | Build command |
|---|---|---|---|---|
| `cell_6p5T` | 6.5T cell | 33 | 8 | `build_project.bat 6.5` |
| `cell_5p5T` | 5.5T cell | 33 | 8 | `build_project.bat 5.5` |
| `cell_4p5T` | 4.5T cell | 33 | 8 | `build_project.bat 4.5` |
| `handshake_4p5T` | 4.5T mid-track handshake with bridge | 56 | 16 | `build_project.bat handshake` |
| `bpr_5T` | 5T cell with buried power rail | 41 | 10 | `build_project.bat bpr` |
| `nanosheet_4p5T` | 4.5T cell with stacked nanosheets | 44 | 8 | `build_project.bat nanosheet` |

Each folder holds:

| File | What it is |
|---|---|
| `<name>.zam` | Model file. Open this in the Manager. |
| `<name>.vproc` | Process file, organised in named sequences. |
| `<name>.gds`, `<name>-layermap.txt` | Layout with one layer per mask, and the layer names. |
| `<name>.vmpd` | Material file, copied from `cellheight/common` when the project is generated. |
| `<name>.json` | The process written as JSON, the input of the tool's JSON Translator. |
| `<name>_check.json` | Comparison of the Manager model with the script model. |

`build_project.bat` regenerates these from `cellheight/cells.py` and builds the model with
ModelBuilder. `check_project.py <variant>` rebuilds a project from scratch and compares it with the
script model. All six agree with their script models to within 1 nm; the buried-rail and nanosheet ones are identical.

Notes:
- The Manager shows a warning about `hdm`, `simultaneousDeposit` and `transport` when it opens a
  file. It comes from the tool's own process library, not from this project.
- The Manager writes an `apps` folder into the folder it is started from. It is safe to delete.
- The Manager models keep mask edges to 0.01 nm; the script models round them to 1 nm. Gaps along
  the fins therefore differ by up to 1 nm between the two (for example 16 nm against 17 nm).
- A process file cannot measure the model, so heights the scripts read back (fin top, gate top)
  are fixed numbers in `make_project.py`.

## Viewing

Open a `.zam` in the **SEMulator3D 13 Viewer**, not the Manager. The models are built by script and
have no process file, so the Manager cannot rebuild them.

- Double-click a step in the Scene list to see the structure after that step.
- For a cross-section: View > Enable Cross Sections, then Options > Align Vertical, then type a
  position into the X field of the Cross-section panel.

## Requirements

- SEMulator3D 13.1 at `C:\Coventor\SEMulator3D13`.
- A reachable license server, named in `COVLD_LICENSE_FILE` (see Setup).

## Where the numbers come from

Cell heights, gate and metal pitch, fin counts, N-to-P separations, the 5 nm contact extension, the
18 nm gate cut and the 31.5 nm buried-rail width are from the papers.

These are assumptions, set at the top of each build script and of `cells.py`:

- Fin height 45 nm (fin width 5 nm and fin pitch 21 nm are confirmed by the paper figures).
- Gate length 16 nm, spacer 6 nm, epi growth 8 nm.
- Stage 3: wing floor 15 nm above the fin top, 20 nm cobalt bridge.
- Stage 4: rail 20 nm into the silicon, 3 nm liner, 15 nm cap, cobalt rail.
- Stage 5: three sheets 6 nm thick with 10 nm between them, 5 nm inner spacers, 10 nm epi.
  The sacrificial layers use `SiGe_Sacrificial`, a renamed copy of SiGe made by
  `make_nanosheet_materials.py` (see Setup).
- Materials are taken from the SEMulator3D 16 nm FinFET example library (not included).

`RESULTS.md` compares every spacing drawn in the papers with what the models measure, and lists the
remaining assumptions and limits. `cellheight/summarize_all.py` prints the measured values, and
`figures/` holds cross-sections of the models.

SEMulator3D is a product of Lam Research (Coventor). This repository is independent of it and contains
none of its files.
