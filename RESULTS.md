# Results: models checked against the papers

Every spacing the two papers draw in their figures is reproduced by the models to within the
1 nm model grid. This file lists each number, where it comes from, and what the model measures.

Sources:
- **[2018]** Sherazi et al., "Track height reduction for standard-cell in below 5nm node: How low
  can you go?", Proc. SPIE 10588.
- **[2019]** Sherazi et al., "Standard-cell design architecture options below 5nm node", Proc. SPIE 10962.

Measured values come from `cellheight/summarize_all.py`, which reads the `.json` written by each build.
The models use 1 nm voxels, so a drawn 10.5 nm appears as 10 or 11 nm.

## 1. Paper versus model

### 6.5T, 5.5T and 4.5T cells (Stage 2)

| Quantity | Paper | Source | Model |
|---|---|---|---|
| Cell height | 136.5 / 115.5 / 94.5 nm | [2018] figs. 2, 3, 6 | 136.5 / 115.5 / 94.5 nm (input) |
| Fin width | 5 nm | [2018] fig. 6 | 5 nm |
| Fin gap across the cell border | 43 / 43 / 44.5 nm | [2018] figs. 2, 4, 7 | 43 / 43 / 44–45 nm |
| N-to-P fin gap, 6.5T | 41.5 nm | [2018] fig. 2 | 41.5 nm (input) |
| N-to-P fin gap, 4.5T | 40 nm | [2018] figs. 6, 7 | 40 nm |
| N-to-P contact tip-to-tip, 5.5T ("tight cut") | 10.5 nm | [2018] fig. 4 | 10–11 nm |
| N-to-P contact tip-to-tip, 4.5T | 30 nm | [2018] fig. 7 | 30 nm |
| Border tip-to-tip, neither contact extended, 5.5T | 33 nm | [2018] fig. 4 | 33 nm |
| Border tip-to-tip, one contact extended to the rail, 5.5T | 16.5 nm | [2018] fig. 4 | 17 nm |
| Border tip-to-tip, one contact extended to the rail, 4.5T | 17.25 nm | [2018] fig. 7 and text | 17–18 nm |
| Gate cut | 18 nm | [2019] fig. 9 | 18 nm |

### 4.5T mid-track handshake (Stage 3)

| Quantity | Paper | Source | Model |
|---|---|---|---|
| Contacts before metal gate, N then P, deep plus shallow pattern, spacer replaced twice | process description | [2018] sec. 4.1, figs. 10–16 | same order of steps |
| N wing and P wing in different contact columns | drawn | [2018] figs. 9, 19 | N wing in column 0, P wing in column 1 |
| Tip-to-tip between a wing and the opposite contact | 10.5 nm | [2018] figs. 9, 19 | 10 nm (two 5.25 nm low-k spacers) |
| Other fin stays covered under the wing | "remains covered by dielectric" | [2018] fig. 18 | oxide from the isolation up to the wing floor |
| Bridge on the centre track joining the two wings | drawn | [2018] fig. 17 | cobalt bridge, vias 10 nm long on each wing |

### 4.5T stacked-nanosheet cell (Stage 5)

| Quantity | Paper | Source | Model |
|---|---|---|---|
| Sheet width | 11 nm | [2018] fig. 8 and text | 11 nm |
| N-to-P sheet gap | 40 nm | [2018] fig. 8 | 40 nm |
| Sheet gap across the cell border | 32.5 nm | [2018] fig. 8 | 32–33 nm |
| N-to-P contact tip-to-tip | 30 nm | [2018] fig. 8 | 30 nm |
| Border tip-to-tip beside a power tap | 10.5 nm | [2018] fig. 8 and text | 10–11 nm |
| Gate wraps each sheet on all four sides | gate all around | [2018] sec. 4 | 2 nm high-k and 3 nm TiN around every sheet, top, bottom and sides |

The power-tap contact is drawn 17 nm past the sheet, which is what leaves the paper's 10.5 nm. That is
0.75 nm beyond the cell boundary; extending only to the boundary would leave 11.25 nm.

### 5T buried-power-rail cell (Stage 4)

| Quantity | Paper | Source | Model |
|---|---|---|---|
| Gate pitch / metal pitch / cell height | 45 / 21 / 105 nm | [2019] sec. 1, fig. 6 | 45 / 21 / 105 nm (input) |
| N-to-P fin gap | 37 nm | [2019] fig. 6 | 37 nm |
| Fin gap across the cell border | 58 nm | [2019] fig. 6 | 58 nm |
| Rail width | 3 CD = 31.5 nm | [2019] sec. 2, fig. 6 | 32 nm trench |
| Fin to rail | 13.25 nm ("above 10 nm") | [2019] fig. 6 and text | 13 nm |
| Contact tip-to-tip at the cell centre | 18.75 nm | [2019] fig. 6 | 19 nm |
| Contact tip-to-tip at the cell border | 39.75 nm | [2019] fig. 6 | 40 nm |
| Power contacts of abutting cells | merge at the rail | [2019] sec. 2.4, fig. 6 | one merged contact landing on the rail |

## 2. What was corrected after reading the figures

The first versions of the models were built from text extracts in which the figures were unreadable.
Reading the PDFs led to three corrections:

1. **Power-tap contacts added to Stage 2.** The paper's 17.25 nm (4.5T) and 16.5 nm (5.5T) border
   spacings are between a contact extended to the cell boundary and the neighbour's short contact.
   The models previously had only short contacts, which shows the relaxed 33–34.5 nm case. One contact
   column now carries the power-tap extension in the middle row; the other column is unchanged.
2. **Handshake spacing in Stage 3 changed from 12 nm to 10.5 nm.** The low-k spacer is now 5.25 nm, so
   two of them give the 10.5 nm the paper draws.
3. **Contact extension under the routing track added to Stage 4.** The 2019 paper uses two extensions:
   5 nm over the fin and 13.25 nm to sit fully under the next routing track. With both, the model
   reproduces 18.75 nm and 39.75 nm. Previously both sides were 5 nm, giving 27 nm.

Also confirmed without change: the 5 nm fin width, the staggered wings, and the 5T fin positions
(37 nm and 58 nm had been derived from the text and match figure 6 exactly).

## 3. Known inconsistency in the 2018 paper

Figure 3 and figure 4 label the 5.5T N-to-P fin gap as 21.5 nm. The model uses 20.5 nm. The paper's
own numbers do not close with 21.5 nm:

- 6.5T fixes the two-fin device width at 26 nm (136.5 − 43 − 41.5 = 52 = 2 × 26).
- At 5.5T, 43 + 2 × 26 + 21.5 = 116.5 nm, which is 1 nm more than the 115.5 nm cell height.
- A 10.5 nm contact tip-to-tip with 5 nm extensions needs a 20.5 nm fin gap, not 21.5 nm.

With 20.5 nm the cell height, the 43 nm border gap and the 10.5 nm tight cut all agree, so the model
keeps 20.5 nm.

## 4. What is still an assumption

The papers are design-rule studies and give no vertical or film dimensions. These values are not from
the papers:

| Assumption | Value | Effect |
|---|---|---|
| Gate length, nitride spacer | 16 nm, 6 nm | Contact opening between gates (14 nm at 42 nm pitch). |
| Epi growth | 8 nm, equal on all faces | N-to-P epi gap (8–9 nm at 5.5T). Real epi is faceted. |
| Fin height, isolation depth | 45 nm, 100 nm | Vertical scale. |
| Dummy gate height, caps | 60 nm, 20–25 nm | Vertical scale. |
| Gate stack, contact liner | 2 nm high-k, 3 nm TiN, W; 2 nm Ti | Film thicknesses. |
| Stage 3 wing floor and thickness | 15 nm above fin top, 20 nm | Staircase height. |
| Stage 3 bridge | 20 nm cobalt, no liner | Bridge thickness. |
| Stage 4 rail depth, liner, cap, metal | 20 nm into silicon, 3 nm, 15 nm, cobalt | Rail cross-section. The 31.5 nm includes the liner. |
| Stage 4 fin pitch option | 21 nm, on grid | The paper's 24 nm option is not modelled. |
| Stage 5 sheet stack | 3 sheets, 6 nm thick, 10 nm apart | Stack height 48 nm; the paper gives the width only. |
| Stage 5 inner spacer and epi | 5 nm, and 10 nm epi | 8 nm of epi does not merge across the 10 nm gaps between sheets, so 10 nm is used. |
| Materials | SEMulator3D 16 nm FinFET example library | Names and colours only. |

## 5. Limits

- The low-k film stays on the sides of the epi in Stage 3, so the contact touches the epi top only.
- Only the handshake track is built above the contacts. Other routing tracks, the power rails above
  the cell in Stages 2 and 3, and gate contacts are not modelled.
- No process variation (overlay, etch bias, film thickness) has been run. All results are nominal.
- Nanosheets are modelled at 4.5T only. The 2019 paper's nanosheet insertion with the buried rail is not.

## 6. Manager projects against script models

Each cell was also built as a Manager project (GDS layout plus process file) with ModelBuilder, then
rebuilt from scratch and measured with the same probes as the script model (`manager/check_project.py`).

| Project | Compared | Result |
|---|---|---|
| 6.5T cell | contact gaps in both columns, epi gaps, profile along the gate direction | within 1 nm; gate-direction profile identical |
| 5.5T cell | same | within 1 nm; gate-direction profile identical |
| 4.5T cell | same | within 1 nm; three of four identical |
| 4.5T handshake | low-k widths, contact spans, bridge, both vias, via stack | within 1 nm; bridge, vias and stack identical |
| 5T buried rail | contact gaps, rail-level profile, three vertical stacks | identical |
| 4.5T nanosheet | contact gaps, sheet gaps, stack under the gate, profiles along and between the sheets, source/drain stack | identical |

The 1 nm differences are along the fin direction only. The script models round mask edges to 1 nm,
while the GDS layouts keep 0.01 nm. Rebuilding the 5.5T script model with fine masks gave exactly
the Manager model's contact gaps (10 / 16 / 33 nm), which confirms the cause.

The 5.5T project was opened in the Manager and its Process Editor; the other four were built and
checked from the command line only.

## 7. Whole-model verification

The comparisons above probe each model along a few lines. `verify/verify_model.py` instead pulls the
full 3D voxel array out of each model and analyses every voxel, at 1 nm and again at 0.5 nm.

### Connectivity: no shorts, no opens

Every connected piece of metal is labelled (scipy, 6-connected) and classified by the materials it
contains. The expected counts follow from the design: three gate lines cut into three segments each,
and one contact per device per column, merged only where the design joins them.

| Cell | Metal pieces found | Expected | Gate-to-contact shorts | Stray slivers | Pieces touching both N and P epi |
|---|---|---|---|---|---|
| 6.5T | 21: 9 gate, 12 contact | 21 | 0 | 0 | 0 |
| 5.5T | 21: 9 gate, 12 contact | 21 | 0 | 0 | 0 |
| 4.5T | 21: 9 gate, 12 contact | 21 | 0 | 0 | 0 |
| 4.5T handshake | 18: 9 gate, 6 contact, 3 bridged N-P pairs | 18 | 0 | 0 | 3, one per bridge |
| 5T buried rail | 19: 9 gate, 6 signal contact, 4 rail groups | 19 | 0 | 0 | 0 |
| 4.5T nanosheet | 21: 9 gate, 12 contact | 21 | 0 | 0 | 0 |

N epi and P epi never touch directly in any cell. The counts are the same at 1 nm and at 0.5 nm.

### Minimum separation: the tight spot is where the papers put it

For every pair of different metal pieces, the smallest dielectric gap along x, y and z is found
over the whole model.

| Cell | Smallest contact-to-contact gap | Paper | Gate metal to contact | Other |
|---|---|---|---|---|
| 6.5T | 16 nm at 1 nm grid, 16.5 nm at 0.5 nm | 16.5 nm (derived) | 8 nm | |
| 5.5T | 10 nm, then 10.5 nm | 10.5 nm | 8 nm | |
| 4.5T | 17 nm, then 17.0 nm | 17.25 nm | 8 nm | |
| 4.5T handshake | 10 nm, then 10.0 nm | 10.5 nm | 7 nm | bridge to gate 25 nm |
| 5T buried rail | 19 nm, then 18.5 nm | 18.75 nm | 7 nm, then 8 nm | rail to gate 17 nm |
| 4.5T nanosheet | 10 nm, then 10.0 nm | 10.5 nm | 8 nm | |

- Nothing tighter than the expected spacing exists anywhere in any model. The smallest gap overall is
  always gate metal to contact: the 6 nm spacer plus the 2 nm high-k (5 nm low-k plus 2 nm in the handshake).
- Earlier sections quote "6 nm gate to contact". That figure counted the high-k as part of the gate;
  metal to metal it is 8 nm.

### Finer grid: the 1 nm differences are grid rounding

At 0.5 nm the measured gaps land on the paper values wherever a 0.5 nm grid can represent them:
10.5 nm at 5.5T and 16.5 nm at 6.5T exactly, and 17.0 and 18.5 nm against 17.25 and 18.75 nm
(both a quarter-nanometre off, the limit of the grid).

One value does not converge: the handshake gap stays at 10.0 nm against the paper's 10.5 nm. It is
made of two 5.25 nm deposited films, and a film thickness of 5.25 nm is rounded to 5.0 nm on both
grids. Reaching 10.5 nm would need a 0.25 nm grid.

### Voids: every fill is complete

The same script searches for empty voxels that are not connected to the open top surface.

| Cell | Voxels in the model | Empty voxels | Enclosed voids |
|---|---|---|---|
| 6.5T | 10,848,600 | 0 | 0 |
| 5.5T | 9,181,620 | 0 | 0 |
| 4.5T | 7,514,640 | 0 | 0 |
| 4.5T handshake | 8,230,320 | 0 | 0 |
| 5T buried rail | 8,930,250 | 0 | 0 |
| 4.5T nanosheet | 9,232,272 | 0 | 0 |

As a cross-check, the voxel counts of the named materials add up to the full array size in each
model, so no voxel is unaccounted for. The contact, gate, low-k spacer, bridge and rail fills all
closed completely.

### Design checks with no process model

`verify/check_design.py` runs 81 checks that involve no build at all, and all 81 pass.

- **Arithmetic (32 checks).** Every paper dimension is recomputed by formula from the 21 nm pitches
  and the 5 nm fin and contact extension, and compared with what `cells.py` produces. This includes
  the closure of each cell height: twice the device width plus the N-to-P gap plus the border gap.
- **Layout read-back (49 checks).** The six GDS files are read with a separate reader and the gaps
  are measured in 2D. They match the paper values exactly, to 0.01 nm: 10.5, 16.5, 17.25, 30, 31.5,
  33 and 34.5 nm for the cells, and 37, 58, 31.5, 13.25, 18.75 and 39.75 nm for the buried-rail cell.

This separates the two possible sources of error. The drawn design matches the papers exactly; the
quarter- to half-nanometre differences in the built models come only from the voxel grid.

### Clean rebuild: the study reproduces from its own files

All existing outputs were moved aside, and everything was rebuilt from the scripts.

| Rebuilt | Builds | Result |
|---|---|---|
| Script models, Stages 1 to 4 | 8, all exit code 0 | All 8 result files identical to the previous ones, value for value |
| Manager projects | 5, all exit code 0 | All 5 comparison files identical to the previous ones |

### Not covered by these checks

- Gaps are measured along the three axes only, not diagonally.
- The void search works on whole voxels, so a pocket smaller than the 1 nm grid would not be seen.
- The checks confirm the models match the design and the design matches the papers' figures as I
  read them. They cannot confirm the dimensions the papers do not give (section 4), and they do not
  replace a second person checking the figure readings.

## 8. Figures

Cross-sections of the corrected models are in `figures/`:

| File | Shows |
|---|---|
| `stage2_5p5T_contact_cut_x42.png` | 5.5T contacts on epi; the middle cell's contacts run out to the cell borders (power tap). |
| `stage3_4p5T_handshake_cut_x42.png` | 4.5T N contacts stepping toward the cell centre, with vias up to the bridge. |
| `stage3_4p5T_bridge_cut_y142.png` | Cut along the centre track: the bridge crosses the gate and lands on both wings. |
| `stage4_5T_BPR_power_column_x42.png` | 5T power contacts landing on the buried rail. |
| `stage4_5T_BPR_signal_column_x90.png` | 5T signal contacts, isolated from the rail. |

## 9. Reproducing

From the `cellheight` folder:

```
run_stage2.bat 6.5
run_stage2.bat 5.5
run_stage2.bat 4.5
run_stage3.bat 4.5
run_stage4.bat
python summarize_all.py
```
