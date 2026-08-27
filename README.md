# il4-trafficking-model

A minimal ODE model of IL-4 receptor trafficking and STAT6 signaling in murine
bone-marrow-derived macrophages (BMDMs), used to interpret the persistence of
phospho-STAT6 signaling after a brief IL-4 pulse.

IL-4 binds its type I receptor (IL-4Rα + common γ-chain, γc) by a **sequential**
mechanism: IL-4 first binds IL-4Rα to form a **binary** complex, which then
recruits γc to form the signaling-competent **ternary complex** ("dimer") that
activates JAK and phosphorylates STAT6. The model has two compartments (plasma
membrane and endosome) connected by constitutive internalization and recycling;
dimerization occurs in both, more efficiently in the endosome (higher local
receptor density; Gandhi et al., 2014). Phospho-STAT6 is a fast, saturable
readout slaved to the total signaling dimer, consistent with STAT6 having no
intrinsic memory (Andrews et al., 2002).

**Key result.** With literature-grounded parameters, canonical (receptor-proximal)
signaling produces **dose-independent** phospho-STAT6 persistence — because
receptor occupancy and receptor-mediated internalization are both saturated at
10 and 100 ng/mL — and slow ligand dissociation clears within ~1–2 h. Neither can
account for the observed dose-dependent, hours-long persistence of a 15-min pulse.

## Requirements

Python ≥ 3.9 with `numpy`, `scipy`, `matplotlib`:

```bash
pip install -r requirements.txt
```

## Usage

```bash
python model.py   # prints Kd, basal receptor numbers, peak pSTAT6 for a 100 ng/mL pulse
python demo.py    # writes the model figures to figures/
```

In code:

```python
from model import default_params, simulate, series
p = default_params()
t, y = simulate(p, dose_ngml=100, pulse_min=15)   # 15-min pulse then washout
pstat6 = series(y, "Sp")
```

## Model structure

State (all in molecules/cell): surface `Ls, Rs, Gs, Bs, Ds`; endosomal
`Le, Re, Ge, Be, De`; phospho-STAT6 `Sp`. `R` = IL-4Rα, `G` = γc, `B` = binary
complex, `D` = signaling dimer. STAT6 readout:
`dSp/dt = kphos·(Ds+De)·(Stot − Sp) − kdeph·Sp`.

Ligand dose is a bath concentration (ng/mL) converted to molecules/cell via an
assumed per-cell volume; because the on-rate is converted with the same volume,
that volume cancels from the dynamics and `kf`/`kr` reproduce the measured Kd.

## Parameters

| Parameter | Symbol | Value | Basis / source |
|---|---|---|---|
| Resting surface IL-4Rα | Rs(0) | ~2,300 /cell | Crawford et al., 1987 |
| IL-4 : IL-4Rα association | kon | 4.4–13.2 × 10⁶ M⁻¹s⁻¹ | Lowenthal 1988; Seipelt 1997; Zhang 2002; Andrews 2006; Kraich 2006 |
| IL-4 : IL-4Rα dissociation | koff | 3–13 × 10⁻⁴ s⁻¹ | (as above) |
| Equilibrium dissociation constant | Kd | 20–300 pM (96 pM used) | (as above) |
| Receptor/complex internalization | ke | 0.025 min⁻¹ (t½ ~13–28 min) | Galizzi 1989; Kuznetsov & Puri 1999 |
| Degradation | kdeg | 0.010 min⁻¹ (t½ 68 min) | Galizzi 1989 |
| STAT6 dephosphorylation | kdeph | 0.006–0.012 min⁻¹ (t½ ~1–2 h) | this study |
| Endosomal concentration factor | f_endo | 10³ | assumed (concept: Gandhi 2014) |
| γc, dimerization, recycling | — | see `model.py` | assumed / derived |

Values without a measured IL-4/macrophage source were set within physiological
ranges or fixed by the basal steady state; the dose-independence result is robust
across the reported ranges (see `demo.py`).

## Reproducing the manuscript figures

`demo.py` regenerates the model curves. The experimental data used for the
data-overlay panels will be added on manuscript submission.

## Citation

This model accompanies:

> [Authors]. [Title]. [Journal / preprint], [year]. [DOI]

(Details to be completed at submission.)

## License

MIT — see [LICENSE](LICENSE).
