"""
Demonstration / reproduction of the model results used in the manuscript.

Produces two figures (saved to figures/):
  1. fig_dose_independence.png - phospho-STAT6 after a 15-min pulse at 10 vs
     100 ng/mL. With canonical (receptor-proximal) signaling the two doses give
     nearly identical persistence, because receptor occupancy and receptor-
     mediated internalization are both saturated at these doses. A parameter
     sweep (koff, endosomal dissociation, degradation) shows this holds across
     the literature ranges.
  2. fig_koff_decay.png - fraction of surface complex remaining after washout for
     the measured range of dissociation rate constants (koff, t1/2 9-38 min),
     showing that slow dissociation clears within ~1-2 h, far short of the ~8 h
     persistence observed experimentally.

Run:  python demo.py
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from model import default_params, simulate, series

OUT = "figures"
os.makedirs(OUT, exist_ok=True)
C10, C100 = "#4a78b5", "#d1495b"


def canonical(kdeph=0.012, kre=0.02, kdeg=0.010):
    p = default_params()
    p["kdeph"] = kdeph
    p["kphos"] = kdeph / 30.0        # keep the dimer->pSTAT6 half-max (K_half) fixed
    p["kre"] = kre
    p["kdeg"] = kdeg
    return p


def fig_dose_independence():
    p = canonical()
    t10, y10 = simulate(p, 10, 15)
    t100, y100 = simulate(p, 100, 15)
    f = 88.0 / max(series(y10, "Sp").max(), series(y100, "Sp").max())

    # envelope across the uncertain parameters, each dose scaled to its own peak
    lo10 = np.full_like(t10, np.inf); hi10 = np.full_like(t10, -np.inf)
    lo100 = np.full_like(t10, np.inf); hi100 = np.full_like(t10, -np.inf)
    for kd in (0.006, 0.012, 0.020):
        for kre in (0.005, 0.02, 0.076):
            for kg in (0.004, 0.010):
                pp = canonical(kd, kre, kg)
                _, a = simulate(pp, 10, 15); _, b = simulate(pp, 100, 15)
                fa = 88.0 / max(series(a, "Sp").max(), series(b, "Sp").max())
                lo10 = np.minimum(lo10, series(a, "Sp") * fa)
                hi10 = np.maximum(hi10, series(a, "Sp") * fa)
                lo100 = np.minimum(lo100, series(b, "Sp") * fa)
                hi100 = np.maximum(hi100, series(b, "Sp") * fa)

    plt.figure(figsize=(7, 4.5))
    plt.fill_between(t10 / 60, lo10, hi10, color=C10, alpha=0.25,
                     label="10 ng/mL (all params)")
    plt.fill_between(t100 / 60, lo100, hi100, color=C100, alpha=0.25,
                     label="100 ng/mL (all params)")
    plt.plot(t10 / 60, series(y10, "Sp") * f, color=C10, lw=2)
    plt.plot(t100 / 60, series(y100, "Sp") * f, color=C100, lw=2)
    plt.xlabel("Time [h]"); plt.ylabel("phospho-STAT6 [a.u.]")
    plt.title("Canonical signaling is dose-independent (15-min pulse)")
    plt.legend(); plt.grid(alpha=0.3); plt.tight_layout()
    plt.savefig(f"{OUT}/fig_dose_independence.png", dpi=150)
    plt.close()


def fig_koff_decay():
    t = np.linspace(0, 8, 200)
    fast, slow = np.log(2) / 9.0, np.log(2) / 38.0        # koff (min^-1), t1/2 9 & 38 min
    plt.figure(figsize=(7, 4.5))
    plt.fill_between(t, np.exp(-fast * t * 60), np.exp(-slow * t * 60),
                     color="#8aa1c1", alpha=0.4, label="measured koff (t½ 9-38 min)")
    plt.axvline(8, color="k", ls="--", lw=1)
    plt.text(6.0, 0.6, "persistence to\nexplain (~8 h)", fontsize=9)
    plt.xlabel("Time after washout [h]"); plt.ylabel("surface complex remaining (frac.)")
    plt.title("Slow dissociation is too fast")
    plt.legend(); plt.grid(alpha=0.3); plt.tight_layout()
    plt.savefig(f"{OUT}/fig_koff_decay.png", dpi=150)
    plt.close()


if __name__ == "__main__":
    fig_dose_independence()
    fig_koff_decay()
    print(f"Wrote {OUT}/fig_dose_independence.png and {OUT}/fig_koff_decay.png")
