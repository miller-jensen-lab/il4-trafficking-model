"""
Basic ODE model of IL-4 receptor trafficking and STAT6 signaling.

IL-4 binds its type I receptor (IL-4Ralpha + common gamma-chain, gamma_c) by a
sequential mechanism: IL-4 first binds IL-4Ralpha to form a binary complex, which
then recruits gamma_c to form the signaling-competent ternary complex ("dimer")
that activates JAK and phosphorylates STAT6. The model has two compartments
(plasma membrane and endosome) connected by constitutive internalization and
recycling; dimerization occurs in both, more efficiently in the endosome
(higher local receptor density; Gandhi et al., 2014). STAT6 phosphorylation is a
fast, saturable readout slaved to the total signaling dimer (STAT6 has no
intrinsic memory; Andrews et al., 2002).

State vector (all species in molecules/cell):
    [Ls, Rs, Gs, Bs, Ds, Le, Re, Ge, Be, De, Sp]
    surface : Ls free IL-4, Rs IL-4Ralpha, Gs gamma_c, Bs binary, Ds dimer
    endosome: Le, Re, Ge, Be, De  (same species)
    Sp      : phospho-STAT6

Units: molecules/cell and min^-1. Ligand dose is a bath concentration (ng/mL)
converted to molecules/cell via an assumed per-cell volume; because the on-rate
kf is converted with the same volume, that volume cancels from the dynamics and
kf/kr reproduce the measured Kd.

See README.md for parameter values and literature sources.
"""

import numpy as np
from scipy.integrate import solve_ivp

SPECIES = ["Ls", "Rs", "Gs", "Bs", "Ds", "Le", "Re", "Ge", "Be", "De", "Sp"]
IDX = {s: i for i, s in enumerate(SPECIES)}

# --- unit bookkeeping (cancels in the dynamics) ------------------------------
NA = 6.022e23           # molecules / mol
V_CELL = 1e-9           # L per cell (any value gives identical dynamics)
MW_IL4 = 14000.0        # g/mol, mouse IL-4
MOLAR_TO_PERCELL = NA * V_CELL
KON_M_S = 1.32e7        # M^-1 s^-1 (Kraich et al., 2006, WT IL-4)
KF_PERCELL = KON_M_S * 60.0 / MOLAR_TO_PERCELL   # (molecules/cell)^-1 min^-1


def ngml_to_molecules(ng_per_mL):
    """IL-4 bath concentration (ng/mL) -> molecules/cell."""
    return ng_per_mL * 1e-6 / MW_IL4 * MOLAR_TO_PERCELL


def default_params():
    """Literature-grounded parameter set (see README for sources)."""
    p = dict(
        # binding: L + R <-> B
        kf=KF_PERCELL, kr=0.076,            # kr = koff (Kraich 2006); Kd = kr/kf ~ 96 pM
        kfe=None, kre=0.02,                 # endosomal binding (kfe set from f_endo below)
        # dimerization: B + gamma_c <-> D
        kdim=1e-7, f_endo=1000.0,           # endosomal concentration factor (Gandhi 2014)
        kundim=0.1,
        # constitutive trafficking
        keR=0.025, keG=0.025, keC=0.025, keD=0.025,
        krecR=0.02, krecG=0.02, krecC=0.02, krecD=0.02, krecL=0.002,
        kdeg=0.010,                         # first-order degradation of endosomal species
        # receptor synthesis (set below to hit resting surface numbers)
        vsR=0.0, vsG=0.0, target_Rs=2300.0, target_Gs=5000.0,
        # STAT6 readout: dSp/dt = kphos*D*(Stot-Sp) - kdeph*Sp
        Stot=1000.0, kphos=4e-4, kdeph=0.012,
        clamp_L=False,
    )
    p["kfe"] = p["kf"] * p["f_endo"]        # endosomal binding is concentration-enhanced
    # ligand-free basal steady state fixes the synthesis rates
    p["vsR"] = p["keR"] * p["target_Rs"] * p["kdeg"] / (p["krecR"] + p["kdeg"])
    p["vsG"] = p["keG"] * p["target_Gs"] * p["kdeg"] / (p["krecG"] + p["kdeg"])
    return p


def basal_state(p):
    """Resting (ligand-free) initial condition."""
    y = np.zeros(len(SPECIES))
    y[IDX["Rs"]] = p["target_Rs"]
    y[IDX["Re"]] = p["keR"] * p["target_Rs"] / (p["krecR"] + p["kdeg"])
    y[IDX["Gs"]] = p["target_Gs"]
    y[IDX["Ge"]] = p["keG"] * p["target_Gs"] / (p["krecG"] + p["kdeg"])
    return y


def rhs(t, y, p):
    Ls, Rs, Gs, Bs, Ds, Le, Re, Ge, Be, De, Sp = y
    kdimE = p["kdim"] * p["f_endo"]

    bind_s = p["kf"] * Ls * Rs - p["kr"] * Bs        # surface  L + R <-> B
    dim_s = p["kdim"] * Bs * Gs - p["kundim"] * Ds   # surface  B + G <-> D
    bind_e = p["kfe"] * Le * Re - p["kre"] * Be      # endosome L + R <-> B
    dim_e = kdimE * Be * Ge - p["kundim"] * De       # endosome B + G <-> D

    dLs = 0.0 if p.get("clamp_L") else (-bind_s + p["krecL"] * Le)
    dRs = p["vsR"] - bind_s - p["keR"] * Rs + p["krecR"] * Re
    dGs = p["vsG"] - dim_s - p["keG"] * Gs + p["krecG"] * Ge
    dBs = bind_s - dim_s - p["keC"] * Bs + p["krecC"] * Be
    dDs = dim_s - p["keD"] * Ds + p["krecD"] * De

    dLe = -bind_e - p["krecL"] * Le - p["kdeg"] * Le
    dRe = p["keR"] * Rs - p["krecR"] * Re - bind_e - p["kdeg"] * Re
    dGe = p["keG"] * Gs - p["krecG"] * Ge - dim_e - p["kdeg"] * Ge
    dBe = p["keC"] * Bs - p["krecC"] * Be + bind_e - dim_e - p["kdeg"] * Be
    dDe = p["keD"] * Ds - p["krecD"] * De + dim_e - p["kdeg"] * De

    J = Ds + De                                       # total signaling dimer ~ active JAK
    dSp = p["kphos"] * J * (p["Stot"] - Sp) - p["kdeph"] * Sp
    return [dLs, dRs, dGs, dBs, dDs, dLe, dRe, dGe, dBe, dDe, dSp]


def integrate_segments(p, y0, segments):
    """Integrate through a list of segments, carrying state across boundaries and
    resetting surface free ligand (Ls) at each. Each segment is a dict:
        {"Ls": <molecules/cell>, "clamp_L": <bool>, "tspan": [t0, t1]}
    """
    times, states = [], []
    y = np.array(y0, dtype=float)
    for seg in segments:
        ps = dict(p)
        ps["clamp_L"] = seg.get("clamp_L", False)
        y = y.copy()
        y[IDX["Ls"]] = seg["Ls"]
        t0, t1 = seg["tspan"]
        t_eval = np.linspace(t0, t1, 400)
        sol = solve_ivp(lambda t, yy: rhs(t, yy, ps), (t0, t1), y,
                        t_eval=t_eval, method="LSODA", rtol=1e-6, atol=1e-3)
        times.append(sol.t)
        states.append(sol.y)
        y = sol.y[:, -1]
    return np.concatenate(times), np.concatenate(states, axis=1)


def equilibrate(p, t_eq=20000.0):
    """Return the ligand-free basal steady state as an initial condition."""
    _, y = integrate_segments(p, basal_state(p), [{"Ls": 0.0, "tspan": [0, t_eq]}])
    return y[:, -1]


def simulate(p, dose_ngml, pulse_min, total_min=600.0):
    """Simulate stimulation with `dose_ngml` for `pulse_min` (a bath held constant,
    then washed out), or continuous stimulation if pulse_min >= total_min.
    Returns (time_minutes, states) where states[IDX[name]] is a species trace."""
    y0 = equilibrate(p)
    dose = ngml_to_molecules(dose_ngml)
    if pulse_min >= total_min:
        segs = [{"Ls": dose, "clamp_L": True, "tspan": [0, total_min]}]
    else:
        segs = [{"Ls": dose, "clamp_L": True, "tspan": [0, pulse_min]},
                {"Ls": 0.0, "clamp_L": True, "tspan": [pulse_min, total_min]}]
    return integrate_segments(p, y0, segs)


def series(states, name):
    return states[IDX[name], :]


if __name__ == "__main__":
    p = default_params()
    yb = equilibrate(p)
    Kd_pM = p["kr"] / p["kf"] / MOLAR_TO_PERCELL * 1e12
    print(f"Kd = {Kd_pM:.0f} pM")
    print("Basal (molecules/cell): "
          + ", ".join(f"{s}={yb[IDX[s]]:.0f}" for s in SPECIES if yb[IDX[s]] > 1))
    t, y = simulate(p, dose_ngml=100, pulse_min=15)
    print(f"peak pSTAT6 (a.u.) = {series(y, 'Sp').max():.1f}")
