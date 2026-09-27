# IHIS
_Implicit Hypersonics In-house Solver_

---

**IHIS** is an implicit, density-based compressible flow solver for **external hypersonic aerodynamics** under **thermal and chemical non-equilibrium**. It is built on **HiSA** (High Speed Aerodynamic solver, CSIR) for OpenFOAM-v2312, with HiSA's OpenFOAM thermophysical library replaced by the in-house **neTCLib** (non-equilibrium ThermoChemical Library). IHIS is the implicit companion of the explicit solver **nexaFoam**: both link against the same neTCLib, so physical models are shared and the two solvers differ only in numerics.

Key capabilities:
- Two-temperature (**TTR** for translational–rotational, **Tᵥ** for vibrational) model
- **Finite-rate reacting flows** of neutral molecular gases
- **Implicit time integration**, steady (local pseudo-time) and unsteady (dual-time BDF2), with no explicit stability limit from stiff chemistry or vibrational relaxation

---

- **Physics**
  - Two-temperature formulation without electronic energy and without charged species
  - Per-molecule vibrational energy equations
  - Vibrational–translational relaxation (Millikan–White with Park's high-temperature correction, Landau–Teller form), vibrational–vibrational exchange, and chemistry–vibration coupling (preferential and non-preferential models)
  - Finite-rate chemistry (e.g. Park 5-species air: N₂, O₂, NO, N, O). Heats of formation are carried inside the total energy, so chemistry is energy-neutral by construction
  - Species diffusion by Fick's law with constant Lewis number
  - Temperature-dependent properties (cₚ(T), e(T), μ(T), κ(T), …) provided by **neTCLib**
  - Laminar flow only (HiSA's turbulence modelling is removed)
- **Numerics**
  - HiSA's coupled implicit solution of the mean-flow block (ρ, ρ**u**, ρE) with a matrix-free GMRES solver and an approximate (Lax–Friedrichs) flux Jacobian as the preconditioning operator
  - Upwind flux schemes: AUSM⁺-up, HLLC, Lax–Friedrichs
  - Limited reconstruction (Minmod)
  - Species and vibrational energies solved as segregated implicit scalars after the mean-flow block, with point-implicit source terms built from neTCLib's analytic chemistry Jacobian
  - Characteristic (Riemann-invariant) far-field, inlet/outlet and wall boundary conditions
  - Smoluchowski temperature-jump and Maxwell velocity-slip wall conditions for rarefied regimes
  - Dual-time stepping (BDF2 in physical time, local pseudo-time with an adaptive pseudo-CFL) for time-accurate unsteady runs
- **Use cases**
  - Blunt and streamlined hypersonic bodies: cylinders, spheres, blunted cones, wedges and ramps
  - Shock–shock and shock–boundary-layer interactions (e.g. Edney Type IV), in both steady and time-accurate modes
  - Shock layers with strong vibrational excitation and dissociation
  - High-enthalpy flows where chemistry matters but **plasma effects are out of scope**

---

```
IHIS/
├─ src/
│  ├─ modules/
│  │  └─ hisa/                 # solver module: coupled outer iteration, field update, residuals
│  └─ finiteVolume/
│     ├─ fluxSchemes/          # AUSM+up, HLLC, Lax-Friedrichs
│     ├─ jacobians/            # inviscid and viscous flux Jacobians
│     └─ boundaryConditions/   # characteristic BCs
├─ tutorials/                  # example/tutorial cases
├─ tests/                      # verification scripts
├─ README.md
```

neTCLib and its interface layer (`TCLibraryBase`) are shared with nexaFoam and linked as external libraries.

---

## Requirements

- **OpenFOAM-v2312** (ESI). The solver relies on headers and ABIs of this release.
- **neTCLib**, compiled against the same OpenFOAM release
- **C++17** toolchain compatible with v2312
- Linux environment

---

## Development status

| Phase | Content | Status |
|---|---|---|
| 1 | Turbulence removal; laminar regression against stock HiSA | ✔ complete |
| 2 | Thermodynamic library swap (`fluidThermo` → neTCLib), 1T frozen N₂ | ✔ complete |
| 3 | Species and vibrational energy equations, finite-rate chemistry, 2T relaxation | in progress |
| 4 | Validation: 0-D heat bath → 1-D normal shock → blunted cone → Edney Type IV | planned |

---

## Limitations

- No electrons or ions: no plasma transport or conduction, no electron energy equation
- Reaction mechanisms restricted to **neutral** species
- No electronic energy mode
- Laminar only; no turbulence or transition modelling
- Non-catalytic walls only; no catalytic or ablating wall conditions
- Radiation not modelled
- HiSA's conjugate-heat-transfer (`thermalSolid`) and modal fluid–structure modules are not supported with neTCLib

---

## Acknowledgements

- Built on **HiSA** (J. A. Heyns, O. F. Oxtoby, CSIR, South Africa), released under the GNU GPL v3.
- Built on **OpenFOAM-v2312** (OpenCFD Ltd.).
- neTCLib: in-house non-equilibrium thermochemical library for molecular gases, shared with nexaFoam.

**Reference:** Heyns, J. A., Oxtoby, O. F., & Steenkamp, A. (2014). *Modelling high-speed flow using a matrix-free coupled solver.* 9th OpenFOAM Workshop, Zagreb, Croatia.
