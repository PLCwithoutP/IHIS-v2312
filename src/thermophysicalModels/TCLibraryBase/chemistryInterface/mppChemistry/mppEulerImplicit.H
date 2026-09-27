/*---------------------------------------------------------------------------*\
    mppEulerImplicit.H

    EulerImplicit chemistry ODE solver for the Mutation++ backend.

    Analogue of neEulerImplicit, but:
      - Not templated — operates directly on mutationppWrapper fields.
      - Uses adapter_.netProductionRates(wdot) instead of per-reaction
        omegaI() decomposition.  Mutation++ handles rate temperature
        selection (T_Park = sqrt(T_tr * T_v)) internally.
      - Uses a diagonal implicit scheme (no full ns×ns Jacobian) because
        Mutation++ does not expose a public Jacobian API.
      - Does NOT update temperatures — operator splitting with the energy
        equation handles T_tr and T_v in the main solver loop.

    Diagonal EulerImplicit scheme (per species, per sub-step):
        rhoi^(n+1) = (rhoi^n + dt * wdot_i^+)
                     / (1 + dt * wdot_i^- / max(rhoi^n, epsilon))

    where wdot_i^+ = max(wdot_i, 0)  (production)
          wdot_i^- = max(-wdot_i, 0) (destruction)

    Sub-step size is controlled by cTauChem * (rhoi / |wdot_i|) to keep
    the correction per step below a user-specified fraction of rhoi.

    RR_[i] = time-averaged rate = (rhoi_final - rhoi_initial) / deltaT
    This accounts for stiffness without requiring a small global timestep.

    Qdot is computed per cell during the solve loop as:
        Qdot = -SUM_i Hc_i * RR_i    [J/m3-s]
    where Hc_i is the formation enthalpy from adapter_.Hc(i).
\*---------------------------------------------------------------------------*/

#ifndef mppEulerImplicit_H
#define mppEulerImplicit_H

#ifdef WITH_MUTATION_PP

#include "volFields.H"
#include "PtrList.H"
#include "dictionary.H"

namespace Foam
{

// Forward declarations
class mutationppWrapper;
class mutationppAdapter;

class mppEulerImplicit
{
public:

    // -------------------------------------------------------------------------
    // Constructor
    // -------------------------------------------------------------------------

    //- Construct from wrapper and solver settings dictionary.
    //  dict should be a subdictionary of mppDict containing at least:
    //    cTauChem   scalar  (fraction of characteristic timescale per sub-step)
    mppEulerImplicit
    (
        mutationppWrapper&  wrapper,
        const dictionary&   dict
    );

    ~mppEulerImplicit() = default;


    // -------------------------------------------------------------------------
    // Main solver
    // -------------------------------------------------------------------------

    //- Advance chemistry by deltaT over all cells.
    //  Sub-steps per cell to control stiffness.
    //  Returns the minimum characteristic chemical timescale [s].
    scalar solve(scalar deltaT);


    // -------------------------------------------------------------------------
    // Result accessors  (valid after solve() has been called)
    // -------------------------------------------------------------------------

    //- Time-averaged species mass production rate for species i [kg/m3-s]
    const volScalarField::Internal& RR(label i) const;

    //- All species production rate fields
    const PtrList<volScalarField::Internal>& RR() const;

    //- Heat release rate [J/m3-s]
    //  Computed during solve() as -SUM_i Hc_i * RR_i
    tmp<volScalarField> Qdot() const;


private:

    // -------------------------------------------------------------------------
    // References to wrapper state (not owned here)
    // -------------------------------------------------------------------------

    mutationppWrapper&  wrapper_;
    mutationppAdapter&  adapter_;   // alias: wrapper_.adapter_


    // -------------------------------------------------------------------------
    // Result fields
    // -------------------------------------------------------------------------

    //- Species production rates [kg/m3-s]
    PtrList<volScalarField::Internal> RR_;

    //- Heat release rate [J/m3-s]
    volScalarField Qdot_;


    // -------------------------------------------------------------------------
    // Solver settings
    // -------------------------------------------------------------------------

    //- Fraction of (rhoi / |wdot_i|) used as sub-step size.
    //  Smaller = more accurate, more expensive.  Typical: 0.1 – 0.5.
    scalar cTauChem_;


    // -------------------------------------------------------------------------
    // Private helpers
    // -------------------------------------------------------------------------

    //- Initialise RR_ fields (zero, dimensioned correctly)
    void initRR_();

    //- Perform one diagonal EulerImplicit sub-step for a single cell.
    //  Updates rhoi in-place.  Returns proposed next sub-step size.
    scalar implicitStep_
    (
        std::vector<double>& rhoi,
        const std::vector<double>& wdot,
        scalar dt_sub
    ) const;
};

} // End namespace Foam

#endif // WITH_MUTATION_PP
#endif // mppEulerImplicit_H
