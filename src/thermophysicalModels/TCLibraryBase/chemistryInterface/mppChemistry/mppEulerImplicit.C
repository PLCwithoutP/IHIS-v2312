/*---------------------------------------------------------------------------*\
    mppEulerImplicit.C

    Implementation of mppEulerImplicit.
    See mppEulerImplicit.H for full design notes.
\*---------------------------------------------------------------------------*/

#ifdef WITH_MUTATION_PP

#include "mppEulerImplicit.H"
#include "mutationppWrapper.H"     // full definition needed for field access
#include "mutationppAdapter.H"
#include "fvMesh.H"

namespace Foam
{

// ─────────────────────────────────────────────────────────────────────────────
//  Constructor
// ─────────────────────────────────────────────────────────────────────────────

mppEulerImplicit::mppEulerImplicit
(
    mutationppWrapper& wrapper,
    const dictionary&  dict
)
:
    wrapper_ (wrapper),
    adapter_ (wrapper.adapter()),   // adapter() accessor — see note below
    Qdot_
    (
        IOobject
        (
            "Qdot",
            wrapper.mesh().time().timeName(),
            wrapper.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        wrapper.mesh(),
        dimensionedScalar("Qdot", dimPower/dimVolume, 0.0)
    ),
    cTauChem_(dict.getOrDefault<scalar>("cTauChem", 0.1))
{
    initRR_();
}


// ─────────────────────────────────────────────────────────────────────────────
//  initRR_
// ─────────────────────────────────────────────────────────────────────────────

void mppEulerImplicit::initRR_()
{
    const int ns = adapter_.nSpecies();
    const fvMesh& mesh = wrapper_.mesh();

    RR_.resize(ns);

    for (int i = 0; i < ns; i++)
    {
        RR_.set
        (
            i,
            new volScalarField::Internal
            (
                IOobject
                (
                    "RR." + word(adapter_.speciesName(i)),
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar("RR", dimDensity/dimTime, 0.0)
            )
        );
    }
}


// ─────────────────────────────────────────────────────────────────────────────
//  implicitStep_
// ─────────────────────────────────────────────────────────────────────────────

scalar mppEulerImplicit::implicitStep_
(
    std::vector<double>& rhoi,
    const std::vector<double>& wdot,
    scalar dt_sub
) const
{
    // Diagonal EulerImplicit per species:
    //   rhoi^(n+1) = (rhoi^n + dt * wdot^+) / (1 + dt * wdot^- / rhoi^n)
    //
    // Also compute the next stable sub-step estimate:
    //   dt_next = cTauChem * min_i (rhoi / |wdot_i|)

    const int ns = adapter_.nSpecies();
    const double eps = 1.0e-20;  // floor for density

    scalar dt_next = GREAT;

    for (int i = 0; i < ns; i++)
    {
        const double src_pos = max(wdot[i],  0.0);  // production
        const double src_neg = max(-wdot[i], 0.0);  // destruction

        // Stable timescale estimate for this species
        if (src_neg > eps)
        {
            const double tau = rhoi[i] / src_neg;
            dt_next = min(dt_next, cTauChem_ * tau);
        }
        if (src_pos > eps)
        {
            // Production rate timescale: how fast would rho change?
            // Use a representative mass density (max of current and 1e-10)
            const double tau = max(rhoi[i], 1.0e-10) / src_pos;
            dt_next = min(dt_next, cTauChem_ * tau);
        }

        // Implicit step
        const double denom = 1.0 + dt_sub * src_neg / max(rhoi[i], eps);
        rhoi[i] = (rhoi[i] + dt_sub * src_pos) / denom;
        rhoi[i] = max(rhoi[i], 0.0);
    }

    // Return clamped next sub-step (never larger than dt_sub)
    return min(dt_next, dt_sub);
}


// ─────────────────────────────────────────────────────────────────────────────
//  solve
// ─────────────────────────────────────────────────────────────────────────────

scalar mppEulerImplicit::solve(scalar deltaT)
{
    if (!adapter_.hasChemistry())
    {
        // No mechanism loaded — zero all RR fields and return
        forAll(RR_, i)
            RR_[i].field() = 0.0;
        Qdot_.primitiveFieldRef() = 0.0;
        return GREAT;
    }

    const int ns = adapter_.nSpecies();

    // Working arrays — allocated once per solve() call (not per cell)
    std::vector<double> Yi(ns);
    std::vector<double> rhoi(ns);
    std::vector<double> rhoi0(ns);
    std::vector<double> wdot(ns);

    // Access state fields from wrapper
    const scalarField& TtrCells  = wrapper_.TTR().internalField();
    const scalarField& TVibCells = wrapper_.TVib().internalField();
    const scalarField& rhoCells  = wrapper_.rho()().internalField();
    const PtrList<volScalarField>& Y = wrapper_.composition().Y();

    // Output references
    scalarField& QdotCells = Qdot_.primitiveFieldRef();

    // Zero RR_ and Qdot_
    forAll(RR_, i)
        RR_[i].field() = 0.0;
    QdotCells = 0.0;

    scalar tauMin = GREAT;   // minimum chemical timescale across all cells

    forAll(TtrCells, cellI)
    {
        const scalar Ttr = TtrCells[cellI];
        const scalar Tv  = TVibCells[cellI];
        const scalar rho = rhoCells[cellI];

        if (rho < SMALL) continue;

        // ── Extract initial state ─────────────────────────────────────────

        for (int i = 0; i < ns; i++)
        {
            Yi[i]    = Y[i][cellI];
            rhoi[i]  = rho * Yi[i];
            rhoi0[i] = rhoi[i];
        }

        // ── Sub-step through deltaT ───────────────────────────────────────

        scalar timeLeft = deltaT;
        scalar dt_sub   = deltaT;   // first guess: full step

        while (timeLeft > SMALL * deltaT)
        {
            // Clamp sub-step to remaining time
            dt_sub = min(dt_sub, timeLeft);

            // Compute current mass fractions from rhoi
            double rhoTot = 0.0;
            for (int i = 0; i < ns; i++) rhoTot += rhoi[i];
            for (int i = 0; i < ns; i++) Yi[i] = rhoi[i] / max(rhoTot, 1.0e-20);

            // Set Mutation++ state and get production rates
            adapter_.setState(Yi.data(), rhoTot, Ttr, Tv);
            adapter_.netProductionRates(wdot.data());

            // Implicit step — updates rhoi, returns next stable dt_sub
            dt_sub = implicitStep_(rhoi, wdot, dt_sub);
            dt_sub = max(dt_sub, SMALL * deltaT);

            timeLeft -= min(dt_sub, timeLeft);

            // Track minimum timescale
            for (int i = 0; i < ns; i++)
            {
                if (std::abs(wdot[i]) > 1.0e-20 && rhoi[i] > 1.0e-20)
                    tauMin = min(tauMin, static_cast<scalar>(rhoi[i] / std::abs(wdot[i])));
            }
        }

        // ── Compute time-averaged production rate ─────────────────────────
        //  RR_[i][cellI] = (rhoi_final - rhoi_initial) / deltaT

        // Update Yi from final rhoi for adapter state (needed for Hc)
        double rhoFinal = 0.0;
        for (int i = 0; i < ns; i++) rhoFinal += rhoi[i];
        for (int i = 0; i < ns; i++) Yi[i] = rhoi[i] / max(rhoFinal, 1.0e-20);

        adapter_.setState(Yi.data(), rhoFinal, Ttr, Tv);

        for (int i = 0; i < ns; i++)
            RR_[i][cellI] = (rhoi[i] - rhoi0[i]) / deltaT;

        // ── Heat release rate ─────────────────────────────────────────────
        //  Qdot = -SUM_i Hc_i * RR_i    [J/m3-s]
        //  adapter_ is now set to final cell state → Hc(i) is valid.

        scalar qdot = 0.0;
        for (int i = 0; i < ns; i++)
            qdot -= adapter_.Hc(i) * RR_[i][cellI];

        QdotCells[cellI] = qdot;
    }

    Qdot_.correctBoundaryConditions();

    return tauMin;
}


// ─────────────────────────────────────────────────────────────────────────────
//  Result accessors
// ─────────────────────────────────────────────────────────────────────────────

const volScalarField::Internal& mppEulerImplicit::RR(label i) const
{
    return RR_[i];
}

const PtrList<volScalarField::Internal>& mppEulerImplicit::RR() const
{
    return RR_;
}

tmp<volScalarField> mppEulerImplicit::Qdot() const
{
    return tmp<volScalarField>(Qdot_);
}

} // End namespace Foam

#endif // WITH_MUTATION_PP
