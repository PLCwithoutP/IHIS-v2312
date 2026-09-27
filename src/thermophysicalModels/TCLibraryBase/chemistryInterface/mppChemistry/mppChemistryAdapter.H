/*---------------------------------------------------------------------------*\
    mppChemistryAdapter.H

    Concrete implementation of chemistryInterface backed by mppEulerImplicit.

    Analogue of neTCChemistryAdapter for the Mutation++ backend.

    This header is included ONLY from chemistryInterface.C — it is an
    implementation detail, not part of the public interface.

    The adapter:
      - Owns mppEulerImplicit (the ODE solver)
      - Reads solver settings from constant/mppDict / mppEulerImplicitCoeffs
      - Delegates all chemistryInterface calls to the solver

    Unlike neTCChemistryAdapter which wraps a templated model,
    mppChemistryAdapter is non-templated because Mutation++ is not
    part of the OpenFOAM chemistry model template hierarchy.
\*---------------------------------------------------------------------------*/

#ifndef mppChemistryAdapter_H
#define mppChemistryAdapter_H

#ifdef WITH_MUTATION_PP

#include "chemistryInterface.H"
#include "mppEulerImplicit.H"
#include "mutationppWrapper.H"
#include "IOdictionary.H"

namespace Foam
{

class mppChemistryAdapter
:
    public chemistryInterface
{
    // -------------------------------------------------------------------------
    // Owned solver
    // -------------------------------------------------------------------------

    mppEulerImplicit solver_;


public:

    // -------------------------------------------------------------------------
    // Constructor
    // -------------------------------------------------------------------------

    //- Construct from Mutation++ wrapper.
    //  Reads mppEulerImplicitCoeffs from constant/mppDict.
    explicit mppChemistryAdapter(mutationppWrapper& wrapper)
    :
        solver_
        (
            wrapper,
            // Read solver settings from mppDict if present,
            // or use an empty dictionary (defaults apply).
            wrapper.mppDict().found("mppEulerImplicitCoeffs")
            ? wrapper.mppDict().subDict("mppEulerImplicitCoeffs")
            : dictionary::null
        )
    {}


    // -------------------------------------------------------------------------
    // chemistryInterface implementation
    // -------------------------------------------------------------------------

    scalar solve(scalar deltaT) override
    {
        return solver_.solve(deltaT);
    }

    const volScalarField::Internal& RR(label i) const override
    {
        return solver_.RR(i);
    }

    const PtrList<volScalarField::Internal>& RR() const override
    {
        return solver_.RR();
    }

    tmp<volScalarField> Qdot() const override
    {
        return solver_.Qdot();
    }
};

} // End namespace Foam

#endif // WITH_MUTATION_PP
#endif // mppChemistryAdapter_H
