/*---------------------------------------------------------------------------*\
    chemistryInterface.H

    Abstract base class for all chemistry backends.
    Derived by:
      neTCChemistryAdapter  — wraps BasicNeChemistryModel<ne2TReactionThermo>
      mutationPPChemistry   — future Mutation++ chemistry backend

    Exposes only the solver-facing API:
      solve()   — advance ODE one timestep
      RR(i)     — species source term field i
      RR()      — all species source term fields
      Qdot()    — heat release rate field

    Factory reads the same 'library' keyword used by tcLibraryInterface.
    The dynamic_cast from tcLibraryInterface to ne2TReactionThermo lives
    exclusively inside the factory — the solver never sees it.
\*---------------------------------------------------------------------------*/

#ifndef chemistryInterface_H
#define chemistryInterface_H

#include "volFields.H"
#include "PtrList.H"
#include "autoPtr.H"

namespace Foam
{

class tcLibraryInterface;

class chemistryInterface
{
public:

    virtual ~chemistryInterface() {}

    // -----------------------------------------------------------------------
    // Factory
    // -----------------------------------------------------------------------

    //- Construct from the abstract thermo interface.
    //  Internally casts to the concrete backend type.
    static autoPtr<chemistryInterface> New(tcLibraryInterface& thermo);


    // -----------------------------------------------------------------------
    // Solver-facing API
    // -----------------------------------------------------------------------

    //- Advance chemistry ODE by deltaT, return characteristic timescale
    virtual scalar solve(scalar deltaT) = 0;

    //- Species mass source term field i [kg/m3/s]
    virtual const volScalarField::Internal& RR(label i) const = 0;

    //- All species mass source term fields
    virtual const PtrList<volScalarField::Internal>& RR() const = 0;

    //- Heat release rate [kg/m/s3]
    virtual tmp<volScalarField> Qdot() const = 0;
};

} // End namespace Foam

#endif
