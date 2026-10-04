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
      calculateWithJacobianDiag() — RR plus Jacobian diagonal at given rho (IHIS)


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

    //- Rates (into RR) and the diagonal of the chemical Jacobian at a
    //  caller-supplied density. No ODE time advance. Used by IHIS for
    //  point-implicit species sources (Phase 4, D6).
    //    dRRdRhoY[i] = d(omega_i)/d(rho Y_i)  [1/s]
    //    dRRdT[i]    = d(omega_i)/dT          [kg/m3/s/K]
    //  Lists are allocated on the first call and overwritten afterwards.
    virtual void calculateWithJacobianDiag
    (
        const volScalarField::Internal& rho,
        PtrList<volScalarField::Internal>& dRRdRhoY,
        PtrList<volScalarField::Internal>& dRRdT
    ) = 0;
};

} // End namespace Foam

#endif
