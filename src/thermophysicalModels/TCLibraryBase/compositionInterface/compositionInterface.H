/*---------------------------------------------------------------------------*\
    compositionInterface.H

    Abstract base for species composition access.
    Location: TCLibraryBase/compositionInterface/

    Derived by:
      basicSpecie2TMixture  — multi-component path (all methods implemented)
      basic2TMixture        — pure mixture path (FatalError stubs for all
                              species-specific methods since mixtureCheck==false)

    Contains the complete per-species interface matching basicSpecie2TMixture.
    All methods added upfront to avoid incremental discovery during compilation.
\*---------------------------------------------------------------------------*/

#ifndef compositionInterface_H
#define compositionInterface_H

#include "volFields.H"
#include "PtrList.H"
#include "speciesTable.H"

namespace Foam
{

class compositionInterface
{
public:

    virtual ~compositionInterface() {}

    // -----------------------------------------------------------------------
    // Species mass fraction fields
    // -----------------------------------------------------------------------

    virtual PtrList<volScalarField>& Y() = 0;
    virtual const PtrList<volScalarField>& Y() const = 0;

    // -----------------------------------------------------------------------
    // Species table
    // -----------------------------------------------------------------------

    virtual const speciesTable& species() const = 0;

    // -----------------------------------------------------------------------
    // Per-species scalar properties
    // Only called when mixtureCheck == true (multi-component path).
    // basic2TMixture provides FatalError stubs for all of these.
    // -----------------------------------------------------------------------

    virtual scalar W(label speciei) const = 0;

    virtual scalar R(label speciei) const = 0;

    virtual scalar Hc(label speciei) const = 0;

    virtual scalar CpTR
    (
        label  speciei,
        scalar p,
        scalar TTR
    ) const = 0;

    virtual scalar CvT
    (
        label  speciei,
        scalar p,
        scalar TTR
    ) const = 0;

    virtual scalar H
    (
        label  speciei,
        scalar p,
        scalar TTR
    ) const = 0;

    virtual scalar Ha
    (
        label  speciei,
        scalar p,
        scalar TTR
    ) const = 0;

    // Called in totalEnergyEquation.H for species diffusive enthalpy flux
    virtual scalar Hs
    (
        label  speciei,
        scalar p,
        scalar TTR
    ) const = 0;

    virtual scalar S
    (
        label  speciei,
        scalar p,
        scalar TTR
    ) const = 0;

    virtual scalar EsT
    (
        label  speciei,
        scalar p,
        scalar TTR
    ) const = 0;

    virtual scalar EsR
    (
        label  speciei,
        scalar p,
        scalar TTR
    ) const = 0;

    // Called in createFields.H
    virtual scalar EsVib
    (
        label  speciei,
        scalar p,
        scalar TTR,
        scalar TVib,
        scalar thetaVib
    ) const = 0;

    virtual scalar G
    (
        label  speciei,
        scalar p,
        scalar TTR
    ) const = 0;

    virtual scalar mu
    (
        label  speciei,
        scalar TTR
    ) const = 0;

    virtual scalar kappaTR
    (
        label  speciei,
        scalar p,
        scalar TTR
    ) const = 0;

    virtual scalar kappaVib
    (
        label  speciei,
        scalar p,
        scalar TTR,
        scalar TVib,
        scalar thetaVib
    ) const = 0;

    virtual scalar alphah
    (
        label  speciei,
        scalar p,
        scalar TTR
    ) const = 0;

    virtual scalar rho
    (
        label  speciei,
        scalar p,
        scalar TTR
    ) const = 0;

    // Called in createFields.H
    virtual scalar thetaVib(label speciei) const = 0;

    // Called in createFields.H and kurganov.H
    virtual bool isSpecieMolecular(label speciei) const = 0;

    // Called in updateVibTemperature.H
    virtual scalar TEVib
    (
        label  speciei,
        scalar esVib,
        scalar p,
        scalar TTR,
        scalar TVib0
    ) const = 0;
};

} // End namespace Foam

#endif
