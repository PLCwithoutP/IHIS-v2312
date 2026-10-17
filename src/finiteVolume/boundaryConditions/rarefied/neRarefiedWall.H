/*---------------------------------------------------------------------------*\
  IHIS - Implicit Hypersonics In-house Solver (OpenFOAM-v2312)
-------------------------------------------------------------------------------
Namespace
    Foam::neRarefiedWall

Description
    Shared helpers for the IHIS slip-regime wall boundary conditions
    (Phase 7, D9): neMaxwellSlipU, neSmoluchowskiJumpTR, neSmoluchowskiJumpTv.

    All functions work on patch (face) values; the BCs evaluate them from the
    current patch state, so the coefficients are lagged by one evaluation.

    Wall state. At a characteristicWallPressure wall IHIS sets the rho patch
    to the face-cell value (updateFields.H), so p/(rho T) is NOT the gas
    constant there. Everything below therefore uses the thermo's own patch
    p and psi = 1/(R T), which are consistent with the patch temperature:
        rho_w = p psi,   R = 1/(psi T),   R T = 1/psi.

    Mean free path: variable hard sphere (VHS), Bird (1994):
        lambda = 2(7 - 2w)(5 - 2w)/15 * mu/(rho sqrt(2 pi R T))
               = 2(7 - 2w)(5 - 2w)/15 * mu/(p sqrt(2 pi psi))
    N2, w = 0.74: prefactor 1.0336 (stock OpenFOAM: sqrt(pi/2) = 1.2533).

SourceFiles
    header only

\*---------------------------------------------------------------------------*/

#ifndef neRarefiedWall_H
#define neRarefiedWall_H

#include "scalarField.H"
#include "dictionary.H"
#include "objectRegistry.H"
#include "volFields.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace neRarefiedWall
{

//- VHS mean free path [m] from patch mu [Pa s], p [Pa], psi [s2/m2]
inline tmp<scalarField> meanFreePathVHS
(
    const scalarField& mu,
    const scalarField& p,
    const scalarField& psi,
    const scalar omega
)
{
    const scalar C = 2.0*(7.0 - 2.0*omega)*(5.0 - 2.0*omega)/15.0;

    return C*mu/(p*sqrt(constant::mathematical::twoPi*psi));
}


//- Mixture gas constant R = 1/(psi T) [J/(kg K)]
inline tmp<scalarField> gasConstant
(
    const scalarField& psi,
    const scalarField& T
)
{
    return 1.0/(psi*T);
}


//- Frozen trans-rotational ratio of specific heats, gamma = 1 + R/cvTR
inline tmp<scalarField> gammaTR
(
    const scalarField& cvTR,
    const scalarField& R
)
{
    return 1.0 + R/cvTR;
}


//- Read an accommodation coefficient (required), check 0 < a <= 1
inline scalar readAccommodation(const dictionary& dict, const word& key)
{
    const scalar a = dict.get<scalar>(key);

    if (a <= 0 || a > 1)
    {
        FatalIOErrorInFunction(dict)
            << "Unphysical " << key << " = " << a
            << " (required: 0 < " << key << " <= 1)"
            << exit(FatalIOError);
    }

    return a;
}


//- (2 - a)/a
inline scalar accommodationFactor(const scalar a)
{
    return (2.0 - a)/a;
}


//- Read the VHS viscosity exponent (required, D8), check 0.5 <= omega <= 1
inline scalar readOmega(const dictionary& dict)
{
    const scalar w = dict.get<scalar>("omega");

    if (w < 0.5 || w > 1.0)
    {
        FatalIOErrorInFunction(dict)
            << "VHS omega = " << w << " outside [0.5, 1.0]"
            << " (N2: 0.74, O2: 0.77, Bird 1994)"
            << exit(FatalIOError);
    }

    return w;
}


//- True once the coefficients can be evaluated safely.
//  The thermo constructs (and may evaluate) TTR before it is complete, and
//  mixed2TEnergy::updateCoeffs() calls TTR.evaluate(). foundObject<ThermoType>
//  is false while the thermo is still a base-class object under construction;
//  rhoE is created by IHIS only after the thermo exists. Until then the BCs
//  keep their initial coefficients (valueFraction 1: U_b = U_w, T_b = T_w).
template<class ThermoType>
inline bool ready(const objectRegistry& db, const word& thermoName)
{
    return
        db.foundObject<ThermoType>(thermoName)
     && db.foundObject<volScalarField>("rhoE");
}


} // End namespace neRarefiedWall
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
