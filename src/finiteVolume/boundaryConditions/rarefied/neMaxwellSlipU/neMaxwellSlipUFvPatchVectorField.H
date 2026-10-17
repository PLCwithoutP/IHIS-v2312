/*---------------------------------------------------------------------------*\
  IHIS - Implicit Hypersonics In-house Solver (OpenFOAM-v2312)
  Derived from rhoCentralFoam/BCs/U/maxwellSlipUFvPatchVectorField
  (Copyright (C) 2011-2016 OpenFOAM Foundation, 2020 OpenCFD Ltd., GPL v3)
-------------------------------------------------------------------------------
Class
    Foam::neMaxwellSlipUFvPatchVectorField

Description
    Maxwell (1879) first-order velocity slip for IHIS / neTCLib (Phase 7).

        U_s = f U_ref + (1 - f) U_c,     f = 1/(1 + deltaCoeffs lambda (2-s)/s)

        U_ref = U_w
              [+ 3/4 nu/T (I - nn).grad(T)]              thermalCreep
              [- lambda (2-s)/s / mu (I - nn).(n.tauMC)] curvature

    lambda: VHS mean free path (neRarefiedWall.H), from the thermo patch
    mu, p and psi. Thermal creep uses the physical sign (gas creeps toward
    higher T); stock OpenFOAM and hy2Foam subtract it. Both optional terms
    are off by default (Phase 7 D7, as the hy2Foam blunted-cone tutorial).

Usage
    \verbatim
    wall
    {
        type                neMaxwellSlipU;
        accommodationCoeff  1;          // tangential momentum, 0 < s <= 1
        omega               0.74;       // VHS exponent, required (N2)
        Uwall               uniform (0 0 0);    // optional, default zero
        thermalCreep        false;      // optional, default false
        curvature           false;      // optional, default false
        T                   TTR;        // optional
        tauMC               tauMC;      // optional
        value               uniform (0 0 0);
    }
    \endverbatim

SourceFiles
    neMaxwellSlipUFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef neMaxwellSlipUFvPatchVectorField_H
#define neMaxwellSlipUFvPatchVectorField_H

#include "partialSlipFvPatchFields.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class neMaxwellSlipUFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class neMaxwellSlipUFvPatchVectorField
:
    public partialSlipFvPatchVectorField
{
    // Private Data

        //- Trans-rotational temperature field name (thermal creep), "TTR"
        word TName_;

        //- tauMC field name (curvature term), "tauMC"
        word tauMCName_;

        //- Tangential momentum accommodation coefficient
        scalar accommodationCoeff_;

        //- VHS viscosity exponent
        scalar omega_;

        //- Wall velocity
        vectorField Uwall_;

        //- Include the thermal creep term (default off)
        Switch thermalCreep_;

        //- Include the boundary curvature term (default off)
        Switch curvature_;


public:

    //- Runtime type information
    TypeName("neMaxwellSlipU");


    // Constructors

        //- Construct from patch and internal field
        neMaxwellSlipUFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        neMaxwellSlipUFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping onto a new patch
        neMaxwellSlipUFvPatchVectorField
        (
            const neMaxwellSlipUFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new neMaxwellSlipUFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        neMaxwellSlipUFvPatchVectorField
        (
            const neMaxwellSlipUFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new neMaxwellSlipUFvPatchVectorField(*this, iF)
            );
        }


    // Member Functions

        // Mapping

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap(const fvPatchFieldMapper&);

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap(const fvPatchVectorField&, const labelList&);


        // Evaluation

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
