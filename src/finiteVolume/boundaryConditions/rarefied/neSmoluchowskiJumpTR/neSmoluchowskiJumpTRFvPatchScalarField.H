/*---------------------------------------------------------------------------*\
  IHIS - Implicit Hypersonics In-house Solver (OpenFOAM-v2312)
  Derived from rhoCentralFoam/BCs/T/smoluchowskiJumpTFvPatchScalarField
  (Copyright (C) 2011-2015 OpenFOAM Foundation, 2020 OpenCFD Ltd., GPL v3)
-------------------------------------------------------------------------------
Class
    Foam::neSmoluchowskiJumpTRFvPatchScalarField

Description
    Smoluchowski (1898) first-order temperature jump of the trans-rotational
    temperature TTR, for IHIS / neTCLib (Phase 7).

        T_s = f T_w + (1 - f) T_c,     f = 1/(1 + deltaCoeffs C_T)

        C_T = lambda 2 gamma/(gamma + 1) kappaTR/(mu cpTR) (2 - a)/a

    gamma = 1 + R/cvTR, cpTR = cvTR + R, R = 1/(psi T): frozen TR values from
    the thermo patch state (no Pr or gamma entry: both are local). kappaTR is
    the Wilke-mixed conductivity IHIS conducts with (kappaTRWilke).
    lambda: VHS mean free path (neRarefiedWall.H).

    Being a mixed patch field, the eTR patch becomes mixed2TEnergy
    automatically (basic2TThermo::hBoundaryTypes).

Usage
    \verbatim
    wall
    {
        type                neSmoluchowskiJumpTR;
        accommodationCoeff  1;          // thermal, 0 < a <= 1
        omega               0.74;       // VHS exponent, required (N2)
        Twall               uniform 297.2;
        value               uniform 297.2;
    }
    \endverbatim

SourceFiles
    neSmoluchowskiJumpTRFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef neSmoluchowskiJumpTRFvPatchScalarField_H
#define neSmoluchowskiJumpTRFvPatchScalarField_H

#include "mixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class neSmoluchowskiJumpTRFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class neSmoluchowskiJumpTRFvPatchScalarField
:
    public mixedFvPatchScalarField
{
    // Private Data

        //- Thermal accommodation coefficient
        scalar accommodationCoeff_;

        //- VHS viscosity exponent
        scalar omega_;

        //- Wall surface temperature
        scalarField Twall_;


public:

    //- Runtime type information
    TypeName("neSmoluchowskiJumpTR");


    // Constructors

        //- Construct from patch and internal field
        neSmoluchowskiJumpTRFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        neSmoluchowskiJumpTRFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping onto a new patch
        neSmoluchowskiJumpTRFvPatchScalarField
        (
            const neSmoluchowskiJumpTRFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new neSmoluchowskiJumpTRFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        neSmoluchowskiJumpTRFvPatchScalarField
        (
            const neSmoluchowskiJumpTRFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new neSmoluchowskiJumpTRFvPatchScalarField(*this, iF)
            );
        }


    // Member Functions

        // Mapping

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap(const fvPatchFieldMapper&);

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap(const fvPatchScalarField&, const labelList&);


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
