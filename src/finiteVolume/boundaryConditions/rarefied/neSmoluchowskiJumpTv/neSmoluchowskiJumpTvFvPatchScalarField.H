/*---------------------------------------------------------------------------*\
  IHIS - Implicit Hypersonics In-house Solver (OpenFOAM-v2312)
  Derived from rhoCentralFoam/BCs/T/smoluchowskiJumpTFvPatchScalarField
  (Copyright (C) 2011-2015 OpenFOAM Foundation, 2020 OpenCFD Ltd., GPL v3)
-------------------------------------------------------------------------------
Class
    Foam::neSmoluchowskiJumpTvFvPatchScalarField

Description
    Smoluchowski first-order jump of one molecule's vibrational temperature
    T_v,m, for the IHIS per-molecule fields TVib_<m> (Phase 7, D3).

        T_v,s = f T_w + (1 - f) T_v,c,     f = 1/(1 + deltaCoeffs C_v)

        C_v = lambda 2 gamma/(gamma + 1) Gamma_v,m/mu (2 - a)/a

    The hy2Foam form (Casseau et al. 2016; hy2Foam nonEqSmoluchowskiJumpTvMix)
    with the trans-rotational kappa/cp replaced by the vibrational
    conductance on e_v. Gamma_v,m = X_m G_m is the field GammaVib_<m>
    registered by IHIS (residualsUpdate.H), i.e. exactly the coefficient the
    vibrational heat flux uses. For pure N2, Gamma_v = kappa_v/cv_v
    (hy2Foam's alphave). gamma is the frozen TR value at the TTR patch.
    (hy2Foam's per-species class nonEqSmoluchowskiJumpTv has alpha in place
    of gamma in the numerator, which is dimensionally inconsistent.)

    The molecule is taken from the field name: TVib_N2 -> N2.
    Until GammaVib_<m> exists (first outer iteration) the coefficients are
    kept (initially valueFraction 1, T_v,b = T_w; Phase 7 D13).

Usage
    \verbatim
    wall
    {
        type                neSmoluchowskiJumpTv;
        accommodationCoeff  1;          // 0 < a <= 1
        omega               0.74;       // VHS exponent, required (N2)
        Twall               uniform 297.2;
        T                   TTR;        // optional
        value               uniform 297.2;
    }
    \endverbatim

SourceFiles
    neSmoluchowskiJumpTvFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef neSmoluchowskiJumpTvFvPatchScalarField_H
#define neSmoluchowskiJumpTvFvPatchScalarField_H

#include "mixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class neSmoluchowskiJumpTvFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class neSmoluchowskiJumpTvFvPatchScalarField
:
    public mixedFvPatchScalarField
{
    // Private Data

        //- Trans-rotational temperature field name, "TTR"
        word TName_;

        //- Molecule name (from the field name TVib_<m>)
        word specie_;

        //- Accommodation coefficient
        scalar accommodationCoeff_;

        //- VHS viscosity exponent
        scalar omega_;

        //- Wall surface temperature
        scalarField Twall_;


    // Private Member Functions

        //- Molecule name from the internal field name TVib_<m>
        static word specieName(const word& fieldName);


public:

    //- Runtime type information
    TypeName("neSmoluchowskiJumpTv");


    // Constructors

        //- Construct from patch and internal field
        neSmoluchowskiJumpTvFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        neSmoluchowskiJumpTvFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping onto a new patch
        neSmoluchowskiJumpTvFvPatchScalarField
        (
            const neSmoluchowskiJumpTvFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new neSmoluchowskiJumpTvFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        neSmoluchowskiJumpTvFvPatchScalarField
        (
            const neSmoluchowskiJumpTvFvPatchScalarField&,
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
                new neSmoluchowskiJumpTvFvPatchScalarField(*this, iF)
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
