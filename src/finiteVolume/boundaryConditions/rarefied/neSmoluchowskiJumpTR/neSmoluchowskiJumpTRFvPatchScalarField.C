/*---------------------------------------------------------------------------*\
  IHIS - Implicit Hypersonics In-house Solver (OpenFOAM-v2312)
  Derived from rhoCentralFoam/BCs/T/smoluchowskiJumpTFvPatchScalarField
  (Copyright (C) 2011-2015 OpenFOAM Foundation, 2020 OpenCFD Ltd., GPL v3)
\*---------------------------------------------------------------------------*/

#include "neSmoluchowskiJumpTRFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "fluid2TThermo.H"
#include "neRarefiedWall.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::neSmoluchowskiJumpTRFvPatchScalarField::
neSmoluchowskiJumpTRFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    accommodationCoeff_(1.0),
    omega_(0.74),       // only without a dictionary; the dictionary requires omega
    Twall_(p.size(), Zero)
{
    refValue() = Zero;
    refGrad() = Zero;
    valueFraction() = Zero;
}


Foam::neSmoluchowskiJumpTRFvPatchScalarField::
neSmoluchowskiJumpTRFvPatchScalarField
(
    const neSmoluchowskiJumpTRFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    accommodationCoeff_(ptf.accommodationCoeff_),
    omega_(ptf.omega_),
    Twall_(ptf.Twall_, mapper)          // mapped (stock copies it unmapped)
{}


Foam::neSmoluchowskiJumpTRFvPatchScalarField::
neSmoluchowskiJumpTRFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    accommodationCoeff_
    (
        neRarefiedWall::readAccommodation(dict, "accommodationCoeff")
    ),
    omega_(neRarefiedWall::readOmega(dict)),
    Twall_("Twall", dict, p.size())
{
    // Until the coefficients are first available: T_b = T_w
    refValue() = Twall_;
    refGrad() = Zero;
    valueFraction() = scalar(1);

    // Restart: keep the written valueFraction
    const auto* hasFrac = dict.findEntry("valueFraction", keyType::LITERAL);
    if (hasFrac)
    {
        valueFraction().assign(*hasFrac, p.size());
    }

    if (!this->readValueEntry(dict))
    {
        fvPatchScalarField::operator=(Twall_);
    }
}


Foam::neSmoluchowskiJumpTRFvPatchScalarField::
neSmoluchowskiJumpTRFvPatchScalarField
(
    const neSmoluchowskiJumpTRFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptf, iF),
    accommodationCoeff_(ptf.accommodationCoeff_),
    omega_(ptf.omega_),
    Twall_(ptf.Twall_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::neSmoluchowskiJumpTRFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    mixedFvPatchScalarField::autoMap(m);
    Twall_.autoMap(m);
}


void Foam::neSmoluchowskiJumpTRFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    mixedFvPatchScalarField::rmap(ptf, addr);

    const auto& tptf =
        refCast<const neSmoluchowskiJumpTRFvPatchScalarField>(ptf);
    Twall_.rmap(tptf.Twall_, addr);
}


void Foam::neSmoluchowskiJumpTRFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    if (!neRarefiedWall::ready<fluid2TThermo>(db(), basic2TThermo::dictName))
    {
        mixedFvPatchScalarField::updateCoeffs();
        return;
    }

    const fluid2TThermo& thermo =
        db().lookupObject<fluid2TThermo>(basic2TThermo::dictName);

    const label patchi = patch().index();

    // Thermo patch state, consistent with the current T_TR face value
    const scalarField& pp = thermo.p().boundaryField()[patchi];
    const scalarField& ppsi = thermo.psi().boundaryField()[patchi];
    const scalarField& pT = *this;

    const scalarField pmu(thermo.mu(patchi));
    const scalarField pkappa(thermo.kappaTRWilke(patchi));
    const scalarField pcv(thermo.CvTR(pp, pT, patchi));

    const scalarField R(neRarefiedWall::gasConstant(ppsi, pT));
    const scalarField gamma(neRarefiedWall::gammaTR(pcv, R));

    // C_T = lambda 2 gamma/(gamma + 1) kappa/(mu cp) (2 - a)/a  [m]
    const scalarField CT
    (
        neRarefiedWall::meanFreePathVHS(pmu, pp, ppsi, omega_)
       *2.0*gamma/(gamma + 1.0)
       *pkappa/(pmu*(pcv + R))
       *neRarefiedWall::accommodationFactor(accommodationCoeff_)
    );

    valueFraction() = 1.0/(1.0 + patch().deltaCoeffs()*CT);
    refValue() = Twall_;
    refGrad() = Zero;

    mixedFvPatchScalarField::updateCoeffs();
}


void Foam::neSmoluchowskiJumpTRFvPatchScalarField::write(Ostream& os) const
{
    // type, refValue, refGradient, valueFraction, value
    mixedFvPatchScalarField::write(os);

    os.writeEntry("accommodationCoeff", accommodationCoeff_);
    os.writeEntry("omega", omega_);
    Twall_.writeEntry("Twall", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        neSmoluchowskiJumpTRFvPatchScalarField
    );
}

// ************************************************************************* //
