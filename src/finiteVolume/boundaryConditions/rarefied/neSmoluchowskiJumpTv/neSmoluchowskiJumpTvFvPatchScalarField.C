/*---------------------------------------------------------------------------*\
  IHIS - Implicit Hypersonics In-house Solver (OpenFOAM-v2312)
  Derived from rhoCentralFoam/BCs/T/smoluchowskiJumpTFvPatchScalarField
  (Copyright (C) 2011-2015 OpenFOAM Foundation, 2020 OpenCFD Ltd., GPL v3)
\*---------------------------------------------------------------------------*/

#include "neSmoluchowskiJumpTvFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "fluid2TThermo.H"
#include "neRarefiedWall.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::word Foam::neSmoluchowskiJumpTvFvPatchScalarField::specieName
(
    const word& fieldName
)
{
    const std::string prefix("TVib_");

    if
    (
        fieldName.size() <= prefix.size()
     || fieldName.compare(0, prefix.size(), prefix) != 0
    )
    {
        FatalErrorInFunction
            << "neSmoluchowskiJumpTv is for the per-molecule fields TVib_<m>,"
            << " not for '" << fieldName << "'." << nl
            << "    Mixture TVib walls take 'calculated' (Phase 7 D12)."
            << exit(FatalError);
    }

    return word(fieldName.substr(prefix.size()));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::neSmoluchowskiJumpTvFvPatchScalarField::
neSmoluchowskiJumpTvFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    TName_("TTR"),
    specie_(),
    accommodationCoeff_(1.0),
    omega_(0.74),       // only without a dictionary; the dictionary requires omega
    Twall_(p.size(), Zero)
{
    refValue() = Zero;
    refGrad() = Zero;
    valueFraction() = Zero;
}


Foam::neSmoluchowskiJumpTvFvPatchScalarField::
neSmoluchowskiJumpTvFvPatchScalarField
(
    const neSmoluchowskiJumpTvFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    TName_(ptf.TName_),
    specie_(ptf.specie_),
    accommodationCoeff_(ptf.accommodationCoeff_),
    omega_(ptf.omega_),
    Twall_(ptf.Twall_, mapper)
{}


Foam::neSmoluchowskiJumpTvFvPatchScalarField::
neSmoluchowskiJumpTvFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    TName_(dict.getOrDefault<word>("T", "TTR")),
    specie_(specieName(iF.name())),
    accommodationCoeff_
    (
        neRarefiedWall::readAccommodation(dict, "accommodationCoeff")
    ),
    omega_(neRarefiedWall::readOmega(dict)),
    Twall_("Twall", dict, p.size())
{
    // Until the coefficients are first available: T_v,b = T_w (D13)
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


Foam::neSmoluchowskiJumpTvFvPatchScalarField::
neSmoluchowskiJumpTvFvPatchScalarField
(
    const neSmoluchowskiJumpTvFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptf, iF),
    TName_(ptf.TName_),
    specie_(ptf.specie_),
    accommodationCoeff_(ptf.accommodationCoeff_),
    omega_(ptf.omega_),
    Twall_(ptf.Twall_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::neSmoluchowskiJumpTvFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    mixedFvPatchScalarField::autoMap(m);
    Twall_.autoMap(m);
}


void Foam::neSmoluchowskiJumpTvFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    mixedFvPatchScalarField::rmap(ptf, addr);

    const auto& tptf =
        refCast<const neSmoluchowskiJumpTvFvPatchScalarField>(ptf);
    Twall_.rmap(tptf.Twall_, addr);
}


void Foam::neSmoluchowskiJumpTvFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const word GammaName("GammaVib_" + specie_);

    if
    (
        !neRarefiedWall::ready<fluid2TThermo>(db(), basic2TThermo::dictName)
     || !db().foundObject<volScalarField>(GammaName)
    )
    {
        // D13: keep the current coefficients
        mixedFvPatchScalarField::updateCoeffs();
        return;
    }

    const fluid2TThermo& thermo =
        db().lookupObject<fluid2TThermo>(basic2TThermo::dictName);

    const label patchi = patch().index();

    const scalarField& pp = thermo.p().boundaryField()[patchi];
    const scalarField& ppsi = thermo.psi().boundaryField()[patchi];
    const scalarField& pTTR =
        patch().lookupPatchField<volScalarField>(TName_);
    const scalarField& pGamma =
        patch().lookupPatchField<volScalarField>(GammaName);

    const scalarField pmu(thermo.mu(patchi));
    const scalarField pcv(thermo.CvTR(pp, pTTR, patchi));

    const scalarField R(neRarefiedWall::gasConstant(ppsi, pTTR));
    const scalarField gamma(neRarefiedWall::gammaTR(pcv, R));

    // C_v = lambda 2 gamma/(gamma + 1) Gamma_v/mu (2 - a)/a  [m]
    const scalarField Cv
    (
        neRarefiedWall::meanFreePathVHS(pmu, pp, ppsi, omega_)
       *2.0*gamma/(gamma + 1.0)
       *pGamma/pmu
       *neRarefiedWall::accommodationFactor(accommodationCoeff_)
    );

    valueFraction() = 1.0/(1.0 + patch().deltaCoeffs()*Cv);
    refValue() = Twall_;
    refGrad() = Zero;

    mixedFvPatchScalarField::updateCoeffs();
}


void Foam::neSmoluchowskiJumpTvFvPatchScalarField::write(Ostream& os) const
{
    // type, refValue, refGradient, valueFraction, value
    mixedFvPatchScalarField::write(os);

    os.writeEntryIfDifferent<word>("T", "TTR", TName_);
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
        neSmoluchowskiJumpTvFvPatchScalarField
    );
}

// ************************************************************************* //
