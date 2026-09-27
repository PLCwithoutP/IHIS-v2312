/*---------------------------------------------------------------------------*\

    IHIS: In-house Hypersonics Implicit solver

    Copyright (C) 2014-2018 Oliver Oxtoby - CSIR, South Africa
    Copyright (C) 2014-2018 Johan Heyns - CSIR, South Africa
    Copyright (C) 2020 Oliver Oxtoby
    Copyright (C) 2011-2013 OpenFOAM Foundation

-------------------------------------------------------------------------------
License
    This file is part of IHIS.

    IHIS is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    IHIS is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with IHIS.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "characteristicWallPressureFvPatchScalarField.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "addToRunTimeSelectionTable.H"
#include "fluid2TThermo.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::characteristicWallPressureFvPatchScalarField::characteristicWallPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF),
    BlockCoupledBoundary(dynamic_cast<fvPatchField<scalar>&>(*this)),
    UCoeff_(p.size(), scalar(0))
{}


Foam::characteristicWallPressureFvPatchScalarField::characteristicWallPressureFvPatchScalarField
(
    const characteristicWallPressureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper),
    BlockCoupledBoundary(dynamic_cast<fvPatchField<scalar>&>(*this)),
    UCoeff_(ptf.UCoeff_, mapper)
{
    patchType() = ptf.patchType();
    // Enforce mapping of values so we have a valid starting value. This
    // constructor is used when reconstructing fields
    this->map(ptf, mapper);
}


Foam::characteristicWallPressureFvPatchScalarField::characteristicWallPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF),
    BlockCoupledBoundary(dynamic_cast<fvPatchField<scalar>&>(*this)),
    UCoeff_(scalarField(p.size(), scalar(0)))
{
    if (dict.found("value") && dict.found("gradient"))
    {
        // Full restart
        fvPatchField<scalar>::operator=
        (
            scalarField("value", dict, p.size())
        );
        gradient() = scalarField("gradient", dict, p.size());
        UCoeff_ = scalarField("UCoeff", dict, p.size());
    }
    else
    {
        // Default to zeroGradient
        fvPatchField<scalar>::operator=(patchInternalField());
        gradient() = 0.0;
    }
}


Foam::characteristicWallPressureFvPatchScalarField::characteristicWallPressureFvPatchScalarField
(
    const characteristicWallPressureFvPatchScalarField& wbppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(wbppsf, iF),
    BlockCoupledBoundary(dynamic_cast<fvPatchField<scalar>&>(*this)),
    UCoeff_(wbppsf.UCoeff_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::characteristicWallPressureFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedGradientFvPatchScalarField::autoMap(m);
    UCoeff_.autoMap(m);
}


void Foam::characteristicWallPressureFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    fixedGradientFvPatchScalarField::rmap(ptf, addr);

    const characteristicWallPressureFvPatchScalarField& mptf =
        refCast<const characteristicWallPressureFvPatchScalarField>(ptf);

    UCoeff_.rmap(mptf.UCoeff_, addr);
}


void Foam::characteristicWallPressureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvMesh& mesh = patch().boundaryMesh().mesh();

    const fluidThermo& thermo =
        mesh.lookupObject<fluidThermo>(basicThermo::dictName);

    const fvPatchVectorField& Up =
        patch().lookupPatchField<volVectorField, vector>("U");
    vectorField Uif
    (
        Up.patchInternalField()
    );
    scalarField pif(this->patchInternalField());

    scalarField psiif
    (
        thermo.psi().boundaryField()[patch().index()].patchInternalField()
    );
    scalarField gammaif
    (
        thermo.gamma()->boundaryField()[patch().index()].patchInternalField()
    );
    scalarField rhoif
    (
        thermo.rho()->boundaryField()[patch().index()].patchInternalField()
    );

    // Reference values (Blazek suggests using internal values at cell centres)
    const scalarField cif(sqrt(gammaif/psiif));

    UCoeff_ = rhoif*cif;
    gradient() = rhoif*cif*((Uif-Up)&patch().nf())*patch().deltaCoeffs();

    // Limit so that p would not go negative based on current values
    scalarField minGrad((pif-SMALL)*patch().deltaCoeffs());
    scalarField limit
    (
        max(min(minGrad/stabilise(gradient(), SMALL), scalar(1)), scalar(0))
    );
    UCoeff_ *= limit;
    gradient() *= limit;

    fixedGradientFvPatchScalarField::updateCoeffs();
}

Foam::tmp<Foam::scalarField>
Foam::characteristicWallPressureFvPatchScalarField::valueInternalCoeffs
(
    const volScalarField& field
) const
{
    if (field.name() == internalField().name())
    {
        return tmp<scalarField>(new scalarField(patch().size(), scalar(1)));
    }
    else if (field.member() == "T")
    {
        return tmp<scalarField>(new scalarField(patch().size(), scalar(0)));
    }
    else
    {
        NotImplemented;
        return tmp<scalarField>(new scalarField());
    }
}

Foam::tmp<Foam::scalarField>
Foam::characteristicWallPressureFvPatchScalarField::gradientInternalCoeffs
(
    const volScalarField& field
) const
{
    if (field.name() == internalField().name())
    {
        return -patch().deltaCoeffs();
    }
    else if (field.member() == "T")
    {
        return tmp<scalarField>(new scalarField(patch().size(), scalar(0)));
    }
    else
    {
        NotImplemented;
        return tmp<scalarField>(new scalarField());
    }
}


Foam::tmp<Foam::vectorField>
Foam::characteristicWallPressureFvPatchScalarField::valueInternalCoeffs
(
    const volVectorField& field
) const
{
    if (field.member() == "U")
    {
        return UCoeff_*patch().nf();
    }
    else
    {
        NotImplemented;
        return tmp<vectorField>(new vectorField());
    }
}

Foam::tmp<Foam::vectorField>
Foam::characteristicWallPressureFvPatchScalarField::gradientInternalCoeffs
(
    const volVectorField& field
) const
{
    if (field.member() == "U")
    {
        return -UCoeff_*patch().nf()*patch().deltaCoeffs();
    }
    else
    {
        NotImplemented;
        return tmp<vectorField>(new vectorField());
    }
}


void Foam::characteristicWallPressureFvPatchScalarField::write(Ostream& os) const
{
    fixedGradientFvPatchScalarField::write(os);
    UCoeff_.writeEntry("UCoeff", os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        characteristicWallPressureFvPatchScalarField
    );
}


// ************************************************************************* //
