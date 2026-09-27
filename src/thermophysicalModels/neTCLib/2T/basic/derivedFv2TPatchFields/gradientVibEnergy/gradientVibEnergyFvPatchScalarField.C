/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2012 OpenFOAM Foundation
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "gradientVibEnergyFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "basic2TThermo.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::gradientVibEnergyFvPatchScalarField::
gradientVibEnergyFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF)
{}


Foam::gradientVibEnergyFvPatchScalarField::
gradientVibEnergyFvPatchScalarField
(
    const gradientVibEnergyFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper)
{}


Foam::gradientVibEnergyFvPatchScalarField::
gradientVibEnergyFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF, dict)
{}


Foam::gradientVibEnergyFvPatchScalarField::
gradientVibEnergyFvPatchScalarField
(
    const gradientVibEnergyFvPatchScalarField& tppsf
)
:
    fixedGradientFvPatchScalarField(tppsf)
{}


Foam::gradientVibEnergyFvPatchScalarField::
gradientVibEnergyFvPatchScalarField
(
    const gradientVibEnergyFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(tppsf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::gradientVibEnergyFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const basic2TThermo& thermo2T = basic2TThermo::lookupThermo(*this);
    const label patchi = patch().index();

    const scalarField& pw   = thermo2T.p().boundaryField()[patchi];
    const scalarField& TTRw = thermo2T.TTR().boundaryField()[patchi];
    fvPatchScalarField& TVw =
        const_cast<fvPatchScalarField&>(thermo2T.TVib().boundaryField()[patchi]);
    TVw.evaluate();

    const scalar theta = thermo2T.thetaVib();
    gradient() = thermo2T.CvVib(pw, TTRw, TVw, theta, patchi)*TVw.snGrad()
    + patch().deltaCoeffs()*
        (
            thermo2T.eVib(pw, TTRw, TVw, theta, patchi)
        - thermo2T.eVib(pw, TTRw, TVw, theta, patch().faceCells())
        );
    fixedGradientFvPatchScalarField::updateCoeffs();
}


void Foam::gradientVibEnergyFvPatchScalarField::write(Ostream& os) const
{
    fixedGradientFvPatchField<scalar>::write(os);
    fvPatchField<scalar>::writeValueEntry(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        gradientVibEnergyFvPatchScalarField
    );
}

// ************************************************************************* //
