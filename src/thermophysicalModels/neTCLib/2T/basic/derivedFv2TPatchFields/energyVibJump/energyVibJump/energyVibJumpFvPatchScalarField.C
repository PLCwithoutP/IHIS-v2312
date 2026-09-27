/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2017 OpenFOAM Foundation
    Copyright (C) 2022 OpenCFD Ltd.
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

#include "addToRunTimeSelectionTable.H"
#include "energyVibJumpFvPatchScalarField.H"
#include "fixedJumpFvPatchFields.H"
#include "basic2TThermo.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::energyVibJumpFvPatchScalarField::energyVibJumpFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedJumpFvPatchField<scalar>(p, iF)
{}


Foam::energyVibJumpFvPatchScalarField::energyVibJumpFvPatchScalarField
(
    const energyVibJumpFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedJumpFvPatchField<scalar>(ptf, p, iF, mapper)
{}


Foam::energyVibJumpFvPatchScalarField::energyVibJumpFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedJumpFvPatchField<scalar>(p, iF)
{
    if (!this->readValueEntry(dict))
    {
        evaluate(Pstream::commsTypes::blocking);
    }
}


Foam::energyVibJumpFvPatchScalarField::energyVibJumpFvPatchScalarField
(
    const energyVibJumpFvPatchScalarField& ptf
)
:
    fixedJumpFvPatchField<scalar>(ptf)
{}


Foam::energyVibJumpFvPatchScalarField::energyVibJumpFvPatchScalarField
(
    const energyVibJumpFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedJumpFvPatchField<scalar>(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::energyVibJumpFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (this->cyclicPatch().owner())
    {
        const basic2TThermo& thermo2T = basic2TThermo::lookupThermo(*this);
        label patchID = patch().index();

        const scalarField& pp   = thermo2T.p().boundaryField()[patchID];
        const scalarField& TTRp = thermo2T.TTR().boundaryField()[patchID];
        const fixedJumpFvPatchScalarField& TbPatch = refCast<const fixedJumpFvPatchScalarField>
            (thermo2T.TVib().boundaryField()[patchID]);
        fixedJumpFvPatchScalarField& Tbp = const_cast<fixedJumpFvPatchScalarField&>(TbPatch);
        Tbp.evaluate(Pstream::commsTypes::blocking);

        const scalar theta = thermo2T.thetaVib();
        const labelUList& faceCells = this->patch().faceCells();
        setJump
        (
            thermo2T.eVib(pp, TTRp, Tbp+Tbp.jump(), theta, faceCells)
        - thermo2T.eVib(pp, TTRp, Tbp,            theta, faceCells)
        );
    }

    fixedJumpFvPatchField<scalar>::updateCoeffs();
}


void Foam::energyVibJumpFvPatchScalarField::write(Ostream& os) const
{
    fixedJumpFvPatchField<scalar>::write(os);
    fvPatchField<scalar>::writeValueEntry(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
   makePatchTypeField
   (
       fvPatchScalarField,
       energyVibJumpFvPatchScalarField
   );
}

// ************************************************************************* //
