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
#include "energy2TJumpAMIFvPatchScalarField.H"
#include "fixedJumpAMIFvPatchFields.H"
#include "basic2TThermo.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::energy2TJumpAMIFvPatchScalarField::energy2TJumpAMIFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedJumpAMIFvPatchField<scalar>(p, iF)
{}


Foam::energy2TJumpAMIFvPatchScalarField::energy2TJumpAMIFvPatchScalarField
(
    const energy2TJumpAMIFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedJumpAMIFvPatchField<scalar>(ptf, p, iF, mapper)
{}


Foam::energy2TJumpAMIFvPatchScalarField::energy2TJumpAMIFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedJumpAMIFvPatchField<scalar>(p, iF)
{
    if (!this->readValueEntry(dict))
    {
        evaluate(Pstream::commsTypes::blocking);
    }
}


Foam::energy2TJumpAMIFvPatchScalarField::energy2TJumpAMIFvPatchScalarField
(
    const energy2TJumpAMIFvPatchScalarField& ptf
)
:
    fixedJumpAMIFvPatchField<scalar>(ptf)
{}


Foam::energy2TJumpAMIFvPatchScalarField::energy2TJumpAMIFvPatchScalarField
(
    const energy2TJumpAMIFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedJumpAMIFvPatchField<scalar>(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::energy2TJumpAMIFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (this->cyclicAMIPatch().owner())
    {
        const basic2TThermo& thermo2T = basic2TThermo::lookupThermo(*this);
        label patchID = patch().index();

        const scalarField& pp = thermo2T.p().boundaryField()[patchID];
        const fixedJumpAMIFvPatchScalarField& TbPatch =
            refCast<const fixedJumpAMIFvPatchScalarField>
            (
                thermo2T.TTR().boundaryField()[patchID]
            );

        fixedJumpAMIFvPatchScalarField& Tbp =
            const_cast<fixedJumpAMIFvPatchScalarField&>(TbPatch);

        // force update of jump
        Tbp.evaluate(Pstream::commsTypes::blocking);

        const labelUList& faceCells = this->patch().faceCells();

        jump_ =
            thermo2T.h(pp, Tbp+Tbp.jump(), faceCells)
          - thermo2T.h(pp, Tbp, faceCells);
    }

    fixedJumpAMIFvPatchField<scalar>::updateCoeffs();
}


void Foam::energy2TJumpAMIFvPatchScalarField::write(Ostream& os) const
{
    fixedJumpAMIFvPatchField<scalar>::write(os);
    fvPatchField<scalar>::writeValueEntry(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
   makePatchTypeField
   (
       fvPatchScalarField,
       energy2TJumpAMIFvPatchScalarField
   );
}

// ************************************************************************* //
