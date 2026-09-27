/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2016 OpenFOAM Foundation
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

Class
    Foam::energy2TJumpAMIFvPatchScalarField

Group
    grpThermoBoundaryConditions grpCoupledBoundaryConditions

Description
    This boundary condition provides an energy jump condition across a pair
    of coupled patches with an arbitrary mesh interface (AMI).  It is not
    applied directly, but is employed on-the-fly when converting temperature
    boundary conditions into energy.

See also
    Foam::fixedJumpAMIFvPatchField

SourceFiles
    energy2TJumpAMIFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef energy2TJumpAMIFvPatchScalarField_H
#define energy2TJumpAMIFvPatchScalarField_H

#include "fixedJumpAMIFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class energy2TJumpAMIFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class energy2TJumpAMIFvPatchScalarField
:
    public fixedJumpAMIFvPatchField<scalar>
{

public:

    //- Runtime type information
    TypeName("energy2TJumpAMI");

    // Constructors

        //- Construct from patch and internal field
        energy2TJumpAMIFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        energy2TJumpAMIFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given energy2TJumpAMIFvPatchScalarField onto a
        //  new patch
        energy2TJumpAMIFvPatchScalarField
        (
            const energy2TJumpAMIFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        energy2TJumpAMIFvPatchScalarField
        (
            const energy2TJumpAMIFvPatchScalarField&
        );

        //- Construct as copy setting internal field reference
        energy2TJumpAMIFvPatchScalarField
        (
            const energy2TJumpAMIFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Return a clone
        virtual tmp<fvPatchField<scalar>> clone() const
        {
            return tmp<fvPatchScalarField>(new energy2TJumpAMIFvPatchScalarField(*this));
        }

        //- Clone with an internal field reference
        virtual tmp<fvPatchField<scalar>> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>(new energy2TJumpAMIFvPatchScalarField(*this, iF));
        }


    // Member functions

        // Evaluation functions

            //- Update the coefficients
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
