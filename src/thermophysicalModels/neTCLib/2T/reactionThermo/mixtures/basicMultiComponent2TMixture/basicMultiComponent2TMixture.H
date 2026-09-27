/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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
    Foam::basicMultiComponent2TMixture

Group
    grpReactionThermophysicalMixtures

Description
    Multi-component mixture.

    Provides a list of mass fraction fields and helper functions to
    query mixture composition.

SourceFiles
    basicMultiComponent2TMixture.C

\*---------------------------------------------------------------------------*/

#ifndef basicMultiComponent2TMixture_H
#define basicMultiComponent2TMixture_H

#include "volFields.H"
#include "PtrList.H"
#include "basic2TMixture.H"
#include "speciesTable.H"
#include "typeInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class basicMultiComponent2TMixture Declaration
\*---------------------------------------------------------------------------*/

class basicMultiComponent2TMixture
:
    public basic2TMixture
{

protected:

    // Protected data

        //- Table of specie names
        speciesTable species_;

        //- List of specie active flags
        List<bool> active_;

        //- Species mass fractions
        PtrList<volScalarField> Y_;


public:

    //- Run time type information
    TypeName("basicMultiComponent2TMixture");


    //- The base class of the mixture
    typedef basicMultiComponent2TMixture basic2TMixtureType;


    // Constructors

        //- Construct from dictionary, species names, mesh and phase name
        basicMultiComponent2TMixture
        (
            const dictionary& thermoDict,
            const wordList& specieNames,
            const fvMesh& mesh,
            const word& phaseName
        );


    //- Destructor
    virtual ~basicMultiComponent2TMixture() = default;


    // Member functions

        //- Return the table of species
        inline const speciesTable& species() const;

        //- Does the mixture include this specie?
        inline bool contains(const word& specieName) const;

        //- Return true for active species
        inline bool active(label speciei) const;

        //- Return the bool list of active species
        inline const List<bool>& active() const;

        //- Set speciei active
        inline void setActive(label speciei);

        //- Set speciei inactive
        inline void setInactive(label speciei);

        //- Return the mass-fraction fields
        inline PtrList<volScalarField>& Y();

        //- Return the const mass-fraction fields
        inline const PtrList<volScalarField>& Y() const;

        //- Return the mass-fraction field for a specie given by index
        inline volScalarField& Y(const label i);

        //- Return the const mass-fraction field for a specie given by index
        inline const volScalarField& Y(const label i) const;

        //- Return the mass-fraction field for a specie given by name
        inline volScalarField& Y(const word& specieName);

        //- Return the const mass-fraction field for a specie given by name
        inline const volScalarField& Y(const word& specieName) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "basicMultiComponent2TMixtureI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
