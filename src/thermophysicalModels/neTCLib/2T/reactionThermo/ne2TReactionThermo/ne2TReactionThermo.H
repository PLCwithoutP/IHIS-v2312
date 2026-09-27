/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2017-2020 OpenCFD Ltd.
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
    Foam::ne2TReactionThermo

Group
    grpPsiThermo

Description
    Foam::ne2TReactionThermo

SourceFiles
    ne2TReactionThermo.C

\*---------------------------------------------------------------------------*/

#ifndef ne2TReactionThermo_H
#define ne2TReactionThermo_H

#include "ne2TThermo.H"
#include "basicSpecie2TMixture.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "specieElement.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

typedef HashTable<List<specieElement>> speciesCompositionTable;

/*---------------------------------------------------------------------------*\
                     Class ne2TReactionThermo Declaration
\*---------------------------------------------------------------------------*/

class ne2TReactionThermo
:
    public ne2TThermo
{

public:

    //- Runtime type information
    TypeName("ne2TReactionThermo");


    //- Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        ne2TReactionThermo,
        fvMesh,
        (const fvMesh& mesh, const word& phaseName),
        (mesh, phaseName)
    );


    //- Declare run-time constructor selection table for dictionary based
    declareRunTimeSelectionTable
    (
        autoPtr,
        ne2TReactionThermo,
        fvMeshDictPhase,
        (const fvMesh& mesh, const word& phaseName, const word& dictName),
        (mesh, phaseName, dictName)
    );


    // Constructors

        //- Construct from mesh and phase name
        ne2TReactionThermo
        (
            const fvMesh&,
            const word& phaseName
        );


        //- Construct from mesh,dictionary,phase name with a single temperature
        ne2TReactionThermo
        (
            const fvMesh&,
            const word& phaseName,
            const word& dictName
        );


    // Selectors

        //- Standard selection based on fvMesh
        static autoPtr<ne2TReactionThermo> New
        (
            const fvMesh&,
            const word& phaseName=word::null
        );


        //- Selector based on mesh, phase name and dictionary
        static autoPtr<ne2TReactionThermo> New
        (
            const fvMesh&,
            const word& phaseName,
            const word& dictName
        );


    //- Destructor
    virtual ~ne2TReactionThermo();


    // Member functions

        //- Return the composition of the multi-component mixture
        virtual compositionInterface& composition() = 0;

        //- Return the composition of the multi-component mixture
        virtual const compositionInterface& composition() const = 0;

        //- Table of species composition
        autoPtr<speciesCompositionTable> specieComposition() const
        {
            return
                autoPtr<speciesCompositionTable>
                (
                    static_cast<const basicSpecie2TMixture&>
                    (
                        composition()
                    ).specieComposition()
                );
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
