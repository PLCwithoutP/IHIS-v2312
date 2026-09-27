/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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
    Foam::neReactingMixture

Group
    grpReactionThermophysicalMixtures

Description
    Foam::neReactingMixture

SourceFiles
    neReactingMixture.C

\*---------------------------------------------------------------------------*/

#ifndef neReactingMixture_H
#define neReactingMixture_H

#include "speciesTable.H"
#include "neChemistryReader.H"
#include "multiComponent2TMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class neReactingMixture Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class neReactingMixture
:
    public speciesTable,
    private autoPtr<neChemistryReader<ThermoType>>,
    public multiComponent2TMixture<ThermoType>,
    public PtrList<neReaction<ThermoType>>
{
    // Private Member Data

        //- Table of species composition
        speciesCompositionTable speciesComposition_;


    // Private Member Functions

        //- The chemistry reader
        autoPtr<neChemistryReader<ThermoType>>& reader();

        //- No copy construct
        neReactingMixture(const neReactingMixture&) = delete;

        //- No copy assignment
        void operator=(const neReactingMixture&) = delete;


public:

    //- The type of thermo package this mixture is instantiated for
    typedef ThermoType thermoType;


    // Constructors

        //- Construct from dictionary, mesh and phase name
        neReactingMixture
        (
            const dictionary& thermoDict,
            const fvMesh& mesh,
            const word& phaseName
        );


    //- Destructor
    virtual ~neReactingMixture() = default;


    // Member functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "neReactingMixture<" + ThermoType::typeName() + '>';
        }

        //- Read dictionary
        void read(const dictionary&);

        using PtrList<neReaction<ThermoType>>::size;
        using PtrList<neReaction<ThermoType>>::operator[];

        //- Table of species composition
        virtual autoPtr<speciesCompositionTable> specieComposition() const
        {
            return
                autoPtr<speciesCompositionTable>
                (
                    new speciesCompositionTable(speciesComposition_)
                );
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "neReactingMixture.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
