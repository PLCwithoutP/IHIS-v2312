/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2021 OpenCFD Ltd.
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
    Foam::foamNeChemistryReader

Group
    grpReactionThermophysicalChemistryReaders

Description
    Chemistry reader for OpenFOAM format

SourceFiles
    foamNeChemistryReader.C

\*---------------------------------------------------------------------------*/

#ifndef foamNeChemistryReader_H
#define foamNeChemistryReader_H

#include "neChemistryReader.H"
#include "fileName.H"
#include "labelList.H"
#include "speciesTable.H"
#include "atomicWeights.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class foamChemistry Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class foamNeChemistryReader
:
    public neChemistryReader<ThermoType>
{
        //- Chemistry/reactions dictionary
        dictionary chemDict_;

        //- Thermo properties dictionary
        dictionary thermoDict_;

        //- List of elements
        DynamicList<word> elementNames_;

        //- Element indices
        HashTable<label> elementIndices_;

        //- Table of species
        speciesTable& speciesTable_;

        //- Table of species composition
        speciesCompositionTable speciesComposition_;

        //- Table of the thermodynamic data given in the foamChemistry file
        ReactionTable<ThermoType> speciesThermo_;

        //- List of the reactions
        neReactionList<ThermoType> reactions_;

        //- Global two-temperature flag read from thermophysicalProperties
        bool twoTemperature_;

    // Private Member Functions

        //- Set the species list
        speciesTable& setSpecies(const dictionary& dict, speciesTable& species);

        //- Read the species composition
        void readSpeciesComposition();

        //- No copy construct
        foamNeChemistryReader(const foamNeChemistryReader&) = delete;

        //- No copy assignment
        void operator=(const foamNeChemistryReader&) = delete;

public:

    //- Runtime type information
    TypeName("foamNeChemistryReader");


    // Constructors

        //- Construct from foamChemistry and thermodynamics file names
        foamNeChemistryReader
        (
            const fileName& reactionsFileName,
            speciesTable& species,
            const fileName& thermoFileName
        );

        //- Construct by getting foamChemistry and thermodynamics file names
        //- from dictionary
        foamNeChemistryReader
        (
            const dictionary& thermoDict,
            speciesTable& species
        );


    //- Destructor
    virtual ~foamNeChemistryReader() = default;


    // Member Functions

        //- List of elements
        const wordList& elementNames() const
        {
            return elementNames_;
        }

        //- Element indices
        const HashTable<label>& elementIndices() const
        {
            return elementIndices_;
        }

        //- Table of species
        const speciesTable& species() const
        {
            return speciesTable_;
        }

        //- Table of species composition
        const speciesCompositionTable& specieComposition() const
        {
            return speciesComposition_;
        }

        //- Table of the thermodynamic data given in the foamChemistry file
        const ReactionTable<ThermoType>& speciesThermo() const
        {
            return speciesThermo_;
        }

        //- List of the reactions
        const neReactionList<ThermoType>& reactions() const
        {
            return reactions_;
        }

        bool twoTemperature() const noexcept
        {
            return twoTemperature_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "foamNeChemistryReader.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
