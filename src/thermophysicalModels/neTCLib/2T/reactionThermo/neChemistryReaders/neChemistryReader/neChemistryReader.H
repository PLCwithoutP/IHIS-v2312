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
    Foam::neChemistryReader

Group
    grpReactionThermophysicalChemistryReaders

Description
    Abstract class for reading chemistry

SourceFiles
    neChemistryReader.C

\*---------------------------------------------------------------------------*/

#ifndef nechemistryReader_H
#define nechemistryReader_H

#include "typeInfo.H"
#include "specieElement.H"
#include "neReaction.H"
#include "neReactionList.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

typedef HashTable<List<specieElement>> speciesCompositionTable;

/*---------------------------------------------------------------------------*\
                      Class neChemistryReader Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class neChemistryReader
{
public:

    //- Runtime type information
    TypeName("neChemistryReader");

    //- The type of thermo package the reader was instantiated for
    typedef ThermoType thermoType;


    // Constructors

        //- Default construct
        neChemistryReader()
        {}


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            neChemistryReader,
            dictionary,
            (
                const dictionary& thermoDict,
                speciesTable& species
            ),
            (thermoDict, species)
        );


    // Selectors

        //- Select constructed from dictionary
        static autoPtr<neChemistryReader> New
        (
            const dictionary& thermoDict,
            speciesTable& species
        );


    //- Destructor
    virtual ~neChemistryReader() = default;


    // Member Functions

        //- Return access to the list of species
        virtual const speciesTable& species() const = 0;

        //- Table of species composition
        virtual const speciesCompositionTable& specieComposition() const = 0;

        //- Return access to the thermo packages
        virtual const ReactionTable<ThermoType>& speciesThermo() const = 0;

        //- Return access to the list of reactions
        virtual const neReactionList<ThermoType>& reactions() const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "neChemistryReader.C"
#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#define makeChemistryReader(Thermo)                                            \
    defineTemplateTypeNameAndDebug(neChemistryReader<Thermo>, 0);                \
    defineTemplateRunTimeSelectionTable(neChemistryReader<Thermo>, dictionary)


#define makeChemistryReaderType(Reader, Thermo)                                \
    defineNamedTemplateTypeNameAndDebug(Reader<Thermo>, 0);                    \
    neChemistryReader<Thermo>::adddictionaryConstructorToTable<Reader<Thermo>>   \
        add##Reader##Thermo##ConstructorToTable_


// for non-templated chemistry readers
#define addChemistryReaderType(Reader, Thermo)                                 \
    defineTypeNameAndDebug(Reader, 0);                                         \
    neChemistryReader<Thermo>::adddictionaryConstructorToTable<Reader>           \
        add##Reader##Thermo##ConstructorToTable_


// for templated chemistry readers
#define addTemplateChemistryReaderType(Reader, Thermo)                         \
    defineNamedTemplateTypeNameAndDebug(Reader, 0);                            \
    neChemistryReader<Thermo>::adddictionaryConstructorToTable<Reader>           \
        add##Reader##Thermo##ConstructorToTable_


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
