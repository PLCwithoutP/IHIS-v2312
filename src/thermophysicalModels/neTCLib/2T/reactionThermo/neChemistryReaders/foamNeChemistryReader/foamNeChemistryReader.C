/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019-2021 OpenCFD Ltd.
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

#include "foamNeChemistryReader.H"
#include "IFstream.H"
#include "stringOps.H"
#include "addToRunTimeSelectionTable.H"
#include "IOdictionary.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

namespace
{
    Foam::fileName expandFileName(const Foam::fileName& raw)
    {
        Foam::fileName f(raw);
        f.expand();
        return f;
    }

    Foam::dictionary readChemistryDictWithTwoTemperature
    (
        const Foam::fileName& reactionsFileName,
        const bool twoTemperature
    )
    {
        Foam::dictionary chemDict
        (
            Foam::IFstream
            (
                expandFileName(reactionsFileName)
            )()
        );

        if (!chemDict.found("reactions"))
        {
            return chemDict;
        }

        Foam::dictionary& reactionsDict = chemDict.subDict("reactions");

        forAllIter(Foam::dictionary, reactionsDict, iter)
        {
            if (!iter().isDict())
            {
                continue;
            }

            Foam::dictionary& reactionDict =
                reactionsDict.subDict(iter().keyword());

            reactionDict.add("twoTemperature", twoTemperature, true);
        }

        return chemDict;
    }
}

template<class ThermoType>
Foam::speciesTable& Foam::foamNeChemistryReader<ThermoType>::setSpecies
(
    const dictionary& dict,
    speciesTable& species
)
{
    wordList s(dict.get<wordList>("species"));
    species.transfer(s);
    return species;
}


template<class ThermoType>
void Foam::foamNeChemistryReader<ThermoType>::readSpeciesComposition()
{
    wordList elems;

    if (!chemDict_.readIfPresent("elements", elems))
    {
        Info<< "    elements not defined in " << chemDict_.name() << endl;
        return;
    }

    elementNames_.clear();
    elementIndices_.clear();

    for (const word& elemName : elems)
    {
        if (elementIndices_.insert(elemName, elementNames_.size()))
        {
            elementNames_.append(elemName);
        }
        else
        {
            IOWarningInFunction(chemDict_)
                << "element " << elemName << " already in table." << endl;
        }
    }

    for (const word& specieName : speciesTable_)
    {
        const dictionary* elemsDict =
            thermoDict_.subDict(specieName).findDict("elements");

        if (!elemsDict)
        {
            FatalIOErrorInFunction(thermoDict_)
                << "Specie " << specieName
                << " does not contain \"elements\" description."
                << exit(FatalIOError);
        }

        wordList elemNames(elemsDict->toc());
        List<specieElement> currentComposition(elemNames.size());

        forAll(elemNames, eni)
        {
            currentComposition[eni].name() = elemNames[eni];
            currentComposition[eni].nAtoms() =
                elemsDict->getOrDefault<label>(elemNames[eni], 0);
        }

        speciesComposition_.erase(specieName);
        speciesComposition_.set(specieName, currentComposition);
    }
}


// * * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::foamNeChemistryReader<ThermoType>::foamNeChemistryReader
(
    const fileName& reactionsFileName,
    speciesTable& species,
    const fileName& thermoFileName
)
:
    neChemistryReader<ThermoType>(),
    chemDict_
    (
        readChemistryDictWithTwoTemperature
        (
            reactionsFileName, 
            false
        )
    ),
    thermoDict_
    (
        IFstream
        (
            expandFileName(thermoFileName)
        )()
    ),
    twoTemperature_(false),
    speciesTable_(setSpecies(chemDict_, species)),
    speciesThermo_(thermoDict_),
    reactions_(speciesTable_, speciesThermo_, chemDict_)
{
    readSpeciesComposition();
}

template<class ThermoType>
Foam::foamNeChemistryReader<ThermoType>::foamNeChemistryReader
(
    const dictionary& thermoDict,
    speciesTable& species
)
:
    neChemistryReader<ThermoType>(),
    chemDict_
    (
        readChemistryDictWithTwoTemperature
        (
            thermoDict.get<fileName>("reactionsList"),
            thermoDict.getOrDefault<bool>("twoTemperature", false)
        )
    ),
    thermoDict_
    (
        IFstream
        (
            expandFileName(thermoDict.get<fileName>("speciesList"))
        )()
    ),
    twoTemperature_(thermoDict.getOrDefault<bool>("twoTemperature", false)),
    speciesTable_(setSpecies(chemDict_, species)),
    speciesThermo_(thermoDict_),
    reactions_(speciesTable_, speciesThermo_, chemDict_)
{
    readSpeciesComposition();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
