/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2017 OpenFOAM Foundation
    Copyright (C) 2017-2021 OpenCFD Ltd.
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

#include "basic2TThermo.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class TThermo, class ThermoConstructTable>
typename ThermoConstructTable::mapped_type
Foam::basic2TThermo::getThermoOrDie
(
    const dictionary& thermoTypeDict,
    ThermoConstructTable& thermoTable,
    const word& thermoTypeName,
    const wordList& cmptNames
)
{
    // Lookup the thermo package

    auto ctorIter = thermoTable.cfind(thermoTypeName);

    // Print error message if package not found in the table
    if (!ctorIter.good())
    {
        FatalIOErrorInLookup
        (
            thermoTypeDict,
            TThermo::typeName,
            word::null, // Suppress long name? Just output dictionary (above)
            thermoTable
        );

        basic2TThermo::printThermoNames
        (
            FatalIOError,
            cmptNames,
            thermoTable.sortedToc()
        ) << exit(FatalIOError);

        // return nullptr;
    }

    return ctorIter.val();
}


template<class TThermo, class ThermoConstructTable>
typename ThermoConstructTable::mapped_type
Foam::basic2TThermo::getThermoOrDie
(
    const dictionary& thermoDict,
    ThermoConstructTable& thermoTable
)
{
    const dictionary* dictptr = thermoDict.findDict("thermoType");
    if (dictptr)
    {
        const auto& thermoTypeDict = *dictptr;

        const wordList* cmptHeaderPtr = &(wordList::null());

        // TThermo package name, constructed from components
        const word thermoTypeName
        (
            basic2TThermo::makeThermoName(thermoTypeDict, cmptHeaderPtr)
        );

        Info<< "Selecting thermodynamics package " << thermoTypeDict << endl;

        return getThermoOrDie<TThermo, ThermoConstructTable>
        (
            thermoTypeDict,
            thermoTable,
            thermoTypeName,
            *cmptHeaderPtr
        );
    }
    else
    {
        const word thermoTypeName(thermoDict.get<word>("thermoType"));

        Info<< "Selecting thermodynamics package " << thermoTypeName << endl;

        auto ctorIter = thermoTable.cfind(thermoTypeName);

        if (!ctorIter.good())
        {
            FatalIOErrorInLookup
            (
                thermoDict,
                TThermo::typeName,
                thermoTypeName,
                thermoTable
            ) << exit(FatalIOError);
        }

        return ctorIter.val();
    }
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

template<class TThermo>
Foam::autoPtr<TThermo> Foam::basic2TThermo::New
(
    const fvMesh& mesh,
    const word& phaseName
)
{
    IOdictionary thermoDict
    (
        IOobject
        (
            phasePropertyName(dictName, phaseName),
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        )
    );

    auto* ctorPtr = getThermoOrDie<TThermo>
    (
        thermoDict,
        *(TThermo::fvMeshConstructorTablePtr_)
    );

    return autoPtr<TThermo>(ctorPtr(mesh, phaseName));
}


template<class TThermo>
Foam::autoPtr<TThermo> Foam::basic2TThermo::New
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& phaseName
)
{
    auto* ctorPtr = getThermoOrDie<TThermo>
    (
        dict,
        *(TThermo::dictionaryConstructorTablePtr_)
    );

    return autoPtr<TThermo>(ctorPtr(mesh, dict, phaseName));
}


template<class TThermo>
Foam::autoPtr<TThermo> Foam::basic2TThermo::New
(
    const fvMesh& mesh,
    const word& phaseName,
    const word& dictName
)
{
    IOdictionary thermoDict
    (
        IOobject
        (
            dictName,
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        )
    );

    auto* ctorPtr = getThermoOrDie<TThermo>
    (
        thermoDict,
        *(TThermo::fvMeshDictPhaseConstructorTablePtr_)
    );

    return autoPtr<TThermo>(ctorPtr(mesh, phaseName, dictName));
}


// ************************************************************************* //
