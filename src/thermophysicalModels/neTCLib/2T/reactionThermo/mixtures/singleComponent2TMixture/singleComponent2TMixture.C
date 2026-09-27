/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017 OpenFOAM Foundation
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

#include "singleComponent2TMixture.H"
#include "IFstream.H"
#include "stringOps.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::fileName Foam::singleComponent2TMixture<ThermoType>::resolveFile
(
    const dictionary& thermoDict,
    const fvMesh& mesh,
    const word& key
)
{
    if (!thermoDict.found(key))
    {
        FatalIOErrorInFunction(thermoDict)
            << "Missing required entry '" << key << "' in "
            << thermoDict.name()
            << Foam::exit(Foam::FatalIOError);
    }

    const fileName raw(thermoDict.get<fileName>(key));
    const string rawStr(raw);
    const string constantPrefix("<constant>/");

    if
    (
        rawStr.size() >= constantPrefix.size()
     && rawStr.substr(0, constantPrefix.size()) == constantPrefix
    )
    {
        return mesh.time().constant()/fileName(rawStr.substr(constantPrefix.size()));
    }

    fileName expanded(raw);
    expanded.expand();
    return expanded;
}

template<class ThermoType>
Foam::word Foam::singleComponent2TMixture<ThermoType>::readSingleSpeciesName
(
    const dictionary& thermoDict,
    const fvMesh& mesh
)
{
    if (thermoDict.get<bool>("chemistry"))
    {
        FatalIOErrorInFunction(thermoDict)
            << "pureMixture/singleComponent2TMixture does not support chemistry=true."
            << exit(FatalIOError);
    }

    // Legacy fallback
    if (thermoDict.found("mixture"))
    {
        return "mixture";
    }

    const fileName speciesFile = resolveFile(thermoDict, mesh, "speciesList");
    IFstream speciesStream(speciesFile);

    if (!speciesStream.good())
    {
        FatalIOErrorInFunction(speciesStream)
            << "Cannot open dictionary file " << speciesFile
            << Foam::exit(Foam::FatalIOError);
    }

    const dictionary speciesDict(speciesStream);

    if (!speciesDict.found("species"))
    {
        FatalIOErrorInFunction(speciesDict)
            << "Entry 'species' not found in dictionary " << speciesFile
            << exit(FatalIOError);
    }

    const wordList speciesNames(speciesDict.get<wordList>("species"));

    if (speciesNames.size() != 1)
    {
        FatalIOErrorInFunction(speciesDict)
            << "singleComponent2TMixture expects exactly one species in "
            << speciesFile << ", but found " << speciesNames
            << exit(FatalIOError);
    }

    return speciesNames[0];
}

template<class ThermoType>
ThermoType Foam::singleComponent2TMixture<ThermoType>::constructThermo
(
    const dictionary& thermoDict
) const
{
    if (thermoDict.found("mixture"))
    {
        return ThermoType(thermoDict.subDict("mixture"));
    }

    const word specieName = readSingleSpeciesName(thermoDict, mesh_);
    const fileName speciesFile = resolveFile(thermoDict, mesh_, "speciesList");
    IFstream speciesStream(speciesFile);

    if (!speciesStream.good())
    {
        FatalIOErrorInFunction(speciesStream)
            << "Cannot open dictionary file " << speciesFile
            << Foam::exit(Foam::FatalIOError);
    }

    const dictionary speciesDict(speciesStream);

    if (!speciesDict.found(specieName))
    {
        FatalIOErrorInFunction(speciesDict)
            << "Species '" << specieName
            << "' not found in species dictionary " << speciesFile
            << exit(FatalIOError);
    }

    return ThermoType(speciesDict.subDict(specieName));
}

template<class ThermoType>
Foam::singleComponent2TMixture<ThermoType>::singleComponent2TMixture
(
    const dictionary& thermoDict,
    const fvMesh& mesh,
    const word& phaseName
)
:
    basicSpecie2TMixture(thermoDict, wordList::null(), mesh, phaseName),
    mesh_(mesh),
    thermo_(constructThermo(thermoDict))
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class ThermoType>
void Foam::singleComponent2TMixture<ThermoType>::read
(
    const dictionary& thermoDict
)
{
    thermo_ = constructThermo(thermoDict);
}


// ************************************************************************* //
