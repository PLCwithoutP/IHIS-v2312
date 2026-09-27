/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2021-2022 OpenCFD Ltd.
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

#include "multiComponent2TMixture.H"
#include "IFstream.H"
#include "stringOps.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class ThermoType>
Foam::fileName Foam::multiComponent2TMixture<ThermoType>::resolveFile
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
Foam::wordList Foam::multiComponent2TMixture<ThermoType>::readSpeciesNames
(
    const dictionary& thermoDict,
    const fvMesh& mesh
)
{
    // Legacy fallback: old style species list directly in thermophysicalProperties
    if (thermoDict.found("species"))
    {
        return thermoDict.get<wordList>("species");
    }

    const bool chemistry = thermoDict.get<bool>("chemistry");

    const fileName listFile =
        chemistry
      ? resolveFile(thermoDict, mesh, "reactionsList")
      : resolveFile(thermoDict, mesh, "speciesList");

    IFstream listStream(listFile);

    if (!listStream.good())
    {
        FatalIOErrorInFunction(listStream)
            << "Cannot open dictionary file " << listFile
            << Foam::exit(Foam::FatalIOError);
    }

    const dictionary listDict(listStream);

    if (!listDict.found("species"))
    {
        FatalIOErrorInFunction(listDict)
            << "Entry 'species' not found in dictionary " << listFile
            << exit(FatalIOError);
    }

    return listDict.get<wordList>("species");
}
template<class ThermoType>
const ThermoType& Foam::multiComponent2TMixture<ThermoType>::constructSpeciesData
(
    const dictionary& thermoDict
)
{
    // Legacy fallback: old style inline species definitions
    if (!thermoDict.found("speciesList"))
    {
        forAll(species_, i)
        {
            speciesData_.set
            (
                i,
                new ThermoType(thermoDict.subDict(species_[i]))
            );
        }

        return speciesData_[0];
    }

    const fileName speciesFile = resolveFile(thermoDict, mesh_, "speciesList");
    IFstream speciesStream(speciesFile);

    if (!speciesStream.good())
    {
        FatalIOErrorInFunction(speciesStream)
            << "Cannot open dictionary file " << speciesFile
            << Foam::exit(Foam::FatalIOError);
    }

    const dictionary speciesDict(speciesStream);

    forAll(species_, i)
    {
        if (!speciesDict.found(species_[i]))
        {
            FatalIOErrorInFunction(speciesDict)
                << "Species '" << species_[i]
                << "' not found in species dictionary " << speciesFile
                << exit(FatalIOError);
        }

        speciesData_.set
        (
            i,
            new ThermoType(speciesDict.subDict(species_[i]))
        );
    }

    return speciesData_[0];
}

template<class ThermoType>
void Foam::multiComponent2TMixture<ThermoType>::correctMassFractions()
{
    // Multiplication by 1.0 changes Yt patches to "calculated"
    volScalarField Yt("Yt", 1.0*Y_[0]);

    for (label n=1; n<Y_.size(); n++)
    {
        Yt += Y_[n];
    }

    if (mag(min(Yt).value()) < ROOTVSMALL)
    {
        FatalErrorInFunction
            << "Sum of mass fractions is zero for species " << this->species()
            << nl
            << incrIndent << indent
            << "Min of mass fraction sum " << min(Yt).value()
            << decrIndent
            << exit(FatalError);
    }

    const scalar diff(mag(max(Yt).value()) - scalar(1));

    if (diff > ROOTVSMALL)
    {
        WarningInFunction
            << "Sum of mass fractions is different from one for species "
            << this->species()
            << nl << incrIndent << indent
            << "Max of mass fraction sum differs from 1 by " << diff
            << decrIndent << nl;
    }

    forAll(Y_, n)
    {
        Y_[n] /= Yt;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::multiComponent2TMixture<ThermoType>::multiComponent2TMixture
(
    const dictionary& thermoDict,
    const wordList& specieNames,
    const ReactionTable<ThermoType>& thermoData,
    const fvMesh& mesh,
    const word& phaseName
)
:
    basicSpecie2TMixture(thermoDict, specieNames, mesh, phaseName),
    mesh_(mesh),
    speciesData_(species_.size()),
    mixture_("mixture", *thermoData[specieNames[0]]),
    mixtureVol_("volMixture", *thermoData[specieNames[0]])
{
    forAll(species_, i)
    {
        speciesData_.set
        (
            i,
            new ThermoType(*thermoData[species_[i]])
        );
    }

    correctMassFractions();
}


template<class ThermoType>
Foam::multiComponent2TMixture<ThermoType>::multiComponent2TMixture
(
    const dictionary& thermoDict,
    const fvMesh& mesh,
    const word& phaseName
)
:
    basicSpecie2TMixture
    (
        thermoDict,
        readSpeciesNames(thermoDict, mesh),
        mesh,
        phaseName
    ),
    mesh_(mesh),
    speciesData_(species_.size()),
    mixture_("mixture", constructSpeciesData(thermoDict)),
    mixtureVol_("volMixture", speciesData_[0])
{
    correctMassFractions();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ThermoType>
const ThermoType& Foam::multiComponent2TMixture<ThermoType>::cellMixture
(
    const label celli
) const
{
    mixture_ = Y_[0][celli]*speciesData_[0];

    for (label n=1; n<Y_.size(); n++)
    {
        mixture_ += Y_[n][celli]*speciesData_[n];
    }

    return mixture_;
}


template<class ThermoType>
const ThermoType& Foam::multiComponent2TMixture<ThermoType>::patchFaceMixture
(
    const label patchi,
    const label facei
) const
{
    mixture_ = Y_[0].boundaryField()[patchi][facei]*speciesData_[0];

    for (label n=1; n<Y_.size(); n++)
    {
        mixture_ += Y_[n].boundaryField()[patchi][facei]*speciesData_[n];
    }

    return mixture_;
}


template<class ThermoType>
const ThermoType& Foam::multiComponent2TMixture<ThermoType>::cellVolMixture
(
    const scalar p,
    const scalar T,
    const label celli
) const
{
    scalar rhoInv = 0.0;
    forAll(speciesData_, i)
    {
        rhoInv += Y_[i][celli]/speciesData_[i].rho(p, T);
    }

    mixtureVol_ =
        Y_[0][celli]/speciesData_[0].rho(p, T)/rhoInv*speciesData_[0];

    for (label n=1; n<Y_.size(); n++)
    {
        mixtureVol_ +=
            Y_[n][celli]/speciesData_[n].rho(p, T)/rhoInv*speciesData_[n];
    }

    return mixtureVol_;
}


template<class ThermoType>
const ThermoType& Foam::multiComponent2TMixture<ThermoType>::
patchFaceVolMixture
(
    const scalar p,
    const scalar T,
    const label patchi,
    const label facei
) const
{
    scalar rhoInv = 0.0;
    forAll(speciesData_, i)
    {
        rhoInv +=
            Y_[i].boundaryField()[patchi][facei]/speciesData_[i].rho(p, T);
    }

    mixtureVol_ =
        Y_[0].boundaryField()[patchi][facei]/speciesData_[0].rho(p, T)/rhoInv
      * speciesData_[0];

    for (label n=1; n<Y_.size(); n++)
    {
        mixtureVol_ +=
            Y_[n].boundaryField()[patchi][facei]/speciesData_[n].rho(p,T)
          / rhoInv*speciesData_[n];
    }

    return mixtureVol_;
}

template<class ThermoType>
void Foam::multiComponent2TMixture<ThermoType>::read
(
    const dictionary& thermoDict
)
{
    const wordList newSpecies = readSpeciesNames(thermoDict, mesh_);

    if (newSpecies.size() != species_.size())
    {
        FatalIOErrorInFunction(thermoDict)
            << "Species count changed during read(). "
            << "Old species list: " << species_ << nl
            << "New species list: " << newSpecies
            << exit(FatalIOError);
    }

    forAll(species_, i)
    {
        if (species_[i] != newSpecies[i])
        {
            FatalIOErrorInFunction(thermoDict)
                << "Species ordering changed during read(). "
                << "Old species list: " << species_ << nl
                << "New species list: " << newSpecies
                << exit(FatalIOError);
        }
    }

    if (!thermoDict.found("speciesList"))
    {
        forAll(species_, i)
        {
            speciesData_[i] = ThermoType(thermoDict.subDict(species_[i]));
        }

        return;
    }

    const fileName speciesFile = resolveFile(thermoDict, mesh_, "speciesList");
    IFstream speciesStream(speciesFile);

    if (!speciesStream.good())
    {
        FatalIOErrorInFunction(speciesStream)
            << "Cannot open dictionary file " << speciesFile
            << Foam::exit(Foam::FatalIOError);
    }

    const dictionary speciesDict(speciesStream);

    forAll(species_, i)
    {
        if (!speciesDict.found(species_[i]))
        {
            FatalIOErrorInFunction(speciesDict)
                << "Species '" << species_[i]
                << "' not found in species dictionary " << speciesFile
                << exit(FatalIOError);
        }

        speciesData_[i] = ThermoType(speciesDict.subDict(species_[i]));
    }
}

// ************************************************************************* //
