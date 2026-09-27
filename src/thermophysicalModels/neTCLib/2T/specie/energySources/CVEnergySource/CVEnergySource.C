/*---------------------------------------------------------------------------*\\
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

\*---------------------------------------------------------------------------*/

#include "CVEnergySource.H"
#include "fvMesh.H"
#include "error.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class MixtureType, class MixingRule>
Foam::CVEnergySource<MixtureType, MixingRule>::CVEnergySource
(
    MixtureType& mixture,
    MixingRule& MR
)
:
    baseEnergySource<MixtureType, MixingRule>(mixture, MR),
    mix_(mixture),
    mr_(MR),
    names_(mixture.species()),
    Q_CVList_(),
    CVTreatment_("non-preferential"),
    alphaDissociation_(0.3),
    dissociationPotential_(0),
    modelDataLoaded_(false)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class MixtureType, class MixingRule>
void Foam::CVEnergySource<MixtureType, MixingRule>::readModelData
(
    const fvMesh& mesh
)
{
    if (modelDataLoaded_)
    {
        return;
    }

    IOdictionary thermoDict
    (
        IOobject
        (
            "thermophysicalProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    IOdictionary speciesDict
    (
        IOobject
        (
            "speciesDict",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    if (!thermoDict.found("CVTreatment"))
    {
        FatalIOErrorInFunction(thermoDict)
            << "Entry 'CVTreatment' was not found in thermophysicalProperties." << nl
            << "Valid values are 'preferential' and 'non-preferential'."
            << exit(FatalIOError);
    }

    CVTreatment_ = word(thermoDict.lookup("CVTreatment"));

    if
    (
        CVTreatment_ != "preferential"
     && CVTreatment_ != "non-preferential"
    )
    {
        FatalIOErrorInFunction(thermoDict)
            << "Invalid CVTreatment = '" << CVTreatment_ << "'." << nl
            << "Valid values are 'preferential' and 'non-preferential'."
            << exit(FatalIOError);
    }

    thermoDict.readIfPresent("alphaDissociation", alphaDissociation_);

    dissociationPotential_.setSize(names_.size());

    forAll(names_, i)
    {
        dissociationPotential_[i] = 0.0;

        if (!mix_.isSpecieMolecular(i))
        {
            continue;
        }

        if (!speciesDict.found(names_[i]))
        {
            FatalIOErrorInFunction(speciesDict)
                << "Species sub-dictionary '" << names_[i] << "' was not found in speciesDict."
                << exit(FatalIOError);
        }

        const dictionary& specieDict = speciesDict.subDict(names_[i]);

        if (!specieDict.found("dissociationPotential"))
        {
            FatalIOErrorInFunction(speciesDict)
                << "Entry 'dissociationPotential' was not found for species '"
                << names_[i] << "' in speciesDict."
                << exit(FatalIOError);
        }

        dissociationPotential_[i] =
            readScalar(specieDict.lookup("dissociationPotential"));
    }

    modelDataLoaded_ = true;
}


template<class MixtureType, class MixingRule>
void Foam::CVEnergySource<MixtureType, MixingRule>::makeQCVSourceFields
(
    const fvMesh& mesh
)
{
    if (Q_CVList_.size())
    {
        return;
    }

    const PtrList<volScalarField>& Y_species = mix_.Y();
    Q_CVList_.setSize(Y_species.size());

    forAll(Y_species, i)
    {
        const volScalarField& Yi = Y_species[i];
        const word fieldName("Q_CV_" + Yi.name());

        Q_CVList_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    fieldName,
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar
                (
                    "zero",
                    dimEnergy/dimVolume/dimTime,
                    0.0
                )
            )
        );
    }
}


template<class MixtureType, class MixingRule>
Foam::scalar Foam::CVEnergySource<MixtureType, MixingRule>::p_s
(
    const scalar p,
    const scalar,
    const label s,
    const label celli
)
{
    const scalar Xs = mr_.Xi(s, celli); 
    
    return Xs*p;
}


template<class MixtureType, class MixingRule>
Foam::scalar Foam::CVEnergySource<MixtureType, MixingRule>::eVib_s
(
    const scalar p,
    const scalar TTR,
    const scalar TVib,
    const label s,
    const label celli
)
{
    const scalar ps = p_s(p, TTR, s, celli);

    return mix_.EsVib
    (
        s,
        ps,
        TTR,
        TVib,
        thetai(s)
    );
}


template<class MixtureType, class MixingRule>
Foam::scalar Foam::CVEnergySource<MixtureType, MixingRule>::eElec_s
(
    const label,
    const label
) const
{
    // Current uploaded thermo stack has no separate electronic-energy model.
    return 0.0;
}


template<class MixtureType, class MixingRule>
Foam::scalar Foam::CVEnergySource<MixtureType, MixingRule>::Dm
(
    const label s
) const
{
    return dissociationPotential_[s];
}


template<class MixtureType, class MixingRule>
Foam::scalar Foam::CVEnergySource<MixtureType, MixingRule>::DmPrime
(
    const scalar p,
    const scalar TTR,
    const scalar TVib,
    const label s,
    const label celli
)
{
    if (!mix_.isSpecieMolecular(s))
    {
        return 0.0;
    }

    if (CVTreatment_ == "preferential")
    {
        return alphaDissociation_*Dm(s);
    }

    return eVib_s(p, TTR, TVib, s, celli);
}


template<class MixtureType, class MixingRule>
Foam::scalar Foam::CVEnergySource<MixtureType, MixingRule>::Q_CV_s
(
    const scalar omegaDot_s,
    const scalar p,
    const scalar TTR,
    const scalar TVib,
    const label s,
    const label celli
)
{
    if (!mix_.isSpecieMolecular(s))
    {
        return 0.0;
    }

    return omegaDot_s
         *
        (
            DmPrime(p, TTR, TVib, s, celli)
          + eElec_s(s, celli)
        );
}


template<class MixtureType, class MixingRule>
Foam::PtrList<Foam::volScalarField>&
Foam::CVEnergySource<MixtureType, MixingRule>::correctCVSource
(
    const PtrList<volScalarField::Internal>& RR,
    const volScalarField& p,
    const volScalarField& TTR,
    const volScalarField& TVib
)
{
    const fvMesh& mesh = TTR.mesh();

    readModelData(mesh);
    makeQCVSourceFields(mesh);
    mr_.precomputeXi();
    
    const PtrList<volScalarField>& Y_species = mix_.Y();

    const scalarField& pCells    = p.primitiveField();
    const scalarField& TTRCells  = TTR.primitiveField();
    const scalarField& TVibCells = TVib.primitiveField();

    forAll(Y_species, s)
    {
        scalarField& Q_CV_Cells = Q_CVList_[s].primitiveFieldRef();
        const scalarField& omegaDotCells = RR[s].field();

        forAll(TTRCells, celli)
        {
            Q_CV_Cells[celli] = Q_CV_s
            (
                omegaDotCells[celli],
                pCells[celli],
                TTRCells[celli],
                TVibCells[celli],
                s,
                celli
            );
        }

        Q_CVList_[s].correctBoundaryConditions();
    }

    return Q_CVList_;
}


// ************************************************************************* //