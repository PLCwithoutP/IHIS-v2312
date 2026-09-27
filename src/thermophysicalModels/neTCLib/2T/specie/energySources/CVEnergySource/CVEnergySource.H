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

Class
    Foam::CVEnergySource

Description
    Chemical-vibrational source term model.

    Implements
        Q_CV,m = omegaDot_m * (D'_m + e_el,m)

    with
        non-preferential: D'_m = e_v,m
        preferential    : D'_m = alpha_m * D_m

    The chemistry source term omegaDot_m is taken from the chemistry model
    species source RR(m), which StandardNeChemistryModel stores as
    [kg/m^3/s].

    In the current 2T stack, e_el,m is taken as zero because no separate
    electronic-energy model/field exists in the uploaded code.

SourceFiles
    CVEnergySource.C

\*---------------------------------------------------------------------------*/

#ifndef CVEnergySource_H
#define CVEnergySource_H

#include "baseEnergySource.H"
#include "speciesTable.H"
#include "wordList.H"
#include "dimensionedScalar.H"
#include "volFields.H"
#include "IOobject.H"
#include "IOdictionary.H"
#include "PtrList.H"

namespace Foam
{

template<class MixtureType, class MixingRule>
class CVEnergySource
:
    public baseEnergySource<MixtureType, MixingRule>
{
        //- Construct as copy (not implemented)
        CVEnergySource(const CVEnergySource&);

private:

    MixtureType& mix_;
    MixingRule& mr_;

    const speciesTable& names_;
    PtrList<volScalarField> Q_CVList_;

    // Model controls read lazily from dictionaries
    word CVTreatment_;
    scalar alphaDissociation_;
    List<scalar> dissociationPotential_;
    bool modelDataLoaded_;

    void readModelData
    (
        const fvMesh& mesh
    );

public:

    CVEnergySource
    (
        MixtureType& mix,
        MixingRule& MR
    );

    static word typeName()
    {
        return "CVEnergySource";
    }

    void makeQCVSourceFields
    (
        const fvMesh& mesh
    );

    label nSpecies() const
    {
        return names_.size();
    }

    const word& specieName(const label i) const
    {
        return names_[i];
    }

    const wordList& specieNames() const
    {
        return names_;
    }

    const MixtureType& mixture() const
    {
        return mix_;
    }

    const volScalarField& Yi(const label i) const
    {
        return mix_.Y(i);
    }

    scalar Wi(const label i) const
    {
        return mix_.W(i);
    }

    scalar thetai(const label i) const
    {
        return mix_.thetaVib(i);
    }

    scalar p_s
    (
        const scalar p,
        const scalar TTR,
        const label s,
        const label celli
    );

    scalar eVib_s
    (
        const scalar p,
        const scalar TTR,
        const scalar TVib,
        const label s,
        const label celli
    );

    scalar eElec_s
    (
        const label s,
        const label celli
    ) const;

    scalar Dm
    (
        const label s
    ) const;

    scalar DmPrime
    (
        const scalar p,
        const scalar TTR,
        const scalar TVib,
        const label s,
        const label celli
    );

    scalar Q_CV_s
    (
        const scalar omegaDot_s,
        const scalar p,
        const scalar TTR,
        const scalar TVib,
        const label s,
        const label celli
    );

    PtrList<volScalarField>& correctCVSource
    (
        const PtrList<volScalarField::Internal>& RR,
        const volScalarField& p,
        const volScalarField& TTR,
        const volScalarField& TVib
    );

    PtrList<volScalarField>& fields()
    {
        return Q_CVList_;
    }

    const PtrList<volScalarField>& fields() const
    {
        return Q_CVList_;
    }
};

} // End namespace Foam

#ifdef NoRepository
    #include "CVEnergySource.C"
#endif

#endif

// ************************************************************************* //