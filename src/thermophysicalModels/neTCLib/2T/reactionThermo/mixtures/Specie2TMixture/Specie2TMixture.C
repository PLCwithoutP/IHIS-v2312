/*---------------------------------------------------------------------------*\
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

#include "Specie2TMixture.H"
#include "fvMesh.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Mixture2TType>
Foam::Specie2TMixture<Mixture2TType>::Specie2TMixture
(
    const dictionary& thermoDict,
    const fvMesh& mesh,
    const word& phaseName
)
:
    Mixture2TType
    (
        thermoDict,
        mesh,
        phaseName
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Mixture2TType>
Foam::scalar Foam::Specie2TMixture<Mixture2TType>::R
(
    const label speciei
) const
{
    return this->getLocalThermo(speciei).R();
}

template<class Mixture2TType>
Foam::scalar Foam::Specie2TMixture<Mixture2TType>::W
(
    const label speciei
) const
{
    return this->getLocalThermo(speciei).W();
}


template<class Mixture2TType>
Foam::scalar Foam::Specie2TMixture<Mixture2TType>::Hc(const label speciei) const
{
    return this->getLocalThermo(speciei).Hc();
}


template<class Mixture2TType>
Foam::scalar Foam::Specie2TMixture<Mixture2TType>::CpTR
(
    const label speciei,
    const scalar p,
    const scalar TTR
) const
{
    return this->getLocalThermo(speciei).CpTR(p, TTR);
}


template<class Mixture2TType>
Foam::scalar Foam::Specie2TMixture<Mixture2TType>::CvT
(
    const label speciei,
    const scalar p,
    const scalar TTR
) const
{
    return this->getLocalThermo(speciei).CvT(p, TTR);
}


template<class Mixture2TType>
Foam::scalar Foam::Specie2TMixture<Mixture2TType>::H
(
    const label speciei,
    const scalar p,
    const scalar TTR
) const
{
    return this->getLocalThermo(speciei).H(p, TTR);
}


template<class Mixture2TType>
Foam::scalar Foam::Specie2TMixture<Mixture2TType>::Ha
(
    const label speciei,
    const scalar p,
    const scalar TTR
) const
{
    return this->getLocalThermo(speciei).Ha(p, TTR);
}


template<class Mixture2TType>
Foam::scalar Foam::Specie2TMixture<Mixture2TType>::Hs
(
    const label speciei,
    const scalar p,
    const scalar TTR
) const
{
    return this->getLocalThermo(speciei).Hs(p, TTR);
}


template<class Mixture2TType>
Foam::scalar Foam::Specie2TMixture<Mixture2TType>::S
(
    const label speciei,
    const scalar p,
    const scalar TTR
) const
{
    return this->getLocalThermo(speciei).S(p, TTR);
}


template<class Mixture2TType>
Foam::scalar Foam::Specie2TMixture<Mixture2TType>::EsT
(
    const label speciei,
    const scalar p,
    const scalar TTR
) const
{
    return this->getLocalThermo(speciei).EsT(p, TTR);
}

template<class Mixture2TType>
Foam::scalar Foam::Specie2TMixture<Mixture2TType>::EsR
(
    const label speciei,
    const scalar p,
    const scalar TTR
) const
{
    return this->getLocalThermo(speciei).EsR(p, TTR);
}

template<class Mixture2TType>
Foam::scalar Foam::Specie2TMixture<Mixture2TType>::EsVib
(
    const label speciei,
    const scalar p,
    const scalar TTR,
    const scalar TVib,
    const scalar ThetaVib
) const
{
    return this->getLocalThermo(speciei).EsVib(p, TTR, TVib, ThetaVib);
}

template<class Mixture2TType>
Foam::scalar Foam::Specie2TMixture<Mixture2TType>::G
(
    const label speciei,
    const scalar p,
    const scalar TTR
) const
{
    return this->getLocalThermo(speciei).G(p, TTR);
}

template<class Mixture2TType>
Foam::scalar Foam::Specie2TMixture<Mixture2TType>::mu
(
    const label speciei,
    const scalar TTR
) const
{
    return this->getLocalThermo(speciei).mu(TTR);
}


template<class Mixture2TType>
Foam::scalar Foam::Specie2TMixture<Mixture2TType>::kappaTR
(
    const label speciei,
    const scalar p,
    const scalar TTR
) const
{
    return this->getLocalThermo(speciei).kappaTR(p, TTR);
}

template<class Mixture2TType>
Foam::scalar Foam::Specie2TMixture<Mixture2TType>::kappaVib
(
    const label speciei,
    const scalar p,
    const scalar TTR,
    const scalar TVib,
    const scalar ThetaVib
) const
{
    return this->getLocalThermo(speciei).kappaVib(p, TTR, TVib, ThetaVib);
}


template<class Mixture2TType>
Foam::scalar Foam::Specie2TMixture<Mixture2TType>::alphah
(
    const label speciei,
    const scalar p,
    const scalar TTR
) const
{
    return this->getLocalThermo(speciei).alphah(p, TTR);
}


template<class Mixture2TType>
Foam::scalar Foam::Specie2TMixture<Mixture2TType>::rho
(
    const label speciei,
    const scalar p,
    const scalar TTR
) const
{
    return this->getLocalThermo(speciei).rho(p, TTR);
}

template<class Mixture2TType>
Foam::scalar Foam::Specie2TMixture<Mixture2TType>::thetaVib
(
    const label speciei
) const
{
    return this->getLocalThermo(speciei).ThetaVib();
}

template<class Mixture2TType>
bool Foam::Specie2TMixture<Mixture2TType>::isSpecieMolecular
(
    const label speciei
) const
{
    return this->getLocalThermo(speciei).isMolecular();
}

template<class Mixture2TType>
Foam::scalar Foam::Specie2TMixture<Mixture2TType>::TEVib
(
    const label speciei,
    const scalar esVib,
    const scalar p,
    const scalar TTR,
    const scalar TVib0
) const
{
    return this->getLocalThermo(speciei).TE_Vib
    (
        esVib,
        p,
        TTR,
        TVib0
    );
}

// ************************************************************************* //
