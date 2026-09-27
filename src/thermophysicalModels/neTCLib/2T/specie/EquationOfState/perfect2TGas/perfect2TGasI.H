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

#include "perfect2TGas.H"
#include "thermodynamicConstants.H"

using namespace Foam::constant::thermodynamic;

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Specie2T>
inline Foam::perfect2TGas<Specie2T>::perfect2TGas(const Specie2T& sp)
:
    Specie2T(sp)
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie2T>
inline Foam::perfect2TGas<Specie2T>::perfect2TGas
(
    const word& name,
    const perfect2TGas<Specie2T>& pg
)
:
    Specie2T(name, pg)
{}


template<class Specie2T>
inline Foam::autoPtr<Foam::perfect2TGas<Specie2T>>
Foam::perfect2TGas<Specie2T>::clone() const
{
    return autoPtr<perfect2TGas<Specie2T>>::New(*this);
}


template<class Specie2T>
inline Foam::autoPtr<Foam::perfect2TGas<Specie2T>> Foam::perfect2TGas<Specie2T>::New
(
    const dictionary& dict
)
{
    return autoPtr<perfect2TGas<Specie2T>>::New(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Specie2T>
inline Foam::scalar Foam::perfect2TGas<Specie2T>::rho(scalar p, scalar T) const
{
    return p/(this->R()*T);
}


template<class Specie2T>
inline Foam::scalar Foam::perfect2TGas<Specie2T>::HT(scalar p, scalar T) const
{
    return 0;
}


template<class Specie2T>
inline Foam::scalar Foam::perfect2TGas<Specie2T>::CpT(scalar p, scalar T) const
{
    return 0;
}


template<class Specie2T>
inline Foam::scalar Foam::perfect2TGas<Specie2T>::ET(scalar p, scalar T) const
{
    return 0;
}

template<class Specie2T>
inline Foam::scalar Foam::perfect2TGas<Specie2T>::ER(scalar p, scalar T) const
{
    return 0;
}

template<class Specie2T>
inline Foam::scalar Foam::perfect2TGas<Specie2T>::EV(scalar p, scalar T) const
{
    return 0;
}

template<class Specie2T>
inline Foam::scalar Foam::perfect2TGas<Specie2T>::CvT(scalar p, scalar T) const
{
    return 0;
}

template<class Specie2T>
inline Foam::scalar Foam::perfect2TGas<Specie2T>::CvR(scalar p, scalar T) const
{
    return 0;
}

template<class Specie2T>
inline Foam::scalar Foam::perfect2TGas<Specie2T>::CvVib(scalar p, scalar T) const
{
    return 0;
}

template<class Specie2T>
inline Foam::scalar Foam::perfect2TGas<Specie2T>::S(scalar p, scalar T) const
{
    return -this->R()*log(p/Pstd);
}


template<class Specie2T>
inline Foam::scalar Foam::perfect2TGas<Specie2T>::psi(scalar p, scalar T) const
{
    return 1.0/(this->R()*T);
}


template<class Specie2T>
inline Foam::scalar Foam::perfect2TGas<Specie2T>::Z(scalar p, scalar T) const
{
    return 1;
}


template<class Specie2T>
inline Foam::scalar Foam::perfect2TGas<Specie2T>::CpMCv(scalar p, scalar T) const
{
    return this->R();
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Specie2T>
inline void Foam::perfect2TGas<Specie2T>::operator+=(const perfect2TGas<Specie2T>& pg)
{
    Specie2T::operator+=(pg);
}


template<class Specie2T>
inline void Foam::perfect2TGas<Specie2T>::operator*=(const scalar s)
{
    Specie2T::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Specie2T>
inline Foam::perfect2TGas<Specie2T> Foam::operator+
(
    const perfect2TGas<Specie2T>& pg1,
    const perfect2TGas<Specie2T>& pg2
)
{
    return perfect2TGas<Specie2T>
    (
        static_cast<const Specie2T&>(pg1) + static_cast<const Specie2T&>(pg2)
    );
}


template<class Specie2T>
inline Foam::perfect2TGas<Specie2T> Foam::operator*
(
    const scalar s,
    const perfect2TGas<Specie2T>& pg
)
{
    return perfect2TGas<Specie2T>(s*static_cast<const Specie2T&>(pg));
}


template<class Specie2T>
inline Foam::perfect2TGas<Specie2T> Foam::operator==
(
    const perfect2TGas<Specie2T>& pg1,
    const perfect2TGas<Specie2T>& pg2
)
{
    return perfect2TGas<Specie2T>
    (
        static_cast<const Specie2T&>(pg1) == static_cast<const Specie2T&>(pg2)
    );
}


// ************************************************************************* //
