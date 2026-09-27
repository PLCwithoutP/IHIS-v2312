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

#include "blottnerEucken2TTransport.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Thermo2T>
Foam::scalar Foam::blottnerEucken2TTransport<Thermo2T>::readCoeff
(
    const word& coeffName,
    const dictionary& dict
)
{
    return dict.subDict("transport").get<scalar>(coeffName);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo2T>
Foam::blottnerEucken2TTransport<Thermo2T>::blottnerEucken2TTransport(const dictionary& dict)
:
    Thermo2T(dict),
    AB_(readCoeff("AB", dict)),
    BB_(readCoeff("BB", dict)),
    CB_(readCoeff("CB", dict)),
    Sc_(dict.subDict("transport").getOrDefault<scalar>("Sc", 0.833))
{}


template<class Thermo2T>
Foam::blottnerEucken2TTransport<Thermo2T>::blottnerEucken2TTransport
(
    const Thermo2T& t,
    const dictionary& dict
)
:
    Thermo2T(t),
    AB_(readCoeff("AB", dict)),
    BB_(readCoeff("BB", dict)),
    CB_(readCoeff("CB", dict)),
    Sc_(dict.subDict("transport").getOrDefault<scalar>("Sc", 0.833))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo2T>
void Foam::blottnerEucken2TTransport<Thermo2T>::write(Ostream& os) const
{
    os.beginBlock(this->specie2T::name());

    Thermo2T::write(os);

    // Entries in dictionary format
    {
        os.beginBlock("transport");
        os.writeEntry("AB", AB_);
        os.writeEntry("BB", BB_);
        os.writeEntry("CB", CB_);
        os.writeEntry("Sc", Sc_);
        os.endBlock();
    }

    os.endBlock();
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Thermo2T>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const blottnerEucken2TTransport<Thermo2T>& st
)
{
    st.write(os);
    return os;
}


// ************************************************************************* //
