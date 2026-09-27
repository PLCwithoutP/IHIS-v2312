/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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

#include "specie2T.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

inline Foam::specie2T::specie2T
(
    const word& name,
    const scalar Y,
    const scalar molWeight,
    const speciesType kind
)
:
    name_(name),
    Y_(Y),
    molWeight_(molWeight),
    speciesType_(kind)
{}


inline Foam::specie2T::specie2T
(
    const scalar Y,
    const scalar molWeight,
    const speciesType kind
)
:
    Y_(Y),
    molWeight_(molWeight),
    speciesType_(kind)
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::specie2T::specie2T(const word& name, const specie2T& st)
:
    name_(name),
    Y_(st.Y_),
    molWeight_(st.molWeight_),
    speciesType_(st.speciesType_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::word& Foam::specie2T::name() const
{
    return name_;
}


inline Foam::scalar Foam::specie2T::W() const
{
    return molWeight_;
}


inline Foam::scalar Foam::specie2T::Y() const
{
    return Y_;
}


inline Foam::scalar Foam::specie2T::R() const
{
    return RR/molWeight_;
}


inline Foam::specie2T::speciesType Foam::specie2T::speciesKind() const
{
    return speciesType_;
}


inline bool Foam::specie2T::isMolecular() const
{
    return speciesType_ == tMolecular;
}


inline bool Foam::specie2T::isAtomic() const
{
    return speciesType_ == tAtomic;
}
// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline void Foam::specie2T::operator=(const specie2T& st)
{
    // Preserve original name
    Y_ = st.Y_;
    molWeight_ = st.molWeight_;
    speciesType_ = st.speciesType_;
}


inline void Foam::specie2T::operator+=(const specie2T& st)
{
    const scalar sumY = Y_ + st.Y_;
    if (mag(sumY) > SMALL)
    {
        molWeight_ = sumY/(Y_/molWeight_ + st.Y_/st.molWeight_);
    }

    Y_ = sumY;
}


inline void Foam::specie2T::operator*=(const scalar s)
{
    Y_ *= s;
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

inline Foam::specie2T Foam::operator+(const specie2T& st1, const specie2T& st2)
{
    const scalar sumY = max(st1.Y_ + st2.Y_, SMALL);

    if (mag(sumY) > SMALL)
    {
        return specie2T
        (
            sumY,
            sumY/(st1.Y_/st1.molWeight_ + st2.Y_/st2.molWeight_),
            st1.speciesType_
        );
    }
    else
    {
        return st1;
    }
}


inline Foam::specie2T Foam::operator*(const scalar s, const specie2T& st)
{
    return specie2T
    (
        s*st.Y_,
        st.molWeight_,
        st.speciesType_
    );
}


inline Foam::specie2T Foam::operator==(const specie2T& st1, const specie2T& st2)
{
    scalar diffY = st2.Y_ - st1.Y_;
    if (mag(diffY) < SMALL)
    {
        diffY = SMALL;
    }

    const scalar diffRW = st2.Y_/st2.molWeight_ - st1.Y_/st1.molWeight_;
    const Foam::specie2T::speciesType kind = st1.speciesType_;
    #ifdef __clang__
    // Using intermediate volatile bool to prevent compiler optimising out the
    // if block (above) - CLANG 3.7.1
    volatile const bool valid = (mag(diffRW) > SMALL);
    const scalar molWeight = valid ? diffY/diffRW : GREAT;
    #else
    scalar molWeight = GREAT;
    if (mag(diffRW) > SMALL)
    {
         molWeight = diffY/diffRW;
    }
    #endif

    return specie2T(diffY, molWeight, kind);
}


// ************************************************************************* //
