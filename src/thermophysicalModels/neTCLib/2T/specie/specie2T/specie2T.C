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
#include "constants.H"
#include "Enum.H"

/* * * * * * * * * * * * * * * public constants  * * * * * * * * * * * * * * */

namespace Foam
{
    defineTypeNameAndDebug(specie2T, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::specie2T::specie2T(const dictionary& dict)
:
    name_(dict.dictName()),
    Y_(dict.subDict("specie").getOrDefault<scalar>("massFraction", 1)),
    molWeight_(dict.subDict("specie").get<scalar>("molWeight")),
    speciesType_(speciesTypeNames.get(dict.subDict("specie").lookupOrDefault<word>("type", "molecular")))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::specie2T::write(Ostream& os) const
{
    // Entries in dictionary format
    {
        os.beginBlock("specie");
        os.writeEntryIfDifferent<scalar>("massFraction", 1, Y_);
        os.writeEntry("molWeight", molWeight_);
        os.writeEntry("type", speciesType_);
        os.endBlock();
    }
}

const Foam::Enum<Foam::specie2T::speciesType> Foam::specie2T::speciesTypeNames
{
    
    {Foam::specie2T::tMolecular, "molecular"},
    {Foam::specie2T::tAtomic, "atomic"}

};

// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const specie2T& st)
{
    st.write(os);
    os.check(FUNCTION_NAME);
    return os;
}


// ************************************************************************* //
