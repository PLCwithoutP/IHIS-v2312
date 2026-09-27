/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2021 OpenCFD Ltd.
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

#include "neIrreversibleReaction.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template
<
    template<class> class ReactionType,
    class ReactionThermo,
    class ReactionRate
>
Foam::neIrreversibleReaction<ReactionType, ReactionThermo, ReactionRate>::
neIrreversibleReaction
(
    const ReactionType<ReactionThermo>& reaction,
    const ReactionRate& k
)
:
    ReactionType<ReactionThermo>(reaction),
    k_(k),
    alphaPark_(-1.0),
    twoTemperature_(false)
{}


template
<
    template<class> class ReactionType,
    class ReactionThermo,
    class ReactionRate
>
Foam::neIrreversibleReaction<ReactionType, ReactionThermo, ReactionRate>::
neIrreversibleReaction
(
    const speciesTable& species,
    const ReactionTable<ReactionThermo>& thermoDatabase,
    const dictionary& dict
)
:
    ReactionType<ReactionThermo>(species, thermoDatabase, dict),
    k_(species, dict),
    alphaPark_(dict.getOrDefault<scalar>("alphaPark", -1.0)),
    twoTemperature_(dict.getOrDefault<bool>("twoTemperature", false))
{}


template
<
    template<class> class ReactionType,
    class ReactionThermo,
    class ReactionRate
>
Foam::neIrreversibleReaction<ReactionType, ReactionThermo, ReactionRate>::
neIrreversibleReaction
(
    const neIrreversibleReaction<ReactionType, ReactionThermo,ReactionRate>& irr,
    const speciesTable& species
)
:
    ReactionType<ReactionThermo>(irr, species),
    k_(irr.k_),
    alphaPark_(irr.alphaPark_),
    twoTemperature_(irr.twoTemperature_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template
<
    template<class> class ReactionType,
    class ReactionThermo,
    class ReactionRate
>
Foam::scalar Foam::neIrreversibleReaction
<
    ReactionType,
    ReactionThermo,
    ReactionRate
>::kf
(
    const scalar p,
    const scalar T,
    const scalarField& c
) const
{
    return k_(p, T, c);
}


template
<
    template<class> class ReactionType,
    class ReactionThermo,
    class ReactionRate
>
void Foam::neIrreversibleReaction<ReactionType, ReactionThermo, ReactionRate>::
write
(
    Ostream& os
) const
{
    ReactionType<ReactionThermo>::write(os);
    k_.write(os);
}


// ************************************************************************* //
