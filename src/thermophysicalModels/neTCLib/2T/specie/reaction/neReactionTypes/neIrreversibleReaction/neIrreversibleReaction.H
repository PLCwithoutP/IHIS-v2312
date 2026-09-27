/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019-2021 OpenCFD Ltd.
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
    Foam::neIrreversibleReaction

Group
    grpSpecieReactions

Description
    Simple extension of neReaction to handle reversible reactions using
    equilibrium thermodynamics.

SourceFiles
    neIrreversibleReaction.C

\*---------------------------------------------------------------------------*/

#ifndef neIrreversibleReaction_H
#define neIrreversibleReaction_H

#include "neReaction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class neIrreversibleReaction Declaration
\*---------------------------------------------------------------------------*/

template
<
    template<class> class ReactionType,
    class ReactionThermo,
    class ReactionRate
>
class neIrreversibleReaction
:
    public ReactionType<ReactionThermo>
{
    // Private Data

        ReactionRate    k_;
        scalar          alphaPark_;
        bool            twoTemperature_;

    // Private Member Functions

        //- No copy assignment
        void operator=
        (
            const neIrreversibleReaction
            <
                ReactionType,
                ReactionThermo,
                ReactionRate
            >&
        ) = delete;


public:

    //- Runtime type information
    TypeName("irreversible");


    // Constructors

        //- Construct from components
        neIrreversibleReaction
        (
            const ReactionType<ReactionThermo>& reaction,
            const ReactionRate& reactionRate
        );

        //- Construct as copy given new speciesTable
        neIrreversibleReaction
        (
            const neIrreversibleReaction
            <
                ReactionType,
                ReactionThermo,
                ReactionRate
            >&,
            const speciesTable& species
        );

        //- Construct from dictionary
        neIrreversibleReaction
        (
            const speciesTable& species,
            const ReactionTable<ReactionThermo>& thermoDatabase,
            const dictionary& dict
        );

        //- Construct and return a clone
        virtual autoPtr<neReaction<ReactionThermo>> clone() const
        {
            return autoPtr<neReaction<ReactionThermo>>
            (
                new neIrreversibleReaction
                <
                    ReactionType,
                    ReactionThermo,
                    ReactionRate
                >(*this)
            );
        }

        //- Construct and return a clone with new speciesTable
        virtual autoPtr<neReaction<ReactionThermo>> clone
        (
            const speciesTable& species
        ) const
        {
            return autoPtr<neReaction<ReactionThermo>>
            (
                new neIrreversibleReaction
                <
                    ReactionType,
                    ReactionThermo,
                    ReactionRate
                >
                (
                    *this,
                    species
                )
            );
        }


    //- Destructor
    virtual ~neIrreversibleReaction() = default;


    // Member Functions

        // neIrreversibleReaction rate coefficients

            //- Forward rate constant
            virtual scalar kf
            (
                const scalar p,
                const scalar T,
                const scalarField& c
            ) const;

        //- Write
        virtual void write(Ostream&) const;

        virtual scalar alphaPark() const
        {
            return alphaPark_;
        }

        virtual bool twoTemperature() const
        {
            return twoTemperature_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "neIrreversibleReaction.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
