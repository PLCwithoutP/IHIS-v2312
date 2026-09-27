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

Class
    Foam::pure2TMixture

Description
    Foam::pure2TMixture

SourceFiles
    pure2TMixture.C

\*---------------------------------------------------------------------------*/

#ifndef pure2TMixture_H
#define pure2TMixture_H

#include "basic2TMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class pure2TMixture Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class pure2TMixture
:
    public basic2TMixture
{
    // Private data

        ThermoType mixture_;

        //- Construct as copy (not implemented)
        pure2TMixture(const pure2TMixture<ThermoType>&);


public:

    //- The type of thermodynamics this mixture is instantiated for
    typedef ThermoType thermoType;


    // Constructors

        //- Construct from dictionary, mesh and phase name
        pure2TMixture(const dictionary&, const fvMesh&, const word&);


    // Member functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "pure2TMixture<" + ThermoType::typeName() + '>';
        }

        const ThermoType& mixture() const
        {
            return mixture_;
        }

        const ThermoType& cellMixture(const label) const
        {
            return mixture_;
        }

        const ThermoType& patchFaceMixture
        (
            const label,
            const label
        ) const
        {
            return mixture_;
        }

        const ThermoType& cellVolMixture
        (
            const scalar,
            const scalar,
            const label
        ) const
        {
            return mixture_;
        }

        const ThermoType& patchFaceVolMixture
        (
            const scalar,
            const scalar,
            const label,
            const label
        ) const
        {
            return mixture_;
        }

        //- Read dictionary
        void read(const dictionary&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

#ifdef NoRepository
    #include "pure2TMixture.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
