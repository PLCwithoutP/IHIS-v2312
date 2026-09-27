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
    Foam::neEulerImplicit

Description
    An Euler implicit solver for chemistry

SourceFiles
    neEulerImplicit.C

\*---------------------------------------------------------------------------*/

#ifndef neEulerImplicit_H
#define neEulerImplicit_H

#include "neChemistrySolver.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template <class Type>
class simpleMatrix;

/*---------------------------------------------------------------------------*\
                       Class neEulerImplicit Declaration
\*---------------------------------------------------------------------------*/

template<class ChemistryModel>
class neEulerImplicit
:
    public neChemistrySolver<ChemistryModel>
{
    // Private data

        //- Coefficients dictionary
        dictionary coeffsDict_;


        // Model constants

            //- Chemistry timescale
            scalar cTauChem_;

            //- Equilibrium rate limiter flag (on/off)
            Switch eqRateLimiter_;

        // Solver data
        mutable scalarField cTp_;


public:

    //- Runtime type information
    TypeName("neEulerImplicit");


    // Constructors

        //- Construct from thermo
        neEulerImplicit(typename ChemistryModel::ne2TReactionThermo& thermo);


    //- Destructor
    virtual ~neEulerImplicit();


    // Member Functions

    void updateRRInReactionI
    (
        const label index,
        const scalar pr,
        const scalar pf,
        const scalar corr,
        const label lRef,
        const label rRef,
        const scalar p,
        const scalar T,
        simpleMatrix<scalar>& RR
    ) const;

    //- Update the concentrations and return the chemical time
    virtual void solve
    (
        scalarField& c,
        scalar& T,
        scalar& p,
        scalar& deltaT,
        scalar& subDeltaT
    ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "neEulerImplicit.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
