/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2017 OpenFOAM Foundation
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
    Foam::sensibleCalculation2T

Group
    grpSpecieThermo

Description
    Thermodynamics mapping class to expose the sensible enthalpy functions.

\*---------------------------------------------------------------------------*/

#ifndef sensibleCalculation2T_H
#define sensibleCalculation2T_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class sensibleCalculation2T Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo2T>
class sensibleCalculation2T
{

public:

    // Constructors

        //- Construct
        sensibleCalculation2T()
        {}


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "sensibleCalculation2T";
        }

        // Fundamental properties

            static word energyName()
            {
                return "h";
            }

            //- Heat capacity at constant pressure [J/(kg K)]
            scalar CpvTR
            (
                const Thermo2T& thermo2T,
                const scalar p,
                const scalar TTR
            ) const
            {
                #ifdef __clang__
                // Using volatile to prevent compiler optimisations leading to
                // a sigfpe
                volatile const scalar cp = thermo2T.CpTR(p, TTR);
                return cp;
                #else
                return thermo2T.CpTR(p, TTR);
                #endif
            }

            //- CpTR/CpTR []
            scalar CpByCpvTR
            (
                const Thermo2T& thermo2T,
                const scalar p,
                const scalar TTR
            ) const
            {
                return 1;
            }

            //- Sensible enthalpy [J/kg]
            scalar H
            (
                const Thermo2T& thermo2T,
                const scalar p,
                const scalar TTR
            ) const
            {
                #ifdef __clang__
                // Using volatile to prevent compiler optimisations leading to
                // a sigfpe
                volatile const scalar hs = thermo2T.Hs(p, TTR);
                return hs;
                #else
                return thermo2T.Hs(p, TTR);
                #endif
            }

            //- Sensible translational internal energy [J/kg]
            scalar ET
            (
                const Thermo2T& thermo2T,
                const scalar p,
                const scalar TTR
            ) const
            {
                #ifdef __clang__
                // Using volatile to prevent compiler optimisations leading to
                // a sigfpe
                volatile const scalar esT = thermo2T.EsT(p, TTR);
                return esT;
                #else
                return thermo2T.EsT(p, TTR);
                #endif
            }

            //- Sensible rotational internal energy [J/kg]
            scalar ER
            (
                const Thermo2T& thermo2T,
                const scalar p,
                const scalar TTR
            ) const
            {
                #ifdef __clang__
                // Using volatile to prevent compiler optimisations leading to
                // a sigfpe
                volatile const scalar esR = thermo2T.EsR(p, TTR);
                return esR;
                #else
                return thermo2T.EsR(p, TTR);
                #endif
            }

            //- Sensible vibrational internal energy [J/kg]
            scalar EV
            (
                const Thermo2T& thermo2T,
                const scalar p,
                const scalar TTR,
                const scalar TVib,
                const scalar ThetaVib
            ) const
            {
                #ifdef __clang__
                // Using volatile to prevent compiler optimisations leading to
                // a sigfpe
                volatile const scalar esVib = thermo2T.EsVib(p, TTR, TVib, ThetaVib);
                return esVib;
                #else
                return thermo2T.EsVib(p, TTR, TVib, ThetaVib);
                #endif
            }

            //- Temperature from sensible enthalpy
            //  given an initial temperature T0
            scalar TH_TR
            (
                const Thermo2T& thermo2T,
                const scalar h,
                const scalar p,
                const scalar TTR0
            ) const
            {
                #ifdef __clang__
                // Using volatile to prevent compiler optimisations leading to
                // a sigfpe
                volatile const scalar ths = thermo2T.THs_TR(h, p, TTR0);
                return ths;
                #else
                return thermo2T.THs_TR(h, p, TTR0);
                #endif
            }

            //- Temperature from sensible internal energy
            //  given an initial temperature T0
            scalar TE_TR
            (
                const Thermo2T& thermo2T,
                const scalar eTR,
                const scalar p,
                const scalar TTR0
            ) const
            {
                #ifdef __clang__
                // Using volatile to prevent compiler optimisations leading to
                // a sigfpe
                volatile const scalar ths = thermo2T.TEs_TR(eTR, p, TTR0);
                return ths;
                #else
                return thermo2T.TEs_TR(eTR, p, TTR0);
                #endif
            }

            //- Temperature from sensible internal energy
            //  given an initial temperature T0
            scalar TE_Vib
            (
                const Thermo2T& thermo2T,
                const scalar eVib,
                const scalar p,
                const scalar TTR0,
                const scalar TVib0
            ) const
            {
                #ifdef __clang__
                // Using volatile to prevent compiler optimisations leading to
                // a sigfpe
                volatile const scalar tesVib = thermo2T.TEs_Vib(eVib, p, TVib0);
                return tesVib;
                #else
                return thermo2T.TEs_Vib(eVib, p, TTR0, TVib0);
                #endif
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
