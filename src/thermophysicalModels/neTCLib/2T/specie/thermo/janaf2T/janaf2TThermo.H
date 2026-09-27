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

Class
    Foam::janaf2TThermo

Group
    grpSpecieThermo

Description
    JANAF tables based thermodynamics package templated
    into the equation of state.

SourceFiles
    janafTThermoI.H
    janaf2TThermo.C

\*---------------------------------------------------------------------------*/

#ifndef janaf2TThermo_H
#define janaf2TThermo_H

#include "scalar.H"
#include "FixedList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations

template<class EquationOfState> class janaf2TThermo;

template<class EquationOfState>
inline janaf2TThermo<EquationOfState> operator+
(
    const janaf2TThermo<EquationOfState>&,
    const janaf2TThermo<EquationOfState>&
);

template<class EquationOfState>
inline janaf2TThermo<EquationOfState> operator*
(
    const scalar,
    const janaf2TThermo<EquationOfState>&
);

template<class EquationOfState>
inline janaf2TThermo<EquationOfState> operator==
(
    const janaf2TThermo<EquationOfState>&,
    const janaf2TThermo<EquationOfState>&
);

template<class EquationOfState>
Ostream& operator<<
(
    Ostream&,
    const janaf2TThermo<EquationOfState>&
);


/*---------------------------------------------------------------------------*\
                           Class janaf2TThermo Declaration
\*---------------------------------------------------------------------------*/

template<class EquationOfState>
class janaf2TThermo
:
    public EquationOfState
{
public:

    // Public Data

        static constexpr int nCoeffs_ = 9;
        typedef FixedList<scalar, nCoeffs_> coeffArray;


private:

    // Private Data

        // Temperature limits of applicability of functions
        scalar TAbsLow_, TAbsHigh_, TFirstInt_, TSecondInt_; 
        scalar ThetaVib_;

        coeffArray highCpCoeffs_;
        coeffArray midCpCoeffs_;
        coeffArray lowCpCoeffs_;

        bool convertCoeffs_;
        bool warnOutOfRange_;   // suppress JANAF out-of-range warnings if false

    // Private Member Functions

        //- Check that input data is valid
        void checkInputData() const;

        //- Return the coefficients corresponding to the given temperature
        inline const coeffArray& coeffs(const scalar TTR) const;


public:

    // Constructors

        //- Construct from components
        inline janaf2TThermo
        (
            const EquationOfState& st,
            const scalar TAbsLow,
            const scalar TAbsHigh,
            const scalar TFirstInt,
            const scalar TSecondInt,
            const scalar ThetaVib,
            const coeffArray& highCpCoeffs,
            const coeffArray& midCpCoeffs,
            const coeffArray& lowCpCoeffs,
            const bool convertCoeffs = false,
            const bool warnOutOfRange = true
        );

        //- Construct from dictionary
        janaf2TThermo(const dictionary& dict);

        //- Construct as a named copy
        inline janaf2TThermo(const word&, const janaf2TThermo&);


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "janaf2T<" + EquationOfState::typeName() + '>';
        }

        //- Limit the temperature to be in the range TAbsLow_ to TAbsHigh_
        inline scalar limit(const scalar TTR) const;


        // Access

            //- Return const access to the lowest possible temperature limit
            inline scalar TAbsLow() const;

            //- Return const access to the highest possible temperature limit
            inline scalar TAbsHigh() const;
            
            //- Return const access to the first interval temperature
            inline scalar TFirstInt() const;

            //- Return const access to the second interval temperature
            inline scalar TSecondInt() const;

            //- Return const access to the vibrational temperature coefficient
            inline scalar ThetaVib() const;

            //- Return const access to the high temperature poly coefficients
            inline const coeffArray& highCpCoeffs() const;

            //- Return const access to the mid temperature poly coefficients
            inline const coeffArray& midCpCoeffs() const;

            //- Return const access to the low temperature poly coefficients
            inline const coeffArray& lowCpCoeffs() const;


        // Fundamental properties

            //- Heat capacity at constant pressure [J/(kg K)]
            inline scalar CpTR(const scalar p, const scalar TTR) const;

            //- Absolute Enthalpy [J/kg]
            inline scalar Ha(const scalar p, const scalar TTR) const;

            //- Sensible enthalpy [J/kg]
            inline scalar Hs(const scalar p, const scalar TTR) const;

            //- Chemical enthalpy [J/kg]
            inline scalar Hc() const;

            //- Entropy [J/(kg K)]
            inline scalar S(const scalar p, const scalar TTR) const;

            //- Gibbs free energy of the mixture in the standard state [J/kg]
            inline scalar Gstd(const scalar p, const scalar TTR) const;

            #include "HtoEthermo2T.H"


        // Derivative term used for Jacobian


            //- Temperature derivative of heat capacity at constant pressure
            inline scalar dCpdT(const scalar p, const scalar TTR) const;


        // I-O

            //- Write to Ostream
            void write(Ostream& os) const;


    // Member operators

        inline void operator+=(const janaf2TThermo&);


    // Friend operators

        friend janaf2TThermo operator+ <EquationOfState>
        (
            const janaf2TThermo&,
            const janaf2TThermo&
        );

        friend janaf2TThermo operator* <EquationOfState>
        (
            const scalar,
            const janaf2TThermo&
        );

        friend janaf2TThermo operator== <EquationOfState>
        (
            const janaf2TThermo&,
            const janaf2TThermo&
        );


    // Ostream Operator

        friend Ostream& operator<< <EquationOfState>
        (
            Ostream&,
            const janaf2TThermo&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "janaf2TThermoI.H"

#ifdef NoRepository
    #include "janaf2TThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
