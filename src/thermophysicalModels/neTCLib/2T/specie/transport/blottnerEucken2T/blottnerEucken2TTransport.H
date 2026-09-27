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
    Foam::blottnerEucken2TTransport

Group
    grpSpecieTransport

Description
    Transport package using Blottner's viscosity formula and Eucken's
    thermal conductivity formula.

    Templated into a given thermodynamics package (needed for thermal
    conductivity).

    Dynamic viscosity [kg/m.s]
    \f[
        \mu_s = 0.1 * exp((A_B*ln(T) + B_B)*ln(T) + C_B)
    \f]
    T is translational temperature.

    Thermal conductivity [W/m.K]
    \f[
        \kappa_tr = \mu_s (\frac{5}{2}*c_{v,t,s} + c_{v,r,s})
        \kappa_v = \mu_s * c_{v,vib,s}(T_v) / Sc_s
    \f]
    Holds until 6000 K.

SourceFiles
    blottnerEuckenTTransportI.H
    blottnerEucken2TTransport.C

\*---------------------------------------------------------------------------*/

#ifndef blottnerEuckenTTransport_H
#define blottnerEuckenTTransport_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations

template<class Thermo2T> class blottnerEucken2TTransport;

template<class Thermo2T>
inline blottnerEucken2TTransport<Thermo2T> operator+
(
    const blottnerEucken2TTransport<Thermo2T>&,
    const blottnerEucken2TTransport<Thermo2T>&
);

template<class Thermo2T>
inline blottnerEucken2TTransport<Thermo2T> operator*
(
    const scalar,
    const blottnerEucken2TTransport<Thermo2T>&
);

template<class Thermo2T>
Ostream& operator<<
(
    Ostream&,
    const blottnerEucken2TTransport<Thermo2T>&
);


/*---------------------------------------------------------------------------*\
                    Class blottnerEucken2TTransport Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo2T>
class blottnerEucken2TTransport
:
    public Thermo2T
{
    // Private Data

        // Blottner's coefficients
        scalar AB_, BB_, CB_;

        scalar Sc_;
    // Private Member Functions

        //- Read coefficient from dictionary
        scalar readCoeff(const word& coeffName, const dictionary& dict);


public:

    // Generated Methods: copy construct, copy assignment


    // Constructors

        //- Construct from components
        inline blottnerEucken2TTransport
        (
            const Thermo2T& t,
            const scalar AB,
            const scalar BB,
            const scalar CB,
            const scalar Sc
        );

        //- Construct as named copy
        inline blottnerEucken2TTransport(const word&, const blottnerEucken2TTransport&);

        //- Construct from dictionary
        explicit blottnerEucken2TTransport(const dictionary& dict);

        //- Construct from base thermo and dictionary
        blottnerEucken2TTransport(const Thermo2T& t, const dictionary& dict);

        //- Construct and return a clone
        inline autoPtr<blottnerEucken2TTransport> clone() const;

        // Selector from dictionary
        inline static autoPtr<blottnerEucken2TTransport> New(const dictionary& dict);


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "blottnerEucken2T<" + Thermo2T::typeName() + '>';
        }

        //- Dynamic viscosity [kg/ms]
        inline scalar mu(const scalar TTR) const;

        //- Thermal translational-rotational conductivity [W/mK]
        inline scalar kappaTR(const scalar p, const scalar TTR) const;

        //- Thermal vibrational conductivity [W/mK]
        inline scalar kappaVib
        (
            const scalar p,
            const scalar TTR,
            const scalar TVib,
            const scalar ThetaVib
        ) const;

        //- Thermal diffusivity of enthalpy [kg/ms]
        inline scalar alphah(const scalar p, const scalar TTR) const;
        
        //- Write to Ostream
        void write(Ostream& os) const;


    // Member Operators

        inline void operator+=(const blottnerEucken2TTransport&);

        inline void operator*=(const scalar);


    // Friend Operators

        friend blottnerEucken2TTransport operator+ <Thermo2T>
        (
            const blottnerEucken2TTransport&,
            const blottnerEucken2TTransport&
        );

        friend blottnerEucken2TTransport operator* <Thermo2T>
        (
            const scalar,
            const blottnerEucken2TTransport&
        );


    // IOstream Operators

        friend Ostream& operator<< <Thermo2T>
        (
            Ostream&,
            const blottnerEucken2TTransport&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "blottnerEucken2TTransportI.H"

#ifdef NoRepository
    #include "blottnerEucken2TTransport.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
