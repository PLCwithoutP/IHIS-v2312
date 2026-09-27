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
    Foam::Specie2TMixture

Group
    grpReactionThermophysicalMixtures

Description
    Foam::Specie2TMixture

SourceFiles
    Specie2TMixture.C

\*---------------------------------------------------------------------------*/

#ifndef Specie2TMixture_H
#define Specie2TMixture_H

#include "scalar.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class dictionary;
class fvMesh;


/*---------------------------------------------------------------------------*\
                     Class Specie2TMixture Declaration
\*---------------------------------------------------------------------------*/

template<class Mixture2TType>
class Specie2TMixture
:
    public Mixture2TType
{

public:

    // Constructors

        //- Construct from dictionary, mesh and phase name
        Specie2TMixture
        (
            const dictionary& thermoDict,
            const fvMesh& mesh,
            const word& phaseName
        );


    //- Destructor
    virtual ~Specie2TMixture() = default;


    // Member functions

        // Per specie properties

            //- Molecular weight of the given specie [kg/kmol]
            virtual scalar W(const label speciei) const override;

            //- Specific gas constant of the given specie [kg/kmol]
            virtual scalar R(const label speciei) const override;

            //- Chemical enthalpy [J/kg]
            virtual scalar Hc(const label speciei) const override;


        // Per specie thermo properties

            //- Heat capacity at constant pressure [J/(kg K)]
            virtual scalar CpTR
            (
                const label speciei,
                const scalar p,
                const scalar TTR
            ) const override;

            //- Heat capacity at constant volume [J/(kg K)]
            virtual scalar CvT
            (
                const label speciei,
                const scalar p,
                const scalar TTR
            ) const override;

            //- Enthalpy/Internal energy [J/kg]
            virtual scalar H
            (
                const label speciei,
                const scalar p,
                const scalar TTR
            ) const override;

            //- Absolute enthalpy [J/kg]
            virtual scalar Ha
            (
                const label speciei,
                const scalar p,
                const scalar TTR
            ) const override;

            //- Sensible enthalpy [J/kg]
            virtual scalar Hs
            (
                const label speciei,
                const scalar p,
                const scalar TTR
            ) const override;


            //- Entropy [J/(kg K)]
            virtual scalar S
            (
                const label speciei,
                const scalar p,
                const scalar TTR
            ) const override;

            //- Sensible translational internal energy [J/kg]
            virtual scalar EsT
            (
                const label speciei,
                const scalar p,
                const scalar TTR
            ) const override;

            //- Sensible rotational internal energy [J/kg]
            virtual scalar EsR
            (
                const label speciei,
                const scalar p,
                const scalar TTR
            ) const override;

            //- Sensible vibrational internal energy [J/kg]
            virtual scalar EsVib
            (
                const label speciei,
                const scalar p,
                const scalar TTR,
                const scalar TVib,
                const scalar ThetaVib
            ) const override;

            //- Gibbs free energy [J/kg]
            virtual scalar G
            (
                const label speciei,
                const scalar p,
                const scalar TTR
            ) const override;

        // Per specie transport properties

            //- Dynamic viscosity [kg/m/s]
            virtual scalar mu
            (
                const label speciei,
                const scalar TTR
            ) const override;

            //- Translational-Rotational thermal conductivity [W/m/K]
            virtual scalar kappaTR
            (
                const label speciei,
                const scalar p,
                const scalar TTR
            ) const override;
            
            //- Vibrational thermal conductivity [W/m/K]
            virtual scalar kappaVib
            (
                const label speciei,
                const scalar p,
                const scalar TTR,
                const scalar TVib,
                const scalar ThetaVib
            ) const override;

            //- Thermal diffusivity of enthalpy [kg/m/s]
            virtual scalar alphah
            (
                const label speciei,
                const scalar p,
                const scalar TTR
            ) const override;

            //- Density [kg/m3]
            virtual scalar rho
            (
                const label speciei,
                const scalar p,
                const scalar TTR
            ) const override;

            //- Vibrational theta constant
            virtual scalar thetaVib
            (
                const label speciei
            ) const override;

            virtual bool isSpecieMolecular
            (
                const label speciei
            ) const override;

            virtual scalar TEVib
            (
                const label speciei,
                const scalar esVib,
                const scalar p,
                const scalar TTR,
                const scalar TVib0
            ) const override;

        // Field properties abstract functions
        // provided here to avoid issues with inheritance of virtual functions
        // in heThermo

            //- Heat capacity at constant pressure for patch [J/kg/K]
            virtual tmp<scalarField> CpTR
            (
                const scalarField& p,
                const scalarField& TTR,
                const label patchi
            ) const = 0;

            //- Heat capacity at constant pressure for patch [J/kg/K]
            virtual tmp<volScalarField> CpTR() const = 0;

            //- Heat capacity at constant pressure at p and TTR  [J/kg/K]
            virtual tmp<scalarField> CpTR
            (
                const scalarField& p,
                const scalarField& TTR,
                const labelList& cells
            ) const = 0;

            //- Heat capacity at constant volume for patch [J/kg/K]
            virtual tmp<scalarField> CvT
            (
                const scalarField& p,
                const scalarField& TTR,
                const label patchi
            ) const = 0;

            //- Heat capacity at constant volume [J/kg/K]
            virtual tmp<volScalarField> CvT() const = 0;

            //- Heat capacity at constant volume for patch [J/kg/K]
            virtual tmp<scalarField> CvR
            (
                const scalarField& p,
                const scalarField& TTR,
                const label patchi
            ) const = 0;

            //- Heat capacity at constant volume [J/kg/K]
            virtual tmp<volScalarField> CvR() const = 0;

            //- Heat capacity at constant volume for patch [J/kg/K]
            virtual tmp<scalarField> CvVib
            (
                const scalarField& p,
                const scalarField& TTR,
                const scalarField& TVib,
                const scalar ThetaVib,
                const label patchi
            ) const = 0;

            //- Heat capacity at constant volume [J/kg/K]
            virtual tmp<volScalarField> CvVib() const = 0;

            //- Molecular weight of the given specie [kg/kmol]
            virtual tmp<volScalarField> W() const = 0;

            //- Thermal diffusivity for temperature
            //  of mixture for patch [J/m/s/K]
            virtual tmp<scalarField> kappaTR
            (
                const label patchi
            ) const = 0;

            //- Thermal diffusivity for temperature of mixture [J/m/s/K]
            virtual tmp<volScalarField> kappaTR() const = 0;

            //- Thermal diffusivity for temperature
            //  of mixture for patch [J/m/s/K]
            virtual tmp<scalarField> kappaVib
            (
                const label patchi
            ) const = 0;

            //- Thermal diffusivity for temperature of mixture [J/m/s/K]
            virtual tmp<volScalarField> kappaVib() const = 0;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "Specie2TMixture.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
