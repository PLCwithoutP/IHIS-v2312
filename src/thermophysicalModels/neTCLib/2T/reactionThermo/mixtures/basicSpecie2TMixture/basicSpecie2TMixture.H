/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2014-2017 OpenFOAM Foundation
    Copyright (C) 2020-2021 OpenCFD Ltd.
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
    Foam::basicSpecie2TMixture

Group
    grpReactionThermophysicalMixtures

Description
    Specialization of basicMultiComponentMixture for a mixture consisting
    of a number for molecular species.

SourceFiles
    basicSpecie2TMixture.C

\*---------------------------------------------------------------------------*/

#ifndef basicSpecie2TMixture_H
#define basicSpecie2TMixture_H

#include "basicMultiComponent2TMixture.H"
#include "specieElement.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

typedef HashTable<List<specieElement>> speciesCompositionTable;

/*---------------------------------------------------------------------------*\
                 Class basicSpecie2TMixture Declaration
\*---------------------------------------------------------------------------*/

class basicSpecie2TMixture
:
    public basicMultiComponent2TMixture
{

public:

    //- Run time type information
    TypeName("basicSpecie2TMixture");

    //- The base class of the mixture
    typedef basicSpecie2TMixture basic2TMixtureType;


    // Constructors

        //- Construct from dictionary, species names, mesh and phase name
        basicSpecie2TMixture
        (
            const dictionary& thermoDict,
            const wordList& specieNames,
            const fvMesh& mesh,
            const word& phaseName
        );


    //- Destructor
    virtual ~basicSpecie2TMixture() = default;


    // Member Functions

        // Per specie properties

            //- Molecular weight of the given specie [kg/kmol]
            virtual scalar W(const label speciei) const = 0;

            //- Specific gas constant of the given specie [kg/kmol]
            virtual scalar R(const label speciei) const = 0;

            //- Chemical enthalpy [J/kg]
            virtual scalar Hc(const label speciei) const = 0;


        // Per specie thermo properties

            //- Heat capacity at constant pressure [J/(kg K)]
            virtual scalar CpTR
            (
                const label speciei,
                const scalar p,
                const scalar TTR
            ) const = 0;

            //- Heat capacity at constant volume [J/(kg K)]
            virtual scalar CvT
            (
                const label speciei,
                const scalar p,
                const scalar TTR
            ) const = 0;

            //- Enthalpy/Internal energy [J/kg]
            virtual scalar H
            (
                const label speciei,
                const scalar p,
                const scalar TTR
            ) const = 0;

            //- Absolute enthalpy [J/kg]
            virtual scalar Ha
            (
                const label speciei,
                const scalar p,
                const scalar TTR
            ) const = 0;

            //- Sensible enthalpy [J/kg]
            virtual scalar Hs
            (
                const label speciei,
                const scalar p,
                const scalar TTR
            ) const = 0;


            //- Entropy [J/(kg K)]
            virtual scalar S
            (
                const label speciei,
                const scalar p,
                const scalar TTR
            ) const = 0;

            //- Sensible translational internal energy [J/kg]
            virtual scalar EsT
            (
                const label speciei,
                const scalar p,
                const scalar TTR
            ) const = 0;

            //- Sensible rotational internal energy [J/kg]
            virtual scalar EsR
            (
                const label speciei,
                const scalar p,
                const scalar TTR
            ) const = 0;

            //- Sensible vibrational internal energy [J/kg]
            virtual scalar EsVib
            (
                const label speciei,
                const scalar p,
                const scalar TTR,
                const scalar TVib,
                const scalar ThetaVib
            ) const = 0;

            //- Gibbs free energy [J/kg]
            virtual scalar G
            (
                const label speciei,
                const scalar p,
                const scalar TTR
            ) const = 0;

        // Per specie transport properties

            //- Dynamic viscosity [kg/m/s]
            virtual scalar mu
            (
                const label speciei,
                const scalar TTR
            ) const = 0;

            //- Translational-Rotational thermal conductivity [W/m/K]
            virtual scalar kappaTR
            (
                const label speciei,
                const scalar p,
                const scalar TTR
            ) const = 0;
            
            //- Vibrational thermal conductivity [W/m/K]
            virtual scalar kappaVib
            (
                const label speciei,
                const scalar p,
                const scalar TTR,
                const scalar TVib,
                const scalar ThetaVib
            ) const = 0;

            //- Thermal diffusivity of enthalpy [kg/m/s]
            virtual scalar alphah
            (
                const label speciei,
                const scalar p,
                const scalar TTR
            ) const = 0;

            //- Density [kg/m3]
            virtual scalar rho
            (
                const label speciei,
                const scalar p,
                const scalar TTR
            ) const = 0;

            //- Vibrational theta constant
            virtual scalar thetaVib
            (
                const label speciei
            ) const override = 0;
            
            virtual bool isSpecieMolecular
            (
                const label speciei
            ) const override = 0;

            virtual scalar TEVib
            (
                const label speciei,
                const scalar esVib,
                const scalar p,
                const scalar TTR,
                const scalar TVib0
            ) const override = 0;
            
            //- Species composition
            virtual autoPtr<speciesCompositionTable> specieComposition() const
            {
                return autoPtr<speciesCompositionTable>::New();
            }

            const speciesTable& species() const override
            {
                return basicMultiComponent2TMixture::species();
            }
            PtrList<volScalarField>& Y() override
            {
                return basicMultiComponent2TMixture::Y();
            }
            const PtrList<volScalarField>& Y() const override
            {
                return basicMultiComponent2TMixture::Y();
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
