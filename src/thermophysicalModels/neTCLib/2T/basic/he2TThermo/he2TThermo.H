/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2015-2021 OpenCFD Ltd.
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
    Foam::he2TThermo

Description
    Enthalpy for a mixture

SourceFiles
    he2TThermo.C

\*---------------------------------------------------------------------------*/

#ifndef he2TThermo_H
#define he2TThermo_H

#include "basic2TMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class he2TThermo Declaration
\*---------------------------------------------------------------------------*/

template<class ne2TThermoType, class MixtureType>
class he2TThermo
:
    public ne2TThermoType,
    public MixtureType
{
protected:

    // Protected data

        //- Enthalpy field
        volScalarField h_;

        //- Translational Energy field
        volScalarField eTR_;

        //- Translational Energy field
        volScalarField eT_;

        //- Rotational Energy field
        volScalarField eR_;

        //- Vibrational Energy field
        volScalarField eVib_;

    // Protected Member Functions

        // Enthalpy

            //- Correct the enthalpy field boundaries
            void hBoundaryCorrection(volScalarField& h);

            void eVibBoundaryCorrection(volScalarField& eVib);


private:

    // Private Member Functions

        //- Construct as copy (not implemented)
        he2TThermo(const he2TThermo<ne2TThermoType, MixtureType>&);


        //- Initialize he2TThermo
        void init
        (
            const volScalarField& p,
            const volScalarField& TTR,
            const volScalarField& TVib,
            volScalarField& h,
            volScalarField& eTR,
            volScalarField& eT,
            volScalarField& eR,
            volScalarField& eVib
        );


public:

    // Constructors

        //- Construct from mesh
        he2TThermo
        (
            const fvMesh&,
            const word& phaseName
        );

        //- Construct from mesh and dictionary
        he2TThermo
        (
            const fvMesh&,
            const dictionary&,
            const word& phaseName
        );

        //- Construct from mesh,dictionary,phase name with a single temperature
        he2TThermo
        (
            const fvMesh&,
            const word& phaseName,
            const word& dictionaryName
        );


    //- Destructor
    virtual ~he2TThermo();


    // Member functions

        //- Return the composition of the mixture
        virtual typename MixtureType::basic2TMixtureType&
        composition()
        {
            return *this;
        }

        //- Return the composition of the mixture
        virtual const typename MixtureType::basic2TMixtureType&
        composition() const
        {
            return *this;
        }

        //- Return the name of the thermo physics
        virtual word thermoName() const
        {
            return MixtureType::thermoType::typeName();
        }

        //- Return true if the equation of state is incompressible
        //  i.e. rho != f(p)
        virtual bool incompressible() const
        {
            return MixtureType::thermoType::incompressible;
        }

        //- Return true if the equation of state is isochoric
        //  i.e. rho = const
        virtual bool isochoric() const
        {
            return MixtureType::thermoType::isochoric;
        }


        // Access to thermodynamic state variables

            //- Enthalpy [J/kg]
            //  Non-const access allowed for transport equations
            virtual volScalarField& h()
            {
                return h_;
            }

            //- Enthalpy [J/kg]
            virtual const volScalarField& h() const
            {
                return h_;
            }

            //- Translational-Rotational Internal Energy [J/kg]
            //  Non-const access allowed for transport equations
            virtual volScalarField& eTR()
            {
                return eTR_;
            }

            //- Translational-Rotational Internal Energy [J/kg]
            virtual const volScalarField& eTR() const
            {
                return eTR_;
            }

            //- Translational Internal Energy [J/kg]
            //  Non-const access allowed for transport equations
            virtual volScalarField& eT()
            {
                return eT_;
            }

            //- Translational Internal Energy [J/kg]
            virtual const volScalarField& eT() const
            {
                return eT_;
            }

            //- Rotational Internal Energy [J/kg]
            //  Non-const access allowed for transport equations
            virtual volScalarField& eR()
            {
                return eR_;
            }

            //- Rotational Internal Energy [J/kg]
            virtual const volScalarField& eR() const
            {
                return eR_;
            }

            //- Vibrational Internal Energy [J/kg]
            //  Non-const access allowed for transport equations
            virtual volScalarField& eVib()
            {
                return eVib_;
            }

            //- Vibrational Internal Energy [J/kg]
            virtual const volScalarField& eVib() const
            {
                return eVib_;
            }

            //- Mixture characteristic vibrational temperature [K]
            virtual scalar thetaVib() const;

        // Fields derived from thermodynamic state variables

            //- Enthalpy
            //  for given pressure and temperature [J/kg]
            virtual tmp<volScalarField> h
            (
                const volScalarField& p,
                const volScalarField& TTR
            ) const;

            //- Enthalpy for cell-set [J/kg]
            virtual tmp<scalarField> h
            (
                const scalarField& p,
                const scalarField& TTR,
                const labelList& cells
            ) const;

            //- Enthalpy for patch [J/kg]
            virtual tmp<scalarField> h
            (
                const scalarField& p,
                const scalarField& TTR,
                const label patchi
            ) const;

            //- Internal Energy
            //  for given pressure and temperature [J/kg]
            virtual tmp<volScalarField> eTR
            (
                const volScalarField& p,
                const volScalarField& TTR
            ) const;

            //- Internal Energy for cell-set [J/kg]
            virtual tmp<scalarField> eTR
            (
                const scalarField& p,
                const scalarField& TTR,
                const labelList& cells
            ) const;

            //- Internal Energy for patch [J/kg]
            virtual tmp<scalarField> eTR
            (
                const scalarField& p,
                const scalarField& TTR,
                const label patchi
            ) const;

            //- Translational Internal Energy
            //  for given pressure and temperature [J/kg]
            virtual tmp<volScalarField> eT
            (
                const volScalarField& p,
                const volScalarField& TTR
            ) const;

            //- Translational Internal Energy for cell-set [J/kg]
            virtual tmp<scalarField> eT
            (
                const scalarField& p,
                const scalarField& TTR,
                const labelList& cells
            ) const;

            //- Translational Internal Energy for patch [J/kg]
            virtual tmp<scalarField> eT
            (
                const scalarField& p,
                const scalarField& TTR,
                const label patchi
            ) const;

            //- Rotational Internal Energy
            //  for given pressure and temperature [J/kg]
            virtual tmp<volScalarField> eR
            (
                const volScalarField& p,
                const volScalarField& TTR
            ) const;

            //- Rotational Internal Energy for cell-set [J/kg]
            virtual tmp<scalarField> eR
            (
                const scalarField& p,
                const scalarField& TTR,
                const labelList& cells
            ) const;

            //- Rotational Internal Energy for patch [J/kg]
            virtual tmp<scalarField> eR
            (
                const scalarField& p,
                const scalarField& TTR,
                const label patchi
            ) const;

            //- Vibrational Internal Energy
            //  for given pressure and temperature [J/kg]
            virtual tmp<volScalarField> eVib
            (
                const volScalarField& p,
                const volScalarField& TTR,
                const volScalarField& TVib,
                const scalar ThetaVib
            ) const;

            //- Vibrational Internal Energy for cell-set [J/kg]
            virtual tmp<scalarField> eVib
            (
                const scalarField& p,
                const scalarField& TTR,
                const scalarField& TVib,
                const scalar ThetaVib,
                const labelList& cells
            ) const;

            //- Vibrational Internal Energy for patch [J/kg]
            virtual tmp<scalarField> eVib
            (
                const scalarField& p,
                const scalarField& TTR,
                const scalarField& TVib,
                const scalar ThetaVib,
                const label patchi
            ) const;

            //- Chemical enthalpy [J/kg]
            virtual tmp<volScalarField> hc() const;

            //- Temperature from translational-rotational energy for cell-set
            virtual tmp<scalarField> TE_TR
            (
                const scalarField& eTR,
                const scalarField& p,
                const scalarField& TTR0,      // starting translational-rotational temperature
                const labelList& cells
            ) const;

            //- Temperature from translational-rotational energy for patch
            virtual tmp<scalarField> TE_TR
            (
                const scalarField& eTR,
                const scalarField& p,
                const scalarField& TTR0,      // starting translational-rotational temperature
                const label patchi
            ) const;

            //- Temperature from enthalpy for cell-set
            virtual tmp<scalarField> TH_TR
            (
                const scalarField& h,
                const scalarField& p,
                const scalarField& TTR0,      // starting translational-rotational temperature
                const labelList& cells
            ) const;

            //- Temperature from enthalpy for patch
            virtual tmp<scalarField> TH_TR
            (
                const scalarField& h,
                const scalarField& p,
                const scalarField& TTR0,      // starting translational-rotational temperature
                const label patchi
            ) const;

            //- Temperature from enthalpy for cell-set
            virtual tmp<scalarField> TE_Vib
            (
                const scalarField& h,
                const scalarField& p,
                const scalarField& TTR0,      // starting translational-rotational temperature
                const scalarField& TVib0,      // starting vibrational temperature
                const labelList& cells
            ) const;

            //- Temperature from enthalpy for patch
            virtual tmp<scalarField> TE_Vib
            (
                const scalarField& h,
                const scalarField& p,
                const scalarField& TTR0,      // starting translational-rotational temperature
                const scalarField& TVib0,      // starting vibrational temperature
                const label patchi
            ) const;

            //- Heat capacity at constant pressure for patch [J/kg/K]
            virtual tmp<scalarField> CpTR
            (
                const scalarField& p,
                const scalarField& TTR,
                const label patchi
            ) const;

            //- Heat capacity using pressure and temperature [J/kg/K]
            virtual tmp<scalarField> CpTR
            (
                const scalarField& p,
                const scalarField& TTR,
                const labelList& cells
            ) const;

            //- Translational-Rotational heat capacity at constant volume [J/kg/K]
            virtual tmp<volScalarField> CvTR() const;

            //- Translational-Rotational heat capacity at constant volume [J/kg/K]
            virtual tmp<scalarField> CvTR
            (
                const scalarField& p,
                const scalarField& T,
                const label patchi
            ) const;

            //- Translational-Rotational heat capacity at constant volume [J/kg/K]
            virtual tmp<scalarField> CvTR
            (
                const scalarField& p,
                const scalarField& T,
                const labelList& cells
            ) const;

            //- Heat capacity at constant pressure [J/kg/K]
            virtual tmp<volScalarField> CpTR() const;

            //- Translational heat capacity at constant volume [J/kg/K]
            virtual tmp<volScalarField> CvT() const;

            //- Translational heat capacity at constant volume for patch [J/kg/K]
            virtual tmp<scalarField> CvT
            (
                const scalarField& p,
                const scalarField& TTR,
                const label patchi
            ) const;

            //- Rotational heat capacity at constant volume [J/kg/K]
            virtual tmp<volScalarField> CvR() const;

            //- Rotational heat capacity at constant volume for patch [J/kg/K]
            virtual tmp<scalarField> CvR
            (
                const scalarField& p,
                const scalarField& TTR,
                const label patchi
            ) const;

            //- Vibrational heat capacity at constant volume [J/kg/K]
            virtual tmp<volScalarField> CvVib() const;

            //- Vibrational heat capacity at constant volume for patch [J/kg/K]
            virtual tmp<scalarField> CvVib
            (
                const scalarField& p,
                const scalarField& TTR,
                const scalarField& TVib,
                const scalar ThetaVib,
                const label patchi
            ) const;

            //- Density from pressure and temperature
            virtual tmp<scalarField> rhoEoS
            (
                const scalarField& p,
                const scalarField& TTR,
                const labelList& cells
            ) const;

            //- Gamma = CpTR/CvT []
            virtual tmp<volScalarField> gamma() const;

            //- Gamma = CpTR/CvT for patch []
            virtual tmp<scalarField> gamma
            (
                const scalarField& p,
                const scalarField& TTR,
                const label patchi
            ) const;

            //- Translational-Rotational heat capacity at constant pressure/volume for patch [J/kg/K]
            virtual tmp<scalarField> CpvTR
            (
                const scalarField& p,
                const scalarField& TTR,
                const label patchi
            ) const;

            //- Translational-Rotational heat capacity at constant pressure/volume [J/kg/K]
            virtual tmp<volScalarField> CpvTR() const;

            //- Translational-Rotational heat capacity ratio []
            virtual tmp<volScalarField> CpByCpvTR() const;

            //- Translational-Rotational heat capacity ratio for patch []
            virtual tmp<scalarField> CpByCpvTR
            (
                const scalarField& p,
                const scalarField& TTR,
                const label patchi
            ) const;

            //- Molecular weight [kg/kmol]
            virtual tmp<volScalarField> W() const;


        // Fields derived from transport state variables

            //- Translational-rotational thermal diffusivity for temperature of mixture [J/m/s/K]
            virtual tmp<volScalarField> kappaTR() const;

            //- Translational-rotational thermal diffusivity for temperature
            //  of mixture for patch [J/m/s/K]
            virtual tmp<scalarField> kappaTR
            (
                const label patchi
            ) const;

            //- Vibrational thermal diffusivity for temperature of mixture [J/m/s/K]
            virtual tmp<volScalarField> kappaVib() const;

            //- Vibrational thermal diffusivity for temperature
            //  of mixture for patch [J/m/s/K]
            virtual tmp<scalarField> kappaVib
            (
                const label patchi
            ) const;

             //- Translational-rotational thermal diffusivity for energy of mixture [kg/m/s]
            virtual tmp<volScalarField> alphaheTR() const;

            //- Translational-rotational thermal diffusivity for energy of mixture for patch [kg/m/s]
            virtual tmp<scalarField> alphaheTR(const label patchi) const;


        //- Read thermophysical properties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

#ifdef NoRepository
    #include "he2TThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
