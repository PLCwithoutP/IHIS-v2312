/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2017-2023 OpenCFD Ltd.
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
    Foam::basic2TThermo

Description
    Abstract base-class for fluid and solid thermodynamic properties

    An important concept is that the thermo can share an existing T
    (similar to p) in which case it will not try to update it. This gets
    triggered purely on construction order - the first one to register
    T is responsible for updating it. Note that the mechanism also means
    that if multiple thermos are running on the same mesh, only the first one
    will update the temperature. The behaviour can be overridden using the
    'updateT' dictionary entry.

SourceFiles
    basic2TThermo.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_basic2TThermo_H
#define Foam_basic2TThermo_H

#include "volFields.H"
#include "typeInfo.H"
#include "IOdictionary.H"
#include "autoPtr.H"
#include "wordIOList.H"
#include "tcLibraryInterface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class basic2TThermo Declaration
\*---------------------------------------------------------------------------*/

class basic2TThermo
:
    public IOdictionary,
    public tcLibraryInterface
{
    // Private Data

        //- Components names/order
        static const wordList componentHeader4;

        //- Components names/order
        static const wordList componentHeader7;


    // Private Member Functions

        //- Construct name of thermo package from dictionary components
        static word makeThermoName
        (
            const dictionary& dict,
            const wordList*& cmptHeaderPtr
        );

        //- Look up field from registry or construct and store
        static volScalarField& lookupOrConstruct
        (
            const fvMesh& mesh,
            const word& fieldName,
            bool& isOwner  //!< Stored to registry by this instance
        );

        //- Generic lookup for thermodynamics package thermoTypeName
        //  \return constructor pointer, or FatalError
        template<class Thermo, class ThermoConstructTable>
        static typename ThermoConstructTable::mapped_type
        getThermoOrDie
        (
            const dictionary& thermoTypeDict,
            ThermoConstructTable& thermoTable,
            const word& thermoTypeName,
            const wordList& cmptNames
        );

        //- Generic lookup for each of the related thermodynamics packages
        //  \return constructor pointer, or FatalError
        template<class Thermo, class ThermoConstructTable>
        static typename ThermoConstructTable::mapped_type
        getThermoOrDie
        (
            const dictionary& thermoDict,
            ThermoConstructTable& thermoTable
        );


protected:

    // Protected Data

        //- Phase-name
        const word& phaseName_;

        //- Pressure created and stored by this instance
        bool pOwner_;

        //- Temperature created and stored by this instance
        bool TOwnerTR_;

        //- Temperature created and stored by this instance
        bool TOwnerVib_;

        //- Include dpdt term in the enthalpy equation?
        bool dpdt_;

        //- True if separate vibrational temperature is physically enabled
        bool twoTemperature_;

        //- Chemistry enabled?
        bool chemistry_;

        //- Chemistry reader type name
        word chemistryReader_;

        //- External species dictionary file
        fileName speciesListFile_;

        //- External reactions dictionary file
        fileName reactionsListFile_;
    // Fields

        //- Pressure [Pa]
        volScalarField& p_;

        //- Translational-Rotational Temperature [K]
        volScalarField& TTR_;

        //- Vibrational Temperature [K]
        volScalarField& TVib_;

        //- Laminar thermal diffusivity [kg/m/s]
        volScalarField alpha_;


    // Protected Member Functions

        //- Return the enthalpy/internal energy field boundary types
        //- by interrogating the temperature field boundary types
        wordList hBoundaryTypes();

        //- Return the enthalpy/internal energy field boundary base types
        //- by interrogating the temperature field boundary types
        wordList hBoundaryBaseTypes();

        //- Return the enthalpy/internal energy field boundary types
        //- by interrogating the temperature field boundary types
        wordList eVibBoundaryTypes();

        //- Return the enthalpy/internal energy field boundary base types
        //- by interrogating the temperature field boundary types
        wordList eVibBoundaryBaseTypes();

        //- No copy construct
        basic2TThermo(const basic2TThermo&) = delete;

        //- Read and validate twoTemperature and chemistry entries
        static bool readTwoTemperatureEntry(const dictionary& dict);
        static bool readChemistryEntry(const dictionary& dict);

        static bool readBoolEntry(const dictionary& dict, const word& key);

        static word readRequiredWordEntry(const dictionary& dict, const word& key);

        static fileName readRequiredFileEntry(const dictionary& dict, const word& key);

        static fileName resolveCaseFile(const fileName& raw, const fvMesh& mesh);
public:

    //- Runtime type information
    TypeName("basic2TThermo");


    //- Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        basic2TThermo,
        fvMesh,
        (const fvMesh& mesh, const word& phaseName),
        (mesh, phaseName)
    );


    //- Declare run-time constructor selection table for dictionary based
    declareRunTimeSelectionTable
    (
        autoPtr,
        basic2TThermo,
        fvMeshDictPhase,
        (const fvMesh& mesh, const word& phaseName, const word& dictName),
        (mesh, phaseName, dictName)
    );


    // Constructors

        //- Construct from mesh and phase name
        basic2TThermo
        (
            const fvMesh&,
            const word& phaseName
        );

        //- Construct from mesh, dictionary and phase name
        basic2TThermo
        (
            const fvMesh&,
            const dictionary&,
            const word& phaseName
        );

        //- Construct from mesh, phase name and explicit naming of the
        //- dictionary (so it can be shared amongst phases).
        basic2TThermo
        (
            const fvMesh&,
            const word& phaseName,
            const word& dictionaryName
        );


    // Selectors

        //- Generic New for each of the related thermodynamics packages
        template<class Thermo>
        static autoPtr<Thermo> New
        (
            const fvMesh&,
            const word& phaseName = word::null
        );

        //- Generic New for each of the related thermodynamics packages
        template<class Thermo>
        static autoPtr<Thermo> New
        (
            const fvMesh&,
            const dictionary&,
            const word& phaseName = word::null
        );

        //- Generic New for each of the related thermodynamics packages
        template<class Thermo>
        static autoPtr<Thermo> New
        (
            const fvMesh&,
            const word& phaseName,
            const word& dictName
        );

        //- Specialisation of the Generic New for basic2TThermo
        static autoPtr<basic2TThermo> New
        (
            const fvMesh&,
            const word& phaseName = word::null
        );


    //- Destructor
    virtual ~basic2TThermo();


    // Member Functions

        //- The dictionary name ("thermophysicalProperties")
        static const word dictName;

        //- The phase property name as property.phase
        //- \sa IOobject::groupName
        static word phasePropertyName(const word& name, const word& phaseName)
        {
            return IOobject::groupName(name, phaseName);
        }

        //- The phase property name as property.phase (using phaseName_)
        //- \sa IOobject::groupName
        word phasePropertyName(const word& name) const
        {
            return IOobject::groupName(name, phaseName_);
        }

        //- The phase property name scoped as scope:property.phase
        //- \sa IOobject::groupName and \sa IOobject::scopedName
        word phaseScopedName(const std::string& scope, const word& name) const
        {
            return IOobject::groupName
            (
                IOobject::scopedName(scope, name),
                phaseName_
            );
        }

        static const basic2TThermo& lookupThermo(const fvPatchScalarField& pf);

        //- Print (filtered) table of thermo names, splits on \c " ,<>"
        static Ostream& printThermoNames
        (
            Ostream& os,
            const wordList& cmptNames,
            const wordList& thermoNames
        );

        //- Check that the thermodynamics package is consistent
        //  with energy forms supported by the application
        void validate
        (
            const string& app,
            const word&
        ) const;

        //- Check that the thermodynamics package is consistent
        //  with energy forms supported by the application
        void validate
        (
            const string& app,
            const word&,
            const word&
        ) const;

        //- Check that the thermodynamics package is consistent
        //  with energy forms supported by the application
        void validate
        (
            const string& app,
            const word&,
            const word&,
            const word&
        ) const;

        //- Check that the thermodynamics package is consistent
        //  with energy forms supported by the application
        void validate
        (
            const string& app,
            const word&,
            const word&,
            const word&,
            const word&
        ) const;

        //- Split thermo package name into a list of components names
        //  Splits on \c " ,<>"
        //  \return empty list if the split name does not have the
        //  expected number of components (non-zero).
        static wordList splitThermoName
        (
            const std::string& thermoName,
            const int nExpectedCmpts
        );

        //- Update properties
        virtual void correct() = 0;

        //- Update properties
        virtual void correctTEnergy() = 0;

        //- Update properties
        virtual void correctREnergy() = 0;

        //- Update properties
        virtual void correctVibEnergy() = 0;

        virtual void correctTVib() = 0;
        virtual PtrList<volScalarField>& correctVibSource() = 0;
        virtual PtrList<volScalarField>& correctVibSource
        (
            const PtrList<volScalarField>& TVibSpecies
        ) = 0;
        virtual PtrList<volScalarField>& correctVibVibSource() = 0;
        virtual PtrList<volScalarField>& correctVibVibSource
        (
            const PtrList<volScalarField>& TVibSpecies
        ) = 0;
        virtual PtrList<volScalarField>& correctVTRelaxationTime() = 0;


        virtual PtrList<volScalarField>& correctCVSource
        (
            const PtrList<volScalarField::Internal>& RR
        ) = 0;

        virtual scalar wilkeKappaAverage() const = 0;
        virtual scalar wilkeMuAverage() const = 0;
        virtual Foam::volScalarField fickDiffusionCoeff() const = 0;

        //- Return the name of the thermo physics
        virtual word thermoName() const = 0;

        //- Return true if the equation of state is incompressible
        //  i.e. rho != f(p)
        virtual bool incompressible() const = 0;

        //- Return true if the equation of state is isochoric
        //  i.e. rho = const
        virtual bool isochoric() const = 0;

        //- True if dpdt term should be included in enthalpy equation
        bool dpdt() const noexcept
        {
            return dpdt_;
        }

        //- True if T_TR should be updated
        bool updateTTR() const noexcept
        {
            return TOwnerTR_;
        }

        //- True only when vibrational temperature is both enabled
        //- physically and owned/updated by this thermo object
        bool updateTVib() const noexcept
        {
            return TOwnerVib_ && twoTemperature_;;
        }

        //- True if two-temperature physics is enabled
        virtual bool twoTemperature() const noexcept
        {
            return twoTemperature_;
        }

        //- Just alias, same functionality as twoTemperature
        bool vibrationalCheck() const noexcept
        {
            return twoTemperature_;
        }

        virtual bool chemistry() const noexcept
        {
            return chemistry_;
        }

        virtual bool chemistryCapable() const noexcept
        {
            return false;
        }
        
        virtual const word& chemistryReader() const noexcept
        {
            return chemistryReader_;
        }

        virtual const fileName& speciesListFile() const noexcept
        {
            return speciesListFile_;
        }

        virtual const fileName& reactionsListFile() const noexcept
        {
            return reactionsListFile_;
        }

        // Access to thermodynamic state variables

            //- Pressure [Pa]
            //  Non-const access allowed for transport equations
            virtual volScalarField& p();

            //- Pressure [Pa]
            virtual const volScalarField& p() const;

            //- Density [kg/m^3]
            virtual tmp<volScalarField> rho() const = 0;

            //- Density for patch [kg/m^3]
            virtual tmp<scalarField> rho(const label patchi) const = 0;

            //- Enthalpy [J/kg]
            //  Non-const access allowed for transport equations
            virtual volScalarField& h() = 0;

            //- Enthalpy [J/kg]
            virtual const volScalarField& h() const = 0;

            //- Translational-Rotational Internal Energy [J/kg]
            //  Non-const access allowed for transport equations
            virtual volScalarField& eTR() = 0;

            //- Translational-Rotational Internal Energy [J/kg]
            virtual const volScalarField& eTR() const = 0;

            //- Translational Internal Energy [J/kg]
            //  Non-const access allowed for transport equations
            virtual volScalarField& eT() = 0;

            //- Translational Internal Energy [J/kg]
            virtual const volScalarField& eT() const = 0;

            //- Rotational Internal Energy [J/kg]
            //  Non-const access allowed for transport equations
            virtual volScalarField& eR() = 0;

            //- Rotational Internal Energy [J/kg]
            virtual const volScalarField& eR() const = 0;

            //- Vibrational Internal Energy [J/kg]
            //  Non-const access allowed for transport equations
            virtual volScalarField& eVib() = 0;

            //- Vibrational Internal Energy [J/kg]
            virtual const volScalarField& eVib() const = 0;
            
            //- Mixture characteristic vibrational temperature [K]
            virtual scalar thetaVib() const = 0;

            //- Enthalpy
            //  for given pressure and temperature [J/kg]
            virtual tmp<volScalarField> h
            (
                const volScalarField& p,
                const volScalarField& TTR
            ) const = 0;

            //- Enthalpy for cell-set [J/kg]
            virtual tmp<scalarField> h
            (
                const scalarField& p,
                const scalarField& TTR,
                const labelList& cells
            ) const = 0;

            //- Enthalpy for patch [J/kg]
            virtual tmp<scalarField> h
            (
                const scalarField& p,
                const scalarField& TTR,
                const label patchi
            ) const = 0;

            //- Internal Energy
            //  for given pressure and temperature [J/kg]
            virtual tmp<volScalarField> eTR
            (
                const volScalarField& p,
                const volScalarField& TTR
            ) const = 0;

            //- Internal Energy for cell-set [J/kg]
            virtual tmp<scalarField> eTR
            (
                const scalarField& p,
                const scalarField& TTR,
                const labelList& cells
            ) const = 0;

            //- Internal Energy for patch [J/kg]
            virtual tmp<scalarField> eTR
            (
                const scalarField& p,
                const scalarField& TTR,
                const label patchi
            ) const = 0;

            //- Translational Internal Energy
            //  for given pressure and temperature [J/kg]
            virtual tmp<volScalarField> eT
            (
                const volScalarField& p,
                const volScalarField& TTR
            ) const = 0;

            //- Translational Internal Energy for cell-set [J/kg]
            virtual tmp<scalarField> eT
            (
                const scalarField& p,
                const scalarField& TTR,
                const labelList& cells
            ) const = 0;

            //- Translational Internal Energy for patch [J/kg]
            virtual tmp<scalarField> eT
            (
                const scalarField& p,
                const scalarField& TTR,
                const label patchi
            ) const = 0;

            //- Rotational Internal Energy
            //  for given pressure and temperature [J/kg]
            virtual tmp<volScalarField> eR
            (
                const volScalarField& p,
                const volScalarField& TTR
            ) const = 0;

            //- Rotational Internal Energy for cell-set [J/kg]
            virtual tmp<scalarField> eR
            (
                const scalarField& p,
                const scalarField& TTR,
                const labelList& cells
            ) const = 0;

            //- Rotational Internal Energy for patch [J/kg]
            virtual tmp<scalarField> eR
            (
                const scalarField& p,
                const scalarField& TTR,
                const label patchi
            ) const = 0;

            //- Vibrational Internal Energy
            //  for given pressure and temperature [J/kg]
            virtual tmp<volScalarField> eVib
            (
                const volScalarField& p,
                const volScalarField& TTR,
                const volScalarField& TVib,
                const scalar ThetaVib
            ) const = 0;

            //- Vibrational Internal Energy for cell-set [J/kg]
            virtual tmp<scalarField> eVib
            (
                const scalarField& p,
                const scalarField& TTR,
                const scalarField& TVib,
                const scalar ThetaVib,
                const labelList& cells
            ) const = 0;

            //- Vibrational Internal Energy for patch [J/kg]
            virtual tmp<scalarField> eVib
            (
                const scalarField& p,
                const scalarField& TTR,
                const scalarField& TVib,
                const scalar ThetaVib,
                const label patchi
            ) const = 0;

            //- Chemical enthalpy [J/kg]
            virtual tmp<volScalarField> hc() const = 0;

            //- Temperature from enthalpy for cell-set
            virtual tmp<scalarField> TH_TR
            (
                const scalarField& h,
                const scalarField& p,
                const scalarField& TTR0,      // starting translational-rotational temperature
                const labelList& cells
            ) const = 0;

            //- Temperature from enthalpy for patch
            virtual tmp<scalarField> TH_TR
            (
                const scalarField& h,
                const scalarField& p,
                const scalarField& TTR0,      // starting translational-rotational temperature
                const label patchi
            ) const = 0;

            //- Temperature from translational-rotational energy for cell-set
            virtual tmp<scalarField> TE_TR
            (
                const scalarField& eTR,
                const scalarField& p,
                const scalarField& TTR0,      // starting translational-rotational temperature
                const labelList& cells
            ) const = 0;

            //- Temperature from translational-rotational energy for patch
            virtual tmp<scalarField> TE_TR
            (
                const scalarField& eTR,
                const scalarField& p,
                const scalarField& TTR0,      // starting translational-rotational temperature
                const label patchi
            ) const = 0;

            //- Temperature from vibrational internal energy for cell-set
            virtual tmp<scalarField> TE_Vib
            (
                const scalarField& esVib,
                const scalarField& p,
                const scalarField& TTR0,      // starting translational-rotational temperature
                const scalarField& TVib0,      // starting vibrational temperature
                const labelList& cells
            ) const = 0;

            //- Temperature from vibrational internal energy for patch
            virtual tmp<scalarField> TE_Vib
            (
                const scalarField& esVib,
                const scalarField& p,
                const scalarField& TTR0,      // starting translational-rotational temperature
                const scalarField& TVib0,      // starting vibrational temperature
                const label patchi
            ) const = 0;

        // Fields derived from thermodynamic state variables

            //- Translational-Rotational Temperature [K]
            virtual const volScalarField& TTR() const;

            //- Vibrational Temperature [K]
            virtual const volScalarField& TVib() const;

            //- Translational-Rotational Temperature [K]
            //  Non-const access allowed for transport equations
            virtual volScalarField& TTR();

            //- Vibrational Temperature [K]
            //  Non-const access allowed for transport equations
            virtual volScalarField& TVib();

            //- Heat capacity at constant pressure [J/kg/K]
            virtual tmp<volScalarField> CpTR() const = 0;

            //- Heat capacity at constant pressure for patch [J/kg/K]
            virtual tmp<scalarField> CpTR
            (
                const scalarField& p,
                const scalarField& TTR,
                const label patchi
            ) const = 0;

            //- Heat capacity using pressure and temperature [J/kg/K]
            virtual tmp<scalarField> CpTR
            (
                const scalarField& p,
                const scalarField& TTR,
                const labelList& cells
            ) const = 0;
            
            //- Translational-Rotational heat capacity at constant volume [J/kg/K]
            virtual tmp<volScalarField> CvTR() const = 0;

            //- Translational-Rotational heat capacity at constant volume [J/kg/K]
            virtual tmp<scalarField> CvTR
            (
                const scalarField& p,
                const scalarField& T,
                const label patchi
            ) const = 0;

            //- Translational-Rotational heat capacity at constant volume [J/kg/K]
            virtual tmp<scalarField> CvTR
            (
                const scalarField& p,
                const scalarField& T,
                const labelList& cells
            ) const = 0;

            //- Translational heat capacity at constant volume [J/kg/K]
            virtual tmp<volScalarField> CvT() const = 0;

            //- Translational heat capacity at constant volume for patch [J/kg/K]
            virtual tmp<scalarField> CvT
            (
                const scalarField& p,
                const scalarField& TTR,
                const label patchi
            ) const = 0;

            //- Rotational heat capacity at constant volume [J/kg/K]
            virtual tmp<volScalarField> CvR() const = 0;

            //- Rotational heat capacity at constant volume for patch [J/kg/K]
            virtual tmp<scalarField> CvR
            (
                const scalarField& p,
                const scalarField& TTR,
                const label patchi
            ) const = 0;

            //- Vibrational heat capacity at constant volume [J/kg/K]
            virtual tmp<volScalarField> CvVib() const = 0;

            //- Vibrational heat capacity at constant volume for patch [J/kg/K]
            virtual tmp<scalarField> CvVib
            (
                const scalarField& p,
                const scalarField& TTR,
                const scalarField& Tvib,
                const scalar ThetaVib,
                const label patchi
            ) const = 0;

            //- Density from pressure and temperature from EoS
            virtual tmp<scalarField> rhoEoS
            (
                const scalarField& p,
                const scalarField& TTR,
                const labelList& cells
            ) const = 0;

            //- Gamma = CpTR/CvT []
            virtual tmp<volScalarField> gamma() const = 0;

            //- Gamma = CpTR/CvT for patch []
            virtual tmp<scalarField> gamma
            (
                const scalarField& p,
                const scalarField& TTR,
                const label patchi
            ) const = 0;

            //- Translational-Rotational heat capacity at constant pressure/volume [J/kg/K]
            virtual tmp<volScalarField> CpvTR() const = 0;

            //- Translational-Rotational heat capacity at constant pressure/volume for patch [J/kg/K]
            virtual tmp<scalarField> CpvTR
            (
                const scalarField& p,
                const scalarField& TTR,
                const label patchi
            ) const = 0;

            //- Translational-Rotational heat capacity ratio []
            virtual tmp<volScalarField> CpByCpvTR() const = 0;

            //- Translational-Rotational heat capacity ratio for patch []
            virtual tmp<scalarField> CpByCpvTR
            (
                const scalarField& p,
                const scalarField& TTR,
                const label patchi
            ) const = 0;

            //- Molecular weight [kg/kmol]
            virtual tmp<volScalarField> W() const = 0;


        // Access to transport state variables

            //- Thermal diffusivity for enthalpy of mixture [kg/m/s]
            virtual const volScalarField& alpha() const;

            //- Thermal diffusivity for enthalpy of mixture for patch [kg/m/s]
            virtual const scalarField& alpha
            (
                const label patchi
            ) const;


        // Fields derived from transport state variables

            //- Translational-rotational thermal diffusivity for temperature of mixture [J/m/s/K]
            virtual tmp<volScalarField> kappaTR() const = 0;

            //- Translational-rotational thermal diffusivity for temperature
            //  of mixture for patch [J/m/s/K]
            virtual tmp<scalarField> kappaTR
            (
                const label patchi
            ) const = 0;

            //- Vibrational thermal diffusivity for temperature of mixture [J/m/s/K]
            virtual tmp<volScalarField> kappaVib() const = 0;

            //- Vibrational thermal diffusivity for temperature
            //  of mixture for patch [J/m/s/K]
            virtual tmp<scalarField> kappaVib
            (
                const label patchi
            ) const = 0;

             //- Translational-rotational thermal diffusivity for energy of mixture [kg/m/s]
            virtual tmp<volScalarField> alphaheTR() const = 0;

            //- Translational-rotational thermal diffusivity for energy of mixture for patch [kg/m/s]
            virtual tmp<scalarField> alphaheTR(const label patchi) const = 0;


        //- Read thermophysical properties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "basic2TThermoTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
