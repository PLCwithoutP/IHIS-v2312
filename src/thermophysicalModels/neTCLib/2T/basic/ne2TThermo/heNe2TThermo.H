/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2015-2017 OpenCFD Ltd.
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
    Foam::heNe2TThermo

Description
    Energy for a mixture based on compressibility

SourceFiles
    heNe2TThermo.C

\*---------------------------------------------------------------------------*/

#ifndef heNe2TThermo_H
#define heNe2TThermo_H

#include "ne2TThermo.H"
#include "he2TThermo.H"
#include "PtrList.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class heNe2TThermo Declaration
\*---------------------------------------------------------------------------*/

template<class ne2TThermoType, class MixtureType>
class heNe2TThermo
:
    public he2TThermo<ne2TThermoType, MixtureType>
{
    // Private Member Functions

        //- Calculate the thermo variables
        void calculate
        (
            const volScalarField& p,
            volScalarField& TTR,
            volScalarField& TVib,
            volScalarField& eTR,
            volScalarField& esT,
            volScalarField& esR,
            volScalarField& esVib,
            volScalarField& psi,
            volScalarField& mu,
            volScalarField& alpha,
            const bool doOldTimes
        );
        
        //- Update the vibrational internal energy
        void calculateVibEnergy
        (
            const volScalarField& p,
            const volScalarField& TTR,
            const volScalarField& TVib,
            volScalarField& esVib
        );

        //- Construct as copy (not implemented)
        heNe2TThermo(const heNe2TThermo<ne2TThermoType, MixtureType>&);


public:

    //- Runtime type information
    TypeName("heNe2TThermo");


    // Constructors

        //- Construct from mesh and phase name
        heNe2TThermo
        (
            const fvMesh&,
            const word& phaseName
        );

        //-  Construct from mesh and phase name and dict name
        heNe2TThermo
        (
            const fvMesh&,
            const word& phaseName,
            const word& dictionaryName
        );

    //- Destructor
    virtual ~heNe2TThermo();


    // Member functions

        virtual bool chemistryCapable() const noexcept override
        {
            return false;
        }

        //- Update properties
        virtual void correct();
        virtual void correctTEnergy();
        virtual void correctREnergy();
        virtual void correctVibEnergy();
        virtual void correctTVib();
        
        virtual PtrList<volScalarField>& correctVibSource();
        virtual PtrList<volScalarField>& correctVibSource
        (
            const PtrList<volScalarField>& TVibSpecies
        );
        virtual PtrList<volScalarField>& correctVibVibSource();
        virtual PtrList<volScalarField>& correctVibVibSource
        (
            const PtrList<volScalarField>& TVibSpecies
        );
        virtual PtrList<volScalarField>& correctVTRelaxationTime();
        virtual PtrList<volScalarField>& correctCVSource
        (
            const PtrList<volScalarField::Internal>& RR
        );
        virtual scalar wilkeKappaAverage() const;
        virtual scalar wilkeMuAverage() const;
        virtual volScalarField fickDiffusionCoeff() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "heNe2TThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
