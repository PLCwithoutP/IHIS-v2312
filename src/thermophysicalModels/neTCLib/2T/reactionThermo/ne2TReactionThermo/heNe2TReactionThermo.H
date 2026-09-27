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
    Foam::heNe2TReactionThermo

Description
    Energy for a mixture based on compressibility

SourceFiles
    heNe2TReactionThermo.C

\*---------------------------------------------------------------------------*/

#ifndef heNe2TReactionThermo_H
#define heNe2TReactionThermo_H

#include "ne2TReactionThermo.H"
#include "he2TThermo.H"
#include "WilkeMR.H"
#include "VTEnergySource.H"
#include "VVEnergySource.H"
#include "CVEnergySource.H"
#include "FickDM.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class heNe2TReactionThermo Declaration
\*---------------------------------------------------------------------------*/

template<class basicNe2TThermo, class MixtureType>
class heNe2TReactionThermo
:
    public he2TThermo<basicNe2TThermo, MixtureType>
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
        heNe2TReactionThermo(const heNe2TReactionThermo<basicNe2TThermo, MixtureType>&);

        scalar calculateRhoMixture
        (
            const scalar pCell,
            const scalar TTRCell,
            const label  celli
        ) const;
        
        scalar calculateCpTRMixture
        (
            const scalar pCell,
            const scalar TTRCell,
            const label  celli
        ) const;
public:

    //- Runtime type information
    TypeName("heNe2TReactionThermo");


    // Constructors

        //- Construct from mesh and phase name
        heNe2TReactionThermo
        (
            const fvMesh&,
            const word& phaseName
        );

        //-  Construct from mesh and phase name and dict name
        heNe2TReactionThermo
        (
            const fvMesh&,
            const word& phaseName,
            const word& dictionaryName
        );

    //- Destructor
    virtual ~heNe2TReactionThermo();


    // Member functions
    
        virtual bool chemistryCapable() const noexcept override
        {
            return true;
        }
        //- Update properties
        virtual void correct();
        virtual void correctTEnergy();
        virtual void correctREnergy();
        virtual void correctVibEnergy();
        virtual void correctTVib();
        
        PtrList<volScalarField>& correctVibSource();
        PtrList<volScalarField>& correctVibSource
        (
            const PtrList<volScalarField>& TVibSpecies
        );
        PtrList<volScalarField>& correctVibVibSource();
        PtrList<volScalarField>& correctVibVibSource
        (
            const PtrList<volScalarField>& TVibSpecies
        );
        PtrList<volScalarField>& correctVTRelaxationTime();
        
        PtrList<volScalarField>& correctCVSource
        (
            const PtrList<volScalarField::Internal>& RR
        );
        
        virtual scalar wilkeKappaAverage() const override;
        virtual scalar wilkeMuAverage() const;
        virtual Foam::volScalarField fickDiffusionCoeff() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "heNe2TReactionThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
