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
    Foam::VTEnergySource

Description
    Foam::VTEnergySource

SourceFiles
    VTEnergySource.C

\*---------------------------------------------------------------------------*/

#ifndef VTEnergySource_H
#define VTEnergySource_H

#include "baseEnergySource.H"
#include "speciesTable.H"
#include "wordList.H"
#include "dimensionedScalar.H"
#include "volFields.H"
#include "IOobject.H"
#include "mathematicalConstants.H"
#include "PtrList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class VTEnergySource Declaration
\*---------------------------------------------------------------------------*/
template<class MixtureType, class MixingRule>
class VTEnergySource
:
    public baseEnergySource<MixtureType, MixingRule>
{

        //- Construct as copy (not implemented)
        VTEnergySource(const VTEnergySource&);

private:

    // Private data
    MixtureType& mix_;
    MixingRule& mr_;

    const speciesTable& names_;
    PtrList<volScalarField> Q_VTList_;
    PtrList<volScalarField> tauVTList_;
    scalarSquareMatrix AMW_;
    scalarSquareMatrix BMW_;
    bool mwCoeffsLoaded_;

    const scalar sigma_prime_; // Effective cross section
    const scalar boltzmann_const_; // Boltzmann constant

    void loadMWCoeffs(const fvMesh& mesh);

public:

    // Constructors

        //- Construct from multiComponent2TMixture& 
        VTEnergySource
        (
            MixtureType& mix,
            MixingRule& MR
        );

    // Member functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "VTEnergySource";
        }

        void makeQVibSourceFields
        (
            const fvMesh& mesh
        );

        label nSpecies() const 
        { 
            return names_.size(); 
        }

        const word& specieName(const label i) const 
        { 
            return names_[i]; 
        }

        const wordList& specieNames() const 
        { 
            return names_; 
        }

        // access to bound mixture
        const MixtureType& mixture() const 
        { 
            return mix_; 
        }

        // Y_i field
        const volScalarField& Yi(const label i) const
        {
            return mix_.Y(i);
        }

        // molecular weight W_i
        scalar Wi(const label i) const
        {
            return mix_.W(i);
        }

        // vibrational constant theta_i
        scalar thetai(const label i) const
        {
            return mix_.thetaVib(i);
        }

        // Park cross-section coefficient
        scalar sigmai
        (
            const scalar TTR 
        );

        // Most Probable Speed
        scalar c_bar_s
        (
            const scalar TTR,
            const label s
        );

        // Total number density
        scalar n_total
        (
            const scalar p,
            const scalar TTR
        );

        // Specie number density
        scalar n_s
        (
            const scalar p,
            const scalar TTR,
            const label s,
            const label celli
        );

        // Partial pressure
        scalar p_s
        (
            const scalar p,
            const scalar TTR,
            const label s,
            const label celli
        );
        // Park Relaxation Time
        scalar T_P_sr
        (
            const scalar p,
            const scalar TTR,
            const label s,
            const label r,
            const label celli
        );

        // Millikan-White A coefficient
        scalar A_sr
        (
            const label s,
            const label r
        );

        // Millikan-White B coefficient
        scalar B_sr
        (
            const label s,
            const label r
        );

        // Millikan-White Relaxation Time
        scalar T_MW_sr
        (
            const scalar p,
            const scalar TTR,
            const label s,
            const label r
        );

        // Total Vibrational-Translational Relaxation Time specie pair s-r
        scalar T_sr
        (
            const scalar p,
            const scalar TTR,
            const label s,
            const label r,
            const label celli
        );

        // Total Vibrational-Translational Relaxation Time for specie s
        scalar T_s
        (
            const scalar p,
            const scalar TTR,
            const label s,
            const label celli
        );

        // Classical Landau-Teller formula
        scalar Q_VT_s
        (
            const scalar p,
            const scalar TTR,
            const scalar TVib,
            const label s,
            const label celli
        );

        PtrList<volScalarField>& correctVibSource
        (
            const volScalarField& p,
            const volScalarField& TTR,
            const PtrList<volScalarField>& TVibSpecies
        );

        PtrList<volScalarField>& correctVTRelaxationTime
        (
            const volScalarField& p,
            const volScalarField& TTR
        );
        
        PtrList<volScalarField>& fields()
        {
            return Q_VTList_;
        }

        const PtrList<volScalarField>& fields() const
        {
            return Q_VTList_;
        }

        PtrList<volScalarField>& tauFields()       
        { 
            return tauVTList_; 
        }
    
        const PtrList<volScalarField>& tauFields() const 
        { 
            return tauVTList_; 
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "VTEnergySource.C"
#endif

#endif

// ************************************************************************* //
