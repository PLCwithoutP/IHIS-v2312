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
    Foam::VVEnergySource

Description
    Foam::VVEnergySource

SourceFiles
    VVEnergySource.C

\*---------------------------------------------------------------------------*/

#ifndef VVEnergySource_H
#define VVEnergySource_H

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
                         Class VVEnergySource Declaration
\*---------------------------------------------------------------------------*/
template<class MixtureType, class MixingRule>
class VVEnergySource
:
    public baseEnergySource<MixtureType, MixingRule>
{

        //- Construct as copy (not implemented)
        VVEnergySource(const VVEnergySource&);

private:

    // Private data
    MixtureType& mix_;
    MixingRule& mr_;

    const speciesTable& names_;
    PtrList<volScalarField> Q_VVList_;

    const scalar boltzmann_const_; // Boltzmann constant
    const scalar avagadro_const_; // Avagadro constant
    // Per species-pair Knab coefficients [nSpec x nSpec]
    scalarSquareMatrix sigma_;   // VV cross section [m²]
    scalarSquareMatrix P_;       // Exchange probability [-]
    bool knabCoeffsLoaded_;

    void loadKnabCoeffs(const fvMesh& mesh);

public:

    // Constructors

        //- Construct from multiComponent2TMixture& 
        VVEnergySource
        (
            MixtureType& mix,
            MixingRule& MR
        );

    // Member functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "VVEnergySource";
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

        // Classical Landau-Teller formula
        scalar Q_VV_s
        (
            const scalar p,
            const scalar TTR,
            const PtrList<volScalarField>& TVibSpecies,
            const label s,
            const label celli
        );

        PtrList<volScalarField>& correctVibVibSource
        (
            const volScalarField& p,
            const volScalarField& TTR,
            const PtrList<volScalarField>& TVibSpecies
        );
        
        PtrList<volScalarField>& fields()
        {
            return Q_VVList_;
        }

        const PtrList<volScalarField>& fields() const
        {
            return Q_VVList_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "VVEnergySource.C"
#endif

#endif

// ************************************************************************* //
