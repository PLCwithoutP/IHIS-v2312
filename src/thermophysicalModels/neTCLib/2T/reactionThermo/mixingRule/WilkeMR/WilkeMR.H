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
    Foam::WilkeMR

Description
    Foam::WilkeMR

SourceFiles
    WilkeMR.C

\*---------------------------------------------------------------------------*/

#ifndef WilkeMR_H
#define WilkeMR_H

#include "baseMR.H"
#include "speciesTable.H"
#include "wordList.H"
#include "dimensionedScalar.H"
#include "volFields.H"
#include "IOobject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class WilkeMR Declaration
\*---------------------------------------------------------------------------*/
template<class MixtureType>
class WilkeMR
:
    public baseMR<MixtureType>
{
    

        //- Construct as copy (not implemented)
        WilkeMR(const WilkeMR&);

private:

    // Private data
    const MixtureType& mix_;
    const speciesTable& names_;
    volScalarField Xi_;
    List<scalarField> XiCells_;

public:


    // Constructors

        //- Construct from multiComponent2TMixture& 
        WilkeMR(const MixtureType& mix);

    // Member functions

        void precomputeXi();

        // Returns Xi field for species speciei — O(S*N) cost.
        // NEVER call inside a per-cell loop. Use precomputeXi() + Xi(s,celli) instead.
        volScalarField& computeXiFromYi
        (
            const label speciei
        );

        scalar scalingFactor
        (
            const label speciei,       // specie index for which φ_i is computed
            const label celli,   // cell index
            const scalar T       // local temperature T in that cell
        );

        template<class QGetter>
        scalar QCell
        (
            const label celli,
            const scalar p,
            const scalar T,
            const QGetter& getQ
        );

        scalar muCell
        (
            const label celli,
            const scalar p,
            const scalar TTR
        );

        scalar kappaTRCell
        (
            const label celli,
            const scalar p,
            const scalar TTR
        );

        scalar alphaTRCell
        (
            const label celli,
            const scalar p,
            const scalar TTR
        );

        //- Return the instantiated type name
        static word typeName()
        {
            return "WilkeMR";
        }

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

        // dynamic viscosity μ_i at given T (scalar)
        scalar mui(const label i, const scalar T) const
        {
            return mix_.mu(i, T);
        }

        volScalarField& Xi()
        {
            return Xi_;
        }

        const volScalarField& Xi() const
        {
            return Xi_;
        }

        // Cached mole fraction accessor — valid after precomputeXi() has been called
        inline scalar Xi(const label speciei, const label celli) const
        {
            return XiCells_[speciei][celli];
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "WilkeMR.C"
#endif

#endif

// ************************************************************************* //
