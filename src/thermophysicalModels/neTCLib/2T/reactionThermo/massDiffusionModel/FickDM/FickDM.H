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
    Foam::FickDM

Description
    Foam::FickDM

SourceFiles
    FickDM.C

\*---------------------------------------------------------------------------*/

#ifndef FickDM_H
#define FickDM_H

#include "baseDM.H"
#include "speciesTable.H"
#include "wordList.H"
#include "dimensionedScalar.H"
#include "volFields.H"
#include "IOobject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class FickDM Declaration
\*---------------------------------------------------------------------------*/
template<class MixtureType>
class FickDM
:
    public baseDM<MixtureType>
{
        //- Construct as copy (not implemented)
        FickDM(const FickDM&);

private:

    // Private data
    const MixtureType& mix_;
    const speciesTable& names_;
    const scalar Le_num_;

public:


    // Constructors

        //- Construct from multiComponent2TMixture& 
        FickDM(const MixtureType& mix);

    // Member functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "FickDM";
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

        scalar DeffCell
        (
            const label celli,
            const scalar pCell,
            const scalar TTRCell,
            const scalar rhoCell,
            const scalar kappaTRCell,
            const scalar CpTRCell
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "FickDM.C"
#endif

#endif

// ************************************************************************* //
