/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
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
    Foam::basic2TMixture

Description
    Foam::basic2TMixture

SourceFiles
    basic2TMixture.C

\*---------------------------------------------------------------------------*/

#ifndef basic2TMixture_H
#define basic2TMixture_H

#include "compositionInterface.H"
#include "PtrList.H"
#include "volFields.H"
#include "speciesTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;
class dictionary;

/*---------------------------------------------------------------------------*\
                        Class basic2TMixture Declaration
\*---------------------------------------------------------------------------*/

class basic2TMixture
:
    public compositionInterface
{

public:

    //- The base class of the mixture
    typedef basic2TMixture basic2TMixtureType;

 PtrList<volScalarField>& Y() override
    {
        FatalErrorInFunction
            << "Y() called on pure mixture — mixtureCheck must be false"
            << exit(FatalError);
        static PtrList<volScalarField> dummy;
        return dummy;
    }

    const PtrList<volScalarField>& Y() const override
    {
        FatalErrorInFunction
            << "Y() const called on pure mixture"
            << exit(FatalError);
        static PtrList<volScalarField> dummy;
        return dummy;
    }

    const speciesTable& species() const override
    {
        FatalErrorInFunction
            << "species() called on pure mixture"
            << exit(FatalError);
        static speciesTable dummy;
        return dummy;
    }

    bool isSpecieMolecular(label) const override
    {
        FatalErrorInFunction
            << "isSpecieMolecular() called on pure mixture"
            << exit(FatalError);
        return false;
    }

    scalar thetaVib(label) const override
    {
        FatalErrorInFunction
            << "thetaVib() called on pure mixture"
            << exit(FatalError);
        return 0;
    }

    scalar EsVib(label, scalar, scalar, scalar, scalar) const override
    {
        FatalErrorInFunction
            << "EsVib() called on pure mixture"
            << exit(FatalError);
        return 0;
    }

    scalar Hs(label, scalar, scalar) const override
    {
        FatalErrorInFunction
            << "Hs() called on pure mixture"
            << exit(FatalError);
        return 0;
    }

    scalar TEVib(label, scalar, scalar, scalar, scalar) const override
    {
        FatalErrorInFunction
            << "TEVib() called on pure mixture"
            << exit(FatalError);
        return 0;
    }

    scalar W(label) const override
    { 
        FatalErrorInFunction 
            << "W() called on pure mixture" 
            << exit(FatalError); 
        return 0; 
    }

    scalar R(label) const override
    { 
        FatalErrorInFunction 
            << "R() called on pure mixture" 
            << exit(FatalError); 
        return 0; 
    }

    scalar Hc(label) const override
    { 
        FatalErrorInFunction 
            << "Hc() called on pure mixture" 
            << exit(FatalError); 
        return 0; 
    }

    scalar CpTR(label, scalar, scalar) const override
    { 
        FatalErrorInFunction 
            << "CpTR() called on pure mixture" 
            << exit(FatalError); 
        return 0; 
    }

    scalar CvT(label, scalar, scalar) const override
    { 
        FatalErrorInFunction 
            << "CvT() called on pure mixture" 
            << exit(FatalError); 
            return 0; 
    }

    scalar H(label, scalar, scalar) const override
    { 
        FatalErrorInFunction 
            << "H() called on pure mixture" 
            << exit(FatalError); 
        return 0; 
    }

    scalar Ha(label, scalar, scalar) const override
    { 
        FatalErrorInFunction 
            << "Ha() called on pure mixture" 
            << exit(FatalError); 
        return 0; 
    }

    scalar S(label, scalar, scalar) const override
    { 
        FatalErrorInFunction 
            << "S() called on pure mixture" 
            << exit(FatalError); 
        return 0; 
    }

    scalar EsT(label, scalar, scalar) const override
    { 
        FatalErrorInFunction 
            << "EsT() called on pure mixture" 
            << exit(FatalError); 
        return 0; 
    }

    scalar EsR(label, scalar, scalar) const override
    { 
        FatalErrorInFunction 
            << "EsR() called on pure mixture" 
            << exit(FatalError); 
        return 0; 
    }

    scalar G(label, scalar, scalar) const override
    { 
        FatalErrorInFunction 
            << "G() called on pure mixture" 
            << exit(FatalError); 
        return 0; 
    }

    scalar mu(label, scalar) const override
    { 
        FatalErrorInFunction 
            << "mu() called on pure mixture" 
            << exit(FatalError); 
        return 0; 
    }

    scalar kappaTR(label, scalar, scalar) const override
    { 
        FatalErrorInFunction 
            << "kappaTR() called on pure mixture" 
            << exit(FatalError); 
        return 0; 
    }

    scalar kappaVib(label, scalar, scalar, scalar, scalar) const override
    { 
        FatalErrorInFunction 
            << "kappaVib() called on pure mixture" 
            << exit(FatalError); 
        return 0; 
    }

    scalar alphah(label, scalar, scalar) const override
    { 
        FatalErrorInFunction 
            << "alphah() called on pure mixture" 
            << exit(FatalError); 
        return 0; 
    }

    scalar rho(label, scalar, scalar) const override
    { 
        FatalErrorInFunction 
            << "rho() called on pure mixture" 
            << exit(FatalError); 
        return 0; 
    }
    // Constructors

        //- Construct from dictionary, mesh and phase name
        basic2TMixture(const dictionary&, const fvMesh&, const word&)
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
