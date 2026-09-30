/*---------------------------------------------------------------------------*\

    IHIS: In-house Hypersonics Implicit solver

    Copyright (C) 2015-2017 Oliver Oxtoby - CSIR, South Africa
    Copyright (C) 2015-2017 Johan Heyns - CSIR, South Africa

-------------------------------------------------------------------------------
License
    This file is part of IHIS.

    IHIS is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    IHIS is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with IHIS.  If not, see <http://www.gnu.org/licenses/>.

Class
    discreteLaxFriedrichs

Description
    Lax-Friedrich approximation of Inviscid LHS Jacobian
    which tries to be more faithful to the discrete form of the flux.

Note
    Causes problems for LUSGS preconditioner - use separate centralLaxFriedrichs
    scheme. 
    NOTE (IHIS): not mixture-ready. The energy row is built from
    aStar*rhoE - rho*theta and ef = rhoE/rho - 0.5|U|^2, which assume a zero
    energy datum (p = (gamma-1)*rho*e). With formation energy (Ns > 1) or
    vibrational energy inside rhoE, dp/drho can have the wrong sign. Use
    laxFriedrichs or centralLaxFriedrichs for mixtures. This affects the
    left-hand side only (convergence, not the converged solution).

Authors
    Oliver Oxtoby
    Johan Heyns
        Council for Scientific and Industrial Research, South Africa

\*---------------------------------------------------------------------------*/

#ifndef discreteLaxFriedrichs_H
#define discreteLaxFriedrichs_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "jacobianInviscid.H"
#include "jacobianMatrix.H"
#include "fluid2TThermo.H"
#include "fvjMatrix.H"
#include "laxFriedrichs.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class discreteLaxFriedrichs Declaration
\*---------------------------------------------------------------------------*/

class discreteLaxFriedrichs
:
    public laxFriedrichs
{

    protected:

        const fvMesh& mesh_;
        const volScalarField& rho_;
        const volVectorField& rhoU_;
        const volScalarField& rhoE_;
        const fluid2TThermo& thermo_;
        tmp< surfaceScalarField > lambdaConv_;
        tmp< volScalarField > gamma_;
        autoPtr< tmp<surfaceScalarField> > meshPhi_;


        //- Private member functions

            void addFluxTerms(compressibleJacobianMatrix& jacobian);


    public:
        //- Runtime type information
        TypeName("discreteLaxFriedrichs");

        discreteLaxFriedrichs
        (
            const dictionary& dict,
            const fvMesh& mesh,
            const volScalarField& rho,
            const volVectorField& rhoU,
            const volScalarField& rhoE,
            const fluid2TThermo& thermo
        );

        virtual void addInviscidJacobian(compressibleJacobianMatrix& jacobian);

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
