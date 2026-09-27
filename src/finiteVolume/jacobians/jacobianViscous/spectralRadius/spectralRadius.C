/*---------------------------------------------------------------------------*\

    IHIS: In-house Hypersonics Implicit solver

    Copyright (C) 2018 Johan Heyns - CSIR, South Africa
    Copyright (C) 2018 Oliver Oxtoby - CSIR, South Africa

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

\*---------------------------------------------------------------------------*/

#include "spectralRadius.H"
#include "emptyFvPatch.H"
#include "transformFvPatchFields.H"
#include "fvjOperators.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceInterpolate.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(spectralRadius, 0);
addToRunTimeSelectionTable(jacobianViscous, spectralRadius, dictionary);


// * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * * * * //

spectralRadius::spectralRadius
(
    const dictionary& dict,
    const fvMesh& mesh,
    const volScalarField& rho,
    const volVectorField& rhoU,
    const volScalarField& rhoE,
    const fluid2TThermo& thermo
)
 :
    jacobianViscous(typeName, dict),
    mesh_(mesh),
    rho_(rho),
    rhoU_(rhoU),
    rhoE_(rhoE),
    thermo_(thermo)
{}


// * * * * * * * * *  Private Member Functions  * * * * * * * * * * * * * * * //

// Assemble full matrix with diagonal and off-diagonal contribution from the
// flux terms.
void spectralRadius::addFluxTerms(compressibleJacobianMatrix& jacobian)
{
    fvjMatrix<scalar>& dContByRho = jacobian.dSByS(0,0);
    fvjMatrix<tensor>& dMomByRhoU = jacobian.dVByV(0,0);
    fvjMatrix<scalar>& dEnergyByRhoE = jacobian.dSByS(1,1);

    // Viscous part
    // Add laplacian term similar to Lax-Friedrich stabilisation term
    tmp<fvjMatrix<scalar> > stab = fvj::laplacian(0.5*lambdaVisc_, geometricOneField());

    dContByRho -= stab();
    dMomByRhoU -= stab()*tensor::I;
    dEnergyByRhoE -= stab;
}

void spectralRadius::boundaryJacobian
(
    label patchi,
    tmp<scalarField>& dContFluxdRho, tmp<vectorField>& dContFluxdRhoU, tmp<scalarField>& dContFluxdRhoE,
    tmp<vectorField>& dMomFluxdRho, tmp<tensorField>& dMomFluxdRhoU, tmp<vectorField>& dMomFluxdRhoE,
    tmp<scalarField>& dEnergyFluxdRho, tmp<vectorField>& dEnergyFluxdRhoU, tmp<scalarField>& dEnergyFluxdRhoE
)
{
    dContFluxdRho = tmp<scalarField>(new scalarField(mesh_.boundary()[patchi].patch().size(), Zero));
    dContFluxdRhoU = tmp<vectorField>(new vectorField(mesh_.boundary()[patchi].patch().size(), Zero));
    dContFluxdRhoE = tmp<scalarField>(new scalarField(mesh_.boundary()[patchi].patch().size(), Zero));
    dMomFluxdRho = tmp<vectorField>(new vectorField(mesh_.boundary()[patchi].patch().size(), Zero));
    dMomFluxdRhoU = tmp<tensorField>(new tensorField(mesh_.boundary()[patchi].patch().size(), Zero));
    dMomFluxdRhoE = tmp<vectorField>(new vectorField(mesh_.boundary()[patchi].patch().size(), Zero));
    dEnergyFluxdRho = tmp<scalarField>(new scalarField(mesh_.boundary()[patchi].patch().size(), Zero));
    dEnergyFluxdRhoU = tmp<vectorField>(new vectorField(mesh_.boundary()[patchi].patch().size(), Zero));
    dEnergyFluxdRhoE = tmp<scalarField>(new scalarField(mesh_.boundary()[patchi].patch().size(), Zero));
}


// * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * * //

void spectralRadius::addViscousJacobian(compressibleJacobianMatrix& jacobian)
{
    muLam_    = fvc::interpolate(thermo_.mu());
    alphaLam_ = fvc::interpolate(thermo_.kappaTR()/thermo_.CvTR());

    const surfaceScalarField& muLam    = muLam_();
    const surfaceScalarField& alphaLam = alphaLam_();

    gamma_ = thermo_.gamma();

    // Viscous wave speed (Luo et al, 2001)
    // CHECK: Blazek proposes an alternative approximation of viscous spectrial
    // radius
    // lambdaVisc_ = (muEff+alphaEff)/fvc::interpolate(rho_);
//    lambda_() += (muEff)/fvc::interpolate(rho_); // Luo et al. (2000) do not include alphaEff
    // Viscous spectral radius (Luo et al. 2001); Blazek Sec. 6.1.4 gives an
    // alternative with a max(4/3, gamma/Pr) prefactor
    lambdaVisc_ = (muLam + alphaLam)/fvc::interpolate(rho_);
    addFluxTerms(jacobian);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
