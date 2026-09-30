/*---------------------------------------------------------------------------*\

    IHIS: In-house Hypersonics Implicit solver

    Copyright (C) 2014-2018 Johan Heyns - CSIR, South Africa
    Copyright (C) 2014-2018 Oliver Oxtoby - CSIR, South Africa
    Copyright (C) 2011-2013 OpenFOAM Foundation

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

    fv::options& fvOptions(*fvOptions_);
    volScalarField rhoSource(fvOptions(rho) & rho);
    volVectorField rhoUSource((fvOptions(rho, U) & U) + rhoSource*U);
    volScalarField rhoESource((fvOptions(rho, neThermo.eTR()) & neThermo.eTR()) + (rhoUSource & U) + rhoSource*neThermo.eTR());

    // Construct cont eqn: [A] \Delta W = [R]
    rhoR =
    (
        -fvc::div(phi) + rhoSource
    );

    // Construct mom eqn: [A] \Delta W = [R]
    rhoUR =
    (
        -fvc::div(phiUp) + rhoUSource
    );

    // Construct energy eqn: [A] \Delta W = [R]
    rhoER =
    (
        -fvc::div(phiEp) + rhoESource
    );

    // Viscous flow
    if (!inviscid_)
    {
        // Laminar transport properties (no turbulence model in IHIS)
        volScalarField muLam("muLam", neThermo.mu());
        volScalarField alphaLam("alphaLam", neThermo.kappaTRWilke()/neThermo.CvTR());

        volTensorField& tauMC = tauMC_();
        tauMC = muLam*dev2(Foam::T(fvc::grad(U)));
        rhoUR += fvc::laplacian(muLam,U);
        rhoUR += fvc::div(tauMC);

        surfaceScalarField sigmaDotU
        (
            "sigmaDotU",
            (
              fvc::interpolate(muLam)*fvc::snGrad(U)
              + (mesh.Sf()/mesh.magSf() & fvc::interpolate(tauMC))
            )
            & Up
        );
        rhoER += fvc::div(sigmaDotU);

        volScalarField eCalc("eCalc", rhoE/rho-0.5*magSqr(U)); // prevent e inheriting BC from T in neThermo
        rhoER += fvc::laplacian(alphaLam, eCalc, "laplacian(alphaLam,e)");
    }

    #include "cellDebug.H"

    // Pseudo and real-time contribution

    tmp<fvScalarMatrix> ddtRho = fvm::ddt(rho);
    tmp<fvVectorMatrix> ddtRhoU = fvm::ddt(rhoU);
    tmp<fvScalarMatrix> ddtRhoE = fvm::ddt(rhoE);

#if OPENFOAM >= 1712
    // Work around bug causing timeIndex not to be stored for ddt0 field (CrankNicolson)
    // Seems to be caused by commented out line in GeometricField::storeOldTimes(),
    // commit 44ce956798
    if (mesh.foundObject<volScalarField>("ddt0(rho)"))
    {
        const_cast<volScalarField&>
        (
            mesh.lookupObject<volScalarField>("ddt0(rho)")
        ).timeIndex() = mesh.time().timeIndex();
    }
    if (mesh.foundObject<volVectorField>("ddt0(rhoU)"))
    {
        const_cast<volVectorField&>
        (
            mesh.lookupObject<volVectorField>("ddt0(rhoU)")
        ).timeIndex() = mesh.time().timeIndex();
    }
    if (mesh.foundObject<volScalarField>("ddt0(rhoE)"))
    {
        const_cast<volScalarField&>
        (
            mesh.lookupObject<volScalarField>("ddt0(rhoE)")
        ).timeIndex() = mesh.time().timeIndex();
    }
#endif

    if (!steadyState_)
    {
        // NOTE: Solving the system A \Delta W^\tau = B and therefore the
        // source term needs to be adjusted as below.
        // For steady state analysis this term reduces it to zero.
        //
        // NOTE: The residual is in strong form and therefore the division by V
        rhoR.primitiveFieldRef() -= (ddtRho().diag()*rho-ddtRho().source())/mesh.V();
        rhoUR.primitiveFieldRef() -= (ddtRhoU().diag()*rhoU-ddtRhoU().source())/mesh.V();
        rhoER.primitiveFieldRef() -= (ddtRhoE().diag()*rhoE-ddtRhoE().source())/mesh.V();
    }

    // NOTE: These should be equal if ddtScheme is consistent
    ddtCoeff = max(max(ddtRho().diag(),ddtRhoU().diag()),ddtRhoE().diag())/mesh.V();
