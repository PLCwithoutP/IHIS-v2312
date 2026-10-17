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
        
        // Physical, non-empty patches only (wall-flux storage)
        const auto skipWall = [&mesh](const label patchi)
        {
            const fvPatch& fp = mesh.boundary()[patchi];
            return fp.coupled() || fp.size() == 0;
        };

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
        
        // Wall viscous traction and work as applied
        forAll(mesh.boundary(), patchi)
        {
            if (skipWall(patchi)) continue;
            const vectorField nf
            (
                mesh.Sf().boundaryField()[patchi]
               /mesh.magSf().boundaryField()[patchi]
            );
            wallTau_->boundaryFieldRef()[patchi] ==
                muLam.boundaryField()[patchi]*U.boundaryField()[patchi].snGrad()
              + (nf & tauMC.boundaryField()[patchi]);
            wallSigmaDotU_->boundaryFieldRef()[patchi] ==
                sigmaDotU.boundaryField()[patchi];
        }

        // Vibrational conduction.
        //   GVib_m = (kappaVib_m/CvVib_m)/phi_m at T* = max(Tv_m, theta_m/20), state k
        //   qVib_m = X_m*GVib_m*snGrad(e_v,m)*|Sf|   [W], shared with vibEqns.H
        // Added here once, so it serves both the Ns > 1 and the Ns = 1 branch below.
        if (twoT_)
        {
            const compositionInterface& comp = neThermo.composition();
            const PtrList<volScalarField>& Yv = comp.Y();
            const label Nm = molI_.size();

            if (GVib_.size() != Nm) GVib_.setSize(Nm);
            if (qVib_.size() != Nm) qVib_.setSize(Nm);
            if (GammaVib_.size() != Nm) GammaVib_.setSize(Nm);
            if (wallQVibM_.size() != Nm) wallQVibM_.setSize(Nm);

            // 1/W = sum_s Y_s/W_s  (mixture molar mass, boundary included)
            volScalarField invW("invW", Yv[0]/comp.W(0));
            for (label s = 1; s < Yv.size(); ++s)
            {
                invW += Yv[s]/comp.W(s);
            }

            wallQVib_() == dimensionedScalar(dimPower, Zero);
            forAll(molI_, m)
            {
                const label s = molI_[m];
                const scalar th = comp.thetaVib(s);

                // Evaluation temperature and Wilke-weighted conductance
                const volScalarField TStar
                (
                    "TStar",
                    max(TVibSpecies_[s], dimensionedScalar(dimTemperature, th/20))
                );
                GVib_.set(m, neThermo.kappaVibWilke(s, TStar).ptr());

                // Mole fraction X_m = (Y_m/W_m)/(1/W); exactly 1 for Ns = 1
                const volScalarField Xm("Xm", (Yv[s]/comp.W(s))/invW);

                // Conduction face flux, Gamma_e,m = X_m*GVib_m
                qVib_.set
                (
                    m,
                    new surfaceScalarField
                    (
                        "qVib_" + comp.species()[s],
                        fvc::interpolate(Xm*GVib_[m])
                       *fvc::snGrad(evVib_[m])*mesh.magSf()
                    )
                );

                // Into the total-energy residual, same sign as div(kappa grad T)
                rhoER += fvc::surfaceIntegrate(qVib_[m]);

                // Conductance for neSmoluchowskiJumpTv
                GammaVib_.set
                (
                    m,
                    new volScalarField
                    (
                        IOobject
                        (
                            "GammaVib_" + comp.species()[s],
                            runTime.timeName(),
                            mesh,
                            IOobject::NO_READ,
                            IOobject::NO_WRITE
                        ),
                        Xm*GVib_[m]
                    )
                );

                // Wall vibrational conduction as applied,
                // summed and per molecule
                if (!wallQVibM_.set(m))
                {
                    wallQVibM_.set
                    (
                        m,
                        new surfaceScalarField
                        (
                            IOobject
                            (
                                "neWall:qVib_" + comp.species()[s],
                                runTime.timeName(),
                                mesh,
                                IOobject::NO_READ,
                                IOobject::NO_WRITE
                            ),
                            mesh,
                            dimensionedScalar(dimPower, Zero)
                        )
                    );
                }
                forAll(mesh.boundary(), patchi)
                {
                    if (skipWall(patchi)) continue;
                    wallQVibM_[m].boundaryFieldRef()[patchi] ==
                        qVib_[m].boundaryField()[patchi];
                    wallQVib_->boundaryFieldRef()[patchi] +=
                        qVib_[m].boundaryField()[patchi];
                }
            }
        }

        volScalarField eCalc("eCalc", rhoE/rho-0.5*magSqr(U)); // prevent e inheriting BC from T in neThermo
        if (trhoD.valid())
        {
            // Frozen-mixture energy diffusion (viscous, Ns > 1):
            //   div(kappa grad T) - div(sum_s h_s J_s),   J_s = -rhoD grad Y_s
            // Reuses the "laplacian(alphaLam,e)" scheme label (no case change).
            rhoER += fvc::laplacian
            (
                neThermo.kappaTRWilke(),
                neThermo.TTR(),
                "laplacian(alphaLam,e)"
            );

            const PtrList<volScalarField>& Ys = neThermo.composition().Y();
            const label last = Ys.size() - 1;
            const surfaceScalarField rhoDf(fvc::interpolate(trhoD()));

            // Wall TR conduction as applied (kappa grad T, Ns > 1)
            forAll(mesh.boundary(), patchi)
            {
                if (skipWall(patchi)) continue;
                wallQTR_->boundaryFieldRef()[patchi] ==
                    neThermo.kappaTRWilke(patchi)
                   *neThermo.TTR().boundaryField()[patchi].snGrad()
                   *mesh.magSf().boundaryField()[patchi];
            }

            surfaceScalarField phiJSum
            (
                IOobject("phiJSum", runTime.timeName(), mesh),
                mesh,
                dimensionedScalar(dimMass/dimTime, Zero)
            );
            surfaceScalarField hJ
            (
                IOobject("hJ", runTime.timeName(), mesh),
                mesh,
                dimensionedScalar(dimPower, Zero)
            );
            #if OPENFOAM >= 1712
            phiJSum.setOriented(true);   // face fluxes are oriented (like phi)
            hJ.setOriented(true);
            #endif

            // Keep each species' face flux for the vibrational
            // equations (D5: rhoE and e_v share the same phiJ)
            if (twoT_ && phiJ_.size() != Ys.size())
            {
                phiJ_.setSize(Ys.size());
            }

            // Face diffusive mass fluxes; last species closes sum_s phiJ_s = 0
            for (label s = 0; s < last; ++s)
            {
                const surfaceScalarField phiJs
                (
                    -rhoDf*fvc::snGrad(Ys[s])*mesh.magSf()
                );
                hJ += fvc::interpolate(hs[s])*phiJs;
                phiJSum += phiJs;
                if (twoT_)
                {
                    phiJ_.set(s, new surfaceScalarField("phiJ_" + Ys[s].name(), phiJs));
                }
            }
            hJ -= fvc::interpolate(hs[last])*phiJSum;
            
            if (twoT_)
            {
                // Closure: phiJ_last = -sum_{s<last} phiJ_s, exactly as in hJ
                phiJ_.set(last, new surfaceScalarField("phiJ_" + Ys[last].name(), -phiJSum));
            }

            rhoER -= fvc::surfaceIntegrate(hJ);
            
            // Wall species-enthalpy diffusion flux as applied
            forAll(mesh.boundary(), patchi)
            {
                if (skipWall(patchi)) continue;
                wallHJ_->boundaryFieldRef()[patchi] == hJ.boundaryField()[patchi];
            }
        }
        else
        {
            volScalarField eCalc("eCalc", rhoE/rho-0.5*magSqr(U)); // prevent e inheriting BC from T in neThermo
            if (twoT_)
            {
                // 2T: conduct the TR part only; vibrational conduction is qVib_ above
                eCalc -= eVibMix_();
            }
            rhoER += fvc::laplacian(alphaLam, eCalc, "laplacian(alphaLam,e)");
            
            // Wall TR conduction as applied (e-form, Ns = 1)
            forAll(mesh.boundary(), patchi)
            {
                if (skipWall(patchi)) continue;
                wallQTR_->boundaryFieldRef()[patchi] ==
                    alphaLam.boundaryField()[patchi]
                   *eCalc.boundaryField()[patchi].snGrad()
                   *mesh.magSf().boundaryField()[patchi];
            }
        }
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
