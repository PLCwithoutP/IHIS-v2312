// Species transport (frozen, decoupled, advective form).
// Spliced into ihisModule::outerIteration() before updateFields.H.
{
    PtrList<volScalarField>& Y = neThermo.composition().Y();

    if (Y.size() > 1)
    {
        // Common-limiter weights shared by all species (sum of Y_f = 1)
        tmp<reconstructionScheme<scalar>> tRecon
        (
            reconstructionScheme<scalar>::New
            (
                mesh,
                mesh.interpolationScheme("reconstruct(rho)")
            )
        );

        tmp<surfaceScalarField> wL, wR;
        forAll(Y, i)
        {
            tmp<surfaceScalarField> psiL, psiR;
            tRecon().limiter(Y[i], psiL, psiR);

            if (i == 0)
            {
                wL = psiL;
                wR = psiR;
            }
            else
            {
                wL.ref() = min(wL(), psiL());
                wR.ref() = min(wR(), psiR());
            }
        }

        tRecon().weights(Y[0], mesh.surfaceInterpolation::weights(), wL, wR);

        const dimensionedScalar zeroW(dimless, Zero);
        const dimensionedScalar oneW(dimless, 1.0);
        wL.ref() = max(min(wL(), oneW), zeroW);
        wR.ref() = max(min(wR(), oneW), zeroW);
        
        const volScalarField ddtRho(fvc::ddt(rho));              
        const volScalarField divPhi(fvc::surfaceIntegrate(phi)); 
        upwind<scalar> upwindYf(mesh, phi);                      



        // Reaction rates and Jacobian diagonal, once per
        // outer iteration, at one state (rho^{k+1}, T^k, Y^k).
        // Temperature feedback at fixed rho and rhoE,
        // complement-balanced (datum-free):
        //   dT/dY_s = -(e_s - e)/((1 - Y_s) c_v),
        //   e_s = EsT_s(p,T) + Hc_s (absolute, E1a),  e = sum_j Y_j e_j.
        PtrList<scalarField> dTdY;   // [K per unit mass fraction], cells only

        if (chemistry_)
        {
            chemistry_->calculateWithJacobianDiag
            (
                rho.internalField(),
                dRRdRhoY_,
                dRRdT_
            );

            const compositionInterface& comp = neThermo.composition();
            const scalarField& pI = neThermo.p().primitiveField();
            const scalarField& TI = neThermo.TTR().primitiveField();
            const tmp<volScalarField> tCv(neThermo.CvTR());
            const scalarField& cvI = tCv().primitiveField();
            const label nCells = TI.size();

            // Absolute species energies e_s and mixture energy e (state k)
            PtrList<scalarField> es(Y.size());
            scalarField eMix(nCells, Zero);
            forAll(Y, i)
            {
                es.set(i, new scalarField(nCells));
                scalarField& esI = es[i];
                const scalarField& YI = Y[i].primitiveField();
                const scalar Hci = comp.Hc(i);
                forAll(esI, celli)
                {
                    esI[celli] = comp.EsT(i, pI[celli], TI[celli]) + Hci;
                    eMix[celli] += YI[celli]*esI[celli];
                }
            }

            // dT/dY_s; zero where the species is (almost) pure
            const scalar oneMinusYMin = 1e-6;
            dTdY.setSize(Y.size());
            forAll(Y, i)
            {
                dTdY.set(i, new scalarField(nCells, Zero));
                scalarField& dTdYI = dTdY[i];
                const scalarField& YI = Y[i].primitiveField();
                const scalarField& esI = es[i];
                forAll(dTdYI, celli)
                {
                    const scalar oneMinusY = 1.0 - YI[celli];
                    if (oneMinusY > oneMinusYMin)
                    {
                        dTdYI[celli] =
                           -(esI[celli] - eMix[celli])/(oneMinusY*cvI[celli]);
                    }
                }
            }
        }

        forAll(Y, i)
        {
            volScalarField& Yi = Y[i];

            // Non-owning tmps: reconstruct() clears the weights it receives
            tmp<surfaceScalarField> YL, YR;
            reconstructionScheme<scalar>::reconstruct
            (
                Yi,
                tmp<surfaceScalarField>(wL()),
                tmp<surfaceScalarField>(wR()),
                YL, YR
            );

            
            tmp<surfaceScalarField> Yf = surfaceFieldSelect(YL, YR, phi, 0.0);         
            tmp<surfaceScalarField> YfUpwind = upwindYf.interpolate(Yi);

            // Implicit upwind + second-order deferred correction
            fvScalarMatrix YiEqn
            (
                fvm::ddt(rho, Yi) - fvm::Sp(ddtRho, Yi)          // real time (advective) + pseudo (dualTime)
              + fvm::div(phi, Yi, "div(phi,Yi)")                  // implicit upwind
              - fvm::Sp(divPhi, Yi)                               // advective-form convection pair
              + fvc::surfaceIntegrate(phi*(Yf() - YfUpwind()))    // deferred correction
            );
            
            if (trhoD.valid())
            {
                YiEqn -= fvm::laplacian(trhoD(), Yi, "laplacian(rhoD,Yi)");
            }
            

            // Point-implicit chemical source, Newton-diagonal
            //   S_i(Y) ~ RR_i^k + a_i (Y_i - Y_i^k),
            //   a_i = rho*dRR_i/d(rho Y_i) + dRR_i/dT * dT/dY_i   [kg/m3/s]
            // Written on the LHS: -S_i. SuSp makes -a_i*Y implicit where
            // a_i < 0 (stabilising) and explicit where a_i > 0; the explicit
            // parts then cancel to RR_i^k. Yi is still Y^k here (Jacobi).
            if (chemistry_)
            {
                volScalarField::Internal aS
                (
                    "aS",
                    rho.internalField()*dRRdRhoY_[i]
                );
                aS.field() += dRRdT_[i].field()*dTdY[i];

                YiEqn -= chemistry_->RR(i) - aS*Yi.internalField();
                YiEqn += fvm::SuSp(-aS, Yi);
            }

            YiEqn.solve(mesh.solverDict("Yi"));
        }

        // Clip, renormalise, refresh patches (before updateEZero/correct)
        volScalarField Ysum("Ysum", 0.0*Y[0]);

        scalar minYRaw = GREAT;
        forAll(Y, i)
        {
            minYRaw = min(minYRaw, gMin(Y[i].primitiveField()));
            Y[i] = max(Y[i], dimensionedScalar(dimless, Zero));
            Ysum += Y[i];
        }
        const scalar sumYErr = gMax(mag(Ysum.primitiveField() - 1.0)());

        const dimensionedScalar sumFloor(dimless, SMALL);
        forAll(Y, i)
        {
            Y[i] /= max(Ysum, sumFloor);
            Y[i].correctBoundaryConditions();
        }
        

        // Mass conservation of the chemistry,
        //   eps = max|sum_i RR_i| / max_i max|RR_i|   (round-off if W additive)
        Info<< "Species: min(Y) = " << minYRaw
            << ", max|sum(Y) - 1| = " << sumYErr;

        if (chemistry_)
        {
            const PtrList<volScalarField::Internal>& RR = chemistry_->RR();

            scalarField sumRR(RR[0].field());
            scalar maxRR = gMax(mag(RR[0].field()));
            for (label i = 1; i < RR.size(); ++i)
            {
                sumRR += RR[i].field();
                maxRR = max(maxRR, gMax(mag(RR[i].field())));
            }

            Info<< ", max|sum(RR)|/max|RR| = "
                << gMax(mag(sumRR))/max(maxRR, VSMALL);
        }

        Info<< endl;
    }
}