// Vibrational energy transport with V-T, V-V and C-V sources (Phase 6).
// Spliced into ihisModule::outerIteration() after speciesEqns.H, before
// updateFields.H. Unknown eps_m = Y_m*e_v,m [J/kg of mixture], one equation per
// molecule, Jacobi across molecules. Sources are point-implicit (Q^k + a(eps-eps^k),
// a <= 0); never call the library correctCVSource (mixture TVib). If inner
// iterations stall in multi-molecule cases, the fallback is a cell-local Nm x Nm block.

if (twoT_)
{
    const compositionInterface& comp = neThermo.composition();
    const PtrList<volScalarField>& Yv = comp.Y();

    // (a) Faces from TVib_<m> with the new Y_b (after speciesEqns.H), so the
    //     fixedValue eps patches carry Y_b^{k+1} * e_v(Tv_b)
    vibBoundaries();

    // (b) Common limiter over the e_v,m fields (D3), same scheme label as species
    tmp<reconstructionScheme<scalar>> tRecon
    (
        reconstructionScheme<scalar>::New
        (
            mesh,
            mesh.interpolationScheme("reconstruct(rho)")
        )
    );

    tmp<surfaceScalarField> wL, wR;
    forAll(evVib_, m)
    {
        tmp<surfaceScalarField> psiL, psiR;
        tRecon().limiter(evVib_[m], psiL, psiR);

        if (m == 0)
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
    tRecon().weights(evVib_[0], mesh.surfaceInterpolation::weights(), wL, wR);

    const dimensionedScalar zeroW(dimless, Zero);
    const dimensionedScalar oneW(dimless, 1.0);
    wL.ref() = max(min(wL(), oneW), zeroW);
    wR.ref() = max(min(wR(), oneW), zeroW);

    // (c) Advective-form pair and upwind interpolator (as speciesEqns.H)
    const volScalarField ddtRho(fvc::ddt(rho));
    const volScalarField divPhi(fvc::surfaceIntegrate(phi));
    upwind<scalar> upwindEps(mesh, phi);

    // 1/W = sum_s Y_s/W_s, for Gamma_eps = rhoD + (W/W_m)*GVib_m (D4)
    volScalarField invW("invW", Yv[0]/comp.W(0));
    for (label s = 1; s < Yv.size(); ++s)
    {
        invW += Yv[s]/comp.W(s);
    }
    
    // (c') Source inputs (Phase 6, D2-D4, D6). Evaluated once per outer iteration:
    //      tau at thermo p, T^k and X(Y^{k+1}); V-V at T_v,m^k.
    const PtrList<volScalarField>& tauVT = neThermo.correctVTRelaxationTime();

    const PtrList<volScalarField>* QVVp  = nullptr;   // V-V needs two molecules
    const PtrList<volScalarField>* nuVVp = nullptr;
    if (molI_.size() > 1)
    {
        QVVp  = &neThermo.correctVibVibSource(TVibSpecies_);
        nuVVp = &neThermo.vibVibFrequency();           // filled by the call above
    }

    const tmp<volScalarField> tCvTR = neThermo.CvTR();
    const scalar TvMinSrc = vibTvMin();

    const scalarField& rhoI = rho.primitiveField();
    const scalarField& pI   = neThermo.p().primitiveField();
    const scalarField& TI   = neThermo.TTR().primitiveField();
    const scalarField& cvI  = tCvTR().primitiveField();

    // Per-molecule maxima for the VibSrc log line [W/m3]
    scalarList QmaxVT(molI_.size(), Zero);
    scalarList QmaxVV(molI_.size(), Zero);
    scalarList QmaxCV(molI_.size(), Zero);

    scalarField sumQVV(rhoI.size(), Zero);   // PHASE6-TEMP T7.3

    // (d) One equation per molecule
    forAll(molI_, m)
    {
        const label s = molI_[m];
        volScalarField& eps = epsVib_[m];
        const volScalarField& ev = evVib_[m];

        // Face e_v (own limiter) times stored Y_m,f^k -> eps_f (D3)
        tmp<surfaceScalarField> evL, evR;
        reconstructionScheme<scalar>::reconstruct
        (
            ev,
            tmp<surfaceScalarField>(wL()),
            tmp<surfaceScalarField>(wR()),
            evL, evR
        );
        tmp<surfaceScalarField> evf = surfaceFieldSelect(evL, evR, phi, 0.0);
        tmp<surfaceScalarField> epsUp = upwindEps.interpolate(eps);

        fvScalarMatrix epsEqn
        (
            fvm::ddt(rho, eps) - fvm::Sp(ddtRho, eps)                 // real + pseudo time
          + fvm::div(phi, eps, "div(phi,epsVib)")                     // implicit upwind
          - fvm::Sp(divPhi, eps)                                      // advective-form pair
          + fvc::surfaceIntegrate(phi*(YfMol_[m]*evf() - epsUp()))    // deferred correction
        );

        if (!inviscid_)
        {
            // Implicit diffusivity on eps (stabiliser; cancels at convergence)
            volScalarField GammaEps("GammaEps", GVib_[m]/(comp.W(s)*invW));
            if (trhoD.valid())
            {
                GammaEps += trhoD();
            }

            // Physical face flux F_m = qVib_m - e_v,f * phiJ_s (D4, D5): same
            // qVib_ and phiJ_ as the rhoE residual, linear e_v as interpolate(h_s)
            surfaceScalarField Fm("Fm", qVib_[m]);
            if (phiJ_.size())
            {
                Fm -= fvc::interpolate(ev)*phiJ_[s];
            }

            epsEqn -= fvm::laplacian(GammaEps, eps, "laplacian(GammaEps,epsVib)");
            epsEqn += fvc::laplacian(GammaEps, eps, "laplacian(GammaEps,epsVib)");
            epsEqn -= fvc::surfaceIntegrate(Fm);
        }

        // (d') Vibrational sources, point-implicit. All three are
        //      collected into one pair: S = Q^k + a*(eps - eps^k), a <= 0.
        {
            volScalarField::Internal Qs
            (
                IOobject("Qvib", mesh.time().timeName(), mesh,
                         IOobject::NO_READ, IOobject::NO_WRITE, IOobject::NO_REGISTER),
                mesh,
                dimensionedScalar(dimEnergy/dimVolume/dimTime, Zero)
            );
            volScalarField::Internal as
            (
                IOobject("avib", mesh.time().timeName(), mesh,
                         IOobject::NO_READ, IOobject::NO_WRITE, IOobject::NO_REGISTER),
                mesh,
                dimensionedScalar(dimDensity/dimTime, Zero)
            );
            scalarField& QI = Qs.field();
            scalarField& aI = as.field();

            const scalar th = comp.thetaVib(s);
            const scalarField& YI   = Yv[s].primitiveField();       // Y^{k+1}
            const scalarField& epsI = eps.primitiveField();         // eps^k
            const scalarField& tauI = tauVT[s].primitiveField();

            // V-T: Q = rho (Y e*(T) - eps)/tau ; a = -(rho/tau)(1 + Y c*_vv/c_v,TR)
            forAll(QI, c)
            {
                const scalar eS  = comp.EsVib(s, pI[c], TI[c], TI[c], th);
                const scalar cvS = comp.CvVib(s, pI[c], TI[c], TI[c], th);
                const scalar r   = rhoI[c]/tauI[c];
                const scalar q   = r*(YI[c]*eS - epsI[c]);

                QI[c] = q;
                aI[c] = -r*(1.0 + YI[c]*cvS/cvI[c]);
                QmaxVT[m] = max(QmaxVT[m], mag(q));
            }

            // V-V: library Q explicit, -rho*nu implicit
            if (QVVp)
            {
                const scalarField& QvvI = (*QVVp)[s].primitiveField();
                const scalarField& nuI  = (*nuVVp)[s].primitiveField();
                forAll(QI, c)
                {
                    QI[c] += QvvI[c];
                    sumQVV[c] += QvvI[c];   // PHASE6-TEMP T7.3
                    aI[c] -= rhoI[c]*nuI[c];
                    QmaxVV[m] = max(QmaxVV[m], mag(QvvI[c]));
                }
            }

            // C-V, Patankar: production explicit, destruction implicit
            if (chemistry_)
            {
                const scalarField& RRI = chemistry_->RR(s).field();
                const scalarField& evI = ev.primitiveField();   // trace-guarded e_v,m
                const scalar aD = alphaCV_*DVib_[m];

                forAll(QI, c)
                {
                    const scalar w    = RRI[c];
                    const scalar eHat = cvPreferential_ ? aD : evI[c];
                    const scalar q    = w*eHat;

                    QI[c] += q;
                    if (w < 0)
                    {
                        const scalar epsMin = max
                        (
                            YI[c]*comp.EsVib(s, pI[c], TI[c], TvMinSrc, th),
                            VSMALL
                        );
                        aI[c] += q/max(epsI[c], epsMin);   // q < 0 -> a decreases
                    }
                    QmaxCV[m] = max(QmaxCV[m], mag(q));
                }
            }

            epsEqn -= Qs - as*eps.internalField();   // explicit part Q^k - a eps^k
            epsEqn += fvm::SuSp(-as, eps);           // -a eps on the LHS, implicit (a <= 0)
        }

        epsEqn.solve(mesh.solverDict("epsVib"));
    }

    // (e) Clip, decode, log (D6, D10)
    label nClip = 0;
    forAll(epsVib_, m)
    {
        scalarField& epsI = epsVib_[m].primitiveFieldRef();
        forAll(epsI, celli)
        {
            if (epsI[celli] < 0)
            {
                epsI[celli] = 0;
                ++nClip;
            }
        }
    }
    nClip = returnReduce(nClip, sumOp<label>());

    const label nFloor = vibDecode();

    Info<< "Vib:";
    forAll(molI_, m)
    {
        const scalarField& TvI = TVibSpecies_[molI_[m]].primitiveField();
        Info<< ' ' << comp.species()[molI_[m]]
            << " Tv [" << gMin(TvI) << ' ' << gMax(TvI) << ']';
    }
    Info<< ", clipped = " << nClip << ", Tv floor hits = " << nFloor << endl;

    Info<< "VibSrc [W/m3] max|QVT|,|QVV|,|QCV|:";
    forAll(molI_, m)
    {
        Info<< ' ' << comp.species()[molI_[m]] << ' '
            << returnReduce(QmaxVT[m], maxOp<scalar>()) << ','
            << returnReduce(QmaxVV[m], maxOp<scalar>()) << ','
            << returnReduce(QmaxCV[m], maxOp<scalar>());
    }
    if (QVVp)   // PHASE6-TEMP T7.3: sum_m Q_VV,m -> 0 at the plateau, not identically
    {
        Info<< "T7.3: max|sum_m QVV| = " << gMax(mag(sumQVV))
            << ", max_m max|QVV| = " << gMax(QmaxVV) << endl;
    }
    Info<< endl;
}