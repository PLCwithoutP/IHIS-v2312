// Vibrational energy transport.
// Spliced into ihisModule::outerIteration() after speciesEqns.H, before
// updateFields.H. Unknown eps_m = Y_m*e_v,m [J/kg of mixture], one equation per
// molecule, same advective form as speciesEqns.H, Jacobi across molecules.

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

        fvScalarMatrix epsEqn
        (
            fvm::ddt(rho, eps) - fvm::Sp(ddtRho, eps)                 // real + pseudo time
          + fvm::div(phi, eps, "div(phi,epsVib)")                     // implicit upwind
          - fvm::Sp(divPhi, eps)                                      // advective-form pair
          + fvc::surfaceIntegrate(phi*(YfMol_[m]*evf() - epsUp()))    // deferred correction
          - fvm::laplacian(GammaEps, eps, "laplacian(GammaEps,epsVib)")   // implicit
          + fvc::laplacian(GammaEps, eps, "laplacian(GammaEps,epsVib)")   // same, at eps^k
          - fvc::surfaceIntegrate(Fm)                                 // physical diffusion
        );

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
}