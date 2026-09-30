// IHIS : species transport step (frozen, decoupled).
// Spliced into ihisModule::outerIteration() after the bounding block and the
// phiUp_/phiEp_/Up_ clear() calls, immediately before updateFields.H, so that
// updateEZero() and neThermo.correct() see the new Y.
// Uses phi from this outer iteration's calcFlux and rho from the 5x5 solve.
// Ns <= 1: skipped entirely (Phase 2 path unchanged).
{
    PtrList<volScalarField>& Y = neThermo.composition().Y();

    if (Y.size() > 1)
    {
        // IHIS : common-limiter weights, built once from all species.
        // HiSA convention: limiter psi in, weights out;
        //   w_L = 1 + psi_L*(w_CD - 1),  w_R = psi_R*w_CD   (both linear in psi),
        // so taking min over species of psi, then converting, gives the most
        // restrictive limiter shared by every species. Common weights give
        // sum_s Y_f,s = 1 on every face.
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

        // Clamp to [0,1]: every face value is then a convex combination of its
        // two cells (Y_f >= 0) whatever limiter is selected.
        const dimensionedScalar zeroW(dimless, Zero);
        const dimensionedScalar oneW(dimless, 1.0);
        wL.ref() = max(min(wL(), oneW), zeroW);
        wR.ref() = max(min(wR(), oneW), zeroW);

        // T5.2 check (temporary): sum of face mass fractions
        surfaceScalarField sumYf
        (
            IOobject("sumYf", runTime.timeName(), mesh),
            mesh,
            dimensionedScalar(dimless, Zero)
        );

        // T5.3: species-independent pieces, computed once
        const volScalarField ddtRho(fvc::ddt(rho));              // real-time d(rho)/dt (no pseudo term)
        const volScalarField divPhi(fvc::surfaceIntegrate(phi)); // div(phi) per unit volume
        upwind<scalar> upwindYf(mesh, phi);                      // upwind interpolation on sign(phi)

        forAll(Y, i)
        {
            volScalarField& Yi = Y[i];

            // Per species: same weights; pass non-owning tmps, because
            // reconstruct() may clear the tmps it receives (tmp::clear is const)
            tmp<surfaceScalarField> YL, YR;
            reconstructionScheme<scalar>::reconstruct
            (
                Yi,
                tmp<surfaceScalarField>(wL()),
                tmp<surfaceScalarField>(wR()),
                YL, YR
            );

            // Upwind selection on sign(phi): phi >= 0 takes the owner-side state
            tmp<surfaceScalarField> Yf = surfaceFieldSelect(YL, YR, phi, 0.0);

            sumYf += Yf();

            // T5.3: species equation, advective form
            //   rho dY/dt + div(phi Y_f) - Y div(phi) = 0
            // LHS: implicit upwind (M-matrix); RHS: second-order deferred correction
            tmp<surfaceScalarField> YfUpwind = upwindYf.interpolate(Yi);

            fvScalarMatrix YiEqn
            (
                fvm::ddt(rho, Yi) - fvm::Sp(ddtRho, Yi)          // real time (advective) + pseudo (dualTime)
              + fvm::div(phi, Yi, "div(phi,Yi)")                  // implicit upwind
              - fvm::Sp(divPhi, Yi)                               // advective-form convection pair
              + fvc::surfaceIntegrate(phi*(Yf() - YfUpwind()))    // deferred correction
            );
            // 3b (T6): if (!inviscid_) { YiEqn -= fvm::laplacian(rhoD, Yi); }
            // T5.4: solve; fvMatrix::solve prints "Solving for <species>, Initial residual ..."
            // Species are not part of the pseudo-time stopping criterion (Q8).
            YiEqn.solve(mesh.solverDict("Yi"));
        }

        Info<< "T5.2 check: max|sum(Yf) - 1| (internal faces) = "
            << gMax(mag(sumYf.primitiveField() - 1.0)()) << endl;
        // T5.5: clip, measure drift, renormalise, refresh patches.
        // Must happen before updateEZero()/neThermo.correct() in updateFields.H,
        // which read Y on cells and patches.
        volScalarField Ysum("Ysum", 0.0*Y[0]);

        scalar minYRaw = GREAT;
        forAll(Y, i)
        {
            minYRaw = min(minYRaw, gMin(Y[i].primitiveField()));
            Y[i] = max(Y[i], dimensionedScalar(dimless, Zero));
            Ysum += Y[i];
        }

        // Drift after clipping, before renormalising (printed in T5.6)
        const scalar sumYErr = gMax(mag(Ysum.primitiveField() - 1.0)());

        const dimensionedScalar sumFloor(dimless, SMALL);
        forAll(Y, i)
        {
            Y[i] /= max(Ysum, sumFloor);
            Y[i].correctBoundaryConditions();
        }
        // T5.6: species log line (min Y before clipping; sum(Y) drift before renormalising)
        Info<< "Species: min(Y) = " << minYRaw
            << ", max|sum(Y) - 1| = " << sumYErr << endl;
    }
}