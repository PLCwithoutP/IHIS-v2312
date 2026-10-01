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
        
        Info<< "Species: min(Y) = " << minYRaw
            << ", max|sum(Y) - 1| = " << sumYErr << endl;
    }
}