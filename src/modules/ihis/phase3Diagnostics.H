// IHIS Phase 3 diagnostics: Checkpoint A (containment) and T4.3.
// Temporary. Enable with   phase3Diagnostics true;   in fvSolution/flowSolver.
// Spliced after updateFields.H, so the thermo state (T, p, gamma, psi) is consistent.
if
(
    mesh.solutionDict().subDict("flowSolver")
        .getOrDefault<bool>("phase3Diagnostics", true)
)
{
    const scalarField& rhoI  = rho.primitiveField();
    const scalarField& rhoEI = rhoE.primitiveField();
    const vectorField& UI    = U.primitiveField();
    const scalarField& pI    = neThermo.p().primitiveField();

    const volScalarField gammaD(neThermo.gamma());
    const volScalarField psiD(neThermo.psi());
    const scalarField& gI   = gammaD.primitiveField();
    const scalarField& psiI = psiD.primitiveField();

    // 1. Energy offset used by T3.3-T3.8
    const scalarField eOffD(rhoEI/rhoI - 0.5*magSqr(UI) - pI/((gI - 1.0)*rhoI));

    // 2. c*: stock HiSA (H-based) vs datum-free (T3.1)
    const scalarField HD(rhoEI/rhoI + pI/rhoI);
    const scalarField cOldD(sqrt(2.0*(gI - 1.0)/(gI + 1.0)*HD));
    const scalarField cNewD
    (
        sqrt(2.0/(gI + 1.0)*(gI/psiI + 0.5*(gI - 1.0)*magSqr(UI)))
    );
    const scalar dcRel = gMax(mag(cNewD - cOldD)/cOldD);

    // 3. Wilke vs mass-weighted kappa (T4.3)
    const volScalarField kW(neThermo.kappaTRWilke());
    const volScalarField kM(neThermo.kappaTR());
    const scalar dkRel =
        gMax(mag(kW.primitiveField() - kM.primitiveField()))
       /max(gMax(kM.primitiveField()), VSMALL);

    Info<< "Phase3 diagnostics: Ns = "
        << neThermo.composition().Y().size() << nl
        << "    eOff [J/kg]           min/max = "
        << gMin(eOffD) << " / " << gMax(eOffD) << nl
        << "    c* rel. change (T3.1) max     = " << dcRel << nl
        << "    kappa Wilke vs mass-w max rel = " << dkRel << endl;
}
