// PHASE6-TEMP T7.1 (remove in T8.13): energy-partition identity after the decode.
// E = rhoE/rho must equal sum_s Y_s (EsT_s(p,T) + Hc_s) + sum_m eps_m + |U|^2/2
// to round-off. 2T only, so the 1T Checkpoint A logs stay bit-identical.
if (twoT_)
{
    const compositionInterface& compP = neThermo.composition();
    const PtrList<volScalarField>& YP = compP.Y();

    const scalarField& pP    = neThermo.p().primitiveField();
    const scalarField& TP    = neThermo.TTR().primitiveField();
    const scalarField& rhoP  = rho.primitiveField();
    const scalarField& rhoEP = rhoE.primitiveField();
    const vectorField& UP    = U.primitiveField();
    const scalarField& eVP   = eVibMix_().primitiveField();   // = sum_m eps_m (vibDecode)

    scalar maxRel = 0;
    forAll(rhoP, c)
    {
        scalar eTRc = 0;
        forAll(YP, s)
        {
            eTRc += YP[s][c]*(compP.EsT(s, pP[c], TP[c]) + compP.Hc(s));
        }

        const scalar E    = rhoEP[c]/rhoP[c];
        const scalar Erec = eTRc + eVP[c] + 0.5*magSqr(UP[c]);

        maxRel = max(maxRel, mag(E - Erec)/max(mag(E), SMALL));
    }

    Info<< "T7.1 energy partition: max|E - E(T,Y,eps,U)|/|E| = "
        << returnReduce(maxRel, maxOp<scalar>()) << endl;
}