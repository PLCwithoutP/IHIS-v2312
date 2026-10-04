// TEMPORARY (Phase 4, T1.3 / T6.1). Remove in T8.10.
// Start-up encode -> decode round trip, before any solve.
{
    const scalarField eDec
    (
        rhoE.primitiveField()/rho.primitiveField()
      - 0.5*magSqr(U.primitiveField())
      - eVibMix_().primitiveField()
      + eZero_().primitiveField()
    );
    const scalarField& eTh = e.primitiveField();

    Info<< setprecision(17)
        << "Probe: max|eDec - e|/|e| = "
        << gMax(mag(eDec - eTh)/max(mag(eTh), SMALL)) << nl
        << "Probe: sum/min/max(rhoE)  = " << gSum(rhoE.primitiveField()) << " "
        << gMin(rhoE.primitiveField()) << " " << gMax(rhoE.primitiveField()) << nl
        << "Probe: sum/min/max(eZero) = " << gSum(eZero_().primitiveField()) << " "
        << gMin(eZero_().primitiveField()) << " " << gMax(eZero_().primitiveField())
        << setprecision(IOstream::defaultPrecision()) << endl;

    // Fields for the cell-wise offset check (written to the start time directory)
    volScalarField rhoEProbe
    (
        IOobject("rhoEProbe", runTime.timeName(), mesh,
                 IOobject::NO_READ, IOobject::NO_WRITE, IOobject::NO_REGISTER),
        rhoE
    );
    volScalarField eZeroProbe
    (
        IOobject("eZeroProbe", runTime.timeName(), mesh,
                 IOobject::NO_READ, IOobject::NO_WRITE, IOobject::NO_REGISTER),
        eZero_()
    );
    rhoEProbe.write();
    eZeroProbe.write();
}