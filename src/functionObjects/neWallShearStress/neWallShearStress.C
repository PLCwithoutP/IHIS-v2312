/*---------------------------------------------------------------------------*\
  IHIS - Implicit Hypersonics In-house Solver (OpenFOAM-v2312)
\*---------------------------------------------------------------------------*/

#include "neWallShearStress.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"
#include "volFields.H"
#include "wallPolyPatch.H"
#include "OFstream.H"
#include "OSspecific.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(neWallShearStress, 0);
    addToRunTimeSelectionTable(functionObject, neWallShearStress, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::neWallShearStress::neWallShearStress
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    patchSet_()
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::neWallShearStress::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

    wordRes patchNames;
    if (dict.readIfPresent("patches", patchNames))
    {
        patchSet_ = pbm.patchSet(patchNames);
    }
    else
    {
        patchSet_.clear();
        forAll(pbm, patchi)
        {
            if (isA<wallPolyPatch>(pbm[patchi]))
            {
                patchSet_.insert(patchi);
            }
        }
    }

    Info<< type() << " " << name() << ": patches";
    for (const label patchi : patchSet_.sortedToc())
    {
        Info<< ' ' << pbm[patchi].name();
    }
    Info<< nl << endl;

    return true;
}


bool Foam::functionObjects::neWallShearStress::execute()
{
    return true;
}


bool Foam::functionObjects::neWallShearStress::write()
{
    const surfaceVectorField* tauPtr =
        mesh_.findObject<surfaceVectorField>("neWall:tau");

    if (!tauPtr)
    {
        WarningInFunction
            << "neWall:tau is not registered (inviscid run, or an IHIS build"
            << " without Phase 7 D10). Nothing written." << endl;
        return true;
    }

    volVectorField tauWall
    (
        IOobject
        (
            "neWallShearStress",
            time_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        mesh_,
        dimensionedVector(dimPressure, Zero)
    );

    const fileName outDir
    (
        time_.globalPath()/functionObject::outputPrefix/name()
       /time_.timeName()
    );
    if (Pstream::master())
    {
        mkDir(outDir);
    }

    Info<< type() << " " << name() << " write:" << nl
        << "    patch: viscous force on the wall [N]" << nl;

    for (const label patchi : patchSet_.sortedToc())
    {
        const fvPatch& fp = mesh_.boundary()[patchi];
        const scalarField& magSf = mesh_.magSf().boundaryField()[patchi];
        const vectorField n(fp.nf());

        // Gas on the wall [Pa]
        const vectorField tauW(-tauPtr->boundaryField()[patchi]);
        const vectorField tauT(tauW - n*(n & tauW));
        const scalarField magTauT(mag(tauT));

        tauWall.boundaryFieldRef()[patchi] = tauW;

        Info<< "    " << fp.name() << ": " << gSum(tauW*magSf) << nl;

        // Gather to the master and write one CSV per patch
        List<vectorField> allCf(Pstream::nProcs());
        List<vectorField> allTau(Pstream::nProcs());
        List<scalarField> allTauT(Pstream::nProcs());

        allCf[Pstream::myProcNo()] = fp.Cf();
        allTau[Pstream::myProcNo()] = tauW;
        allTauT[Pstream::myProcNo()] = magTauT;

        Pstream::gatherList(allCf);
        Pstream::gatherList(allTau);
        Pstream::gatherList(allTauT);

        if (Pstream::master())
        {
            OFstream os(outDir/(fp.name() + ".csv"));
            os.precision(12);
            os  << "# " << type() << ", time " << time_.timeName()
                << ", stress of the gas on the wall [Pa]" << nl
                << "x,y,z,tau_x,tau_y,tau_z,tau_t" << nl;

            forAll(allCf, proci)
            {
                const vectorField& Cf = allCf[proci];
                const vectorField& T = allTau[proci];
                const scalarField& Tt = allTauT[proci];

                forAll(Cf, facei)
                {
                    os  << Cf[facei].x() << ',' << Cf[facei].y() << ','
                        << Cf[facei].z() << ','
                        << T[facei].x() << ',' << T[facei].y() << ','
                        << T[facei].z() << ',' << Tt[facei] << nl;
                }
            }
        }
    }

    tauWall.write();
    Info<< endl;

    return true;
}


// ************************************************************************* //
