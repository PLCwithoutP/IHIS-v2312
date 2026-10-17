/*---------------------------------------------------------------------------*\
  IHIS - Implicit Hypersonics In-house Solver (OpenFOAM-v2312)
\*---------------------------------------------------------------------------*/

#include "neWallHeatFlux.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"
#include "volFields.H"
#include "wallPolyPatch.H"
#include "OFstream.H"
#include "OSspecific.H"
#include "ListListOps.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(neWallHeatFlux, 0);
    addToRunTimeSelectionTable(functionObject, neWallHeatFlux, dictionary);
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::functionObjects::neWallHeatFlux::patchFlux
(
    const word& fieldName,
    const label patchi
) const
{
    const surfaceScalarField* fPtr =
        mesh_.findObject<surfaceScalarField>(fieldName);

    if (fPtr)
    {
        return tmp<scalarField>::New(fPtr->boundaryField()[patchi]);
    }

    return tmp<scalarField>::New(mesh_.boundary()[patchi].size(), Zero);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::neWallHeatFlux::neWallHeatFlux
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

bool Foam::functionObjects::neWallHeatFlux::read(const dictionary& dict)
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


bool Foam::functionObjects::neWallHeatFlux::execute()
{
    return true;
}


bool Foam::functionObjects::neWallHeatFlux::write()
{
    if (!mesh_.foundObject<surfaceScalarField>("neWall:qTR"))
    {
        WarningInFunction
            << "neWall:qTR is not registered (inviscid run, or an IHIS build"
            << " without Phase 7 D10). Nothing written." << endl;
        return true;
    }

    volScalarField qWall
    (
        IOobject
        (
            "neWallHeatFlux",
            time_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        mesh_,
        dimensionedScalar(dimPower/dimArea, Zero)
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
        << "    patch: integrals [W] total, cond, slip, vib, diff" << nl;

    for (const label patchi : patchSet_.sortedToc())
    {
        const fvPatch& fp = mesh_.boundary()[patchi];
        const scalarField& magSf = mesh_.magSf().boundaryField()[patchi];

        // Positive into the wall [W/m2]
        const scalarField qCond(-patchFlux("neWall:qTR", patchi)/magSf);
        const scalarField qSlip(-patchFlux("neWall:sigmaDotU", patchi)/magSf);
        const scalarField qVib(-patchFlux("neWall:qVib", patchi)/magSf);
        const scalarField qDiff(patchFlux("neWall:hJ", patchi)/magSf);
        const scalarField qTot(qCond + qSlip + qVib + qDiff);

        qWall.boundaryFieldRef()[patchi] = qTot;

        Info<< "    " << fp.name() << ": "
            << gSum(qTot*magSf) << ' '
            << gSum(qCond*magSf) << ' '
            << gSum(qSlip*magSf) << ' '
            << gSum(qVib*magSf) << ' '
            << gSum(qDiff*magSf) << nl;

        // Gather to the master and write one CSV per patch
        List<vectorField> allCf(Pstream::nProcs());
        List<List<scalarField>> allQ(Pstream::nProcs());

        allCf[Pstream::myProcNo()] = fp.Cf();
        List<scalarField>& myQ = allQ[Pstream::myProcNo()];
        myQ.setSize(5);
        myQ[0] = qTot;
        myQ[1] = qCond;
        myQ[2] = qSlip;
        myQ[3] = qVib;
        myQ[4] = qDiff;

        Pstream::gatherList(allCf);
        Pstream::gatherList(allQ);

        if (Pstream::master())
        {
            OFstream os(outDir/(fp.name() + ".csv"));
            os.precision(12);
            os  << "# " << type() << ", time " << time_.timeName()
                << ", positive into the wall [W/m2]" << nl
                << "x,y,z,q_tot,q_cond,q_slip,q_vib,q_diff" << nl;

            forAll(allCf, proci)
            {
                const vectorField& Cf = allCf[proci];
                const List<scalarField>& Q = allQ[proci];

                forAll(Cf, facei)
                {
                    os  << Cf[facei].x() << ',' << Cf[facei].y() << ','
                        << Cf[facei].z();
                    forAll(Q, k)
                    {
                        os << ',' << Q[k][facei];
                    }
                    os << nl;
                }
            }
        }
    }

    qWall.write();
    Info<< endl;

    return true;
}


// ************************************************************************* //
