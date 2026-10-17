/*---------------------------------------------------------------------------*\
  IHIS - Implicit Hypersonics In-house Solver (OpenFOAM-v2312)
\*---------------------------------------------------------------------------*/

#include "neVibWallHeatFlux.H"
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
    defineTypeNameAndDebug(neVibWallHeatFlux, 0);
    addToRunTimeSelectionTable(functionObject, neVibWallHeatFlux, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::neVibWallHeatFlux::neVibWallHeatFlux
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

bool Foam::functionObjects::neVibWallHeatFlux::read(const dictionary& dict)
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


bool Foam::functionObjects::neVibWallHeatFlux::execute()
{
    return true;
}


bool Foam::functionObjects::neVibWallHeatFlux::write()
{
    // Per-molecule fluxes registered by the solver: neWall:qVib_<m>
    const word prefix("neWall:qVib_");
    wordList fluxNames;
    for
    (
        const word& objName
      : mesh_.sortedNames<surfaceScalarField>()
    )
    {
        if
        (
            objName.size() > prefix.size()
         && objName.compare(0, prefix.size(), prefix) == 0
        )
        {
            fluxNames.append(objName);
        }
    }

    if (fluxNames.empty())
    {
        WarningInFunction
            << "No neWall:qVib_<m> registered (1T or inviscid run)."
            << " Nothing written." << endl;
        return true;
    }

    const label nMol = fluxNames.size();

    wordList molNames(nMol);
    forAll(fluxNames, k)
    {
        molNames[k] = fluxNames[k].substr(prefix.size());
    }

    volScalarField qWall
    (
        IOobject
        (
            "neVibWallHeatFlux",
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
        << "    patch: integrals [W] q_vib";
    for (const word& m : molNames)
    {
        Info<< ", " << m;
    }
    Info<< nl;

    for (const label patchi : patchSet_.sortedToc())
    {
        const fvPatch& fp = mesh_.boundary()[patchi];
        const scalarField& magSf = mesh_.magSf().boundaryField()[patchi];

        // Column 0: total; columns 1..nMol: per molecule. Positive into the wall
        List<scalarField> Q(nMol + 1, scalarField(fp.size(), Zero));

        forAll(fluxNames, k)
        {
            const surfaceScalarField& qm =
                mesh_.lookupObject<surfaceScalarField>(fluxNames[k]);

            Q[k + 1] = -qm.boundaryField()[patchi]/magSf;
            Q[0] += Q[k + 1];
        }

        qWall.boundaryFieldRef()[patchi] = Q[0];

        Info<< "    " << fp.name() << ":";
        forAll(Q, k)
        {
            Info<< ' ' << gSum(Q[k]*magSf);
        }
        Info<< nl;

        // Gather to the master and write one CSV per patch
        List<vectorField> allCf(Pstream::nProcs());
        List<List<scalarField>> allQ(Pstream::nProcs());

        allCf[Pstream::myProcNo()] = fp.Cf();
        allQ[Pstream::myProcNo()] = Q;

        Pstream::gatherList(allCf);
        Pstream::gatherList(allQ);

        if (Pstream::master())
        {
            OFstream os(outDir/(fp.name() + ".csv"));
            os.precision(12);
            os  << "# " << type() << ", time " << time_.timeName()
                << ", positive into the wall [W/m2]" << nl
                << "x,y,z,q_vib";
            for (const word& m : molNames)
            {
                os << ",q_vib_" << m;
            }
            os << nl;

            forAll(allCf, proci)
            {
                const vectorField& Cf = allCf[proci];
                const List<scalarField>& Qp = allQ[proci];

                forAll(Cf, facei)
                {
                    os  << Cf[facei].x() << ',' << Cf[facei].y() << ','
                        << Cf[facei].z();
                    forAll(Qp, k)
                    {
                        os << ',' << Qp[k][facei];
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
