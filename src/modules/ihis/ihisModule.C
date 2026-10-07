/*---------------------------------------------------------------------------*\

    IHIS: In-house Hypersonics Implicit solver

    Copyright (C) 2014-2018 Johan Heyns - CSIR, South Africa
    Copyright (C) 2014-2018 Oliver Oxtoby - CSIR, South Africa
    Copyright (C) 2022 Oliver Oxtoby

-------------------------------------------------------------------------------
License
    This file is part of IHIS.

    IHIS is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    IHIS is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with IHIS.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "ihisModule.H"
#include "addToRunTimeSelectionTable.H"
#include "bound.H"
#include "preconditioner.H"
#include "ihisSolver.H"
#include "fvcSmooth.H"
#include "jacobian.H"
#include "discreteLaxFriedrichs.H"
#include "reconstructionScheme.H"
#include "cellFaceFunctions.H"
#include "characteristicWallPressureFvPatchScalarField.H"
#include "dynamicFvMesh.H"
#include "fvcLaplacian.H"
#include "fvcMeshPhi.H"
#include "surfaceInterpolate.H"   // fvc::interpolate
#include "fvcDiv.H"               // fvc::div
#include "fvcSnGrad.H"            // fvc::snGrad
#include "fvmDdt.H"               // fvm::ddt
#include "fvmDiv.H"
#include "fvmSup.H"
#include "fvcDdt.H"
#include "fvcSurfaceIntegrate.H"
#include "upwind.H"
#include "fvmLaplacian.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(ihisModule, 0);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


void ihisModule::createPreconditioners
(
    PtrList<preconditioner<2,1>>& preconditioners,
    PtrList<jacobian>& jacobians,
    const dictionary& parentDict
)
{
    if (parentDict.found("preconditioner"))
    {
        word preconditionerName(parentDict.lookup("preconditioner"));
        label nj = jacobians.size();
        label np = preconditioners.size();

        if (parentDict.found(preconditionerName))
        {
            const dictionary& dict = parentDict.subDict(preconditionerName);

            // Create preconditioner-specific Jacobian if specified
            if (dict.found("inviscidJacobian"))
            {
                jacobians.append
                (
                    new jacobian
                    (
                        parentDict.subDict(preconditionerName),
                        mesh(),
                        scalarVars_[0],
                        vectorVars_[0],
                        scalarVars_[1],
                        pThermo_(),
                        U_(),
                        ddtCoeff_(),
                        inviscid_
                    )
                );
                nj = jacobians.size();
            }

            // Recursively create further preconditioners if applicable
            createPreconditioners(preconditioners, jacobians, dict);
        }

        // Create this preconditioner; pass newly created Jacobian if created,
        // else previous one, and pointer to child preconditioner if one
        // was created otherwise NULL.
        preconditioners.append
        (
            preconditioner<2,1>::New
            (
                parentDict,
                jacobians[nj-1].matrix(),
                #if OPENFOAM >= 2212
                np > 0 ? preconditioners.get(np-1) : NULL
                #else
                np > 0 ? preconditioners(np-1) : NULL
                #endif
            )
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

ihisModule::ihisModule
(
    const word& name,
    const Time& t
)
:
    solverModule(name),
    time_(t)
{
}

ihisModule::ihisModule
(
    const word& name,
    const Time& t,
    const dictionary& dict
)
:
    solverModule(name),
    time_(t)
{
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void ihisModule::initialise()
{

    const Time& runTime = time_;

    #include "createDynamicFvMesh.H"

    // Detect steady-state analysis
    const dictionary& ddtControls = mesh.schemesDict().subDict("ddtSchemes");

    wordList ddtToc (ddtControls.toc());
    steadyState_ = false;
    forAll(ddtToc,s)
    {
        word ddtScheme(ddtToc[s]);
        word ddtSchemeLastWord;
        const tokenList& tokens = ddtControls.lookup(ddtToc[s]);
        if (tokens.last().isWord() && tokens.last().wordToken() == "steadyState")
        {
            if (ddtToc[s] == "default" || ddtToc[s] == "rhoU")
            {
                steadyState_ = true;
            }
        }
    }
    if (steadyState_)
    {
        Info << "Steady-state analysis detected" << nl << endl;
    }
    else
    {
        Info << "Transient analysis detected" << nl << endl;
    }

    residualIO defaultPseudoTol(2, 1, residualOrdering, 1e-5);
    residualIO defaultPseudoTolRel(2, 1, residualOrdering, 0.0);

    solnControl_.reset
    (
        new pseudotimeControl
        (
            mesh,
            steadyState_,
            2,
            1,
            defaultPseudoTol,
            defaultPseudoTolRel
        )
    );
    if (steadyState_)
    {
        solnControl_->setCorr(runTime.startTimeIndex());
    }

    localTimestepping_ =
        solnControl_->dict().lookupOrDefault<Switch>
        (
            "localTimestepping",
            true
        );
    if (localTimestepping_)
    {
        Info << "Local timestepping selected" << nl << endl;
        localTimesteppingBounding_ =
            solnControl_->dict().lookupOrDefault<Switch>
            (
                "localTimesteppingBounding",
                true
            );

        localTimesteppingLowerBound_ =
            solnControl_->dict().lookupOrDefault<scalar>
            (
                "localTimesteppingLowerBound",
                0.95
            );
        localTimesteppingLowerBound_ =
            (localTimesteppingLowerBound_ > 0 ? localTimesteppingLowerBound_ : 0.0);
        localTimesteppingLowerBound_ =
            (localTimesteppingLowerBound_ < 0.999 ? localTimesteppingLowerBound_ : 0.999);

        localTimesteppingUpperBound_ =
            solnControl_->dict().lookupOrDefault<scalar>
            (
                "localTimesteppingUpperBound",
                1.05
            );
        localTimesteppingUpperBound_ =
            (localTimesteppingUpperBound_ > 1.001 ? localTimesteppingUpperBound_ : 1.001);

        boundingDecrease_ =
            solnControl_->dict().lookupOrDefault<scalar>
            (
                "boundingDecrease",
                0.5
            );
        boundingDecrease_ = min(boundingDecrease_, scalar(1));
        boundingIncrease_ =
            solnControl_->dict().lookupOrDefault<scalar>
            (
                "boundingIncrease",
                2.0
            );
        boundingIncrease_ = max(boundingIncrease_, scalar(1));
    }
    else
    {
        Info << "Global timestepping selected" << nl << endl;
        localTimesteppingBounding_ = false;
    }

    #include "createFields.H"

    // IHIS: discreteLaxFriedrichs assumes a zero energy datum, which the
    // absolute-energy datum (E1a) violates for every Ns. LHS only: affects convergence, not the answer.
    // Walks the same dictionaries as the main Jacobian (flowSolver/<solver>)
    // and createPreconditioners (preconditioner chain).
    {
        const dictionary& flowDict = mesh.solutionDict().subDict("flowSolver");
        const word solverType(flowDict.get<word>("solver"));
        const dictionary solverDict(flowDict.subOrEmptyDict(solverType));

        const dictionary* dPtr = &solverDict;
        while (dPtr)
        {
            const dictionary& d = *dPtr;

            if
            (
                d.getOrDefault<word>("inviscidJacobian", word::null)
             == discreteLaxFriedrichs::typeName
            )
            {
                WarningInFunction
                    << "inviscidJacobian " << discreteLaxFriedrichs::typeName
                    << " in " << d.name()
                    << " assumes a zero energy datum; IHIS stores the absolute energy" << nl
                    << "    (formation energy included), so dp/drho can have the wrong sign." << nl
                    << "    Use laxFriedrichs or centralLaxFriedrichs." << endl;
            }

            dPtr = nullptr;
            if (d.found("preconditioner"))
            {
                const word pName(d.get<word>("preconditioner"));
                if (d.isDict(pName))
                {
                    dPtr = &d.subDict(pName);
                }
            }
        }
    }

    // rhoE is READ_IF_PRESENT but NO_WRITE, so any rhoE file in the start
    // time comes from another run and bypasses the eZero(Y) encoding.
    {
        IOobject rhoEIO
        (
            "rhoE",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        );

        if (rhoEIO.typeHeaderOk<volScalarField>(true))
        {
            WarningInFunction
                << "rhoE was read from " << rhoEIO.objectPath() << nl
                << "    IHIS never writes rhoE, so this file comes from another run"
                << " and bypasses the eZero(Y) encoding." << nl
                << "    Delete it unless you know it was encoded with this case's"
                << " species and thermo coefficients." << endl;
        }
    }
    
    // User config check
    #include "config.H"

    // Read or initialise pseudo Co number
    // NOTE: It is not necessarily read before every outer iteration
    // (see resetPseudo)
    pseudoCoNum_.reset
    (
        new uniformDimensionedScalarField
        (
            IOobject
            (
                "pseudoCoNum",
                runTime.timeName(),
                "uniform",
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            dimensionedScalar
            (
                "pseudoCoNumInit",
                dimless,
                solnControl_().dict().lookupOrDefault<scalar>("pseudoCoNum", 1)
            )
        )
    );

    if (localTimestepping_)
    {
        pseudoCoField_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    "pseudoCoField",
                    runTime.timeName(),
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                pseudoCoNum_()
            )
        );
        // Local damping is part of the pseudo-time state: write it and read it
        // back on restart. Rebuilding it as pseudoCoField/pseudoCoNum (stock
        // HiSA) mixes pseudoCoField from before the last damping/CFL update
        // with pseudoCoNum from after it, so a restart starts one update off.
        IOobject dampIO
        (
            "pseudoCoDamping",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        );

        if (dampIO.typeHeaderOk<volScalarField>(true))
        {
            pseudoCoDamping_.reset(new volScalarField(dampIO, mesh));
        }
        else
        {
            dampIO.readOpt(IOobject::NO_READ);
            pseudoCoDamping_.reset
            (
                new volScalarField(dampIO, pseudoCoField_()/pseudoCoNum_())
            );
        }

        scalar totCells = mesh.globalData().nTotalCells();
        Info<< "Initial pseudo Courant No: "
            << "Min: "
            << min(pseudoCoField_()).value()
            << " Mean: "
            << sum(pseudoCoField_()).value()/totCells
            << " Max: "
            << max(pseudoCoField_()).value() << nl << endl;
    }
    else
    {
        Info << "Initial pseudo Courant No: " << pseudoCoNum_->value() << nl << endl;
    }

    findDebugCell();
}

scalar ihisModule::timeStepScaling(const scalar& maxCoNum)
{

    const Time& runTime = time_;
    fvMesh& mesh = mesh_();
    ne2TReactionThermo& neThermo = pThermo_();
    volVectorField& U = U_();

    Info << "Mesh region: " << name() << endl;

    #include "compressibleCourantNo.H"

    scalar scale = maxCoNum/stabilise(CoNum, SMALL);

    return scale;
}

void ihisModule::beginTimeStep()
{
    // Clear out residuals from previous time step
    if (!steadyState_)
    {
        res_.clear();
        prevRes_.clear();
    }

    pseudotimeControl& solnControl = solnControl_();
    moveMeshOuterCorrectors_ =
        solnControl.dict().lookupOrDefault<label>
        (
            "moveMeshOuterCorrectors",
            0
        );

    // Reset pseudo Co number before every time step
    pseudoCoNumInit_ =
        solnControl.dict().lookupOrDefault<scalar>("pseudoCoNum", 1);
    if (solnControl.dict().lookupOrDefault<bool>("resetPseudo", false) )
    {
        pseudoCoNum_->value() = pseudoCoNumInit_;
    }
    pseudoCoNumMin_ =
        solnControl.dict().lookupOrDefault<scalar>("pseudoCoNumMin", 0.1);
    pseudoCoNumMax_ =
        solnControl.dict().lookupOrDefault<scalar>("pseudoCoNumMax", 25);
    pseudoCoNumMaxIncr_ = solnControl.dict().lookupOrDefault<scalar>
    (
        "pseudoCoNumMaxIncreaseFactor",
        1.2
    );
    pseudoCoNumMinDecr_ = solnControl.dict().lookupOrDefault<scalar>
    (
        "pseudoCoNumMinDecreaseFactor",
        0.1
    );
}

void ihisModule::outerIteration()
{

    #include "ptrRefs.H"

    if
    (
        solnControl.firstIter()
     || steadyState_
     || (moveMeshOuterCorrectors_ && !(solnControl.corr() % moveMeshOuterCorrectors_))
    )
    {
        if (isA<dynamicFvMesh>(mesh))
        {
            // Do any mesh changes

            scalar timeBeforeMeshUpdate = runTime.elapsedCpuTime();

            // Work around issues with mapping of oriented fields
            surfaceVectorField Un(phi*mesh.Sf()/sqr(mesh.magSf()));
            #if OPENFOAM >= 1712
            phi.setOriented(false);
            #endif


            dynamic_cast<dynamicFvMesh&>(mesh).update();


            if (mesh.changing())
            {
                Info<< "Execution time for mesh.update() = "
                    << runTime.elapsedCpuTime() - timeBeforeMeshUpdate
                    << " s" << endl;
            }

            #if OPENFOAM >= 1712
            phi.setOriented(true);
            #endif
            phi = (Un & mesh.Sf());
        }
    }

    // Store value before solve for bounding. Don't use storePrevIter here
    // as it introduces too many problems for mesh refinement/redistribution
    volScalarField rhoPrevIter("rhoPrevIter", rho);
    bounded_ = false;

    phiUp_.reset
    (
        new surfaceVectorField
        (
            IOobject
            (
                "phiUp",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimArea*rhoU.dimensions()*U.dimensions()
        )
    );
    phiEp_.reset
    (
        new surfaceScalarField
        (
            IOobject
            (
                "phiEp",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimArea*rhoE.dimensions()*U.dimensions()
        )
    );
    Up_.reset
    (
        new surfaceVectorField
        (
            IOobject
            (
                "Up",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimArea*U.dimensions()
        )
    );

    surfaceVectorField& phiUp = phiUp_();
    surfaceScalarField& phiEp = phiEp_();
    surfaceVectorField& Up = Up_();

    // Interpolation of primitive fields on faces
    flux.calcFlux(phi, phiUp, phiEp, Up);

    // Species diffusion coefficient, constant Lewis number:
    //   rho*D = Le*kappa_Wilke/Cp_TR   (same for all species -> sum of J_s = 0)
    // Built once per outer iteration from the current thermo state; used by the
    // energy residual and the species Laplacian (speciesEqns.H).
    tmp<volScalarField> trhoD;
    if (!inviscid_ && neThermo.composition().Y().size() > 1)
    {
        const scalar Le = neThermo.getOrDefault<scalar>("Le", 1.4);
        trhoD = Le*neThermo.kappaTRWilke()/neThermo.CpTR();
        trhoD.ref().rename("rhoD");
    }

    // Species enthalpies on IHIS's energy datum (same as updateEZero()):
    //   h_s = Hs(s,p,T) - hOff_s,   hOff_s = -Hc(s),   i.e. h_s = Ha(s,p,T) (absolute)
    // so that sum_s Y_s h_s = (e - eZero + eVibMix) + p/rho. Built only when diffusion is active.
    PtrList<volScalarField> hs;
    if (trhoD.valid())
    {
        const compositionInterface& comp = neThermo.composition();
        const volScalarField& pTh = neThermo.p();
        const volScalarField& TTh = neThermo.TTR();

        hs.setSize(comp.Y().size());

        forAll(hs, s)
        {
            const scalar hOff = -comp.Hc(s);

            hs.set
            (
                s,
                new volScalarField
                (
                    IOobject
                    (
                        "h_" + comp.species()[s],
                        mesh.time().timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE,
                        IOobject::NO_REGISTER
                    ),
                    mesh,
                    dimensionedScalar(dimEnergy/dimMass, Zero)
                )
            );

            // Cells
            const scalarField& pI = pTh.primitiveField();
            const scalarField& TI = TTh.primitiveField();
            scalarField& hI = hs[s].primitiveFieldRef();
            forAll(hI, celli)
            {
                hI[celli] = comp.Hs(s, pI[celli], TI[celli]) - hOff;
            }

            // Patches (all, including processor: h_s depends only on p and T)
            volScalarField::Boundary& hBf = hs[s].boundaryFieldRef();
            forAll(hBf, patchi)
            {
                const scalarField& pP = pTh.boundaryField()[patchi];
                const scalarField& TP = TTh.boundaryField()[patchi];
                scalarField& hP = hBf[patchi];
                forAll(hP, facei)
                {
                    hP[facei] = comp.Hs(s, pP[facei], TP[facei]) - hOff;
                }
            }
        }
        // Each molecule's enthalpy carries its own vibrational
        // energy, h_s += e_v,s(T_v,s), cells and all patches. evVib_ holds
        // neighbour values on processor patches, as p and T do here. Then
        //   sum_s Y_s h_s = (e - eZero + eVibMix) + p/rho,
        // and -div(sum_s h_s J_s) in rhoE contains -div(sum_m e_v,m J_m),
        // the same term each vibrational equation subtracts (D5).
        if (twoT_)
        {
            forAll(molI_, m)
            {
                hs[molI_[m]] += evVib_[m];
            }
        }
    }

    // Update residuals
    #include "residualsUpdate.H"

    // Pseudo Courant number relaxation
    setPseudoCoNum();

    // Update pseudo time step
    setPseudoDeltaT();

    // Set based on latest rPseudoDeltaT
    ddtCoeff = max(max(fvm::ddt(rho)->diag(),fvm::ddt(rhoU)->diag()),fvm::ddt(rhoE)->diag())/mesh.V();

    #include "cellDebug.H"

    // Store previous iteration values. Don't use storePrevIter here as it
    // introduces too many problems for mesh refinement/redistribution
    scalarVarsPrevIter_.clear();
    vectorVarsPrevIter_.clear();
    forAll(scalarVars_, i)
    {
        scalarVarsPrevIter_.append(new volScalarField(scalarVars_[i].name() + "PrevIter", scalarVars_[i]));
    }
    forAll(vectorVars_, i)
    {
        vectorVarsPrevIter_.append(new volVectorField(vectorVars_[i].name() + "PrevIter", vectorVars_[i]));
    }

    // Solver
    {
        residualIO defaultSolverTol(2, 1, residualOrdering, 1e-12);

        // Storage of jacobians and preconditioners
        PtrList<jacobian> jacobians;
        PtrList<preconditioner<2,1>> preconditioners;


        const dictionary& dict = mesh.solutionDict().subDict("flowSolver");

        const word solverType(dict.lookup("solver"));

        // Create main Jacobian
        jacobians.append
        (
            new jacobian
            (
                dict.subOrEmptyDict(solverType),
                mesh,
                rho,
                rhoU,
                rhoE,
                neThermo,
                U,
                ddtCoeff_(),
                inviscid_
            )
        );

        // Recursively create needed preconditioners and their jacobians (if
        // applicable)
        createPreconditioners
        (
            preconditioners,
            jacobians,
            dict.subOrEmptyDict(solverType)
        );

        // Create solver
        autoPtr<ihisSolver<2,1>> sol =
            ihisSolver<2,1>::New
            (
                dict,
                jacobians[0].matrix(),
                #if OPENFOAM >= 2212
                preconditioners.size() ? preconditioners.get(0) : NULL,
                #else
                preconditioners.size() ? preconditioners(0) : NULL,
                #endif
                defaultSolverTol
            );

        // Call solver
        res_.clear();
        sol->solve
        (
            scalarVars_,
            vectorVars_,
            scalarResiduals_,
            vectorResiduals_,
            res_
        );
        solnControl_->setResidual(res_());
    }

    if (localTimesteppingBounding_)
    {
        volScalarField rhoMin
        (
            localTimesteppingLowerBound_*scalarVarsPrevIter_[0]
        );
        volScalarField rhoMax(localTimesteppingUpperBound_*scalarVarsPrevIter_[0]);
        // IHIS: datum-free bounding. Relative bounds and the
        // positivity check act on the thermal (pressure-producing) energy
        // eTh = p/((gamma-1)*rho). eOff = e - eTh (formation energy; later e_v)
        // is frozen across the solve: Y is unchanged until speciesEqns.H, and
        // neThermo still holds the previous-iteration p and gamma here.
        const volScalarField eThPrev
        (
            neThermo.p()/((neThermo.gamma() - 1.0)*scalarVarsPrevIter_[0])
        );
        const volScalarField eOff
        (
            scalarVarsPrevIter_[1]/scalarVarsPrevIter_[0]
          - 0.5*magSqr(vectorVarsPrevIter_[0]/scalarVarsPrevIter_[0])
          - eThPrev
        );
        volScalarField eMin(localTimesteppingLowerBound_*eThPrev);
        volScalarField eMax(localTimesteppingUpperBound_*eThPrev);
        volScalarField eTemp(rhoE/rho - 0.5*magSqr(rhoU/rho) - eOff);

        volScalarField factor
        (
            IOobject
            (
                "factor",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("factor", dimless, 1)
        );
        factor.primitiveFieldRef() = GREAT;

        forAll(mesh.owner(), facei)
        {
            label own = mesh.owner()[facei];
            label nei = mesh.neighbour()[facei];

            if
            (
                (rho[own] < rhoMin[own]) || (rho[own] > rhoMax[own])
             || (eTemp[own] < eMin[own]) || (eTemp[own] > eMax[own])
             || (eTemp[own] < SMALL)
            )
            {
                factor[own] = min(boundingDecrease_, factor[own]);
                factor[nei] = min((1+boundingDecrease_)/2, factor[nei]);
            }
            else
            {
                factor[own] = min(boundingIncrease_, factor[own]);
            }
            if
            (
                (rho[nei] < rhoMin[nei]) || (rho[nei] > rhoMax[nei])
             || (eTemp[nei] < eMin[nei]) || (eTemp[nei] > eMax[nei])
             || (eTemp[nei] < SMALL)
            )
            {
                factor[nei] = min(boundingDecrease_, factor[nei]);
                factor[own] = min((1+boundingDecrease_)/2, factor[own]);
            }
            else
            {
                factor[nei] = min(boundingIncrease_, factor[nei]);
            }
        }

        forAll(mesh.boundary(), patchi)
        {
            if (mesh.boundary()[patchi].coupled())
            {
                scalarField rhoNei
                (
                    rho.boundaryField()[patchi].patchNeighbourField()
                );
                scalarField rhoMinNei
                (
                    rhoMin.boundaryField()[patchi].patchNeighbourField()
                );
                scalarField rhoMaxNei
                (
                    rhoMax.boundaryField()[patchi].patchNeighbourField()
                );
                scalarField eTempNei
                (
                    eTemp.boundaryField()[patchi].patchNeighbourField()
                );
                scalarField eMinNei
                (
                    eMin.boundaryField()[patchi].patchNeighbourField()
                );
                scalarField eMaxNei
                (
                    eMax.boundaryField()[patchi].patchNeighbourField()
                );
                const labelUList& fc = mesh.boundary()[patchi].faceCells();
                forAll(fc, bfacei)
                {
                    if
                    (
                        (rho[fc[bfacei]] < rhoMin[fc[bfacei]])
                     || (rho[fc[bfacei]] > rhoMax[fc[bfacei]])
                     || (eTemp[fc[bfacei]] < eMin[fc[bfacei]])
                     || (eTemp[fc[bfacei]] > eMax[fc[bfacei]])
                     || (eTemp[fc[bfacei]] < SMALL)
                    )
                    {
                        factor[fc[bfacei]] =
                            min(boundingDecrease_, factor[fc[bfacei]]);
                    }
                    else
                    {
                        factor[fc[bfacei]] =
                            min(boundingIncrease_, factor[fc[bfacei]]);
                    }
                    if
                    (
                        (rhoNei[bfacei] < rhoMinNei[bfacei])
                     || (rhoNei[bfacei] > rhoMaxNei[bfacei])
                     || (eTempNei[bfacei] < eMinNei[bfacei])
                     || (eTempNei[bfacei] > eMaxNei[bfacei])
                     || (eTempNei[bfacei] < SMALL)
                    )
                    {
                        factor[fc[bfacei]] =
                            min((1+boundingDecrease_)/2, factor[fc[bfacei]]);
                    }
                }
            }
        }

        // Also clamp the variables themselves
        rho = max(rho, rhoMin);
        rho = min(rho, rhoMax);
        eTemp = max(eTemp, eMin);
        eTemp = min(eTemp, eMax);
        rhoE = rho*(eTemp + eOff + 0.5*magSqr(rhoU/rho));

        pseudoCoDamping_() *= factor;
        pseudoCoDamping_() = min(scalar(1), pseudoCoDamping_());
        pseudoCoDamping_() =
            max(pseudoCoNumMin_/pseudoCoNum_(), pseudoCoDamping_());
    }

    phiUp_.clear();
    phiEp_.clear();
    Up_.clear();

    #include "speciesEqns.H"

    #include "vibEqns.H" 
    
    // Update fields
    #include "updateFields.H"
    Info<< "min(rho) = " << gMin(rho.primitiveField())
        << ", min(TTR) = " << gMin(neThermo.TTR().primitiveField()) << endl;
    #include "phase6Probe.H"   // PHASE6-TEMP T7.1
    
}

void ihisModule::updateEZero()
{
    const ne2TReactionThermo& neThermo = pThermo_();

    eZero_() = -neThermo.hc();

    // eZero must be finite everywhere (moved here from createFields.H)
    {
        bool finite = true;

        for (const scalar v : eZero_().primitiveField())
        {
            if (!std::isfinite(v)) { finite = false; break; }
        }
        forAll(eZero_().boundaryField(), patchi)
        {
            for (const scalar v : eZero_().boundaryField()[patchi])
            {
                if (!std::isfinite(v)) { finite = false; break; }
            }
        }

        if (!returnReduce(finite, andOp<bool>()))
        {
            FatalErrorInFunction
                << "eZero = -hc() is not finite. Check the formation-energy"
                << " data (b1 coefficients / Hc) in the species dictionary."
                << exit(FatalError);
        }
    }
}

void ihisModule::findDebugCell()
{
    // Code adapted from findRefCell.C
    const dictionary& dict = mesh_->time().controlDict();
    if (dict.found("debugCell"))
    {
        cellDebugging_ = true;

        if (!Pstream::parRun() || Pstream::myProcNo() == readLabel(dict.lookup("debugProcessor")))
        {
            debugCell_ = readLabel(dict.lookup("debugCell"));
            if (debugCell_ < 0 || debugCell_ >= mesh_->nCells())
            {
                debugCell_ = -1;
                Pout << "Debug cell " << debugCell_ << " is out of range." << endl;
            }
        }
        else
        {
            debugCell_ = -1;
        }

        label hasRef = (debugCell_ >= 0 ? 1 : 0);
        if (returnReduce<label>(hasRef, sumOp<label>()) != 1)
        {
            WarningInFunction
                << "debugCell not found." << endl;
            cellDebugging_ = false;
        }
    }
    else if (dict.found("debugPoint"))
    {
        point debugPoint(dict.lookup("debugPoint"));

        // Try fast approximate search avoiding octree construction
        debugCell_ = mesh_->findCell(debugPoint, polyMesh::FACE_PLANES);

        label hasRef = (debugCell_ >= 0 ? 1 : 0);
        label sumHasRef = returnReduce<label>(hasRef, sumOp<label>());

        // If reference cell no found use octree search
        // with cell tet-decompositoin
        if (sumHasRef != 1)
        {
            debugCell_ = mesh_->findCell(debugPoint);

            hasRef = (debugCell_ >= 0 ? 1 : 0);
            sumHasRef = returnReduce<label>(hasRef, sumOp<label>());
        }

        if (sumHasRef != 1)
        {
            WarningInFunction
                << "Unable to set debug cell at point " << debugPoint << ":"
                << nl << "Found on " << sumHasRef << " domains (should be one)"
                << endl;
            cellDebugging_ = false;
        }
        else
        {
            cellDebugging_ = true;
        }
    }
    else
    {
        cellDebugging_ = false;
    }
    if (cellDebugging_)
    {
        if (debugCell_ >= 0)
        {
            Pout << "Found debug cell " << debugCell_ << nl << endl;
        }
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    addToRunTimeSelectionTable
    (
        solverModule,
        ihisModule,
        dictionary
    );
}

// ************************************************************************* //

