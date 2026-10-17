/*---------------------------------------------------------------------------*\

    IHIS: In-house Hypersonics Implicit solver

    Copyright (C) 2014-2018 Johan Heyns - CSIR, South Africa
    Copyright (C) 2014-2018 Oliver Oxtoby - CSIR, South Africa
    Copyright (C) 2020-2021 Oliver Oxtoby
    Copyright (C) 2011-2013 OpenFOAM Foundation

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

Info<< "Reading thermophysical properties\n" << endl;
pThermo_ =
(
    ne2TReactionThermo::New(mesh)
);
Info << endl;
ne2TReactionThermo& neThermo = pThermo_();

volScalarField& e = neThermo.eTR();

tmp<volScalarField> tmu = neThermo.mu();
const volScalarField& mu = tmu();

inviscid_ = neThermo.getOrDefault<bool>("inviscid", false);
Info<< (inviscid_ ? "Inviscid" : "Viscous")
    << " analysis detected" << nl << endl;

Info<< "Reading field U" << nl << endl;
U_.reset
(
    new volVectorField
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    )
);
volVectorField& U = U_();

scalarVars_.setSize(2);
vectorVars_.setSize(1);

scalarVars_.set
(
    0,
    new volScalarField
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        neThermo.rho(),
        calculatedFvPatchScalarField::typeName
    )
);
volScalarField& rho = scalarVars_[0];

vectorVars_.set
(
    0,
    new volVectorField
    (
        IOobject
        (
            "rhoU",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        rho*U
    )
);
volVectorField& rhoU = vectorVars_[0];

eZero_.reset
(
    new volScalarField
    (
        IOobject("eZero", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
        mesh,
        dimensionedScalar("eZero", e.dimensions(), 0.0)
    )
);
updateEZero();


// Finite-rate chemistry. Master switch is 'chemistry' in
// thermophysicalProperties (the same key that makes the mixture read
// reactionsList); built only for mixtures. Never use Qdot() in rhoE:
// chemistry is energy-neutral by construction (E1a).
if
(
    neThermo.getOrDefault<bool>("chemistry", false)
 && neThermo.composition().Y().size() > 1
)
{
    chemistry_ = chemistryInterface::New(neThermo);

    // The model registered constant/chemistryProperties during construction;
    // its own switch must agree, otherwise calculateWithJacobianDiag()
    // would refuse to run (and calculate() would silently return).
    const IOdictionary& chemProps =
        mesh.lookupObject<IOdictionary>("chemistryProperties");

    if (!chemProps.get<Switch>("chemistry"))
    {
        FatalErrorInFunction
            << "'chemistry' is true in thermophysicalProperties but false in "
            << chemProps.objectPath() << nl
            << "    Both switches must agree for a reacting run."
            << exit(FatalError);
    }

    Info<< "IHIS: finite-rate chemistry active ("
        << neThermo.composition().Y().size() << " species)" << nl << endl;
}

// Two-temperature (Phase 5): one transported vibrational energy per molecule,
// V–T, V–V and C–V sources are added point-implicitly in vibEqns.H . IHIS owns the vibrational
// state; the thermo only receives the mixture TVib (D1). All cases use a
// multi-component mixture, so Y() is never empty (the thermo enforces this in 2T).
twoT_ = neThermo.twoTemperature();

if (twoT_)
{
    // A 2T case with chemistry on needs the chemistry model; chemistry_ is
    // built only for Ns > 1 (Phase 4 D7), so a pure-gas case would silently
    // run frozen.
    if (neThermo.getOrDefault<bool>("chemistry", false) && !chemistry_)
    {
        FatalErrorInFunction
            << "chemistry is true but the mixture has a single species:"
            << " there is nothing to react." << nl
            << "    Set chemistry false or add the product species."
            << exit(FatalError);
    }
    const compositionInterface& comp = neThermo.composition();
    const PtrList<volScalarField>& Ys = comp.Y();

    // Molecules = species with a vibrational mode; molI_[m] = species index
    DynamicList<label> mols(Ys.size());
    forAll(Ys, s)
    {
        if (comp.isSpecieMolecular(s))
        {
            mols.append(s);
        }
    }
    molI_.transfer(mols);

    if (molI_.empty())
    {
        FatalErrorInFunction
            << "twoTemperature is true but no species is molecular:"
            << " there is no vibrational mode to transport."
            << exit(FatalError);
    }

    vvOn_ = neThermo.getOrDefault<bool>("VVrelaxation", true);
    Info<< "IHIS: V-V relaxation " << ((vvOn_ && molI_.size() > 1) ? "on" : "off") << nl;

    Info<< "IHIS: two-temperature active, " << molI_.size()
        << " vibrational equation(s):";
    forAll(molI_, m)
    {
        Info<< ' ' << comp.species()[molI_[m]];
    }

    Info << nl;

    // C-V settings. Read IHIS-side: the library correctCVSource
    // uses the mixture TVib and is never called by IHIS.
    if (chemistry_)
    {
        const word cvTreatment(neThermo.get<word>("CVTreatment"));

        if (cvTreatment == "preferential")
        {
            cvPreferential_ = true;
        }
        else if (cvTreatment != "non-preferential")
        {
            FatalIOErrorInFunction(neThermo)
                << "CVTreatment '" << cvTreatment << "' unknown." << nl
                << "    Valid: preferential, non-preferential"
                << exit(FatalIOError);
        }

        alphaCV_ = neThermo.getOrDefault<scalar>("alphaDissociation", 0.3);

        DVib_.setSize(molI_.size(), Zero);

        if (cvPreferential_)
        {
            const IOdictionary speciesDict
            (
                IOobject
                (
                    "speciesDict",
                    mesh.time().constant(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE,
                    IOobject::NO_REGISTER
                )
            );

            forAll(molI_, m)
            {
                const word& name = comp.species()[molI_[m]];
                DVib_[m] =
                    speciesDict.subDict(name).get<scalar>("dissociationPotential");

                // Units guard: J/kg expected (N2 ~3.36e7). Catches eV or J/mol data.
                if (DVib_[m] < 1e6 || DVib_[m] > 1e8)
                {
                    FatalErrorInFunction
                        << "dissociationPotential of " << name << " = "
                        << DVib_[m] << " is outside 1e6-1e8 J/kg." << nl
                        << "    Give D/W in J/kg (N2: 3.36e7)."
                        << exit(FatalError);
                }
            }
        }

        Info<< "IHIS: 2T chemistry active, C-V " << cvTreatment;
        if (cvPreferential_)
        {
            Info<< " (alpha = " << alphaCV_ << "): alpha*D =";
            forAll(molI_, m)
            {
                Info<< ' ' << comp.species()[molI_[m]] << ' '
                    << alphaCV_*DVib_[m];
            }
            Info<< " J/kg";
        }
        Info<< nl;
    }

    // Per-species vibrational temperature, Ns entries in species order
    // (the layout the V-T/V-V sources expect).
    // Molecules: TVib_<name>, read if present (restart), otherwise a copy of the
    // mixture TVib (values and patch types). Atoms: hold T_TR, never written.
    TVibSpecies_.setSize(Ys.size());

    forAll(Ys, s)
    {
        const word fieldName("TVib_" + comp.species()[s]);

        if (comp.isSpecieMolecular(s))
        {
            IOobject io
            (
                fieldName,
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            );

            if (io.typeHeaderOk<volScalarField>(true))
            {
                TVibSpecies_.set(s, new volScalarField(io, mesh));
                Info<< "    " << fieldName << ": read" << nl;
            }
            else
            {
                io.readOpt(IOobject::NO_READ);
                TVibSpecies_.set(s, new volScalarField(io, neThermo.TVib()));
                Info<< "    " << fieldName << ": initialised from TVib" << nl;
            }
        }
        else
        {
            TVibSpecies_.set
            (
                s,
                new volScalarField
                (
                    IOobject
                    (
                        fieldName,
                        runTime.timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE,
                        IOobject::NO_REGISTER
                    ),
                    neThermo.TTR()
                )
            );
        }
    }
    Info<< endl;
    
    // A molecule's TVib_<m> must not end up with a 'calculated'
    // wall (copied from a mixture TVib wall when 0/TVib_<m> is missing)
    forAll(molI_, m)
    {
        const volScalarField& Tvm = TVibSpecies_[molI_[m]];
        forAll(Tvm.boundaryField(), patchi)
        {
            const fvPatchScalarField& pf = Tvm.boundaryField()[patchi];
            if
            (
                !polyPatch::constraintType(mesh.boundaryMesh()[patchi].type())
             && isA<calculatedFvPatchScalarField>(pf)
            )
            {
                FatalErrorInFunction
                    << Tvm.name() << " has a 'calculated' patch '"
                    << mesh.boundary()[patchi].name() << "'." << nl
                    << "    Provide 0/" << Tvm.name() << " with a real BC there"
                    << " (e.g. neSmoluchowskiJumpTv); the mixture TVib wall is"
                    << " 'calculated' by design (Phase 7 D12)."
                    << exit(FatalError);
            }
        }
    }

    // eps_m = Y_m*e_v,m (unknowns) and e_v,m (derived), Nm entries.
    // Patch types follow the molecule's own TVib_<m> patch types (D7):
    //   constraint patch (empty, symmetry, wedge, cyclic, processor, ...) -> same type
    //   TVib mixed (jump wall) or fixes the value -> fixedValue (refreshed every
    //                                                iteration by vibBoundaries)
    //   anything else                            -> zeroGradient
    // Values are set by vibEncode() (T3.5); allocated as zero here.
    epsVib_.setSize(molI_.size());
    evVib_.setSize(molI_.size());

    forAll(molI_, m)
    {
        const label s = molI_[m];
        const word& sName = comp.species()[s];
        const volScalarField::Boundary& TvBf = TVibSpecies_[s].boundaryField();

        wordList types(TvBf.size());
        forAll(TvBf, patchi)
        {
            const word& pType = mesh.boundaryMesh()[patchi].type();

            if (polyPatch::constraintType(pType))
            {
                types[patchi] = pType;
            }
            else if (isA<mixedFvPatchScalarField>(TvBf[patchi]))
            {
                types[patchi] = fixedValueFvPatchScalarField::typeName;
            }
            else if (TvBf[patchi].fixesValue())
            {
                types[patchi] = fixedValueFvPatchScalarField::typeName;
            }
            else
            {
                types[patchi] = zeroGradientFvPatchScalarField::typeName;
            }
        }

        evVib_.set
        (
            m,
            new volScalarField
            (
                IOobject
                (
                    "evVib_" + sName,
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar(dimEnergy/dimMass, Zero),
                types
            )
        );

        epsVib_.set
        (
            m,
            new volScalarField
            (
                IOobject
                (
                    "epsVib_" + sName,
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar(dimEnergy/dimMass, Zero),
                types
            )
        );
    }

    // Molecule face mass fractions used by vibEqns.H (D3).
    // Ns > 1: refreshed every outer iteration from the species loop (T6.1).
    // Ns = 1: the species loop never runs and Y_f = 1 exactly, so this
    //         uniform value is final (T3.6).
    YfMol_.setSize(molI_.size());
    forAll(molI_, m)
    {
        YfMol_.set
        (
            m,
            new surfaceScalarField
            (
                IOobject
                (
                    "YfMol_" + comp.species()[molI_[m]],
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    IOobject::NO_REGISTER
                ),
                mesh,
                dimensionedScalar(dimless, 1.0)
            )
        );
    }

}

Info<< nl << endl;

eVibMix_.reset
(
    new volScalarField
    (
        IOobject
        (
            "eVibMix",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", e.dimensions(), 0.0)
    )
);

// Model checks, then Tv_m -> e_v,m, eps_m, mixture TVib and
// eVibMix_ = sum eps_m, before rhoE is encoded below (D8, F33)
if (twoT_)
{
    vibCheckModel();
    vibEncode();
}

scalarVars_.set
(
    1,
    new volScalarField
    (
        IOobject
        (
            "rhoE",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        rho*(e - eZero_() + eVibMix_() + 0.5*magSqr(U))
    )
);
volScalarField& rhoE = scalarVars_[1];

phi_.reset
(
    new surfaceScalarField
    (
        IOobject
        (
            "phi",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh.Sf() & fvc::interpolate(rhoU)
    )
);
surfaceScalarField& phi = phi_();


Info<< "Laminar analysis selected.\n" << endl;

Info<< "Creating flux scheme\n" << endl;
flux_ =
(
    fluxScheme::New(mesh.schemesDict(), neThermo, rho, U, rhoU, rhoE)
);
Info << endl;

rPseudoDeltaT_.reset
(
    new volScalarField
    (
        IOobject
        (
            "rPseudoDeltaT",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("rPseudoDeltaT", dimless/dimTime, 0.0),
        zeroGradientFvPatchScalarField::typeName
    )
);

ddtCoeff_.reset(new scalarField);

scalarResiduals_.setSize(2);
vectorResiduals_.setSize(1);

// Continuity residual
scalarResiduals_.set
(
    0,
    new volScalarField
    (
        IOobject
        (
            "rhoR",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", rho.dimensions()/dimTime, 0.0)
    )
);


// Mom Residual
vectorResiduals_.set
(
    0,
    new volVectorField
    (
        IOobject
        (
            "rhoUR",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector("zero", rhoU.dimensions()/dimTime, vector::zero)
    )
);

// Energy residual
scalarResiduals_.set
(
    1,
    new volScalarField
    (
        IOobject
        (
            "rhoER",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", rhoE.dimensions()/dimTime, 0.0)
    )
);

if (!inviscid_)
{
    // tauMC is stored for the benefit of the maxwellSlip BC - maybe better
    // to modify the BC to recalculate it rather
    tauMC_.reset
    (
        new volTensorField
        (
            IOobject
            (
                "tauMC",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mu*dev2(Foam::T(fvc::grad(U)))
        )
    );
    // Wall-face flux storage for the neWall* function objects
    auto makeWallScalar = [&](const word& name)
    {
        return new surfaceScalarField
        (
            IOobject(name, runTime.timeName(), mesh,
                     IOobject::NO_READ, IOobject::NO_WRITE),
            mesh,
            dimensionedScalar(dimPower, Zero)
        );
    };
    wallSigmaDotU_.reset(makeWallScalar("neWall:sigmaDotU"));
    wallQTR_.reset(makeWallScalar("neWall:qTR"));
    wallQVib_.reset(makeWallScalar("neWall:qVib"));
    wallHJ_.reset(makeWallScalar("neWall:hJ"));
    wallTau_.reset
    (
        new surfaceVectorField
        (
            IOobject("neWall:tau", runTime.timeName(), mesh,
                     IOobject::NO_READ, IOobject::NO_WRITE),
            mesh,
            dimensionedVector(dimPressure, Zero)
        )
    );
}

#include "createFvOptions.H"
fvOptions_ = &fvOptions;
