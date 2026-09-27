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

if (gMax(mu.internalField()) > 0.0)
{
    inviscid_ = false;
    Info << "Viscous analysis detected" << nl << endl;
}
else
{
    inviscid_ = true;
    Info << "Inviscid analysis detected" << nl << endl;
}

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
}

#include "createFvOptions.H"
fvOptions_ = &fvOptions;

