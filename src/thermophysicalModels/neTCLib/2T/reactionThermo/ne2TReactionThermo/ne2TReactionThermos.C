/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2023 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "make2TReactionThermo.H"

#include "ne2TReactionThermo.H"
#include "heNe2TReactionThermo.H"
#include "heNe2TThermo.H"

#include "specie2T.H"
#include "perfect2TGas.H"
#include "janaf2TThermo.H"
#include "thermo2T.H"

#include "multiComponent2TMixture.H"
#include "singleComponent2TMixture.H"
#include "neReactingMixture.H"

#include "thermoPhysics2TTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// Multi-component thermo for sensible calculation

make2TThermoPhysicsReactionThermos
(
    ne2TThermo,
    ne2TReactionThermo,
    heNe2TReactionThermo,
    multiComponent2TMixture,
    blottnerEucken2TGasCThermoPhysics
);

// Single-component thermo for sensible calculation

make2TThermoPhysicsReactionThermo
(
    ne2TReactionThermo,
    heNe2TReactionThermo,
    singleComponent2TMixture,
    blottnerEucken2TGasCThermoPhysics
);

// Reacting thermo for sensible calculation

make2TThermoPhysicsReactionThermos
(
    ne2TThermo,
    ne2TReactionThermo,
    heNe2TReactionThermo,
    neReactingMixture,
    blottnerEucken2TGasCThermoPhysics
);

// Multi-component thermo for sensible calculation

make2TThermoPhysicsReactionThermos
(
    ne2TThermo,
    ne2TReactionThermo,
    heNe2TThermo,
    multiComponent2TMixture,
    blottnerEucken2TGasCThermoPhysics
);

// Single-component thermo for sensible calculation

make2TThermoPhysicsReactionThermo
(
    ne2TReactionThermo,
    heNe2TThermo,
    singleComponent2TMixture,
    blottnerEucken2TGasCThermoPhysics
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
