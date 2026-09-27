/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
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

#include "WilkeMR.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class MixtureType>
Foam::WilkeMR<MixtureType>::WilkeMR
(
    const MixtureType& mixture
)
:
    baseMR<MixtureType>(mixture),
    mix_(mixture),
    names_(mixture.species()),
    Xi_
    (
        IOobject
        (
            "Xi",
            mixture.Y()[0].time().timeName(),  // instance: current time
            mixture.Y()[0].mesh(),             // registry: mesh (fvMesh is an objectRegistry)
            Foam::IOobject::NO_READ,
            Foam::IOobject::NO_WRITE,
            false                           // do not register in DB (scratch field)
        ),
        mixture.Y()[0].mesh(),
        Foam::dimensionedScalar("zero", Foam::dimless, 0.0)
    )
{
    precomputeXi();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class MixtureType>
void Foam::WilkeMR<MixtureType>::precomputeXi()
{
    const PtrList<volScalarField>& Y = mix_.Y();
    const label nSpec = Y.size();
    const label nCells = Y[0].primitiveField().size();

    // Resize storage if needed — happens on first call or after mesh change
    if (XiCells_.size() != nSpec)
    {
        XiCells_.setSize(nSpec);
        for (label i = 0; i < nSpec; ++i)
        {
            XiCells_[i].setSize(nCells, 0.0);
        }
    }

    // Compute denominator: sum_j(Yj/Wj) for all internal cells
    // This is a single O(n) pass over the mesh
    scalarField denom(nCells, SMALL);
    for (label j = 0; j < nSpec; ++j)
    {
        const scalarField& Yj = Y[j].primitiveField();
        const scalar Wj = mix_.W(j);
        forAll(denom, celli)
        {
            denom[celli] += Yj[celli] / Wj;
        }
    }

    // Compute Xi for each species — one O(n) pass per species
    for (label speciei = 0; speciei < nSpec; ++speciei)
    {
        const scalarField& Yi = Y[speciei].primitiveField();
        const scalar Wi = mix_.W(speciei);
        forAll(XiCells_[speciei], celli)
        {
            XiCells_[speciei][celli] = (Yi[celli] / Wi) / denom[celli];
        }
    }
}

template<class MixtureType>
Foam::volScalarField&
Foam::WilkeMR<MixtureType>::computeXiFromYi
(
    const label speciei
)
{
    const PtrList<volScalarField>& Y = mix_.Y();
    const volScalarField& Yi = Y[speciei];
    const scalar Wi = mix_.W(speciei);

    // Denominator: sum_j ( Yj / Wj )
    // Start from a small positive value to avoid divide-by-zero.
    tmp<volScalarField> tDenom
    (
        new volScalarField
        (
            IOobject
            (
                "Xi_denom",
                Yi.time().timeName(),
                Yi.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            Yi.mesh(),
            dimensionedScalar("small", dimless, SMALL)
        )
    );

    volScalarField& denom = tDenom.ref();

    for (label j = 0; j < Y.size(); ++j)
    {
        const scalar Wj = mix_.W(j);
        denom += Y[j]/Wj;            // denom += Y_j / W_j
    }

    // Xi_i = (Yi / Wi) / denom
    Xi_ = (Yi/Wi)/denom;

    return Xi_;
}

template<class MixtureType>
Foam::scalar
Foam::WilkeMR<MixtureType>::scalingFactor
(
    const label speciei,
    const label celli,
    const scalar TTR
)
{
    const label nSpec = mix_.Y().size();

    // Use precomputed scalar lookup — O(1) per cell, no field operation
    const scalar Xs   = XiCells_[speciei][celli];
    const scalar Ms   = mix_.W(speciei);
    const scalar mu_s = mix_.mu(speciei, TTR);

    scalar phi_s = Xs;

    for (label r = 0; r < nSpec; ++r)
    {
        if (r == speciei) continue;

        // O(1) lookup — no computeXiFromYi call
        const scalar Xr   = XiCells_[r][celli];
        const scalar Mr   = mix_.W(r);
        const scalar mu_r = mix_.mu(r, TTR);

        const scalar term1 =
            1.0 + Foam::sqrt(mu_s/mu_r)*Foam::pow(Mr/Ms, 0.25);

        phi_s += Xr * Foam::sqr(term1) / Foam::sqrt(8.0*(1.0 + Ms/Mr));
    }

    return phi_s;
}

template<class MixtureType>
template<class QGetter>
Foam::scalar
Foam::WilkeMR<MixtureType>::QCell
(
    const label celli,
    const scalar p,
    const scalar TTR,
    const QGetter& getQ
)
{
    const label nSpec = mix_.Y().size();

    // For a single species, Wilke mixing rule reduces to the
    // species property directly — no loop needed.
    if (nSpec == 1)
    {
        return getQ(0, p, TTR);
    }

    scalar Qmix = 0.0;

    for (label speciei = 0; speciei < nSpec; ++speciei)
    {
        const scalar Xs    = XiCells_[speciei][celli];
        const scalar phi_s = scalingFactor(speciei, celli, TTR);
        const scalar Qs    = getQ(speciei, p, TTR);

        Qmix += Qs * Xs / phi_s;
    }

    return Qmix;
}

template<class MixtureType>
Foam::scalar
Foam::WilkeMR<MixtureType>::muCell
(
    const label celli,
    const scalar p,
    const scalar TTR
)
{
    auto getMu = [&](const label speciei, const scalar p, const scalar TTRi)
    {
        return mix_.mu(speciei, TTRi);       
    };

    return QCell(celli, p, TTR, getMu);
}

template<class MixtureType>
Foam::scalar
Foam::WilkeMR<MixtureType>::kappaTRCell
(
    const label celli,
    const scalar p,
    const scalar TTR
)
{
    auto getK = [&](const label speciei, const scalar p, const scalar TTRi)
    {
        return mix_.kappaTR(speciei, p, TTRi);    
    };

    return QCell(celli, p, TTR, getK);
}

template<class MixtureType>
Foam::scalar
Foam::WilkeMR<MixtureType>::alphaTRCell
(
    const label celli,
    const scalar p,
    const scalar TTR
)
{
    auto getAlpha = [&](const label speciei, const scalar p, const scalar TTRi)
    {
        return mix_.kappaTR(speciei, p, TTRi)/mix_.CpTR(speciei, p, TTRi);    
    };

    return QCell(celli, p, TTR, getAlpha);
}

// ************************************************************************* //
