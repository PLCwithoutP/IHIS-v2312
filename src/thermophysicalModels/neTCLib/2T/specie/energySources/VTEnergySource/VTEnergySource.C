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

#include "VTEnergySource.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class MixtureType, class MixingRule>
Foam::VTEnergySource<MixtureType, MixingRule>::VTEnergySource
(
    MixtureType& mixture,
    MixingRule& MR
)
:
    baseEnergySource<MixtureType, MixingRule>(mixture,MR),
    mix_(mixture),
    mr_(MR),
    names_(mixture.species()),
    sigma_prime_ (3e-21), // N2, O2, NO
    boltzmann_const_ (1.380649e-23),
    AMW_(0),
    BMW_(0),
    mwCoeffsLoaded_(false),
    Q_VTList_()
{

}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class MixtureType, class MixingRule>
void Foam::VTEnergySource<MixtureType, MixingRule>::loadMWCoeffs
(
    const fvMesh& mesh
)
{
    if (mwCoeffsLoaded_) return;

    const label nSpec = mix_.Y().size();
    AMW_ = scalarSquareMatrix(nSpec);
    BMW_ = scalarSquareMatrix(nSpec);

    // 1) Fill every pair with the Millikan-White correlation
    //    (identical to the previous hard-coded behaviour)
    for (label s = 0; s < nSpec; ++s)
    {
        for (label r = 0; r < nSpec; ++r)
        {
            const scalar mu_sr = (Wi(s)*Wi(r))/(Wi(s) + Wi(r));
            AMW_[s][r] = 1.16e-3*pow(mu_sr, 0.5)*pow(thetai(s), 4.0/3.0);
            BMW_[s][r] = 0.015*pow(mu_sr, 0.25);
        }
    }

    // 2) Override from thermophysicalProperties, where entries exist
    IOdictionary thermophysicalProperties
    (
        IOobject
        (
            "thermophysicalProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    if (!thermophysicalProperties.found("MillikanWhiteCoefficients"))
    {
        WarningInFunction
            << "'MillikanWhiteCoefficients' not found in thermophysicalProperties. "
            << "Using the MW correlation for all pairs." << nl;
        mwCoeffsLoaded_ = true;
        return;
    }

    const dictionary& mwDict =
        thermophysicalProperties.subDict("MillikanWhiteCoefficients");

    // Molecule-first keys ONLY (RELAXER_COLLIDER). No reversed-key
    // fallback: N2_O and O_N2 are physically different statements,
    // and atomic relaxers never use these values anyway.
    for (label s = 0; s < nSpec; ++s)
    {
        if (!mix_.isSpecieMolecular(s)) continue;

        for (label r = 0; r < nSpec; ++r)
        {
            const word key = names_[s] + "_" + names_[r];

            if (mwDict.found(key))
            {
                const dictionary& pairDict = mwDict.subDict(key);
                pairDict.readIfPresent("A", AMW_[s][r]);
                pairDict.readIfPresent("B", BMW_[s][r]);

                Info<< "MW override: " << key
                    << " A=" << AMW_[s][r]
                    << " B=" << BMW_[s][r] << nl;
            }
        }
    }

    mwCoeffsLoaded_ = true;
}

template<class MixtureType, class MixingRuleType>
void Foam::VTEnergySource<MixtureType, MixingRuleType>::makeQVibSourceFields
(
    const fvMesh& mesh
)
{
    if (Q_VTList_.size() || tauVTList_.size())
    {
        // already built
        return;
    }

    const PtrList<volScalarField>& Y_species = mix_.Y();
    Q_VTList_.setSize(Y_species.size());
    tauVTList_.setSize(Y_species.size());

    forAll(Y_species, i)
    {
        const volScalarField& Yi = Y_species[i];

        const word fieldName("Q_VT_" + Yi.name());
        const word fieldNameTau("tau_VT_" + Yi.name());

        Q_VTList_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    fieldName,
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar
                (
                    "zero",
                    dimEnergy/dimVolume/dimTime,
                    0.0
                )
            )
        );
        tauVTList_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    fieldNameTau,
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE   // or AUTO_WRITE if you want to write them
                ),
                mesh,
                dimensionedScalar
                (
                    "zero",
                    dimTime,  // [s] – relaxation time
                    0.0
                )
            )
        );
    }
}

template<class MixtureType, class MixingRule>
Foam::scalar 
Foam::VTEnergySource<MixtureType, MixingRule>::sigmai
(
    const scalar TTR 
)
{
    //Info << "Cross Section: " << sigma_prime_*(50000.0/TTR)*(50000.0/TTR) << nl; 
    return sigma_prime_*(50000.0/TTR)*(50000.0/TTR);
}

template<class MixtureType, class MixingRule>
Foam::scalar 
Foam::VTEnergySource<MixtureType, MixingRule>::c_bar_s
(
    const scalar TTR,
    const label s 
)
{
    //Info << "Most Probable Speed: " << pow((8*mix_.R(s)*TTR)/(constant::mathematical::pi), 0.5) << nl; 
    return pow((8*mix_.R(s)*TTR)/(constant::mathematical::pi) , 0.5);
}

template<class MixtureType, class MixingRule>
Foam::scalar 
Foam::VTEnergySource<MixtureType, MixingRule>::n_total
(
    const scalar p,
    const scalar TTR
)
{
    //Info << "Total Number Density: " << p/(TTR*boltzmann_const_) << nl; 
    return p/(TTR*boltzmann_const_);
}

template<class MixtureType, class MixingRule>
Foam::scalar 
Foam::VTEnergySource<MixtureType, MixingRule>::n_s
(
    const scalar p,
    const scalar TTR,
    const label s,
    const label celli
)
{
    const scalar nTot = n_total(p, TTR);
    const scalar Xr = mr_.Xi(s, celli);                  
    return Xr*nTot;
}

template<class MixtureType, class MixingRule>
Foam::scalar 
Foam::VTEnergySource<MixtureType, MixingRule>::p_s
(
    const scalar p,
    const scalar TTR,   // kept for symmetry; not used
    const label s,
    const label celli
)
{
    const scalar Xs = mr_.Xi(s, celli);
    const scalar pSpecie = Xs*p;
    return pSpecie;
}


template<class MixtureType, class MixingRule>
Foam::scalar 
Foam::VTEnergySource<MixtureType, MixingRule>::T_P_sr
(
    const scalar p,
    const scalar TTR,
    const label s,
    const label r,
    const label celli
)
{
    //Info << "Tau Park: " << 1/(c_bar_s(TTR,s) * sigmai(TTR) * n_s(p, TTR, s, celli)) << nl; 
    return 1/(c_bar_s(TTR,s) * sigmai(TTR) * n_s(p, TTR, r, celli));
}

template<class MixtureType, class MixingRule>
Foam::scalar 
Foam::VTEnergySource<MixtureType, MixingRule>::A_sr
(
    const label s,
    const label r
)
{
    return AMW_[s][r];
}

template<class MixtureType, class MixingRule>
Foam::scalar 
Foam::VTEnergySource<MixtureType, MixingRule>::B_sr
(
    const label s,
    const label r
)
{
    return BMW_[s][r];
}

template<class MixtureType, class MixingRule>
Foam::scalar 
Foam::VTEnergySource<MixtureType, MixingRule>::T_MW_sr
(
    const scalar p,
    const scalar TTR,
    const label s,
    const label r
)
{
    scalar A = A_sr(s,r);
    scalar B = B_sr(s,r);
    //Info << "Tau Millikan-White: " << (1/p)*exp(A * (pow(TTR, -1/3) - B) - 18.42) << nl; 
    const scalar p_atm = p / 101325.0;
    return (1/p_atm)*exp(A * (pow(TTR, -1.0/3.0) - B) - 18.42);
}

template<class MixtureType, class MixingRule>
Foam::scalar 
Foam::VTEnergySource<MixtureType, MixingRule>::T_sr
(
    const scalar p,
    const scalar TTR,
    const label s,
    const label r,
    const label celli
)
{
    //Info << "Tau total: " << T_MW_sr(p, TTR, s, r) + T_P_sr(p, TTR, s) << nl; 
    return T_MW_sr(p, TTR, s, r) + T_P_sr(p, TTR, s, r, celli);
}

template<class MixtureType, class MixingRule>
Foam::scalar 
Foam::VTEnergySource<MixtureType, MixingRule>::T_s
(
    const scalar p,
    const scalar TTR,
    const label s,
    const label celli
)
{
    // Atomic species have no vibrational mode: no V-T relaxation.
    if (!mix_.isSpecieMolecular(s))
    {
        return GREAT;
    }

    scalar sumX            = 0.0;
    scalar sumX_over_tauMW = 0.0;
    const label nSpec = mix_.Y().size();

    // Mixture-averaged Millikan-White time only.
    // The Park limiting-cross-section correction is NOT a pairwise
    // quantity and must not enter this harmonic average.
    for (label r = 0; r < nSpec; ++r)
    {
        const scalar Xr = mr_.Xi(r, celli);

        if (Xr <= SMALL) continue;

        const scalar tau_sr_MW = T_MW_sr(p, TTR, s, r);

        if (tau_sr_MW <= SMALL) continue;

        sumX            += Xr;
        sumX_over_tauMW += Xr/tau_sr_MW;
    }

    // safety: no valid colliders -> avoid 0/0
    if (sumX_over_tauMW <= SMALL)
    {
        return GREAT;
    }

    const scalar tau_MW = sumX/sumX_over_tauMW;

    // Park correction, evaluated once per species using the TOTAL
    // number density (not the collider partial density).
    const scalar nTot = n_total(p, TTR);

    scalar tau_P = GREAT;
    const scalar denom = c_bar_s(TTR, s)*sigmai(TTR)*nTot;

    if (denom > SMALL)
    {
        tau_P = 1.0/denom;
    }

    return tau_MW + tau_P;
}

/*template<class MixtureType, class MixingRule>
Foam::scalar 
Foam::VTEnergySource<MixtureType, MixingRule>::T_s
(
    const scalar p,
    const scalar TTR,
    const label s,
    const label celli
)
{
    scalar sumX          = 0.0;
    scalar sumX_over_tau = 0.0;
    const label nSpec = mix_.Y().size();
    
    for (label r = 0; r < nSpec; ++r)
    {
        const scalar Xr = mr_.Xi(r, celli);

        if (Xr <= SMALL) continue;

        const scalar tau_sr = T_sr(p, TTR, s, r, celli);

        sumX          += Xr;
        sumX_over_tau += Xr/tau_sr;
    }

    // safety: no valid colliders → avoid 0/0
    if (sumX_over_tau <= SMALL)
    {
        return GREAT;
    }

    return sumX/sumX_over_tau;
}*/

// Classical Landau-Teller Formula
template<class MixtureType, class MixingRule>
Foam::scalar 
Foam::VTEnergySource<MixtureType, MixingRule>::Q_VT_s
(
    const scalar p,
    const scalar TTR,
    const scalar TVib,
    const label s,
    const label celli
)
{
    return mix_.rho(s,p_s(p, TTR, s, celli),TTR)*(mix_.EsVib(s,p,TTR,TTR,thetai(s)) 
    - mix_.EsVib(s,p,TTR,TVib,thetai(s)))/(T_s(p, TTR, s, celli));
}

template<class MixtureType, class MixingRule>
Foam::PtrList<Foam::volScalarField>&
Foam::VTEnergySource<MixtureType, MixingRule>::correctVibSource
(
    const volScalarField& p,
    const volScalarField& TTR,
    const PtrList<volScalarField>& TVibSpecies
)
{
    makeQVibSourceFields(TTR.mesh());
    loadMWCoeffs(TTR.mesh());

    mr_.precomputeXi();

    const scalarField& pCells   = p.primitiveField();
    const scalarField& TTRCells = TTR.primitiveField();
    const PtrList<volScalarField>& Y_species = mix_.Y();

    forAll(Y_species, speciei)
    {
        scalarField& Q_VT_Cells = Q_VTList_[speciei].primitiveFieldRef();
        const scalarField& TVibCells = TVibSpecies[speciei].primitiveField();

        forAll(TTRCells, celli)
        {
            Q_VT_Cells[celli] = Q_VT_s
            (
                pCells[celli],
                TTRCells[celli],
                TVibCells[celli],    // species-own Tv
                speciei,
                celli
            );
        }
    }

    return Q_VTList_;
}

template<class MixtureType, class MixingRule>
Foam::PtrList<Foam::volScalarField>&
Foam::VTEnergySource<MixtureType, MixingRule>::correctVTRelaxationTime
(
    const volScalarField& p,
    const volScalarField& TTR
)
{
    makeQVibSourceFields(TTR.mesh());
    loadMWCoeffs(TTR.mesh());
    mr_.precomputeXi(); 
    
    const scalarField& pCells   = p.primitiveField();
    const scalarField& TTRCells = TTR.primitiveField();

    const label nSpec = mix_.Y().size();
    const label nCells = TTRCells.size();

    for (label s = 0; s < nSpec; ++s)
    {
        // If you want to skip non-molecular species entirely:
        if (!mix_.isSpecieMolecular(s))
        {
            scalarField& tauCells = tauVTList_[s].primitiveFieldRef();
            forAll(tauCells, celli)
            {
                tauCells[celli] = GREAT;
            }
            tauVTList_[s].correctBoundaryConditions();
            continue;
        }

        scalarField& tauCells = tauVTList_[s].primitiveFieldRef();

        for (label celli = 0; celli < nCells; ++celli)
        {
            const scalar pCell   = pCells[celli];
            const scalar TTRCell = TTRCells[celli];

            const scalar tau_sVT =
                T_s
                (
                    pCell,
                    TTRCell,
                    s,
                    celli
                );

            tauCells[celli] = tau_sVT;
        }

        tauVTList_[s].correctBoundaryConditions();
    }

    return tauVTList_;
}


// ************************************************************************* //
