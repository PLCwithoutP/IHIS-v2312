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

#include "VVEnergySource.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class MixtureType, class MixingRule>
Foam::VVEnergySource<MixtureType, MixingRule>::VVEnergySource
(
    MixtureType& mixture,
    MixingRule& MR
)
:
    baseEnergySource<MixtureType, MixingRule>(mixture,MR),
    mix_(mixture),
    mr_(MR),
    names_(mixture.species()),
    boltzmann_const_ (1.380649e-23),
    avagadro_const_ (6.022e23),
    sigma_(0),
    P_(0),
    knabCoeffsLoaded_(false),
    Q_VVList_()
{
     
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class MixtureType, class MixingRuleType>
void Foam::VVEnergySource<MixtureType, MixingRuleType>::makeQVibSourceFields
(
    const fvMesh& mesh
)
{
    if (Q_VVList_.size())
    {
        // already built
        return;
    }

    const PtrList<volScalarField>& Y_species = mix_.Y();
    Q_VVList_.setSize(Y_species.size());

    forAll(Y_species, i)
    {
        const volScalarField& Yi = Y_species[i];

        const word fieldName("Q_VV_" + Yi.name());

        Q_VVList_.set
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
    }
}

template<class MixtureType, class MixingRule>
Foam::scalar 
Foam::VVEnergySource<MixtureType, MixingRule>::c_bar_s
(
    const scalar TTR,
    const label s 
)
{
    //Info << "Most Probable Speed: " << pow((8*mix_.R(s)*TTR)/(constant::mathematical::pi), 0.5) << nl; 
    return pow((8.0*mix_.R(s)*TTR)/(constant::mathematical::pi) , 0.5);
}

template<class MixtureType, class MixingRule>
Foam::scalar 
Foam::VVEnergySource<MixtureType, MixingRule>::n_total
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
Foam::VVEnergySource<MixtureType, MixingRule>::n_s
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
Foam::VVEnergySource<MixtureType, MixingRule>::p_s
(
    const scalar p,
    const scalar TTR,   // kept for symmetry; not used
    const label s,
    const label celli
)
{
    const scalar Xs = mr_.Xi(s, celli);
    
    return Xs*p;
}


// Knab's V–V formula: Q_{m,V-V} for one species m (=s)
template<class MixtureType, class MixingRule>
Foam::scalar 
Foam::VVEnergySource<MixtureType, MixingRule>::Q_VV_s
(
    const scalar p,      
    const scalar TTR,
    const PtrList<volScalarField>& TVibSpecies,
    const label  s,      
    const label  celli
)
{
    scalar Qm = 0.0;

    const label nSpec = mix_.Y().size();

    const scalar p_m   = p_s(p, TTR, s, celli);          
    const scalar rho_m = mix_.rho(s, p_m, TTR);          

    const scalar ev_s_Ttr =
        mix_.EsVib(s, p_m, TTR, TTR, thetai(s));         
    const scalar ev_s_Tvm =
        mix_.EsVib(s, p_m, TTR, TVibSpecies[s][celli], thetai(s));        

    const scalar Mm = Wi(s);                             
             
    //Info << "Most probable speed is : " << cbar_m << nl;
    
    for (label r = 0; r < nSpec; ++r)
    {
        if (r == s) continue;                            // r ≠ s
        if (!mix_.isSpecieMolecular(r)) continue;

        const scalar Ml = Wi(r);                         

        const scalar p_r   = p_s(p, TTR, r, celli);      
        const scalar rho_l = mix_.rho(r, p_r, TTR);      

        // Reduced molar mass [kg/kmol]
        const scalar Msr = (Mm * Ml) / (Mm + Ml);

        // Relative mean speed using reduced mass
        const scalar c_rel = Foam::sqrt
        (
            8.0 * mix_.R(s) * (Mm/Msr) * TTR
            / Foam::constant::mathematical::pi
        );

        // Molar concentration of species l [kmol/m³] — linear, NOT sqrt
        const scalar n_l = rho_l / Ml * 1000;

        const scalar velFactor = c_rel * n_l;

        //Info << "Vel factor is : " << velFactor << nl;

        //const scalar sigma_ml = sigmai(TTR);
        const scalar sigma_ml = sigma_[s][r];
        const scalar P_ml     = P_[s][r];           

        const scalar ev_r_Ttr =
            mix_.EsVib(r, p_r, TTR, TTR, thetai(r));     
        const scalar ev_r_Tvl =
            mix_.EsVib(r, p_r, TTR, TVibSpecies[r][celli], thetai(r));  
        const scalar energyBracket =
            ev_s_Ttr*(ev_r_Tvl/ev_r_Ttr) - ev_s_Tvm;

        //Info << "Energy bracket is : " << energyBracket << nl;

        Qm +=
            avagadro_const_   
          * sigma_ml          
          * P_ml              
          * velFactor
          * rho_m
          * energyBracket;
    }

    return Qm;   
}

template<class MixtureType, class MixingRule>
void Foam::VVEnergySource<MixtureType, MixingRule>::loadKnabCoeffs
(
    const fvMesh& mesh
)
{
    if (knabCoeffsLoaded_) return;

    const label nSpec = mix_.Y().size();
    sigma_ = scalarSquareMatrix(nSpec);
    sigma_ = 1.5e-18;

    P_ = scalarSquareMatrix(nSpec);
    P_ = 0.01;

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

    if (!thermophysicalProperties.found("KnabVVCoefficients"))
    {
        WarningInFunction
            << "'KnabVVCoefficients' not found in thermophysicalProperties. "
            << "Using defaults: sigma=1.5e-18, P=0.01 for all pairs." << nl;
        knabCoeffsLoaded_ = true;
        return;
    }

    const dictionary& knabDict = thermophysicalProperties.subDict("KnabVVCoefficients");

    // Read defaults first
    scalar defaultSigma = 1.5e-18;
    scalar defaultP     = 0.01;

    if (knabDict.found("default"))
    {
        const dictionary& defDict = knabDict.subDict("default");
        defDict.readIfPresent("sigma12", defaultSigma);
        defDict.readIfPresent("P21",     defaultP);
    }

    // Fill all pairs with defaults
    for (label s = 0; s < nSpec; s++)
    {
        for (label r = 0; r < nSpec; r++)
        {
            sigma_[s][r] = defaultSigma;
            P_[s][r]     = defaultP;
        }
    }

    // Override with species-pair specific values
    for (label s = 0; s < nSpec; s++)
    {
        for (label r = 0; r < nSpec; r++)
        {
            if (r == s) continue;

            const word keyFwd = names_[s] + "_" + names_[r];
            const word keyRev = names_[r] + "_" + names_[s];

            const dictionary* pairDict = nullptr;

            if (knabDict.found(keyFwd))
            {
                pairDict = &knabDict.subDict(keyFwd);
            }
            else if (knabDict.found(keyRev))
            {
                pairDict = &knabDict.subDict(keyRev);
            }

            if (pairDict)
            {
                pairDict->readIfPresent("sigma12", sigma_[s][r]);
                pairDict->readIfPresent("P21",     P_[s][r]);

                Info << "KnabVV: " << names_[s] << "-" << names_[r]
                     << " sigma=" << sigma_[s][r]
                     << " P=" << P_[s][r] << nl;
            }
        }
    }

    knabCoeffsLoaded_ = true;
}

template<class MixtureType, class MixingRule>
Foam::PtrList<Foam::volScalarField>&
Foam::VVEnergySource<MixtureType, MixingRule>::correctVibVibSource
(
    const volScalarField& p,
    const volScalarField& TTR,
    const PtrList<volScalarField>& TVibSpecies
)
{
    loadKnabCoeffs(TTR.mesh());
    makeQVibSourceFields(TTR.mesh());
    mr_.precomputeXi();
    
    const scalarField& pCells    = p.primitiveField();
    const scalarField& TTRCells  = TTR.primitiveField();

    const PtrList<volScalarField>& Y_species = mix_.Y();

    forAll(Y_species, m)
    {
        scalarField& Q_VV_Cells = Q_VVList_[m].primitiveFieldRef();

        forAll(TTRCells, celli)
        {
            const scalar pCell    = pCells[celli];
            const scalar TTRCell  = TTRCells[celli];

            Q_VV_Cells[celli] = Q_VV_s
            (
                pCell,
                TTRCell,
                TVibSpecies,
                m,
                celli
            );
        }

        Q_VVList_[m].correctBoundaryConditions();
    }

    return Q_VVList_;
}


// ************************************************************************* //
