/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2015-2017 OpenCFD Ltd.
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

#include "heNe2TReactionThermo.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class BasicNe2TThermo, class MixtureType>
void Foam::heNe2TReactionThermo<BasicNe2TThermo, MixtureType>::calculate
(
    const volScalarField& p,
    volScalarField& TTR,
    volScalarField& TVib,
    volScalarField& eTR,        // was: volScalarField& h
    volScalarField& eT,
    volScalarField& eR,
    volScalarField& eVib,
    volScalarField& psi,
    volScalarField& mu,
    volScalarField& alphaTR,
    const bool doOldTimes
)
{
    if (doOldTimes && (p.nOldTimes() || TTR.nOldTimes() || TVib.nOldTimes()))
    {
        calculate
        (
            p.oldTime(),
            TTR.oldTime(),
            TVib.oldTime(),
            eTR.oldTime(),      // was: h.oldTime()
            eT.oldTime(),
            eR.oldTime(),
            eVib.oldTime(),
            psi.oldTime(),
            mu.oldTime(),
            alphaTR.oldTime(),
            true
        );
    }

    (void)TVib;
    (void)eT;
    (void)eR;
    (void)eVib;

    using mixingMixtureType = typename MixtureType::basicSpecie2TMixture;
    mixingMixtureType& mix = this->composition();

    const bool multiSpecies = mix.Y().size() > 1;
    autoPtr<WilkeMR<mixingMixtureType>> wilkeMixPtr(nullptr);

    if (multiSpecies)
    {
        wilkeMixPtr.reset(new WilkeMR<mixingMixtureType>(mix));
    }

    const scalarField& eTRCells = eTR.primitiveField();  // was: hCells
    const scalarField& pCells   = p.primitiveField();

    scalarField& TTRCells     = TTR.primitiveFieldRef();
    scalarField& psiCells     = psi.primitiveFieldRef();
    scalarField& muCells      = mu.primitiveFieldRef();
    scalarField& alphaTRCells = alphaTR.primitiveFieldRef();

    forAll(TTRCells, celli)
    {
        const typename MixtureType::thermoType& cellMixture_ =
            this->cellMixture(celli);

        if (this->updateTTR())
        {
            TTRCells[celli] = cellMixture_.TE_TR    // was: TH_TR
            (
                eTRCells[celli],                    // was: hCells[celli]
                pCells[celli],
                TTRCells[celli]
            );
        }

        psiCells[celli] = cellMixture_.psi(pCells[celli], TTRCells[celli]);

        if (multiSpecies)
        {
            muCells[celli] =
                wilkeMixPtr().muCell(celli, pCells[celli], TTRCells[celli]);

            alphaTRCells[celli] =
                wilkeMixPtr().alphaTRCell(celli, pCells[celli], TTRCells[celli]);
        }
        else
        {
            muCells[celli]      = cellMixture_.mu(TTRCells[celli]);
            alphaTRCells[celli] = cellMixture_.alphah(pCells[celli], TTRCells[celli]);
        }
    }

    const volScalarField::Boundary& pBf = p.boundaryField();

    volScalarField::Boundary& TTRBf   = TTR.boundaryFieldRef();
    volScalarField::Boundary& psiBf   = psi.boundaryFieldRef();
    volScalarField::Boundary& eTRBf   = eTR.boundaryFieldRef();  // was: hBf
    volScalarField::Boundary& muBf    = mu.boundaryFieldRef();
    volScalarField::Boundary& alphaBf = alphaTR.boundaryFieldRef();

    forAll(pBf, patchi)
    {
        const fvPatchScalarField& pp = pBf[patchi];

        fvPatchScalarField& pTTR     = TTRBf[patchi];
        fvPatchScalarField& ppsi     = psiBf[patchi];
        fvPatchScalarField& peTR     = eTRBf[patchi];  // was: ph
        fvPatchScalarField& pmu      = muBf[patchi];
        fvPatchScalarField& palphaTR = alphaBf[patchi];

        const fvPatch& patch = pp.patch();
        const labelList& fc  = patch.faceCells();

        if (pTTR.fixesValue())
        {
            forAll(pTTR, facei)
            {
                const label celli = fc[facei];
                const typename MixtureType::thermoType& cellMixture_ =
                    this->patchFaceMixture(patchi, facei);

                peTR[facei] = cellMixture_.EsT(pp[facei], pTTR[facei]);  // was: H(...)

                ppsi[facei] = cellMixture_.psi(pp[facei], pTTR[facei]);

                if (multiSpecies)
                {
                    pmu[facei] =
                        wilkeMixPtr().muCell(celli, pp[facei], pTTR[facei]);

                    palphaTR[facei] =
                        wilkeMixPtr().alphaTRCell(celli, pp[facei], pTTR[facei]);
                }
                else
                {
                    pmu[facei]      = cellMixture_.mu(pTTR[facei]);
                    palphaTR[facei] = cellMixture_.alphah(pp[facei], pTTR[facei]);
                }
            }
        }
        else
        {
            forAll(pTTR, facei)
            {
                const label celli = fc[facei];
                const typename MixtureType::thermoType& cellMixture_ =
                    this->patchFaceMixture(patchi, facei);

                if (this->updateTTR())
                {
                    pTTR[facei] = cellMixture_.TE_TR    // was: TH_TR
                    (
                        peTR[facei],                    // was: ph[facei]
                        pp[facei],
                        pTTR[facei]
                    );
                }

                ppsi[facei] = cellMixture_.psi(pp[facei], pTTR[facei]);

                if (multiSpecies)
                {
                    pmu[facei] =
                        wilkeMixPtr().muCell(celli, pp[facei], pTTR[facei]);

                    palphaTR[facei] =
                        wilkeMixPtr().alphaTRCell(celli, pp[facei], pTTR[facei]);
                }
                else
                {
                    pmu[facei]      = cellMixture_.mu(pTTR[facei]);
                    palphaTR[facei] = cellMixture_.alphah(pp[facei], pTTR[facei]);
                }
            }
        }
    }
}

template<class BasicNe2TThermo, class MixtureType>
void Foam::heNe2TReactionThermo<BasicNe2TThermo, MixtureType>::calculateVibEnergy
(
    const volScalarField& p,
    const volScalarField& TTR,
    const volScalarField& TVib,
    volScalarField& eVib
)
{
    const bool use2T = this->twoTemperature();
    const scalar thetaVib_ = this->cellMixture(0).ThetaVib(); 

    scalarField& eVibCells = eVib.primitiveFieldRef();

    const scalarField& pCells = p.primitiveField();
    const scalarField& TTRCells = TTR.primitiveField();
    const scalarField& TVibCells = TVib.primitiveField();

    forAll(TTRCells, celli)
    {
        const typename MixtureType::thermoType& cellMixture_ =
            this->cellMixture(celli);

        const scalar TVibUse = use2T ? TVibCells[celli] : TTRCells[celli];

        eVibCells[celli] = cellMixture_.EV
        (
            pCells[celli],
            TTRCells[celli],
            TVibUse,
            thetaVib_
        );
    }

    const volScalarField::Boundary& pBf = p.boundaryField();
    const volScalarField::Boundary& TTRBf = TTR.boundaryField();
    const volScalarField::Boundary& TVibBf = TVib.boundaryField();
    volScalarField::Boundary& eVibBf = eVib.boundaryFieldRef();

    forAll(pBf, patchi)
    {
        const fvPatchScalarField& pp = pBf[patchi];
        const fvPatchScalarField& pTTR = TTRBf[patchi];
        const fvPatchScalarField& pTVib = TVibBf[patchi];
        fvPatchScalarField& pesVib = eVibBf[patchi];

        if (pTTR.fixesValue())
        {
            forAll(pTTR, facei)
            {
                const typename MixtureType::thermoType& cellMixture_ =
                    this->patchFaceMixture(patchi, facei);

                const scalar TVibUse = use2T ? pTVib[facei] : pTTR[facei];

                pesVib[facei] = cellMixture_.EV
                (
                    pp[facei],
                    pTTR[facei],
                    TVibUse,
                    thetaVib_
                );

            }
        }
    }
}

template<class BasicNe2TThermo, class MixtureType>
Foam::scalar
Foam::heNe2TReactionThermo<BasicNe2TThermo, MixtureType>::wilkeKappaAverage
() const
{
    using mixingMixtureType = const typename MixtureType::basicSpecie2TMixture;
    mixingMixtureType& mix = this->composition();

    const bool pureMixture = (mix.Y().size() == 0);

    const scalarField& TTRCells = this->TTR_.primitiveField();
    const scalarField& pCells = this->p_.primitiveField();

    const fvMesh& mesh = this->TTR_.mesh();
    const scalarField& V = mesh.V();
    const scalar Vtot = gSum(V);

    scalar sumKappaV = 0.0;

    if (pureMixture)
    {
        forAll(TTRCells, celli)
        {
            const typename MixtureType::thermoType& cellMixture_ =
                this->cellMixture(celli);

            const scalar kappaCell =
                cellMixture_.kappaTR(pCells[celli], TTRCells[celli]);

            sumKappaV += kappaCell * V[celli];
        }
    }
    else
    {
        WilkeMR<mixingMixtureType> wilkeMix(mix);

        forAll(TTRCells, celli)
        {
            const scalar kappaCell =
                wilkeMix.kappaTRCell(celli, pCells[celli], TTRCells[celli]);

            sumKappaV += kappaCell * V[celli];
        }
    }

    return sumKappaV / Vtot;
}

template<class BasicNe2TThermo, class MixtureType>
Foam::scalar
Foam::heNe2TReactionThermo<BasicNe2TThermo, MixtureType>::wilkeMuAverage
() const
{
    using mixingMixtureType = const typename MixtureType::basicSpecie2TMixture;
    mixingMixtureType& mix = this->composition();

    const bool pureMixture = (mix.Y().size() == 0);

    const scalarField& TTRCells = this->TTR_.primitiveField();
    const scalarField& pCells = this->p_.primitiveField();

    const fvMesh& mesh = this->TTR_.mesh();
    const scalarField& V = mesh.V();
    const scalar Vtot = gSum(V);

    scalar sumMuV = 0.0;

    if (pureMixture)
    {
        forAll(TTRCells, celli)
        {
            const typename MixtureType::thermoType& cellMixture_ =
                this->cellMixture(celli);

            const scalar muCell = cellMixture_.mu(TTRCells[celli]);
            sumMuV += muCell * V[celli];
        }
    }
    else
    {
        WilkeMR<mixingMixtureType> wilkeMix(mix);

        forAll(TTRCells, celli)
        {
            const scalar muCell =
                wilkeMix.muCell(celli, pCells[celli], TTRCells[celli]);

            sumMuV += muCell * V[celli];
        }
    }

    return sumMuV / Vtot;
}

template<class BasicNe2TThermo, class MixtureType>
Foam::scalar
Foam::heNe2TReactionThermo<BasicNe2TThermo, MixtureType>::calculateRhoMixture
(
    const scalar pCell,
    const scalar TTRCell,
    const label  celli
) const
{
    const typename MixtureType::thermoType& mixCell =
        this->cellMixture(celli);

    return mixCell.rho(pCell, TTRCell);
}

template<class BasicNe2TThermo, class MixtureType>
Foam::scalar
Foam::heNe2TReactionThermo<BasicNe2TThermo, MixtureType>::calculateCpTRMixture
(
    const scalar pCell,
    const scalar TTRCell,
    const label  celli
) const
{
    using mixingMixtureType = const typename MixtureType::basicSpecie2TMixture;
    mixingMixtureType& comp = this->composition();

    scalar cpMix = 0.0;

    // PtrList of Y-fields, one per species
    const auto& Y = comp.Y();
    const label nSpec = Y.size();

    for (label i = 0; i < nSpec; ++i)
    {
        const scalar Yi   = Y[i][celli];   // mass fraction of species i in this cell

        const scalar cp_i = comp.CpTR(i, pCell, TTRCell);

        cpMix += Yi*cp_i;
    }

    return cpMix;   // mixture cp_TR [J/(kg·K)] if cp_i is mass-based
}


template<class BasicNe2TThermo, class MixtureType>
Foam::volScalarField
Foam::heNe2TReactionThermo<BasicNe2TThermo, MixtureType>::fickDiffusionCoeff() const
{
    using mixingMixtureType = const typename MixtureType::basicSpecie2TMixture;
    mixingMixtureType& mixingComp = this->composition();

    const fvMesh& mesh = this->TTR_.mesh();

    volScalarField Deff
    (
        IOobject
        (
            "Deff",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        mesh,
        dimensionedScalar("zero", dimArea/dimTime, 0.0)
    );

    // Pure mixture: no Fick diffusion
    if (mixingComp.Y().size() == 0)
    {
        return Deff;
    }

    WilkeMR<mixingMixtureType> wilkeMix(mixingComp);
    FickDM<mixingMixtureType> fickDM(mixingComp);

    forAll(Deff, celli)
    {
        const scalar pCell   = this->p_[celli];
        const scalar TTRCell = this->TTR_[celli];

        const scalar CpTRCell =
            this->calculateCpTRMixture(pCell, TTRCell, celli);

        const scalar kappaTRCell =
            wilkeMix.kappaTRCell(celli, pCell, TTRCell);

        const scalar rhoCell =
            this->calculateRhoMixture(pCell, TTRCell, celli);

        Deff[celli] = fickDM.DeffCell
        (
            celli,
            pCell,
            TTRCell,
            rhoCell,
            kappaTRCell,
            CpTRCell
        );
    }

    return Deff;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicNe2TThermo, class MixtureType>
Foam::heNe2TReactionThermo<BasicNe2TThermo, MixtureType>::heNe2TReactionThermo
(
    const fvMesh& mesh,
    const word& phaseName
)
:
    he2TThermo<BasicNe2TThermo, MixtureType>(mesh, phaseName)
{
    calculate
    (
        this->p_,
        this->TTR_,
        this->TVib_,
        this->eTR_,
        this->eT_,
        this->eR_,
        this->eVib_,
        this->psi_,
        this->mu_,
        this->alpha_,
        true
    );

    if (this->twoTemperature())
    {
        // Mixture eVib from the current mixture TVib.
        // This is acceptable as initialization/compatibility.
        calculateVibEnergy(this->p_, this->TTR_, this->TVib_, this->eVib_);
    }
    else
    {
        // 1T compatibility only
        this->TVib_ = this->TTR_;
        this->TVib_.correctBoundaryConditions();

        this->eVib_ =
            dimensionedScalar("zero", this->eVib_.dimensions(), 0.0);
        this->eVib_.correctBoundaryConditions();
    }
}


template<class BasicNe2TThermo, class MixtureType>
Foam::heNe2TReactionThermo<BasicNe2TThermo, MixtureType>::heNe2TReactionThermo
(
    const fvMesh& mesh,
    const word& phaseName,
    const word& dictionaryName
)
:
    he2TThermo<BasicNe2TThermo, MixtureType>(mesh, phaseName, dictionaryName)
{
    calculate
    (
        this->p_,
        this->TTR_,
        this->TVib_,
        this->eTR_,
        this->eT_,
        this->eR_,
        this->eVib_,
        this->psi_,
        this->mu_,
        this->alpha_,
        true
    );

    if (this->twoTemperature())
    {
        // Mixture eVib from the current mixture TVib.
        // This is acceptable as initialization/compatibility.
        calculateVibEnergy(this->p_, this->TTR_, this->TVib_, this->eVib_);
    }
    else
    {
        // 1T compatibility only
        this->TVib_ = this->TTR_;
        this->TVib_.correctBoundaryConditions();

        this->eVib_ =
            dimensionedScalar("zero", this->eVib_.dimensions(), 0.0);
        this->eVib_.correctBoundaryConditions();
    }
}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasicNe2TThermo, class MixtureType>
Foam::heNe2TReactionThermo<BasicNe2TThermo, MixtureType>::~heNe2TReactionThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicNe2TThermo, class MixtureType>
void Foam::heNe2TReactionThermo<BasicNe2TThermo, MixtureType>::correct()
{
    DebugInFunction << endl;

    calculate
    (
        this->p_,
        this->TTR_,
        this->TVib_,
        this->eTR_,
        this->eT_,
        this->eR_,
        this->eVib_,
        this->psi_,
        this->mu_,
        this->alpha_,
        false           // No need to update old times
    );

    DebugInFunction << "Finished" << endl;
}

template<class BasicNe2TThermo, class MixtureType>
void Foam::heNe2TReactionThermo<BasicNe2TThermo, MixtureType>::correctTEnergy()
{
    // No-op.
    // In the 2T model used here, TTR represents the combined
    // translational-rotational heavy-mode temperature.
    // Rotational energy is therefore already embedded in h_ / TTR.
}

template<class BasicNe2TThermo, class MixtureType>
void Foam::heNe2TReactionThermo<BasicNe2TThermo, MixtureType>::correctREnergy()
{
    // No-op.
    // In the 2T model used here, TTR represents the combined
    // translational-rotational heavy-mode temperature.
    // Rotational energy is therefore already embedded in h_ / TTR.
}

template<class BasicNe2TThermo, class MixtureType>
void Foam::heNe2TReactionThermo<BasicNe2TThermo, MixtureType>::correctVibEnergy()
{
    DebugInFunction << endl;

    calculateVibEnergy
    (
        this->p_,
        this->TTR_,
        this->TVib_,
        this->eVib_
    );

    DebugInFunction << "Finished" << endl;
}

template<class BasicNe2TThermo, class MixtureType>
Foam::PtrList<Foam::volScalarField>&
Foam::heNe2TReactionThermo<BasicNe2TThermo, MixtureType>::correctVibSource()
{
    static Foam::PtrList<Foam::volScalarField> emptyList;
    return emptyList;
}

template<class BasicNe2TThermo, class MixtureType>
Foam::PtrList<Foam::volScalarField>&
Foam::heNe2TReactionThermo<BasicNe2TThermo, MixtureType>::correctVibSource
(
    const PtrList<volScalarField>& TVibSpecies
)
{
    static Foam::PtrList<Foam::volScalarField> emptyList;
    if (!this->twoTemperature()) return emptyList;

    using mixingMixtureType = typename MixtureType::basicSpecie2TMixture;
    mixingMixtureType& mix = this->composition();
    using mixingRuleType = Foam::WilkeMR<mixingMixtureType>;

    static mixingRuleType wilkeMix(mix);
    wilkeMix.precomputeXi(); 
    static VTEnergySource<mixingMixtureType, mixingRuleType> Q_vt_source(mix, wilkeMix);

    return Q_vt_source.correctVibSource(this->p_, this->TTR_, TVibSpecies);
}

template<class BasicNe2TThermo, class MixtureType>
void Foam::heNe2TReactionThermo<BasicNe2TThermo, MixtureType>::correctTVib()
{
    if (!this->twoTemperature())
    {
        return;
    }

    const scalarField& pCells    = this->p_.primitiveField();
    const scalarField& TTRCells  = this->TTR_.primitiveField();
    const scalarField& eVCells   = this->eVib_.primitiveField();
    scalarField& TVibCells       = this->TVib_.primitiveFieldRef();

    forAll(TVibCells, celli)
    {
        const typename MixtureType::thermoType& cellMixture_ =
            this->cellMixture(celli);

        TVibCells[celli] = cellMixture_.TE_Vib
        (
            eVCells[celli],
            pCells[celli],
            TTRCells[celli],
            TVibCells[celli]
        );
    }

    const volScalarField::Boundary& pBf   = this->p_.boundaryField();
    const volScalarField::Boundary& TTRBf = this->TTR_.boundaryField();
    const volScalarField::Boundary& eVBf  = this->eVib_.boundaryField();
    volScalarField::Boundary& TVBf        = this->TVib_.boundaryFieldRef();

    forAll(pBf, patchi)
    {
        const fvPatchScalarField& pp   = pBf[patchi];
        const fvPatchScalarField& pTTR = TTRBf[patchi];
        const fvPatchScalarField& peV  = eVBf[patchi];
        fvPatchScalarField& pTV        = TVBf[patchi];

        forAll(pp, facei)
        {
            const typename MixtureType::thermoType& cellMixture_ =
                this->patchFaceMixture(patchi, facei);

            pTV[facei] = cellMixture_.TE_Vib
            (
                peV[facei],
                pp[facei],
                pTTR[facei],
                pTV[facei]
            );
        }
    }

    this->TVib_.correctBoundaryConditions();
}

template<class BasicNe2TThermo, class MixtureType>
Foam::PtrList<Foam::volScalarField>&
Foam::heNe2TReactionThermo<BasicNe2TThermo, MixtureType>::correctCVSource
(
    const PtrList<volScalarField::Internal>& RR
)
{
    static Foam::PtrList<Foam::volScalarField> emptyList;

    if (!this->twoTemperature())
    {
        return emptyList; // no vibrational energy mode
    }

    DebugInFunction << endl;

    using mixingMixtureType = typename MixtureType::basicSpecie2TMixture;
    mixingMixtureType& mix = this->composition();

    using mixingRuleType = Foam::WilkeMR<mixingMixtureType>;

    static mixingRuleType wilkeMix(mix);
    wilkeMix.precomputeXi(); 
    static CVEnergySource<mixingMixtureType, mixingRuleType>
        Q_cv_source(mix, wilkeMix);

    PtrList<volScalarField>& QCVlist =
        Q_cv_source.correctCVSource
        (
            RR,
            this->p_,
            this->TTR_,
            this->TVib_
        );

    DebugInFunction << "Finished" << endl;

    return QCVlist;
}


template<class BasicNe2TThermo, class MixtureType>
Foam::PtrList<Foam::volScalarField>&
Foam::heNe2TReactionThermo<BasicNe2TThermo, MixtureType>::correctVibVibSource()
{
    static Foam::PtrList<Foam::volScalarField> emptyList;
    return emptyList;
}

template<class BasicNe2TThermo, class MixtureType>
Foam::PtrList<Foam::volScalarField>&
Foam::heNe2TReactionThermo<BasicNe2TThermo, MixtureType>::correctVibVibSource
(
    const PtrList<volScalarField>& TVibSpecies
)
{
    static Foam::PtrList<Foam::volScalarField> emptyList;

    if (!this->twoTemperature())
    {
        return emptyList;
    }

    using mixingMixtureType = typename MixtureType::basicSpecie2TMixture;
    mixingMixtureType& mix = this->composition();

    using mixingRuleType = Foam::WilkeMR<mixingMixtureType>;

    static mixingRuleType wilkeMix(mix);
    wilkeMix.precomputeXi(); 
    static VVEnergySource<mixingMixtureType, mixingRuleType> Q_vv_source(mix, wilkeMix);

    return Q_vv_source.correctVibVibSource(this->p_, this->TTR_, TVibSpecies);
}

template<class BasicNe2TThermo, class MixtureType>
Foam::PtrList<Foam::volScalarField>&
Foam::heNe2TReactionThermo<BasicNe2TThermo, MixtureType>::correctVTRelaxationTime()
{
    static Foam::PtrList<Foam::volScalarField> emptyList;

    if (!this->twoTemperature())
    {
        return emptyList; // no vibrational energy mode
    }

    DebugInFunction << endl;

    using mixingMixtureType = typename MixtureType::basicSpecie2TMixture;
    mixingMixtureType& mix = this->composition();

    using mixingRuleType = Foam::WilkeMR<mixingMixtureType>;

    static mixingRuleType wilkeMix(mix);
    wilkeMix.precomputeXi(); 
    static VTEnergySource<mixingMixtureType, mixingRuleType> Q_vt_source(mix, wilkeMix);

    PtrList<volScalarField>& TauList =
        Q_vt_source.correctVTRelaxationTime(this->p_, this->TTR_);

    DebugInFunction << "Finished" << endl;

    return TauList;
}

// ************************************************************************* //
