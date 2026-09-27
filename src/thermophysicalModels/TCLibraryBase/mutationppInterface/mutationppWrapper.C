/*---------------------------------------------------------------------------*\
    mutationppWrapper.C

    Implementation of mutationppWrapper.
    See mutationppWrapper.H for full design notes.
\*---------------------------------------------------------------------------*/

#ifdef WITH_MUTATION_PP

#include "mutationppWrapper.H"
#include "fvMesh.H"
#include "Time.H"

namespace Foam
{

// ─────────────────────────────────────────────────────────────────────────────
//  Constructor
// ─────────────────────────────────────────────────────────────────────────────

mutationppWrapper::mutationppWrapper(const fvMesh& mesh)
:
    tcLibraryInterface(),
    compositionInterface(),

    mesh_(mesh),

    // Read constant/mppDict at construction
    mppDict_
    (
        IOobject
        (
            "mppDict",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),

    // Build Mutation++ adapter from mppDict
    adapter_(mppDict_),

    // ── Primary state fields ──────────────────────────────────────────────
    // Read from case 0/ directory.  Field names match neTCLib convention
    // so existing case setups work unchanged on the Mutation++ path.

    p_
    (
        IOobject
        (
            "p",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),

    TTR_
    (
        IOobject
        (
            "TTR",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),

    TVib_
    (
        IOobject
        (
            "TVib",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),

    // psi is derived — construct with NO_READ, compute in correct()
    psi_
    (
        IOobject
        (
            "psi",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("psi", dimTime*dimTime/dimLength/dimLength, 0.0)
    ),

    // h is written by solver before correct() — start from NO_READ
    h_
    (
        IOobject
        (
            "h",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("h", dimEnergy/dimMass, 0.0)
    ),

    eTR_
    (
        IOobject
        (
            "eTR",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("eTR", dimEnergy/dimMass, 0.0)
    ),
    
    eVib_
    (
        IOobject
        (
            "eVib",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("eVib", dimEnergy/dimMass, 0.0)
    ),

    pe_
    (
        IOobject
        (
            "pe",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("pe", dimPressure, 0.0)
    ),
    
    // ── Transport fields ──────────────────────────────────────────────────
    // All constructed NO_READ — filled by correct()

    mu_
    (
        IOobject
        (
            "mu",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("mu", dimMass/dimLength/dimTime, 1.0e-5)
    ),

    kappaTR_
    (
        IOobject
        (
            "kappaTR",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("kappaTR", dimPower/dimLength/dimTemperature, 0.0)
    ),

    kappaVib_
    (
        IOobject
        (
            "kappaVib",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("kappaVib", dimPower/dimLength/dimTemperature, 0.0)
    ),

    gamma_
    (
        IOobject
        (
            "gamma",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("gamma", dimless, 1.4)
    ),

    rho_
    (
        IOobject
        (
            "rho",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("rho", dimDensity, 1.0)
    ),

    fickCoeff_
    (
        IOobject
        (
            "fickCoeff",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("fickCoeff", dimArea/dimTime, 0.0)
    ),

    // Configuration strings
    speciesListFile_  (fileName::null),
    reactionsListFile_(fileName::null),
    chemistryReader_  ("mutationpp")
{
    // Initialise species table and Y_ fields
    initSpecies_();

    // Initialise source term lists to zero
    initSourceTerms_();

    const IOdictionary thermoProps
    (
        IOobject
        (
            "thermophysicalProperties",
            mesh_.time().constant(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        )
    );
    chemistry_ = thermoProps.lookupOrDefault<bool>("chemistry", false);

    // Run initial correct() to fill all fields from initial conditions
    initForward_();

    Info<< "mutationppWrapper: initialised with "
        << adapter_.nSpecies() << " species, "
        << adapter_.nReactions() << " reactions, "
        << (adapter_.twoTemperature() ? "2T" : "1T")
        << " state model." << nl << endl;
}


// ─────────────────────────────────────────────────────────────────────────────
//  initSpecies_
// ─────────────────────────────────────────────────────────────────────────────

void mutationppWrapper::initSpecies_()
{
    const int ns = adapter_.nSpecies();

    // Build speciesTable from Mutation++ species names
    wordList speciesNames(ns);
    for (int i = 0; i < ns; i++)
        speciesNames[i] = word(adapter_.speciesName(i));

    species_ = speciesTable(speciesNames);

    // Construct Y_ — read each species mass fraction from 0/ directory
    Y_.resize(ns);
    forAll(species_, i)
    {
        Y_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "Y." + species_[i],
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_
            )
        );
    }
}


// ─────────────────────────────────────────────────────────────────────────────
//  initSourceTerms_
// ─────────────────────────────────────────────────────────────────────────────

void mutationppWrapper::initSourceTerms_()
{
    const int ns = adapter_.nSpecies();

    // TBD_ENERGY_EXCHANGE: All source term lists are zero until
    // Mutation++ energy transfer API (VT, VV) is confirmed.

    auto makeZeroList = [&](PtrList<volScalarField>& list)
    {
        list.resize(ns);
        forAll(species_, i)
        {
            list.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        "Q." + species_[i],
                        mesh_.time().timeName(),
                        mesh_,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh_,
                    dimensionedScalar
                    (
                        "Q",
                        dimEnergy/dimVolume/dimTime,
                        0.0
                    )
                )
            );
        }
    };

    makeZeroList(QVibSource_);
    makeZeroList(QVibVibSource_);
    makeZeroList(QCVSource_);
}


// ─────────────────────────────────────────────────────────────────────────────
//  correct  — main update loop
// ─────────────────────────────────────────────────────────────────────────────
void mutationppWrapper::initForward_()
{
    const int ns = adapter_.nSpecies();
    std::vector<double> Yi(ns);

    const scalarField& TtrCells  = TTR_.internalField();
    const scalarField& TVibCells = TVib_.internalField();
    const scalarField& pCells    = p_.internalField();

    scalarField& hCells         = h_.primitiveFieldRef();
    scalarField& eTRCells       = eTR_.primitiveFieldRef();
    scalarField& eVibCells      = eVib_.primitiveFieldRef();
    scalarField& psiCells       = psi_.primitiveFieldRef();
    scalarField& rhoCells       = rho_.primitiveFieldRef();
    scalarField& muCells        = mu_.primitiveFieldRef();
    scalarField& kappaTRCells   = kappaTR_.primitiveFieldRef();
    scalarField& kappaVibCells  = kappaVib_.primitiveFieldRef();
    scalarField& gammaCells     = gamma_.primitiveFieldRef();
    scalarField& fickCells      = fickCoeff_.primitiveFieldRef();

    forAll(TtrCells, cellI)
    {
        for (int i = 0; i < ns; i++) Yi[i] = Y_[i][cellI];

        double R_mix = 0.0;
        for (int i = 0; i < ns; i++) R_mix += Yi[i] * adapter_.R(i);
        const double rho0 = pCells[cellI]
                          / max(R_mix * TtrCells[cellI], 1.0e-20);

        // Forward: state is known → compute all properties
        adapter_.setState(Yi.data(), rho0, TtrCells[cellI], TVibCells[cellI]);

        hCells[cellI]        = adapter_.hTR();   // h FROM T (not inverted)
        eTRCells[cellI] = adapter_.hTR() - pCells[cellI]/rhoCells[cellI];
        eVibCells[cellI]     = adapter_.eVib();
        rhoCells[cellI]      = adapter_.rho();
        psiCells[cellI]      = adapter_.psi();
        muCells[cellI]       = adapter_.mu();
        kappaTRCells[cellI]  = adapter_.kappaTR();
        kappaVibCells[cellI] = adapter_.kappaVib();
        gammaCells[cellI]    = adapter_.gamma();
        fickCells[cellI]     = (ns > 1) ? averageDi_(cellI) : 0.0;
    }

    h_.correctBoundaryConditions();
    eTR_.correctBoundaryConditions();
    eVib_.correctBoundaryConditions();
    rho_.correctBoundaryConditions();
    psi_.correctBoundaryConditions();
    mu_.correctBoundaryConditions();
    kappaTR_.correctBoundaryConditions();
    kappaVib_.correctBoundaryConditions();
    gamma_.correctBoundaryConditions();
    fickCoeff_.correctBoundaryConditions();
}

void mutationppWrapper::correct()
{
    const int ns = adapter_.nSpecies();
    std::vector<double> Yi(ns);

    // eTR_ is the TARGET TR energy written by the solver
    // TTR_ is what we UPDATE (output)
    // TVib_ is the current vibrational temperature (input, updated by correctTVib)
    const scalarField& eTRCells = eTR_.internalField();
    const scalarField& TVibCells = TVib_.internalField();
    const scalarField& pCells    = p_.internalField();

    scalarField& TtrCells      = TTR_.primitiveFieldRef();  // OUTPUT
    scalarField& psiCells      = psi_.primitiveFieldRef();
    scalarField& eVibCells     = eVib_.primitiveFieldRef();
    scalarField& rhoCells      = rho_.primitiveFieldRef();
    scalarField& muCells       = mu_.primitiveFieldRef();
    scalarField& kappaTRCells  = kappaTR_.primitiveFieldRef();
    scalarField& kappaVibCells = kappaVib_.primitiveFieldRef();
    scalarField& gammaCells    = gamma_.primitiveFieldRef();
    scalarField& fickCells     = fickCoeff_.primitiveFieldRef();

    forAll(eTRCells, cellI)
    {
        for (int i = 0; i < ns; i++)
            Yi[i] = Y_[i][cellI];

        // Initial rho estimate from current (old) T_tr
        double R_mix = 0.0;
        for (int i = 0; i < ns; i++) R_mix += Yi[i] * adapter_.R(i);
        const double rhoGuess = pCells[cellI] / max(R_mix * TtrCells[cellI], 1.0e-20);

        // ── Step 1: Invert hTR → T_tr ─────────────────────────────────
        TtrCells[cellI] = invertETR_
        (
            Yi, 
            rhoGuess, 
            TtrCells[cellI], 
            TVibCells[cellI],
            eTRCells[cellI]     // target
        );

        // ── Step 2: Set final state ────────────────────────────────────
        const double rhoFinal = pCells[cellI]
                              / max(R_mix * TtrCells[cellI], 1.0e-20);

        adapter_.setState(Yi.data(), rhoFinal, TtrCells[cellI], TVibCells[cellI]);

        // ── Step 3: Fill all property fields ──────────────────────────
        rhoCells[cellI]      = adapter_.rho();
        psiCells[cellI]      = adapter_.psi();
        h_[cellI]            = adapter_.hTR();
        eVibCells[cellI]     = adapter_.eVib();
        pe_[cellI]           = adapter_.electronPressure();
        muCells[cellI]       = adapter_.mu();
        kappaTRCells[cellI]  = adapter_.kappaTR();
        kappaVibCells[cellI] = adapter_.kappaVib();
        gammaCells[cellI]    = adapter_.gamma();
        fickCells[cellI]     = (ns > 1) ? averageDi_(cellI) : 0.0;
    }

    // Boundary conditions
    TTR_.correctBoundaryConditions();
    psi_.correctBoundaryConditions();
    h_.correctBoundaryConditions();
    eVib_.correctBoundaryConditions();
    rho_.correctBoundaryConditions();
    mu_.correctBoundaryConditions();
    kappaTR_.correctBoundaryConditions();
    kappaVib_.correctBoundaryConditions();
    gamma_.correctBoundaryConditions();
    fickCoeff_.correctBoundaryConditions();

    // Refresh pressure field from updated psi and rho
    // (mirrors what updateTTRTemperature.H does after calling correct())
    // p is NOT updated here — that is done by the solver after correct()
}


// ─────────────────────────────────────────────────────────────────────────────
//  correctVibEnergy
// ─────────────────────────────────────────────────────────────────────────────

void mutationppWrapper::correctVibEnergy()
{
    // On Mutation++ path: eVib_ stores e_ve = e_vib + e_elec
    // Recompute from current TVib_ after it has been updated.
    // This is only called from initForward_() and applyChemistry.H.
    // During the time loop, eVib_ is updated from rhoEv/rho directly.

    const int ns = adapter_.nSpecies();
    std::vector<double> Yi(ns);

    forAll(TVib_.internalField(), cellI)
    {
        for (int i = 0; i < ns; i++) Yi[i] = Y_[i][cellI];

        const double rho0 = rho_[cellI];
        adapter_.setState(Yi.data(), rho0, TTR_[cellI], TVib_[cellI]);

        eVib_[cellI] = adapter_.eVib();   // now returns e_ve
    }

    eVib_.correctBoundaryConditions();
}


// ─────────────────────────────────────────────────────────────────────────────
//  correctTVib
// ─────────────────────────────────────────────────────────────────────────────

void mutationppWrapper::correctTVib()
{
    const int ns = adapter_.nSpecies();

    forAll(eVib_.internalField(), cellI)
    {
        int iMol = -1;
        for (int i = 0; i < ns; i++)
        {
            if (adapter_.isSpecieMolecular(i)) { iMol = i; break; }
        }

        if (iMol < 0) { TVib_[cellI] = TTR_[cellI]; continue; }

        const double Yi_mol = Y_[iMol][cellI];
        double esVe_i = (Yi_mol > SMALL) ? eVib_[cellI] / Yi_mol : 0.0;

        TVib_[cellI] = adapter_.TEVib(iMol, esVe_i, TVib_[cellI]);
    }

    TVib_.correctBoundaryConditions();
}


// ─────────────────────────────────────────────────────────────────────────────
//  correctVibSource 
// ─────────────────────────────────────────────────────────────────────────────

PtrList<volScalarField>& mutationppWrapper::correctVibSource
(
    const PtrList<volScalarField>& TVibSpecies
)
{
    // OmegaVT is a purely molecular vibrational source.
    // Distribute over molecular species only, weighted by mass fraction.
    // Atomic species receive zero — their electronic energy is handled
    // diagnostically in updateVibTemperature.H as a passive function of T_ve.

    const int ns = adapter_.nSpecies();
    std::vector<double> Yi(ns);

    const scalarField& TtrCells  = TTR_.internalField();
    const scalarField& TVibCells = TVib_.internalField();

    // Zero all species first
    forAll(species_, i)
        QVibSource_[i].primitiveFieldRef() = 0.0;

    std::vector<int> molIdx;
    for (int i = 0; i < ns; i++)
        if (adapter_.isSpecieMolecular(i)) molIdx.push_back(i);

    if (molIdx.empty()) return QVibSource_;

    forAll(TtrCells, cellI)
    {
        for (int i = 0; i < ns; i++) Yi[i] = Y_[i][cellI];

        // ── Per-species VT: evaluate OmegaVT at each species' own TVib ──────
        for (int iMol : molIdx)
        {
            const double TVib_i = TVibSpecies[iMol][cellI];

            // Set state with THIS species' vibrational temperature
            adapter_.setState(Yi.data(), rho_[cellI], TtrCells[cellI], TVib_i);

            // energyTransferSource now sees (T_tr - TVib_i) correctly
            const double Q_i = adapter_.vibEnergySource();  // [J/m3-s]

            QVibSource_[iMol][cellI] = Q_i;
        }
    }

    forAll(species_, i)
        QVibSource_[i].correctBoundaryConditions();

    // TEMPORARY DIAGNOSTIC — remove after confirming Q_VT is non-zero
    {
        std::vector<double> Yi_debug(ns);
        for (int i = 0; i < ns; i++)
            Yi_debug[i] = Y_[i][0];

        adapter_.setState
        (
            Yi_debug.data(),
            rho_[0],
            TTR_.internalField()[0],
            TVib_.internalField()[0]
        );

        Info<< "DEBUG vibEnergySource at cell 0: "
            << "T_tr=" << TTR_.internalField()[0]
            << " T_v=" << TVib_.internalField()[0]
            << " Q_VT=" << adapter_.vibEnergySource() << " J/m3-s" << nl;
    }

    return QVibSource_;
}

// ─────────────────────────────────────────────────────────────────────────────
//  correctVibVibSource  — TBD_ENERGY_EXCHANGE
// ─────────────────────────────────────────────────────────────────────────────

PtrList<volScalarField>& mutationppWrapper::correctVibVibSource
(
    const PtrList<volScalarField>& TVibSpecies
)
{
    // TBD_ENERGY_EXCHANGE:
    // Knab VV relaxation source terms are not in Mutation++'s public API.
    // QVibVibSource_ remains zero.
    return QVibVibSource_;
}


// ─────────────────────────────────────────────────────────────────────────────
//  correctCVSource  — TBD_ENERGY_EXCHANGE
// ─────────────────────────────────────────────────────────────────────────────

PtrList<volScalarField>& mutationppWrapper::correctCVSource
(
    const PtrList<volScalarField::Internal>& RR
)
{
    // OmegaCV is already included inside correctVibSource() via
    // Mutation++'s energyTransferSource() which sums all transfer
    // models registered to the vibrational mode.
    // Returning zero here avoids double-counting.
    return QCVSource_;
}


// ─────────────────────────────────────────────────────────────────────────────
//  tmp<volScalarField> accessors
// ─────────────────────────────────────────────────────────────────────────────

tmp<volScalarField> mutationppWrapper::rho() const
{
    return tmp<volScalarField>(rho_);
}

tmp<volScalarField> mutationppWrapper::mu() const
{
    return tmp<volScalarField>(mu_);
}

tmp<volScalarField> mutationppWrapper::kappaTR() const
{
    return tmp<volScalarField>(kappaTR_);
}

tmp<volScalarField> mutationppWrapper::kappaVib() const
{
    return tmp<volScalarField>(kappaVib_);
}

tmp<volScalarField> mutationppWrapper::gamma() const
{
    return tmp<volScalarField>(gamma_);
}

volScalarField mutationppWrapper::fickDiffusionCoeff() const
{
    return fickCoeff_;
}


// ─────────────────────────────────────────────────────────────────────────────
//  Configuration flags
// ─────────────────────────────────────────────────────────────────────────────

bool mutationppWrapper::twoTemperature()  const
{
    return adapter_.twoTemperature();
}

bool mutationppWrapper::chemistry() const
{
    return chemistry_;
}

bool mutationppWrapper::chemistryCapable() const
{
    return adapter_.hasChemistry();
}

const fileName& mutationppWrapper::speciesListFile() const
{
    return speciesListFile_;
}

const fileName& mutationppWrapper::reactionsListFile() const
{
    return reactionsListFile_;
}

const word& mutationppWrapper::chemistryReader() const
{
    return chemistryReader_;
}


// ─────────────────────────────────────────────────────────────────────────────
//  compositionInterface — per-species scalar properties
//  NOTE: These use the state from the LAST correct() call.
//        Re-setState per (p, TTR) argument is NOT done here because
//        Mutation++ is stateful and re-entrant setState would corrupt
//        the mixture state mid-solve.  The arguments p, TTR in the
//        function signatures exist for neTCLib API compatibility only.
// ─────────────────────────────────────────────────────────────────────────────

scalar mutationppWrapper::W(label i) const
{
    return adapter_.W(i);
}

scalar mutationppWrapper::R(label i) const
{
    return adapter_.R(i);
}

scalar mutationppWrapper::Hc(label i) const
{
    return adapter_.Hc(i);
}

scalar mutationppWrapper::CpTR(label i, scalar /*p*/, scalar /*TTR*/) const
{
    return adapter_.CpTR(i);
}

scalar mutationppWrapper::CvT(label i, scalar /*p*/, scalar /*TTR*/) const
{
    return adapter_.CvTR(i);
}

scalar mutationppWrapper::H(label i, scalar /*p*/, scalar /*TTR*/) const
{
    return adapter_.H(i);
}

scalar mutationppWrapper::Ha(label i, scalar /*p*/, scalar /*TTR*/) const
{
    return adapter_.Ha(i);
}

scalar mutationppWrapper::Hs(label i, scalar /*p*/, scalar /*TTR*/) const
{
    return adapter_.Hs(i);
}

scalar mutationppWrapper::S(label i, scalar /*p*/, scalar /*TTR*/) const
{
    return adapter_.S(i);
}

scalar mutationppWrapper::EsT(label i, scalar /*p*/, scalar TTR) const
{
    // Translational energy: (3/2) * R_i * T_tr
    return 1.5 * adapter_.R(i) * TTR;
}

scalar mutationppWrapper::EsR(label i, scalar /*p*/, scalar TTR) const
{
    // Rotational energy: R_i * T_tr for linear diatomic, 0 for atoms
    return adapter_.isSpecieMolecular(i) ? adapter_.R(i) * TTR : 0.0;
}

scalar mutationppWrapper::EsVib
(
    label i, 
    scalar /*p*/, 
    scalar /*TTR*/, 
    scalar TVib, 
    scalar /*thetaVib*/
) const
{
    // Use TVib argument explicitly so createFields.H and updateVibTemperature.H
    // get the correct energy at the specified second temperature.
    return adapter_.eVibAtTve(i, TVib);
}

tmp<volScalarField> mutationppWrapper::alphaheTR() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "alphaheTR",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            kappaTR_
           *(gamma_ - dimensionedScalar("one", dimless, 1.0))
           * psi_
           * TTR_
           / gamma_
        )
    );
}

scalar mutationppWrapper::G(label i, scalar /*p*/, scalar /*TTR*/) const
{
    return adapter_.G(i);
}

scalar mutationppWrapper::mu(label i, scalar /*TTR*/) const
{
    return adapter_.mu_i(i);
}

scalar mutationppWrapper::kappaTR(label i, scalar /*p*/, scalar /*TTR*/) const
{
    return adapter_.kappaTR_i(i);
}

scalar mutationppWrapper::kappaVib(label i, scalar /*TTR*/) const
{
    return adapter_.kappaVib_i(i);
}

scalar mutationppWrapper::alphah(label i, scalar /*p*/, scalar /*TTR*/) const
{
    return adapter_.alphah_i(i);
}

scalar mutationppWrapper::rho(label i, scalar /*p*/, scalar /*TTR*/) const
{
    return adapter_.rho_i(i);
}

scalar mutationppWrapper::thetaVib(label i) const
{
    return adapter_.thetaVib(i);
}

bool mutationppWrapper::isSpecieMolecular(label i) const
{
    return adapter_.isSpecieMolecular(i);
}

scalar mutationppWrapper::TEVib
(
    label  i,
    scalar esVib,
    scalar /*p*/,
    scalar /*TTR*/,
    scalar TVib0
) const
{
    return adapter_.TEVib(i, esVib, TVib0);
}


// ─────────────────────────────────────────────────────────────────────────────
//  Private helpers
// ─────────────────────────────────────────────────────────────────────────────

scalar mutationppWrapper::averageDi_(label cellI) const
{
    // Mass-fraction-weighted average of per-species diffusion coefficients.
    // Di_ was updated by adapter_.setState() in correct().
    double Di_avg = 0.0;
    for (int i = 0; i < adapter_.nSpecies(); i++)
        Di_avg += Y_[i][cellI] * adapter_.Di(i);
    return Di_avg;
}

scalar mutationppWrapper::esVibHO_(label i, scalar Tv) const
{
    // Harmonic oscillator vibrational energy [J/kg]
    // (R_u/Mw_i) * thetaVib_i / (exp(thetaVib_i / Tv) - 1)
    const scalar tV = adapter_.thetaVib(i);
    if (tV <= 0.0 || Tv <= 0.0) return 0.0;
    return adapter_.R(i) * tV / (std::exp(tV / Tv) - 1.0);
}

scalar mutationppWrapper::invertHTR_
(
    const std::vector<double>& Yi,
    double  rhoGuess,
    scalar  Ttr0,
    scalar  Tv,
    scalar  hTR_target
)
{
    scalar Ttr = max(Ttr0, scalar(200.0));

    for (int iter = 0; iter < 50; iter++)
    {
        // Recompute rho at current Ttr guess
        double R_mix = 0.0;
        for (int i = 0; i < adapter_.nSpecies(); i++)
            R_mix += Yi[i] * adapter_.R(i);
        double rho = max(rhoGuess * Ttr0 / Ttr, 1.0e-20);

        adapter_.setState(Yi.data(), rho, Ttr, Tv);

        const scalar hTR_current = adapter_.hTR();
        const scalar CpTR        = adapter_.CpTRMix();

        if (CpTR < SMALL) break;

        const scalar dT = (hTR_target - hTR_current) / CpTR;
        Ttr += dT;
        Ttr  = max(min(Ttr, scalar(100000.0)), scalar(200.0));

        if (mag(dT) < 1.0e-4 * Ttr) break;
    }

    return Ttr;
}

scalar mutationppWrapper::invertETR_
(
    const std::vector<double>& Yi,
    double  rhoGuess,
    scalar  Ttr0,
    scalar  Tv,
    scalar  eTR_target
)
{
    scalar Ttr = max(Ttr0, scalar(200.0));

    // Mixture gas constant — needed for eTR = hTR - R_mix*Ttr
    double R_mix = 0.0;
    for (int i = 0; i < adapter_.nSpecies(); i++)
        R_mix += Yi[i] * adapter_.R(i);

    for (int iter = 0; iter < 50; iter++)
    {
        const double rho = max(rhoGuess * Ttr0 / Ttr, 1.0e-20);

        adapter_.setState(Yi.data(), rho, Ttr, Tv);

        // eTR = hTR - p/rho = hTR - R_mix * Ttr  (ideal gas)
        const scalar eTR_current = adapter_.hTR() - R_mix * Ttr;
        // CvTR = CpTR - R_mix  (Mayer's relation)
        const scalar CvTR = adapter_.CpTRMix() - R_mix;

        if (CvTR < SMALL) break;

        const scalar dT = (eTR_target - eTR_current) / CvTR;
        Ttr += dT;
        Ttr  = max(min(Ttr, scalar(100000.0)), scalar(200.0));

        if (mag(dT) < 1.0e-4 * Ttr) break;
    }

    return Ttr;
}

} // End namespace Foam

#endif // WITH_MUTATION_PP
