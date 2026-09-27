/*---------------------------------------------------------------------------*\
    mutationppAdapter.C

    Implementation of mutationppAdapter.
    See mutationppAdapter.H for full design notes and unit conventions.
\*---------------------------------------------------------------------------*/

#ifdef WITH_MUTATION_PP

#include "mutationppAdapter.H"
#include <algorithm>
#include <stdexcept>

using namespace Mutation;
using namespace Mutation::Thermodynamics;

namespace Foam
{

// ─────────────────────────────────────────────────────────────────────────────
//  buildMixture_  (static private helper)
// ─────────────────────────────────────────────────────────────────────────────

std::unique_ptr<Mutation::Mixture>
mutationppAdapter::buildMixture_(const IOdictionary& mppDict)
{
    const word     mixtureName = mppDict.get<word>("mixtureName");
    const fileName mppCaseDir  = mppDict.get<fileName>("mppCaseDir");
    const word     viscAlgo    = mppDict.lookupOrDefault<word>("viscosityAlgo", "Wilke");
    const word     diffAlgo    = mppDict.lookupOrDefault<word>("diffusionAlgo", "Ramshaw");

    // ── Step 1: chdir FIRST ───────────────────────────────────────────────
    char caseRoot[PATH_MAX];
    if (::getcwd(caseRoot, sizeof(caseRoot)) == nullptr)
        FatalErrorInFunction << "getcwd() failed." << exit(FatalError);

    const std::string absMppDir =
        std::string(caseRoot) + "/" + std::string(mppCaseDir.c_str());

    if (::chdir(absMppDir.c_str()) != 0)
        FatalErrorInFunction
            << "Cannot chdir to mppCaseDir: " << absMppDir << nl
            << "Ensure the directory exists and contains the mixture XML."
            << exit(FatalError);

    // ── Step 2: MixtureOptions AFTER chdir ───────────────────────────────
    // MixtureOptions reads the XML immediately in its constructor.
    // The current directory must already be constant/mpp/ at this point.
    MixtureOptions opts(mixtureName.c_str());
    opts.setThermodynamicDatabase("RRHO");
    opts.setViscosityAlgorithm(viscAlgo.c_str());

    // ── Step 3: Construct Mixture ─────────────────────────────────────────
    auto mix = std::make_unique<Mutation::Mixture>(opts);

    // ── Step 4: Restore case root ─────────────────────────────────────────
    if (::chdir(caseRoot) != 0)
        FatalErrorInFunction
            << "Cannot chdir back to case root: " << caseRoot
            << exit(FatalError);

    mix->setDiffusionMatrixAlgo(diffAlgo.c_str());

    return mix;
}


// ─────────────────────────────────────────────────────────────────────────────
//  Constructor
// ─────────────────────────────────────────────────────────────────────────────

mutationppAdapter::mutationppAdapter(const IOdictionary& mppDict)
:
    mix_( buildMixture_(mppDict) ),
    ns_( mix_->nSpecies()    ),
    nT_( mix_->nEnergyEqns() ),

    // All dirty flags start true — no state has been set yet
    hDirty_      (true),
    cpDirty_     (true),
    cvDirty_     (true),
    sDirty_      (true),
    gDirty_      (true),
    lambdaDirty_ (true),
    omegaDirty_  (true),
    diDirty_     (true)
{
    allocateArrays_();
    buildSpeciesTable_(mppDict);

    Info<< "mutationppAdapter: loaded mixture '"
        << mix_->speciesName(0);
    for (int i = 1; i < ns_; i++) Info<< " " << mix_->speciesName(i);
    Info<< "'" << nl
        << "  nSpecies    = " << ns_    << nl
        << "  nReactions  = " << mix_->nReactions() << nl
        << "  nTempEqns   = " << nT_    << nl
        << endl;
}


// ─────────────────────────────────────────────────────────────────────────────
//  allocateArrays_
// ─────────────────────────────────────────────────────────────────────────────

void mutationppAdapter::allocateArrays_()
{
    rhoi_  .assign(ns_, 0.0);
    T_     .assign(nT_, 0.0);

    lambda_.assign(nT_, 0.0);
    Di_    .assign(ns_, 0.0);

    h_   .assign(ns_, 0.0);  ht_  .assign(ns_, 0.0);  hr_ .assign(ns_, 0.0);
    hv_  .assign(ns_, 0.0);  hel_ .assign(ns_, 0.0);  hf_ .assign(ns_, 0.0);

    cp_  .assign(ns_, 0.0);  cpt_ .assign(ns_, 0.0);  cpr_.assign(ns_, 0.0);
    cpv_ .assign(ns_, 0.0);  cpel_.assign(ns_, 0.0);

    cv_  .assign(ns_, 0.0);  cvt_ .assign(ns_, 0.0);  cvr_.assign(ns_, 0.0);
    cvv_ .assign(ns_, 0.0);  cvel_.assign(ns_, 0.0);

    s_   .assign(ns_, 0.0);
    g_   .assign(ns_, 0.0);
    wdot_.assign(ns_, 0.0);
    omega_.assign(nT_, 0.0);

    thetaVib_   .assign(ns_, 0.0);
    isMolecular_.assign(ns_, false);
}


// ─────────────────────────────────────────────────────────────────────────────
//  buildSpeciesTable_
// ─────────────────────────────────────────────────────────────────────────────

void mutationppAdapter::buildSpeciesTable_(const IOdictionary& mppDict)
{
    const dictionary& thetaDict = mppDict.subDict("thetaVib");

    for (int i = 0; i < ns_; i++)
    {
        // isSpecieMolecular: nAtoms > 1
        isMolecular_[i] = (mix_->species(i).nAtoms() > 1);

        if (isMolecular_[i])
        {
            const word name(mix_->speciesName(i));

            if (!thetaDict.found(name))
            {
                FatalErrorInFunction
                    << "Molecular species '" << name
                    << "' is in the mixture but has no thetaVib entry"
                    << " in constant/mppDict." << nl
                    << "Add the following line to the thetaVib subdictionary:"
                    << nl << "    " << name << "    <value_in_Kelvin>;" << nl
                    << exit(FatalError);
            }

            thetaVib_[i] = thetaDict.get<scalar>(name);

            Info<< "  thetaVib[" << name << "] = "
                << thetaVib_[i] << " K" << nl;
        }
        else
        {
            // Atoms: no vibrational mode
            thetaVib_[i] = 0.0;
        }
    }
}


// ─────────────────────────────────────────────────────────────────────────────
//  setState
// ─────────────────────────────────────────────────────────────────────────────

void mutationppAdapter::setState
(
    const double* Yi,
    double        rho,
    double        Ttr,
    double        Tv
)
{
    // Convert mass fractions + mixture density → species mass densities
    for (int i = 0; i < ns_; i++)
        rhoi_[i] = rho * Yi[i];

    T_[0] = Ttr;
    if (nT_ > 1) T_[1] = Tv;

    // Commit to Mutation++
    // Flag 1: state vector is {rhoi [kg/m3], T_array [K]}
    mix_->setState(rhoi_.data(), T_.data(), 1);

    markAllDirty_();
}


// ─────────────────────────────────────────────────────────────────────────────
//  markAllDirty_
// ─────────────────────────────────────────────────────────────────────────────

void mutationppAdapter::markAllDirty_()
{
    hDirty_      = true;
    cpDirty_     = true;
    cvDirty_     = true;
    sDirty_      = true;
    gDirty_      = true;
    lambdaDirty_ = true;
    omegaDirty_  = true;
    diDirty_     = true;
}


// ─────────────────────────────────────────────────────────────────────────────
//  Lazy array updaters
// ─────────────────────────────────────────────────────────────────────────────

void mutationppAdapter::ensureH_() const
{
    if (!hDirty_) return;

    // 2T mode-decomposed enthalpy.
    // Temperature mapping for heavy-particle 2T without free electrons:
    //   Th  = T_tr   (translational)
    //   Te  = T_tr   (no free electrons — conservative coupling)
    //   Tr  = T_tr   (rotation couples to translation in 2T)
    //   Tv  = T_[1]  (vibrational)
    //   Tel = T_[1]  (electronic coupled to vibrational — standard 2T)

    const double Th = T_[0];
    const double Tv = (nT_ > 1) ? T_[1] : T_[0];

    mix_->speciesHOverRT
    (
        Th, Th, Th, Tv, Tv,
        h_.data(),
        ht_.data(),
        hr_.data(),
        hv_.data(),
        hel_.data(),
        hf_.data()
    );

    hDirty_ = false;
}

void mutationppAdapter::ensureCp_() const
{
    if (!cpDirty_) return;

    const double Th = T_[0];
    const double Tv = (nT_ > 1) ? T_[1] : T_[0];

    mix_->speciesCpOverR
    (
        Th, Th, Th, Tv, Tv,
        cp_.data(),
        cpt_.data(),
        cpr_.data(),
        cpv_.data(),
        cpel_.data()
    );

    cpDirty_ = false;
}

void mutationppAdapter::ensureCv_() const
{
    if (!cvDirty_) return;

    const double Th = T_[0];
    const double Tv = (nT_ > 1) ? T_[1] : T_[0];

    mix_->speciesCvOverR
    (
        Th, Th, Th, Tv, Tv,
        cv_.data(),
        cvt_.data(),
        cvr_.data(),
        cvv_.data(),
        cvel_.data()
    );

    cvDirty_ = false;
}

void mutationppAdapter::ensureS_() const
{
    if (!sDirty_) return;
    mix_->speciesSOverR(s_.data());
    sDirty_ = false;
}

void mutationppAdapter::ensureG_() const
{
    if (!gDirty_) return;
    mix_->speciesGOverRT(T_[0], mix_->P(), g_.data());
    gDirty_ = false;
}

void mutationppAdapter::ensureLambda_() const
{
    if (!lambdaDirty_) return;
    mix_->frozenThermalConductivityVector(lambda_.data());
    lambdaDirty_ = false;
}

void mutationppAdapter::ensureDi_() const
{
    if (!diDirty_) return;

    if (ns_ < 2)
    {
        std::fill(Di_.begin(), Di_.end(), 0.0);
        diDirty_ = false;
        return;
    }

    mix_->averageDiffusionCoeffs(Di_.data());
    diDirty_ = false;
}


// ─────────────────────────────────────────────────────────────────────────────
//  esVibHO_  — harmonic oscillator vibrational energy [J/kg]
// ─────────────────────────────────────────────────────────────────────────────

double mutationppAdapter::esVibHO_(int i, double Tv) const
{
    const double thetaV = thetaVib_[i];
    if (thetaV <= 0.0 || Tv <= 0.0) return 0.0;
    return (R_UNIV_ / mix_->speciesMw(i)) * thetaV
           / (std::exp(thetaV / Tv) - 1.0);
}


// ─────────────────────────────────────────────────────────────────────────────
//  Model information
// ─────────────────────────────────────────────────────────────────────────────

int  mutationppAdapter::nSpecies()   const { return ns_; }
int  mutationppAdapter::nReactions() const { return mix_->nReactions(); }
bool mutationppAdapter::twoTemperature() const { return nT_ > 1; }
bool mutationppAdapter::hasChemistry()   const { return mix_->nReactions() > 0; }

std::string mutationppAdapter::speciesName(int i) const
{
    return mix_->speciesName(i);
}

int mutationppAdapter::speciesIndex(const std::string& name) const
{
    return mix_->speciesIndex(name);
}


// ─────────────────────────────────────────────────────────────────────────────
//  tcLibraryInterface mapping — mixture-level
// ─────────────────────────────────────────────────────────────────────────────

bool mutationppAdapter::hasElectrons() const
{
    for (int i = 0; i < ns_; i++)
        if (mix_->speciesName(i) == "e-") return true;
    return false;
}

double mutationppAdapter::P() const
{
    return mix_->P();
}

double mutationppAdapter::Ttr() const
{
    return T_[0];
}

double mutationppAdapter::Tv() const
{
    return (nT_ > 1) ? T_[1] : T_[0];
}

double mutationppAdapter::psi() const
{
    // psi = rho/P  (ideal gas: P = rho * R_mix * T_tr)
    return mix_->density() / mix_->P();
}

double mutationppAdapter::eVib() const
{
    // On Mutation++ path: eVib is actually eVe = e_vib + e_elec
    // In 2T model, T_vib = T_elec = T_ve (shared second temperature)
    // Both modes governed by T_[1]
    if (nT_ < 2) return 0.0;
    ensureH_();
    const double* Y = mix_->Y();
    double eVe = 0.0;
    for (int i = 0; i < ns_; i++)
    {
        // hv_[i] = H_vib,i / (R_u * T_tr) → e_vib,i [J/kg]
        // hel_[i] = H_el,i / (R_u * T_tr) → e_el,i [J/kg]
        const double eVib_i =
            hv_[i]  * R_UNIV_ * T_[0] / mix_->speciesMw(i);
        const double eEl_i =
            hel_[i] * R_UNIV_ * T_[0] / mix_->speciesMw(i);
        eVe += Y[i] * (eVib_i + eEl_i);
    }
    return eVe;
}

double mutationppAdapter::hTR() const
{
    // mixtureHMass() = h_tr + h_rot + h_vib + h_elec
    // eVib()         = e_vib + e_elec  (new definition)
    // hTR            = h_tr + h_rot    (correct for nexaFoam)
    return mix_->mixtureHMass() - eVib();
}

double mutationppAdapter::rho() const
{
    return mix_->density();
}

double mutationppAdapter::mu() const
{
    return mix_->viscosity();
}

double mutationppAdapter::kappaTR() const
{
    ensureLambda_();
    return lambda_[0];
}

double mutationppAdapter::kappaVib() const
{
    if (nT_ < 2) return 0.0;
    ensureLambda_();
    return lambda_[1];
}

double mutationppAdapter::Di(int i) const
{
    ensureDi_();
    return Di_[i];
}

double mutationppAdapter::gamma() const
{
    // gamma = Cp / Cv = Cp / (Cp - R_mix)   ideal gas
    const double Cp    = mix_->mixtureFrozenCpMass();
    const double R_mix = mix_->P() / (mix_->density() * T_[0]);
    return Cp / (Cp - R_mix);
}

double mutationppAdapter::CpTRMix() const
{
    ensureCp_();
    const double* Y = mix_->Y();
    double cpTR = 0.0;
    for (int i = 0; i < ns_; i++)
        cpTR += Y[i] * (cpt_[i] + cpr_[i]) * R_UNIV_ / mix_->speciesMw(i);
    return cpTR;
}

// ─────────────────────────────────────────────────────────────────────────────
//  compositionInterface mapping — species-level
// ─────────────────────────────────────────────────────────────────────────────

double mutationppAdapter::W(int i) const
{
    return mix_->speciesMw(i);
}

double mutationppAdapter::R(int i) const
{
    return R_UNIV_ / mix_->speciesMw(i);
}

double mutationppAdapter::thetaVib(int i) const
{
    return thetaVib_[i];
}

bool mutationppAdapter::isSpecieMolecular(int i) const
{
    return isMolecular_[i];
}

double mutationppAdapter::CpTR(int i) const
{
    ensureCp_();
    // CpTR_i [J/kg-K] = (Cp_tr/R_u + Cp_rot/R_u) * R_UNIV / Mw_i
    return (cpt_[i] + cpr_[i]) * R_UNIV_ / mix_->speciesMw(i);
}

double mutationppAdapter::CvTR(int i) const
{
    ensureCv_();
    return (cvt_[i] + cvr_[i]) * R_UNIV_ / mix_->speciesMw(i);
}

double mutationppAdapter::H(int i) const
{
    ensureH_();
    // h_[i] = H_i / (R_u * T_tr)
    // H_i [J/kg] = h_[i] * R_UNIV_ * T_tr / Mw_i
    return h_[i] * R_UNIV_ * T_[0] / mix_->speciesMw(i);
}

double mutationppAdapter::Hc(int i) const
{
    ensureH_();
    // hf_[i] = H_f,i / (R_u * T_tr)
    return hf_[i] * R_UNIV_ * T_[0] / mix_->speciesMw(i);
}

double mutationppAdapter::Ha(int i) const
{
    // Mutation++ H already includes formation enthalpy
    return H(i);
}

double mutationppAdapter::Hs(int i) const
{
    // Sensible = total - formation
    return H(i) - Hc(i);
}

double mutationppAdapter::S(int i) const
{
    ensureS_();
    // s_[i] = S_i / R_u
    // S_i [J/kg-K] = s_[i] * R_UNIV_ / Mw_i
    return s_[i] * R_UNIV_ / mix_->speciesMw(i);
}

double mutationppAdapter::G(int i) const
{
    ensureG_();
    // g_[i] = G_i / (R_u * T_tr)
    // G_i [J/kg] = g_[i] * R_UNIV_ * T_tr / Mw_i
    return g_[i] * R_UNIV_ * T_[0] / mix_->speciesMw(i);
}

double mutationppAdapter::mu_i(int i) const
{
    // Mole-weighted decomposition: mu_i = mu_mix * X[i]
    // Consistent with neTCLib mole-fraction mixing convention
    return mix_->viscosity() * mix_->X()[i];
}

double mutationppAdapter::kappaTR_i(int i) const
{
    ensureLambda_();
    return lambda_[0] * mix_->X()[i];
}

double mutationppAdapter::kappaVib_i(int i) const
{
    if (nT_ < 2) return 0.0;
    ensureLambda_();
    return lambda_[1] * mix_->X()[i];
}

double mutationppAdapter::alphah_i(int i) const
{
    // alphah_i = kappaTR_i / (rho * Cp_mix)
    const double Cp = mix_->mixtureFrozenCpMass();
    return kappaTR_i(i) / (mix_->density() * Cp);
}

double mutationppAdapter::rho_i(int i) const
{
    // Ideal gas per species: rho_i = P * W_i / (R_u * T_tr)
    return mix_->P() * mix_->speciesMw(i) / (R_UNIV_ * T_[0]);
}

double mutationppAdapter::TEVib
(
    int    i,
    double esVib_target, // target = e_vib + e_elec for species i
    double Tv0
) const
{
    // For atoms: only electronic contributes to Ve mode.
    // If target is negligible (T_ve ≈ 0), keep current estimate.
    if (thetaVib_[i] <= 0.0 && esVib_target < 1.0)
        return Tv0;
    // Otherwise fall through to Newton — hv_trial[i]=0, hel_trial[i]≠0
    // → converges to T_ve where H_el(i, T_ve) = esVib_target * Mw
    if (esVib_target <= 0.0)  return TEV_MIN_;

    const double Mw = mix_->speciesMw(i);

    double Tv = std::max(TEV_MIN_, std::min(Tv0, TEV_MAX_));

    for (int iter = 0; iter < TEV_ITER; iter++)
    {
        // Temporarily set Ve temperature to trial Tv to evaluate e_ve and Cv_ve
        // Save current T_tr
        double Ttr_save = T_[0];

        // Use 5-temperature speciesCvOverR(Th, Te, Tr, Tv, Tel, ...)
        // with Tve = Tv for both vibrational and electronic modes
        std::vector<double> cv_dummy(ns_), cvt_d(ns_), cvr_d(ns_),
                            cvv_trial(ns_), cvel_trial(ns_);

        mix_->speciesCvOverR
        (
            Ttr_save, Ttr_save, Ttr_save, Tv, Tv,
            cv_dummy.data(), cvt_d.data(), cvr_d.data(),
            cvv_trial.data(), cvel_trial.data()
        );

        // e_ve_i(Tv) from speciesHOverRT at Tv for vib+elec modes
        std::vector<double> h_d(ns_), ht_d(ns_), hr_d(ns_),
                            hv_trial(ns_), hel_trial(ns_), hf_d(ns_);

        mix_->speciesHOverRT
        (
            Ttr_save, Ttr_save, Ttr_save, Tv, Tv,
            h_d.data(), ht_d.data(), hr_d.data(),
            hv_trial.data(), hel_trial.data(), hf_d.data()
        );

        const double eVe_current =
            hv_trial[i]  * R_UNIV_ * Ttr_save / Mw
          + hel_trial[i] * R_UNIV_ * Ttr_save / Mw;

        const double Cv_ve =
            (cvv_trial[i] + cvel_trial[i]) * R_UNIV_ / Mw;

        if (Cv_ve < 1.0e-30) break;

        const double dTv = (esVib_target - eVe_current) / Cv_ve;
        Tv += dTv;
        Tv  = std::max(TEV_MIN_, std::min(Tv, TEV_MAX_));

        if (std::abs(dTv) < TEV_TOL * Tv) break;
    }

    return Tv;
}

double mutationppAdapter::electronPressure() const
{
    if (nT_ < 2) return 0.0;

    // Find electron species index
    for (int i = 0; i < ns_; i++)
    {
        if (mix_->speciesName(i) == "e-")
        {
            // p_e = rho_e * R_e * T_ve
            // R_e = R_UNIV / M_e
            const double rho_e = mix_->density() * mix_->Y()[i];
            return rho_e * (R_UNIV_ / mix_->speciesMw(i)) * T_[1];
        }
    }
    return 0.0;   // no electrons in mixture
}

double mutationppAdapter::eVibSpecies(int i) const
{
    if (nT_ < 2) return 0.0;
    ensureH_();
    // hv_[i] = H_vib_i(T_ve) / (R_u * T_tr)  → zero for atoms
    // hel_[i] = H_el_i(T_ve)  / (R_u * T_tr)  → non-zero for N, O etc.
    return (hv_[i] + hel_[i]) * R_UNIV_ * T_[0] / mix_->speciesMw(i);
}

double mutationppAdapter::eVibAtTve(int i, double T_ve) const
{
    // Compute e_vib+e_el for species i at T_ve without changing adapter state.
    std::vector<double> h_d(ns_), ht_d(ns_), hr_d(ns_),
                        hv_t(ns_), hel_t(ns_), hf_d(ns_);
    mix_->speciesHOverRT
    (
        T_[0], T_[0], T_[0], T_ve, T_ve,
        h_d.data(), ht_d.data(), hr_d.data(),
        hv_t.data(), hel_t.data(), hf_d.data()
    );
    return (hv_t[i] + hel_t[i]) * R_UNIV_ * T_[0] / mix_->speciesMw(i);
}

void mutationppAdapter::ensureOmega_() const
{
    if (!omegaDirty_) return;
    if (nT_ < 2)
    {
        std::fill(omega_.begin(), omega_.end(), 0.0);
        omegaDirty_ = false;
        return;
    }
    mix_->energyTransferSource(omega_.data());
    omegaDirty_ = false;
}

double mutationppAdapter::vibEnergySource() const
{
    ensureOmega_();
    // omega_[0] = vibrational mode source = OmegaVT + OmegaCV [J/m3-s]
    return (nT_ > 1) ? omega_[0] : 0.0;
}

// ─────────────────────────────────────────────────────────────────────────────
//  Chemistry
// ─────────────────────────────────────────────────────────────────────────────

void mutationppAdapter::netProductionRates(double* wdot) const
{
    if (mix_->nReactions() == 0)
    {
        std::fill(wdot, wdot + ns_, 0.0);
        return;
    }
    mix_->netProductionRates(wdot);
}

} // End namespace Foam

#endif // WITH_MUTATION_PP
