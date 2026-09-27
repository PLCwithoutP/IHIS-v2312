/*---------------------------------------------------------------------------*\
    mutationppAdapter.H

    Translation layer between Mutation++ and nexaFoam/neTCLib.

    Role:
      Owns exactly one Mutation::Mixture object (heavy — constructed once
      at solver startup, never copied).  Exposes one scalar-returning
      function per virtual method in tcLibraryInterface and
      compositionInterface.  All OpenFOAM field extraction (volScalarField,
      PtrList) lives in mutationPPWrapper, which calls this adapter with
      plain scalars.

    Compile-time guard:
      Entire class is inside #ifdef WITH_MUTATION_PP.  Without the flag
      the compiler sees two empty translation units — neTCLib is unchanged.

    Runtime configuration (constant/mppDict):
      mixtureName   — Mutation++ mixture XML name in constant/mpp/
      mppCaseDir    — path to case-local XML files (relative to case root)
      thetaVib      — characteristic vib temperatures per molecular species
      viscosityAlgo — optional transport override
      diffusionAlgo — optional diffusion override

    XML file resolution:
      Mutation++ searches the current working directory first.  The
      constructor temporarily changes directory to constant/mpp/ so that
      case-local mixture and mechanism XMLs are found.  Installation
      databases (RRHO, collisions, transport) are read from
      $MPP_DATA_DIRECTORY as usual.

    Thread safety note:
      The chdir() call is process-level.  Safe for MPI parallelism
      (each process has independent CWD); not safe if threads share
      a process.  OpenFOAM's standard MPI model is unaffected.

    Dirty-flag pattern:
      Species thermo arrays (speciesHOverRT, speciesCpOverR, …) are
      expensive.  Each is computed at most once per setState() call,
      lazily on first property access.  setState() marks all arrays dirty.

    Unit conventions (Mutation++ OverR / OverRT → SI per-mass):
      Cp_i  [J/kg-K] = (Cp_i/R_u)        * R_UNIV / Mw_i
      H_i   [J/kg]   = (H_i / R_u T_tr)  * R_UNIV * T_tr / Mw_i
      S_i   [J/kg-K] = (S_i / R_u)        * R_UNIV / Mw_i
      G_i   [J/kg]   = (G_i / R_u T_tr)  * R_UNIV * T_tr / Mw_i
\*---------------------------------------------------------------------------*/

#ifndef mutationppAdapter_H
#define mutationppAdapter_H

#ifdef WITH_MUTATION_PP

// OpenFOAM (for dict reading and error reporting only — no field types)
#include "IOdictionary.H"
#include "error.H"

// Mutation++ — macro collision guard required.
// OpenFOAM defines Log as a preprocessor macro (messageStream.H).
// Mutation++ defines a class named Log (Functors.h).
// Temporarily suppress the macro around the Mutation++ include.
#pragma push_macro("Log")
#undef Log
#include "mutation++.h"
#pragma pop_macro("Log")

// Standard
#include <vector>
#include <string>
#include <memory>      // std::unique_ptr
#include <cmath>
#include <climits>     // PATH_MAX
#include <unistd.h>    // getcwd, chdir

// ─────────────────────────────────────────────────────────────────────────────

namespace Foam
{

class mutationppAdapter
{
public:

    // -------------------------------------------------------------------------
    // Constructor / Destructor
    // -------------------------------------------------------------------------

    //- Construct from constant/mppDict.
    //  Reads mixtureName, mppCaseDir, thetaVib, and algorithm choices.
    //  Temporarily chdirs to constant/mpp/ so Mutation++ finds case XMLs.
    //  Expensive — construct once at solver startup.
    explicit mutationppAdapter(const IOdictionary& mppDict);

    ~mutationppAdapter() = default;

    // Non-copyable: Mixture is a heavy stateful object
    mutationppAdapter(const mutationppAdapter&)            = delete;
    mutationppAdapter& operator=(const mutationppAdapter&) = delete;


    // -------------------------------------------------------------------------
    // State
    // -------------------------------------------------------------------------

    //- Set thermochemical state for the current cell.
    //  MUST be called before any property accessor.
    //  Marks all cached arrays dirty.
    //
    //  Yi   : species mass fractions [nSpecies], Yi[i] in [0,1]
    //  rho  : mixture density [kg/m3]
    //  Ttr  : translational-rotational temperature [K]
    //  Tv   : vibrational temperature [K]
    //         (pass Ttr for 1T mixtures)
    void setState
    (
        const double* Yi,
        double        rho,
        double        Ttr,
        double        Tv
    );


    // -------------------------------------------------------------------------
    // Model information  (no setState required)
    // -------------------------------------------------------------------------

    int         nSpecies()    const;
    int         nReactions()  const;
    bool        twoTemperature() const;   // nEnergyEqns() > 1
    bool        hasChemistry()   const;   // nReactions()  > 0
    std::string speciesName(int i)          const;
    int         speciesIndex(const std::string& name) const;


    // -------------------------------------------------------------------------
    // tcLibraryInterface mapping — mixture-level scalars
    // (setState must have been called for the current cell)
    // -------------------------------------------------------------------------

    //- Pressure [Pa]                        mix_->P()
    double P() const;

    //- Translational temperature [K]        T_[0]
    double Ttr() const;

    //- Vibrational temperature [K]          T_[1]  (= Ttr for 1T)
    double Tv() const;

    //- Compressibility psi = rho/P  [s2/m2]
    //  mix_->density() / mix_->P()
    double psi() const;

    //- Translational-rotational enthalpy [J/kg]
    //  mix_->mixtureHMass() - eVib()
    double hTR() const;

    //- Vibrational energy of mixture [J/kg]
    //  SUM_i Y_i * (R_u/Mw_i) * thetaV_i / (exp(thetaV_i/Tv) - 1)
    //  Atoms contribute 0.
    double eVib() const;

    //- Mixture density [kg/m3]              mix_->density()
    double rho() const;

    //- Mixture dynamic viscosity [Pa.s]     mix_->viscosity()
    double mu() const;

    //- TR frozen thermal conductivity [W/m-K]
    //  frozenThermalConductivityVector()[0]
    double kappaTR() const;

    //- Vibrational thermal conductivity [W/m-K]
    //  frozenThermalConductivityVector()[1]  (0 for 1T)
    double kappaVib() const;

    //- Mixture-averaged diffusion coeff for species i [m2/s]
    //  averageDiffusionCoeffs(Di_)
    double Di(int i) const;

    //- Frozen ratio of specific heats Cp/Cv  [-]
    //  Cp / (Cp - R_mix)   ideal gas
    double gamma() const;

    //- Mixture TR-only frozen Cp [J/kg-K]
    //  SUM_i Y_i * (Cp_tr_i + Cp_rot_i) * R_UNIV / Mw_i
    //  Used as Newton denominator for hTR → T_tr inversion.
    double CpTRMix() const;

    // -------------------------------------------------------------------------
    // compositionInterface mapping — species-level scalars
    // -------------------------------------------------------------------------

    //- Molecular weight [kg/mol]            mix_->speciesMw(i)
    double W(int i) const;

    //- Specific gas constant [J/kg-K]       R_UNIV_ / Mw_i
    double R(int i) const;

    //- Characteristic vibrational temperature [K]
    //  Read from mppDict/thetaVib at construction. 0 for atoms.
    double thetaVib(int i) const;

    //- True if species i is molecular       mix_->species(i).nAtoms() > 1
    bool isSpecieMolecular(int i) const;

    //- TR specific heat at constant pressure [J/kg-K]
    //  (Cp_tr + Cp_rot) from speciesCpOverR(Th,Te,Tr,Tv,Tel,...)
    //  (cpt_[i] + cpr_[i]) * R_UNIV_ / Mw_i
    double CpTR(int i) const;

    //- TR specific heat at constant volume [J/kg-K]
    //  (Cv_tr + Cv_rot) from speciesCvOverR(...)
    //  (cvt_[i] + cvr_[i]) * R_UNIV_ / Mw_i
    double CvTR(int i) const;

    //- Total species enthalpy [J/kg]
    //  speciesHOverRT(Th,Te,Tr,Tv,Tel, h_,...) → h_[i]*R_UNIV_*Ttr/Mw_i
    double H(int i) const;

    //- Formation enthalpy [J/kg]
    //  hf_ optional output of speciesHOverRT → hf_[i]*R_UNIV_*Ttr/Mw_i
    double Hc(int i) const;

    //- Absolute enthalpy [J/kg]
    //  In Mutation++ H already includes formation → same as H(i)
    double Ha(int i) const;

    //- Sensible enthalpy [J/kg]             H(i) - Hc(i)
    double Hs(int i) const;

    //- Species entropy [J/kg-K]
    //  speciesSOverR(s_) → s_[i] * R_UNIV_ / Mw_i
    double S(int i) const;

    //- Species Gibbs free energy [J/kg]
    //  speciesGOverRT(Ttr, P, g_) → g_[i] * R_UNIV_ * Ttr / Mw_i
    double G(int i) const;

    //- Species dynamic viscosity [Pa.s]
    //  Mole-weighted: mu() * X[i]
    //  (per-species mu not in Mutation++ public API)
    double mu_i(int i) const;

    //- Species TR thermal conductivity [W/m-K]
    //  kappaTR() * X[i]
    double kappaTR_i(int i) const;

    //- Species vibrational thermal conductivity [W/m-K]
    //  kappaVib() * X[i]
    double kappaVib_i(int i) const;

    //- Species thermal diffusivity [m2/s]
    //  kappaTR_i(i) / (rho() * Cp_mix)
    double alphah_i(int i) const;

    //- Species partial density from ideal gas [kg/m3]
    //  P() * W(i) / (R_UNIV_ * Ttr())
    double rho_i(int i) const;

    //- Vibrational temperature from vibrational energy  [K]
    //  Newton inversion of harmonic oscillator:
    //  esVib = (R_u/Mw_i) * thetaV / (exp(thetaV/Tv) - 1)
    //  Tv0 is the initial guess (previous cell value recommended).
    //  Returns Ttr for atoms (no vibrational mode).
    double TEVib
    (
        int    i,
        double esVib_target,
        double Tv0
    ) const;

    //- Energy transfer source for vibrational mode [J/m3-s]
    //  = OmegaVT + OmegaCV (summed internally by Mutation++)
    //  Call after setState(). Returns 0 for 1T mixtures.
    double vibEnergySource() const;

    //- Electron pressure [Pa]
    //  p_e = n_e * k_B * T_ve = rho_e * R_e * T_ve
    //  Needed for p_e * div(U) work term in vibEnergyEquation.H
    double electronPressure() const;

    bool hasElectrons() const;

    //- Per-species vibrational-electronic energy [J/kg] at current state.
    //  Molecules: e_vib(i, T_ve) + e_el(i, T_ve)
    //  Atoms:     e_el(i, T_ve) only  (thetaVib_[i] = 0)
    double eVibSpecies(int i) const;

    //- Per-species Ve energy [J/kg] at a specified T_ve (does not modify adapter state).
    //  Used to diagnostically compute atomic electronic energy from known T_ve.
    double eVibAtTve(int i, double T_ve) const;
    
    // -------------------------------------------------------------------------
    // Chemistry
    // -------------------------------------------------------------------------

    //- Net species mass production rates [kg/m3-s]
    //  mix_->netProductionRates(wdot)
    //  Reactions defined in constant/mpp/<mechanism>.xml.
    //  Returns zero array if no mechanism is loaded.
    void netProductionRates(double* wdot) const;


private:

    // -------------------------------------------------------------------------
    // Mutation++ mixture
    // -------------------------------------------------------------------------

    //- Heap-allocated to avoid copy construction during chdir initialisation
    std::unique_ptr<Mutation::Mixture> mix_;


    // -------------------------------------------------------------------------
    // Cached dimensions  (set after mix_ is constructed)
    // -------------------------------------------------------------------------

    int ns_;   // nSpecies
    int nT_;   // nEnergyEqns (1 for 1T, 2 for 2T)


    // -------------------------------------------------------------------------
    // Working arrays — allocated once at construction, reused per cell
    // -------------------------------------------------------------------------

    // setState inputs
    std::vector<double> rhoi_;    // species mass densities [ns_]   kg/m3
    std::vector<double> T_;       // temperature array [nT_]        K

    // Transport (filled lazily)
    mutable std::vector<double> lambda_;   // thermal conductivities [nT_]
    mutable std::vector<double> Di_;       // diffusion coefficients [ns_]

    // Species thermo — all DIMENSIONLESS (OverR or OverRT)
    mutable std::vector<double> h_;    // H_i  / (R_u * T_tr)  [ns_]
    mutable std::vector<double> ht_;   // Ht_i / (R_u * T_tr)  [ns_]
    mutable std::vector<double> hr_;   // Hr_i / (R_u * T_tr)  [ns_]
    mutable std::vector<double> hv_;   // Hv_i / (R_u * T_tr)  [ns_]
    mutable std::vector<double> hel_;  // Hel_i/ (R_u * T_tr)  [ns_]
    mutable std::vector<double> hf_;   // Hf_i / (R_u * T_tr)  [ns_]

    mutable std::vector<double> cp_;   // Cp_i   / R_u  [ns_]
    mutable std::vector<double> cpt_;  // Cp_tr_i/ R_u  [ns_]
    mutable std::vector<double> cpr_;  // Cp_rot_i/R_u  [ns_]
    mutable std::vector<double> cpv_;  // Cp_vib_i/R_u  [ns_]
    mutable std::vector<double> cpel_; // Cp_el_i/ R_u  [ns_]

    mutable std::vector<double> cv_;   // Cv_i   / R_u  [ns_]
    mutable std::vector<double> cvt_;  // Cv_tr_i/ R_u  [ns_]
    mutable std::vector<double> cvr_;  // Cv_rot_i/R_u  [ns_]
    mutable std::vector<double> cvv_;  // Cv_vib_i/R_u  [ns_]
    mutable std::vector<double> cvel_; // Cv_el_i/ R_u  [ns_]

    mutable std::vector<double> s_;    // S_i / R_u            [ns_]
    mutable std::vector<double> g_;    // G_i / (R_u * T_tr)   [ns_]

    mutable std::vector<double> wdot_; // production rates      [ns_]
    mutable std::vector<double> omega_;   // energy transfer source [nT_]

    // -------------------------------------------------------------------------
    // Species-constant lookup tables (built at construction from mppDict)
    // -------------------------------------------------------------------------

    std::vector<double> thetaVib_;     // characteristic vib temp [K]
    std::vector<bool>   isMolecular_;  // true if nAtoms > 1


    // -------------------------------------------------------------------------
    // Dirty flags
    // Set true by setState(), set false when the array is freshly computed.
    // -------------------------------------------------------------------------

    mutable bool hDirty_;
    mutable bool cpDirty_;
    mutable bool cvDirty_;
    mutable bool sDirty_;
    mutable bool gDirty_;
    mutable bool lambdaDirty_;
    mutable bool diDirty_;
    mutable bool omegaDirty_;

    // -------------------------------------------------------------------------
    // Physical constants
    // -------------------------------------------------------------------------

    static constexpr double R_UNIV_  = 8.314462618;  // J/mol/K
    static constexpr double TEV_MIN_ = 50.0;          // Newton lower bound [K]
    static constexpr double TEV_MAX_ = 60000.0;       // Newton upper bound [K]
    static constexpr int    TEV_ITER = 50;            // max Newton iterations
    static constexpr double TEV_TOL  = 1.0e-8;       // Newton convergence


    // -------------------------------------------------------------------------
    // Private helpers
    // -------------------------------------------------------------------------

    //- Build Mutation::Mixture from mppDict.
    //  Handles chdir to mppCaseDir so case XMLs are found.
    //  Returns heap-allocated Mixture (ownership transferred to mix_).
    static std::unique_ptr<Mutation::Mixture>
    buildMixture_(const IOdictionary& mppDict);

    //- Resize all working arrays to ns_ and nT_.
    void allocateArrays_();

    //- Populate thetaVib_ and isMolecular_ from mppDict/thetaVib.
    //  FatalError if a molecular species has no thetaVib entry.
    void buildSpeciesTable_(const IOdictionary& mppDict);

    //- Mark all cached arrays dirty (called at end of setState).
    void markAllDirty_();

    //- Lazy updaters — no-op if not dirty.
    void ensureH_()      const;
    void ensureCp_()     const;
    void ensureCv_()     const;
    void ensureS_()      const;
    void ensureG_()      const;
    void ensureLambda_() const;
    void ensureOmega_()  const;
    void ensureDi_()     const;

    //- Harmonic oscillator vibrational energy for species i at Tv [J/kg].
    //  Returns 0 for atoms (thetaVib_[i] == 0).
    double esVibHO_(int i, double Tv) const;
};

} // End namespace Foam

#endif // WITH_MUTATION_PP
#endif // mutationppAdapter_H

