/*---------------------------------------------------------------------------*\
    mutationppWrapper.H

    Concrete implementation of tcLibraryInterface and compositionInterface
    using Mutation++ as the thermochemical backend.

    Hierarchy:
      tcLibraryInterface      (abstract — solver-facing API)
          └── mutationppWrapper   (this class — Mutation++ path)
      compositionInterface    (abstract — species composition API)
          └── mutationppWrapper   (this class — also implements this)

    Activation:
      constant/thermophysicalProperties:
          library  mutationpp;

      The factory in tcLibraryInterface::New() constructs this class
      when the keyword is found AND WITH_MUTATION_PP is defined at
      compile time.  Without the compile flag, a FatalError is thrown
      at runtime.

    Configuration (constant/mppDict):
      mixtureName      — name of mixture XML in constant/mpp/
      mppCaseDir       — path to case-local Mutation++ XML files
      thetaVib         — characteristic vibrational temperatures [K]
      viscosityAlgo    — optional: Wilke | Gupta-Yos  (default Wilke)
      diffusionAlgo    — optional: Exact | Ramshaw     (default Ramshaw)

    Energy exchange source terms:
      correctVibSource()    — VT relaxation (Park model)
      correctVibVibSource() — VV relaxation (Knab model)
      correctCVSource()     — C-V coupling
      These are NOT in Mutation++'s public API.  They return zero-filled
      lists until a Mutation++ energy-transfer API is confirmed.
      Search TBD_ENERGY_EXCHANGE to find all affected sites.

    Architectural note:
      Mutation++ is STATEFUL — mix.setState() must be called before any
      property query.  All field updates happen in correct() which calls
      adapter_.setState(cellI) for every cell and fills every stored
      volScalarField.  Per-species compositionInterface methods assume
      correct() has been called for the current timestep.
\*---------------------------------------------------------------------------*/

#ifndef mutationppWrapper_H
#define mutationppWrapper_H

#ifdef WITH_MUTATION_PP

#include "tcLibraryInterface.H"
#include "compositionInterface.H"
#include "mutationppAdapter.H"
#include "volFields.H"
#include "PtrList.H"
#include "speciesTable.H"
#include "IOdictionary.H"
#include "fvMesh.H"

namespace Foam
{

class mutationppWrapper
:
    public tcLibraryInterface,
    public compositionInterface
{
public:

    // -------------------------------------------------------------------------
    // Constructor / Destructor
    // -------------------------------------------------------------------------

    //- Construct from mesh.
    //  Reads constant/mppDict, builds adapter, initialises all fields.
    explicit mutationppWrapper(const fvMesh& mesh);

    virtual ~mutationppWrapper() = default;


    // =========================================================================
    // tcLibraryInterface — pure virtual implementations
    // =========================================================================

    // ── Composition access ───────────────────────────────────────────────────

    virtual compositionInterface& composition() override
    {
        return *this;
    }

    virtual const compositionInterface& composition() const override
    {
        return *this;
    }


    // ── Thermodynamic state fields ───────────────────────────────────────────

    virtual volScalarField& p()       override { return p_; }

    virtual const volScalarField& TTR() const override { return TTR_; }

    virtual volScalarField& TTR() override { return TTR_; }

    virtual volScalarField& TVib()    override { return TVib_; }

    virtual const volScalarField& psi() const override { return psi_; }

    //- TR enthalpy [J/kg]  (solver writes this before calling correct())
    virtual volScalarField& h()       override { return h_; }
    virtual volScalarField& eTR() override { return eTR_; }
    
    //- Vibrational energy [J/kg]
    virtual volScalarField& eVib()    override { return eVib_; }

    //- Electron pressure field [Pa]
    //  Zero for non-ionized mixtures (no e- species).
    virtual tmp<volScalarField> pe() const override
    {
        return tmp<volScalarField>(pe_);
    }

    //- Mixture density [kg/m3]
    virtual tmp<volScalarField> rho() const override;


    // ── Transport properties ─────────────────────────────────────────────────

    //- Dynamic viscosity [Pa.s]
    virtual tmp<volScalarField> mu() const override;

    //- Mixture Fick diffusion coefficient [m2/s]
    virtual volScalarField fickDiffusionCoeff() const override;

    //- TR frozen thermal conductivity [W/m-K]
    virtual tmp<volScalarField> kappaTR() const override;

    //- Vibrational thermal conductivity [W/m-K]
    virtual tmp<volScalarField> kappaVib() const override;

    //- Frozen ratio of specific heats [-]
    virtual tmp<volScalarField> gamma() const override;

    virtual tmp<volScalarField> alphaheTR() const override;

    // ── Lifecycle ────────────────────────────────────────────────────────────

    //- Main update: calls adapter_.setState() per cell, fills all fields.
    virtual void correct() override;

    //- Update vibrational energy field from current TVib_
    virtual void correctVibEnergy() override;

    //- Update TVib_ field by inverting eVib_ per species
    virtual void correctTVib() override;

    //- VT energy source terms per species [J/m3-s]
    //  TBD_ENERGY_EXCHANGE: Park VT not in Mutation++ API.
    //  Returns zero-filled list until energy transfer API is confirmed.
    virtual PtrList<volScalarField>& correctVibSource
    (
        const PtrList<volScalarField>& TVibSpecies
    ) override;

    //- VV energy source terms per species [J/m3-s]
    //  TBD_ENERGY_EXCHANGE: Knab VV not in Mutation++ API.
    virtual PtrList<volScalarField>& correctVibVibSource
    (
        const PtrList<volScalarField>& TVibSpecies
    ) override;

    //- Chemistry-vibrational coupling source terms [J/m3-s]
    //  TBD_ENERGY_EXCHANGE: C-V coupling not in Mutation++ API.
    virtual PtrList<volScalarField>& correctCVSource
    (
        const PtrList<volScalarField::Internal>& RR
    ) override;


    // ── Configuration flags ──────────────────────────────────────────────────

    virtual bool twoTemperature()  const override;
    virtual bool chemistry()       const override;
    virtual bool chemistryCapable() const override;

    //- Not applicable on Mutation++ path — XML files own the species list.
    //  Returns empty fileName.
    virtual const fileName& speciesListFile()   const override;

    //- Not applicable on Mutation++ path — mechanism XML owns reactions.
    //  Returns empty fileName.
    virtual const fileName& reactionsListFile() const override;

    //- Not applicable on Mutation++ path.  Returns "mutationpp".
    virtual const word& chemistryReader() const override;


    // =========================================================================
    // compositionInterface — pure virtual implementations
    // =========================================================================

    // ── Species mass fraction fields ─────────────────────────────────────────

    virtual PtrList<volScalarField>& Y() override       { return Y_; }
    virtual const PtrList<volScalarField>& Y() const override { return Y_; }

    // ── Species table ────────────────────────────────────────────────────────

    virtual const speciesTable& species() const override { return species_; }

    // ── Per-species scalar properties ────────────────────────────────────────
    //  NOTE: These query adapter_ which uses the state from the LAST
    //  correct() call.  Do not call before correct() has run at least once.

    virtual scalar W(label i) const override;
    virtual scalar R(label i) const override;
    virtual scalar Hc(label i) const override;

    virtual scalar CpTR
    (
        label i, scalar p, scalar TTR
    ) const override;

    virtual scalar CvT
    (
        label i, scalar p, scalar TTR
    ) const override;

    virtual scalar H
    (
        label i, scalar p, scalar TTR
    ) const override;

    virtual scalar Ha
    (
        label i, scalar p, scalar TTR
    ) const override;

    virtual scalar Hs
    (
        label i, scalar p, scalar TTR
    ) const override;

    virtual scalar S
    (
        label i, scalar p, scalar TTR
    ) const override;

    virtual scalar EsT
    (
        label i, scalar p, scalar TTR
    ) const override;

    virtual scalar EsR
    (
        label i, scalar p, scalar TTR
    ) const override;

    virtual scalar EsVib
    (
        label  i,
        scalar p,
        scalar TTR,
        scalar TVib,
        scalar thetaVib
    ) const override;

    virtual scalar G
    (
        label i, scalar p, scalar TTR
    ) const override;

    virtual scalar mu
    (
        label i, scalar TTR
    ) const override;

    virtual scalar kappaTR
    (
        label i, scalar p, scalar TTR
    ) const override;

    virtual scalar kappaVib
    (
        label i, scalar TTR
    ) const override;

    virtual scalar alphah
    (
        label i, scalar p, scalar TTR
    ) const override;

    virtual scalar rho
    (
        label i, scalar p, scalar TTR
    ) const override;

    virtual scalar thetaVib(label i) const override;

    virtual bool isSpecieMolecular(label i) const override;

    virtual scalar TEVib
    (
        label  i,
        scalar esVib,
        scalar p,
        scalar TTR,
        scalar TVib0
    ) const override;

    //- Newton inversion: find T_tr such that hTR(T_tr, Tv) = hTR_target [J/kg]
    //  For RRHO, hTR is linear in T_tr → converges in 1–2 iterations.
    scalar invertHTR_
    (
        const std::vector<double>& Yi,
        double  rhoGuess,
        scalar  Ttr0,
        scalar  Tv,
        scalar  hTR_target
    );

    scalar invertETR_
    (
        const std::vector<double>& Yi,
        double  rhoGuess,
        scalar  Ttr0,
        scalar  Tv,
        scalar  eTR_target
    );
    
    // Forward initialisation: TTR_, TVib_ → h_, eVib_, transport fields.
    // Called ONCE from constructor. correct() is for solver calls only.
    void initForward_();
    
    // ── Accessors for mppEulerImplicit and mppChemistryAdapter ──────────────

    const fvMesh&        mesh()    const { return mesh_;    }
    const IOdictionary&  mppDict() const { return mppDict_; }

    mutationppAdapter&       adapter()       { return adapter_; }
    const mutationppAdapter& adapter() const { return adapter_; }

private:

    // -------------------------------------------------------------------------
    // Mesh reference
    // -------------------------------------------------------------------------

    const fvMesh& mesh_;


    // -------------------------------------------------------------------------
    // Runtime dictionary
    // -------------------------------------------------------------------------

    IOdictionary mppDict_;


    // -------------------------------------------------------------------------
    // Mutation++ translation layer
    // -------------------------------------------------------------------------

    mutationppAdapter adapter_;


    // -------------------------------------------------------------------------
    // Primary state fields (read from case 0/ directory)
    // -------------------------------------------------------------------------

    volScalarField p_;      // pressure [Pa]
    volScalarField TTR_;    // translational-rotational temperature [K]
    volScalarField TVib_;   // vibrational temperature [K]
    volScalarField psi_;    // compressibility [s2/m2]
    volScalarField h_;      // TR enthalpy [J/kg]
    volScalarField eTR_;   // TR sensible internal energy [J/kg]
    volScalarField eVib_;   // vibrational energy [J/kg]
    volScalarField pe_;    // electron pressure [Pa]  — zero for non-ionized mixtures

    // -------------------------------------------------------------------------
    // Transport fields (updated in correct(), returned as tmp)
    // -------------------------------------------------------------------------

    volScalarField mu_;         // dynamic viscosity [Pa.s]
    volScalarField kappaTR_;    // TR conductivity [W/m-K]
    volScalarField kappaVib_;   // V conductivity [W/m-K]
    volScalarField gamma_;      // ratio of specific heats [-]
    volScalarField rho_;        // mixture density [kg/m3]
    volScalarField fickCoeff_;  // Fick diffusion coefficient [m2/s]


    // -------------------------------------------------------------------------
    // Species composition (read from case 0/ directory)
    // -------------------------------------------------------------------------

    speciesTable          species_;
    PtrList<volScalarField> Y_;


    // -------------------------------------------------------------------------
    // Energy exchange source term lists
    // (TBD_ENERGY_EXCHANGE — zero until API confirmed)
    // -------------------------------------------------------------------------

    PtrList<volScalarField> QVibSource_;     // VT source [J/m3-s]
    PtrList<volScalarField> QVibVibSource_;  // VV source [J/m3-s]
    PtrList<volScalarField> QCVSource_;      // C-V source [J/m3-s]


    // -------------------------------------------------------------------------
    // Configuration strings returned by interface flags
    // -------------------------------------------------------------------------

    fileName speciesListFile_;    // empty — Mutation++ owns species
    fileName reactionsListFile_;  // empty — Mutation++ owns reactions
    word     chemistryReader_;    // "mutationpp"
    bool     chemistry_;

    // -------------------------------------------------------------------------
    // Private helpers
    // -------------------------------------------------------------------------

    //- Initialise species table and Y_ fields from adapter species list
    void initSpecies_();

    //- Initialise all source term PtrLists to zero
    void initSourceTerms_();

    //- Compute mass-fraction-weighted average diffusion coefficient [m2/s]
    //  from adapter_.Di(i) per species
    scalar averageDi_(label cellI) const;

    //- Harmonic oscillator vibrational energy [J/kg]
    //  used by EsVib and correctVibEnergy
    scalar esVibHO_(label i, scalar Tv) const;
};

} // End namespace Foam

#endif // WITH_MUTATION_PP
#endif // mutationppWrapper_H
