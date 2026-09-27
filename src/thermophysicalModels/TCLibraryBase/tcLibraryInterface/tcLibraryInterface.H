/*---------------------------------------------------------------------------*\
    tcLibraryInterface.H

    Abstract base class for all thermochemical library backends.
    Derived by:
      - ne2TReactionThermo  (neTCLib native path)
      - mutationPPWrapper   (Mutation++ path, future)

    Every pure virtual corresponds to at least one call site in the solver.
    composition() is NOT declared here — that belongs to Step 2
    (compositionInterface).

    Factory reads 'library' keyword from thermophysicalProperties:
      library  neTCLib;     -> ne2TReactionThermo  (default)
      library  default;     -> ne2TReactionThermo  (alias)
      library  mutationpp;  -> FatalError (not yet implemented)
\*---------------------------------------------------------------------------*/

#ifndef tcLibraryInterface_H
#define tcLibraryInterface_H

#include "volFields.H"
#include "surfaceFields.H"
#include "fvMesh.H"
#include "autoPtr.H"
#include "PtrList.H"
#include "fileName.H"
#include "compositionInterface.H"

namespace Foam
{

class tcLibraryInterface
{
public:

    TypeName("tcLibraryInterface");

    virtual ~tcLibraryInterface() {}

    // -----------------------------------------------------------------------
    // Factory
    // -----------------------------------------------------------------------

    static autoPtr<tcLibraryInterface> New(const fvMesh& mesh);

    virtual compositionInterface& composition() = 0;
    
    virtual const compositionInterface& composition() const = 0;

    // -----------------------------------------------------------------------
    // Thermodynamic state fields
    // -----------------------------------------------------------------------

    virtual volScalarField& p() = 0;

    virtual const volScalarField& TTR() const = 0;
    virtual volScalarField& TTR() = 0;
    virtual volScalarField& TVib() = 0;

    virtual const volScalarField& psi() const = 0;

    // h() is read/write: solver writes hTR before calling correct()
    virtual volScalarField& h() = 0;
    virtual volScalarField& eTR() = 0;      // new cell-centred eTR field
    virtual volScalarField& eVib() = 0;

    //- Electron pressure [Pa]
    //  Returns zero field for backends without free electrons.
    virtual tmp<volScalarField> pe() const = 0;

    virtual tmp<volScalarField> rho() const = 0;


    // -----------------------------------------------------------------------
    // Transport properties
    // -----------------------------------------------------------------------

    virtual tmp<volScalarField> mu() const = 0;

    virtual volScalarField fickDiffusionCoeff() const = 0;

    virtual tmp<volScalarField> kappaTR() const = 0;

    virtual tmp<volScalarField> kappaVib() const = 0;

    //- Translational-rotational thermal diffusivity for energy [kg/m/s]
    virtual tmp<volScalarField> alphaheTR() const = 0;

    // Called in kurganov.H: c = sqrt(gamma()/psi)
    virtual tmp<volScalarField> gamma() const = 0;


    // -----------------------------------------------------------------------
    // Lifecycle
    // -----------------------------------------------------------------------

    virtual void correct() = 0;

    virtual void correctVibEnergy() = 0;

    virtual void correctTVib() = 0;

    virtual PtrList<volScalarField>& correctVibSource
    (
        const PtrList<volScalarField>& TVibSpecies
    ) = 0;

    virtual PtrList<volScalarField>& correctVibVibSource
    (
        const PtrList<volScalarField>& TVibSpecies
    ) = 0;

    virtual PtrList<volScalarField>& correctCVSource
    (
        const PtrList<volScalarField::Internal>& RR
    ) = 0;


    // -----------------------------------------------------------------------
    // Configuration flags
    // -----------------------------------------------------------------------

    virtual bool twoTemperature() const = 0;

    virtual bool chemistry() const = 0;

    virtual bool chemistryCapable() const = 0;

    virtual const fileName& speciesListFile() const = 0;

    virtual const fileName& reactionsListFile() const = 0;

    virtual const word& chemistryReader() const = 0;

    //- True if backend couples electronic energy into the Ve mode (Mutation++ RRHO).
    //  False for neTCLib (no electronic energy in 2T model).
    virtual bool hasElectronicVibCoupling() const { return false; }

};

} // End namespace Foam

#endif
