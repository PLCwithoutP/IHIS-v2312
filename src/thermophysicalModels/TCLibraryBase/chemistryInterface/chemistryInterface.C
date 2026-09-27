/*---------------------------------------------------------------------------*\
    chemistryInterface.C

    Factory for chemistryInterface::New().

    Selection logic:
      1. If thermo is a mutationppWrapper (WITH_MUTATION_PP compiled in),
         construct mppChemistryAdapter.
      2. Otherwise dynamic_cast to ne2TReactionThermo and construct
         neTCChemistryAdapter.

    The dynamic_cast lives here and nowhere else.
    The solver never sees concrete types — only chemistryInterface*.
\*---------------------------------------------------------------------------*/

#include "chemistryInterface.H"
#include "neTCChemistryAdapter.H"
#include "tcLibraryInterface.H"
#include "ne2TReactionThermo.H"

#ifdef WITH_MUTATION_PP
#include "mppChemistryAdapter.H"
#include "mutationppWrapper.H"
#endif

namespace Foam
{

autoPtr<chemistryInterface> chemistryInterface::New(tcLibraryInterface& thermo)
{
#ifdef WITH_MUTATION_PP
    // ── Mutation++ path ───────────────────────────────────────────────────────
    // Detect wrapper type via dynamic_cast.
    // dynamic_cast returns nullptr (not exception) for non-matching types,
    // so this is safe even when library = neTCLib.

    mutationppWrapper* mppPtr =
        dynamic_cast<mutationppWrapper*>(&thermo);

    if (mppPtr != nullptr)
    {
        Info<< "chemistryInterface::New : using mppChemistryAdapter" << nl
            << endl;

        return autoPtr<chemistryInterface>
        (
            new mppChemistryAdapter(*mppPtr)
        );
    }
#endif

    // ── neTCLib path (default) ────────────────────────────────────────────────
    ne2TReactionThermo& nativeThermo =
        dynamic_cast<ne2TReactionThermo&>(thermo);

    Info<< "chemistryInterface::New : using neTCChemistryAdapter" << nl
        << endl;

    return autoPtr<chemistryInterface>
    (
        new neTCChemistryAdapter(nativeThermo)
    );
}

} // End namespace Foam
