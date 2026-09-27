/*---------------------------------------------------------------------------*\
    tcLibraryInterface.C

    Factory implementation for tcLibraryInterface::New().

    Reads 'library' keyword from constant/thermophysicalProperties:
      library  neTCLib;     -> constructs ne2TReactionThermo
      library  default;     -> constructs ne2TReactionThermo (alias)
      library  mutationpp;  -> constructs mutationppWrapper
                               (requires WITH_MUTATION_PP compile flag)
      (keyword absent)      -> constructs ne2TReactionThermo (safe default)
\*---------------------------------------------------------------------------*/

#include "tcLibraryInterface.H"
#include "ne2TReactionThermo.H"
#include "IOdictionary.H"
#include "Time.H"

// Mutation++ wrapper — included only when compiled with support
#ifdef WITH_MUTATION_PP
#include "mutationppWrapper.H"
#endif

defineTypeName(Foam::tcLibraryInterface);

namespace Foam
{

autoPtr<tcLibraryInterface> tcLibraryInterface::New(const fvMesh& mesh)
{
    

    IOdictionary dict
    (
        IOobject
        (
            "thermophysicalProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        )
    );

    const word libName =
        dict.lookupOrDefault<word>("library", "neTCLib");

    Info<< "tcLibraryInterface::New : library = " << libName << nl << endl;

    if (libName == "mutationpp")
    {
#ifdef WITH_MUTATION_PP
        // Mutation++ path: mutationppWrapper satisfies both
        // tcLibraryInterface and compositionInterface.
        return autoPtr<tcLibraryInterface>
        (
            new mutationppWrapper(mesh)
        );
#else
        FatalErrorInFunction
            << "library 'mutationpp' selected but nexaFoam was compiled "
            << "without Mutation++ support." << nl
            << "Recompile with:" << nl
            << "    export WITH_MUTATION_PP=1" << nl
            << "    ./Allwmake" << nl
            << "Or set 'library' to 'neTCLib' or 'default'."
            << exit(FatalError);

        return autoPtr<tcLibraryInterface>(nullptr);  // unreachable
#endif
    }
    else if
    (
        libName == "neTCLib"
     || libName == "default"
    )
    {
        autoPtr<ne2TReactionThermo> native =
            ne2TReactionThermo::New(mesh);

        return autoPtr<tcLibraryInterface>(native.release());
    }
    else
    {
        FatalErrorInFunction
            << "Unknown library value: '" << libName << "'." << nl
            << "Valid options: 'neTCLib', 'default', 'mutationpp'."
            << exit(FatalError);

        return autoPtr<tcLibraryInterface>(nullptr);  // unreachable
    }
}

} // End namespace Foam
