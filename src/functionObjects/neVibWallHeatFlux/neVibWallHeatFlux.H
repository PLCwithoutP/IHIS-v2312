/*---------------------------------------------------------------------------*\
  IHIS - Implicit Hypersonics In-house Solver (OpenFOAM-v2312)
-------------------------------------------------------------------------------
Class
    Foam::functionObjects::neVibWallHeatFlux

Description
    Vibrational wall heat flux for IHIS / neTCLib (Phase 7), per molecule.
    Counterpart of nexaFoam's neVibWallHeatFlux.

    Reads the per-molecule vibrational conduction fluxes the solver applied
    in residualsUpdate.H (registered as neWall:qVib_<m>). Per wall face,
    positive INTO the wall, [W/m2]:

        q_vib,m = -neWall:qVib_<m>/|Sf|,     q_vib = sum_m q_vib,m

    q_vib equals the q_vib column of neWallHeatFlux. 1T runs register no
    neWall:qVib_<m>: a warning, nothing written.

    Output at write time:
      - volScalarField neVibWallHeatFlux (q_vib on the selected patches)
      - postProcessing/<name>/<time>/<patch>.csv:
            x,y,z,q_vib,q_vib_<m1>,q_vib_<m2>,...
      - log: patch integrals [W], total and per molecule

Usage
    \verbatim
    vibWallHeatFlux
    {
        type            neVibWallHeatFlux;
        libs            (ihisFunctionObjects);
        patches         (wall);         // optional, default: all wall patches
        writeControl    writeTime;
    }
    \endverbatim

SourceFiles
    neVibWallHeatFlux.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_neVibWallHeatFlux_H
#define functionObjects_neVibWallHeatFlux_H

#include "fvMeshFunctionObject.H"
#include "HashSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                      Class neVibWallHeatFlux Declaration
\*---------------------------------------------------------------------------*/

class neVibWallHeatFlux
:
    public fvMeshFunctionObject
{
    // Private Data

        //- Selected patches
        labelHashSet patchSet_;


public:

    //- Runtime type information
    TypeName("neVibWallHeatFlux");


    // Constructors

        //- Construct from Time and dictionary
        neVibWallHeatFlux
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- No copy construct
        neVibWallHeatFlux(const neVibWallHeatFlux&) = delete;

        //- No copy assignment
        void operator=(const neVibWallHeatFlux&) = delete;


    //- Destructor
    virtual ~neVibWallHeatFlux() = default;


    // Member Functions

        //- Read the settings
        virtual bool read(const dictionary&);

        //- Nothing to do: the fluxes are computed by the solver
        virtual bool execute();

        //- Write the vibrational wall heat flux
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
