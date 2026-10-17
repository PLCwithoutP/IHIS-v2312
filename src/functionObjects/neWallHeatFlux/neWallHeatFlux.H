/*---------------------------------------------------------------------------*\
  IHIS - Implicit Hypersonics In-house Solver (OpenFOAM-v2312)
-------------------------------------------------------------------------------
Class
    Foam::functionObjects::neWallHeatFlux

Description
    Wall heat flux for IHIS / neTCLib (Phase 7, D4/D10).

    Reads the wall-face fluxes the solver itself applied in residualsUpdate.H
    (registered as neWall:*), so the result is exactly the energy that left
    the gas through the wall in the last outer iteration. Per wall face,
    positive INTO the wall, [W/m2]:

        q_cond = -neWall:qTR/|Sf|        trans-rotational conduction
        q_slip = -neWall:sigmaDotU/|Sf|  viscous work (non-zero only with slip)
        q_vib  = -neWall:qVib/|Sf|       vibrational conduction (2T)
        q_diff = +neWall:hJ/|Sf|         species enthalpy diffusion (Ns > 1)
        q_tot  = sum of the above

    Missing fields (1T, Ns = 1, inviscid) contribute zero.

    Output at write time:
      - volScalarField neWallHeatFlux (q_tot on the selected patches)
      - postProcessing/<name>/<time>/<patch>.csv:
            x,y,z,q_tot,q_cond,q_slip,q_vib,q_diff
      - log: patch integrals [W] of every component

Usage
    \verbatim
    wallHeatFlux
    {
        type            neWallHeatFlux;
        libs            (ihisFunctionObjects);
        patches         (wall);         // optional, default: all wall patches
        writeControl    writeTime;
    }
    \endverbatim

SourceFiles
    neWallHeatFlux.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_neWallHeatFlux_H
#define functionObjects_neWallHeatFlux_H

#include "fvMeshFunctionObject.H"
#include "HashSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class neWallHeatFlux Declaration
\*---------------------------------------------------------------------------*/

class neWallHeatFlux
:
    public fvMeshFunctionObject
{
    // Private Data

        //- Selected patches
        labelHashSet patchSet_;


    // Private Member Functions

        //- Boundary values of a registered surface field on patchi,
        //  zero if the field is not registered
        tmp<scalarField> patchFlux
        (
            const word& fieldName,
            const label patchi
        ) const;


public:

    //- Runtime type information
    TypeName("neWallHeatFlux");


    // Constructors

        //- Construct from Time and dictionary
        neWallHeatFlux
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- No copy construct
        neWallHeatFlux(const neWallHeatFlux&) = delete;

        //- No copy assignment
        void operator=(const neWallHeatFlux&) = delete;


    //- Destructor
    virtual ~neWallHeatFlux() = default;


    // Member Functions

        //- Read the settings
        virtual bool read(const dictionary&);

        //- Nothing to do: the fluxes are computed by the solver
        virtual bool execute();

        //- Write the wall heat flux
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
