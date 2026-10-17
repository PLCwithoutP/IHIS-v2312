/*---------------------------------------------------------------------------*\
  IHIS - Implicit Hypersonics In-house Solver (OpenFOAM-v2312)
-------------------------------------------------------------------------------
Class
    Foam::functionObjects::neWallShearStress

Description
    Wall shear stress for IHIS / neTCLib (Phase 7, D4/D10).

    Reads the viscous wall traction the solver itself applied in
    residualsUpdate.H (neWall:tau = mu snGrad(U) + n.tauMC, the traction of
    the wall on the gas) and reports the stress the GAS exerts on the WALL:

        tau_w   = -neWall:tau                [Pa]
        tau_t   = (I - n n).tau_w            tangential part
        |tau_t|                              for the skin friction coefficient

    Output at write time:
      - volVectorField neWallShearStress (tau_w on the selected patches)
      - postProcessing/<name>/<time>/<patch>.csv:
            x,y,z,tau_x,tau_y,tau_z,tau_t
      - log: patch force integrals [N] (viscous only)

Usage
    \verbatim
    wallShearStress
    {
        type            neWallShearStress;
        libs            (ihisFunctionObjects);
        patches         (wall);         // optional, default: all wall patches
        writeControl    writeTime;
    }
    \endverbatim

SourceFiles
    neWallShearStress.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_neWallShearStress_H
#define functionObjects_neWallShearStress_H

#include "fvMeshFunctionObject.H"
#include "HashSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class neWallShearStress Declaration
\*---------------------------------------------------------------------------*/

class neWallShearStress
:
    public fvMeshFunctionObject
{
    // Private Data

        //- Selected patches
        labelHashSet patchSet_;


    // Private Member Functions

public:

    //- Runtime type information
    TypeName("neWallShearStress");


    // Constructors

        //- Construct from Time and dictionary
        neWallShearStress
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- No copy construct
        neWallShearStress(const neWallShearStress&) = delete;

        //- No copy assignment
        void operator=(const neWallShearStress&) = delete;


    //- Destructor
    virtual ~neWallShearStress() = default;


    // Member Functions

        //- Read the settings
        virtual bool read(const dictionary&);

        //- Nothing to do: the fluxes are computed by the solver
        virtual bool execute();

        //- Write the wall shear stress
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
