/*---------------------------------------------------------------------------*\

    IHIS: In-house Hypersonics Implicit solver

    Copyright (C) 2014-2018 Johan Heyns - CSIR, South Africa
    Copyright (C) 2014-2018 Oliver Oxtoby - CSIR, South Africa
    Copyright (C) 2020-2022 Oliver Oxtoby

-------------------------------------------------------------------------------
License
    This file is part of IHIS.

    IHIS is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    IHIS is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with IHIS.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::ihisModule

Description
    Coupled, density-based compressible flow solver

See also
    rhoCentralFoam

SourceFiles
    ihisModule.C

Authors
    Johan Heyns
    Oliver Oxtoby
        Council for Scientific and Industrial Research, South Africa

\*---------------------------------------------------------------------------*/

#ifndef ihisModule_H
#define ihisModule_H

#include "solverModule.H"
#include "pimpleControl.H"
#include "ne2TReactionThermo.H"
#include "chemistryInterface.H"
#include "zeroGradientFvPatchFields.H"
#include "fixedValueFvPatchFields.H"
#include "fluxScheme.H"
#include "pseudotimeControl/pseudotimeControl.H"
#include "jacobian.H"
#include "upwind.H"
#include "preconditioner.H"
#include "faceSet.H"
#include "fvOptions.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class ihisModule Declaration
\*---------------------------------------------------------------------------*/

class ihisModule
:
    public solverModule
{

protected:

    // Ordering of residual printed output
    const label residualOrderingArray[3] = {0, 2, 1};
    const labelList residualOrdering =
        labelList(FixedList<label,3>(residualOrderingArray));

    // Private data

        //- Reference to the time database
        const Time& time_;

        autoPtr<fvMesh> mesh_;

        autoPtr< pseudotimeControl > solnControl_;

        autoPtr<ne2TReactionThermo> pThermo_;

        autoPtr<volVectorField> U_;

        // Stored (if viscous) for use in maxwellSlip boundary condition
        autoPtr<volTensorField> tauMC_;

        autoPtr<surfaceScalarField> phi_;
        autoPtr<surfaceVectorField> phiUp_;
        autoPtr<surfaceScalarField> phiEp_;
        autoPtr<surfaceVectorField> Up_;

        // The value of neThermo.eTR() at T = 0 and p = 0, used as offset to get
        // absolute internal energy
        autoPtr<volScalarField> eZero_;
        autoPtr<volScalarField> eVibMix_;

        //- Two-temperature active (= neThermo.twoTemperature(), set in createFields.H)
        bool twoT_ = false;

        //- C-V settings (Phase 6, D4), read when twoT_ && chemistry_
        bool cvPreferential_ = false;
        bool vvOn_ = true;
        scalar alphaCV_ = 0.3;
        scalarList DVib_;     // dissociation potential per molecule [J/kg]
        
        //- Species index of molecule m, m = 0..Nm-1 (from isSpecieMolecular)
        labelList molI_;

        //- Unknowns eps_m = Y_m*e_v,m [J/kg] (Nm entries)
        PtrList<volScalarField> epsVib_;

        //- Vibrational energy per unit mass of molecule m, e_v,m [J/kg] (Nm, derived)
        PtrList<volScalarField> evVib_;

        //- Vibrational temperature per species [K] (Ns entries; atoms hold T_TR)
        PtrList<volScalarField> TVibSpecies_;

        //- (kappaVib_m/CvVib_m)/phi_m from kappaVibWilke [kg/m/s] (Nm),
        //  computed once per outer iteration, used by residualsUpdate.H and vibEqns.H
        PtrList<volScalarField> GVib_;

        //- Species diffusive mass fluxes [kg/s] (Ns), shared by rhoE and vib equations
        PtrList<surfaceScalarField> phiJ_;

        //- Vibrational conduction face fluxes Gamma_e,m*snGrad(e_v,m)*|Sf| [W] (Nm)
        PtrList<surfaceScalarField> qVib_;

        //- Molecule face mass fractions from the species loop, state k (Nm)
        PtrList<surfaceScalarField> YfMol_;

        //- Finite-rate chemistry. Constructed only when
        //  thermophysicalProperties 'chemistry' is true and Ns > 1.
        autoPtr<chemistryInterface> chemistry_;

        //- Chemical Jacobian diagonal, filled each outer
        //  iteration by chemistry_->calculateWithJacobianDiag()
        PtrList<volScalarField::Internal> dRRdRhoY_;
        PtrList<volScalarField::Internal> dRRdT_;

        autoPtr<volScalarField> rPseudoDeltaT_;
        // 1/dt (for Euler), or equivalent for other schemes
        autoPtr<scalarField> ddtCoeff_;
        // Whether variables had to be bounded this iteration
        bool bounded_;

        //rho and rhoE variables
        PtrList<volScalarField> scalarVars_;
        PtrList<volScalarField> scalarVarsPrevIter_;
        //rhoU variable
        PtrList<volVectorField> vectorVars_;
        PtrList<volVectorField> vectorVarsPrevIter_;

        //rho and rhoE residuals
        PtrList<volScalarField> scalarResiduals_;
        //rhoU residual
        PtrList<volVectorField> vectorResiduals_;

        autoPtr<fluxScheme> flux_;

        bool inviscid_;

        //- How frequently mesh is moved during outer corrector loop.
        //  0 means do not move (only at beginning of timestep).
        //  Only applies to transient analysis.
        label moveMeshOuterCorrectors_;

        bool steadyState_;
        bool localTimestepping_;
        bool localTimesteppingBounding_;
        scalar localTimesteppingLowerBound_;
        scalar localTimesteppingUpperBound_;
        scalar boundingDecrease_;
        scalar boundingIncrease_;

        //- Store previous residual for CoNum relaxation
        autoPtr<residualIO> res_;
        autoPtr<residualIO> prevRes_;

        //- Current pseudo-courant number
        autoPtr<uniformDimensionedScalarField> pseudoCoNum_; // global value
        autoPtr<volScalarField> pseudoCoField_;   // locally damped value
        autoPtr<volScalarField> pseudoCoDamping_; // local damping depending on relative change in variables
        //- Initial pseudo-Courant number
        scalar pseudoCoNumInit_;
        //- Minimum pseudo-Courant number
        scalar pseudoCoNumMin_;
        //- Maximum pseudo-Courant number
        scalar pseudoCoNumMax_;
        //- Maximum increase as a fraction of current pseudo-Courant number
        scalar pseudoCoNumMaxIncr_;
        //- Minimum decrease as a fraction of current pseudo-Courant number
        scalar pseudoCoNumMinDecr_;

        //- Print debug info for nearest cell
        bool cellDebugging_;
        label debugCell_;

        fv::options* fvOptions_;


    // Protected member functions

        //- Set timestep size based on current Courant no.
        void setPseudoDeltaT();

        //- Adjust Courant number
        void setPseudoCoNum();

        //- Recursively create preconditioners and their Jacobians if applicable
        void createPreconditioners
        (
            PtrList<preconditioner<2,1>>& preconditioners,
            PtrList<jacobian>& jacobians,
            const dictionary& parentDict
        );

        void redistributePar();

        //- Search for debug point specified in controlDict
        void findDebugCell();

        //- Recompute eZero_ = -hc(); aborts if non-finite.
        //  Called at initialisation and before every decode (energy-neutral chemistry in Phase 3).
        void updateEZero();

        //- Vibrational temperature floor: fvSolution pseudoTime TMin (same as T_TR)
        scalar vibTvMin() const;

        //- Start-up: Tv_m -> e_v,m, eps_m (cells), then vibBoundaries()
        void vibEncode();

        //- After the vibrational solve: eps_m -> e_v,m, Tv_m (cells), then
        //  vibBoundaries(). Returns the number of cells where a Tv_m hit the floor.
        label vibDecode();

        //- Faces from TVib_<m> BCs; atoms <- T_TR; mixture TVib; eVibMix_ = sum eps_m
        void vibBoundaries();

        //- Start-up checks (permanent): closed-form inverse matches the library's
        //  EsVib (FatalError otherwise); report whether kappaVib/CvVib is Eucken-type
        void vibCheckModel() const;
public:

    //- Runtime type information
    TypeName("ihis");


    // Constructors

        //- Construct from components
        ihisModule
        (
            const word& name,
            const Time&
        );

        //- Construct from components
        ihisModule
        (
            const word& name,
            const Time&,
            const dictionary&
        );

    // Member Functions

        //- Returns ratio between desired and current time step
        virtual scalar timeStepScaling(const scalar& maxCoNum);

        //- Called before time loop. Typically create mesh and fields here.
        virtual void initialise();

        //- Called at the start of a time step
        virtual void beginTimeStep();

        //- Called for each outer ('pimple') iteration
        virtual void outerIteration();

        //- Should return the mesh created during initialise
        virtual fvMesh& mesh()
        {
            return mesh_();
        }

        //- Should return the solution control created during initialise
        virtual solutionControl& solnControl()
        {
            return solnControl_();
        }

        //- Returns whether analysis should be steady state - makes time steps
        // and outer iterations one and the same.
        virtual bool steadyState()
        {
            return steadyState_;
        }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

