/*---------------------------------------------------------------------------*\

    IHIS: In-house Hypersonics Implicit solver

    Copyright (C) 2014-2018 Oliver Oxtoby - CSIR, South Africa
    Copyright (C) 2014-2018 Johan Heyns - CSIR, South Africa
    Copyright (C) 1991-2008 OpenCFD Ltd.

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
    Foam::solver

Description
    Base class for matrix solvers

SourceFiles
    solverNew.C

Authors
    Oliver Oxtoby
    Johan Heyns
        Council for Scientific and Industrial Research, South Africa

\*---------------------------------------------------------------------------*/

#ifndef ihisSolver_H
#define ihisSolver_H

#include "IOdictionary.H"
#include "volFields.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "jacobianMatrix.H"
#include "preconditioner.H"
#include "residualIO.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class solver Declaration
\*---------------------------------------------------------------------------*/

template <int nScalar, int nVector>
class ihisSolver
{

protected:
    // Protected data

        const jacobianMatrix<nScalar,nVector>& jacobian_;

        const fvMesh& mesh_;

        const preconditioner<nScalar,nVector>* preconditioner_;

        const dictionary& dict_;

        const residualIO defaultTol_;

private:
    // Private Member Functions

        //- Disallow copy construct
        ihisSolver(const ihisSolver&);

        //- Disallow default bitwise assignment
        void operator=(const ihisSolver&);

public:
    //- Runtime type information
    TypeName("ihisSolver");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            ihisSolver,
            dictionary,
            (
                const dictionary& dict,
                const jacobianMatrix<nScalar,nVector>& jacobian,
                const preconditioner<nScalar,nVector>* preconditioner,
                residualIO& defaultTol
            ),
            (dict, jacobian, preconditioner, defaultTol)
        );


	// Constructors

        //- Construct from components
        ihisSolver
        (
            const word& type,
            const dictionary& dict,
            const jacobianMatrix<nScalar,nVector>& jacobian,
            const preconditioner<nScalar,nVector>* preconditioner,
            residualIO& defaultTol
        )
         :
            jacobian_(jacobian),
            mesh_(jacobian.mesh()),
            preconditioner_(preconditioner),
            dict_(dict),
            defaultTol_(defaultTol)
        {}


    // Selectors

        //- Return a pointer to the selected solver
        static autoPtr<ihisSolver> New
        (
            const dictionary& dict,
            const jacobianMatrix<nScalar,nVector>& jacobian,
            const preconditioner<nScalar,nVector>* preconditioner,
            // Default residual tolerance and ordering of residual I/O:
            residualIO& defaultTol
        );


    // Destructor

        virtual ~ihisSolver()
        {}


    // Member Functions

        //- Solve
        virtual label solve
        (
            // Initial/returned solution variables:
            PtrList<volScalarField>& dsW, PtrList<volVectorField>& dvW,
            // Residuals:
            const PtrList<volScalarField>& sR, const PtrList<volVectorField>& vR,
            // Returns initial residual:
            autoPtr< residualIO >& pInitRes
        ) const = 0;

        //- Call preconditioner if there is one
        virtual void precondition
        (
            PtrList<volScalarField>& sVec,
            PtrList<volVectorField>& vVec
        ) const
        {
            if (preconditioner_)
            {
               preconditioner_->precondition(sVec, vVec);
            }
        }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeSolverType(SS, nScalar, nVector)                                   \
                                                                               \
typedef SS<nScalar, nVector> SS##nScalar##nVector;                             \
defineNamedTemplateTypeNameAndDebug(SS##nScalar##nVector, 0);                  \
                                                                               \
ihisSolver<nScalar, nVector>::adddictionaryConstructorToTable<SS<nScalar, nVector>> \
    add##SS##nScalar##nVector##ConstructorToTable_;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ihisSolverNew.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
