/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::BasicNeChemistryModel

Description
    Basic chemistry model templated on thermodynamics

SourceFiles
    BasicNeChemistryModelI.H
    BasicNeChemistryModel.C

\*---------------------------------------------------------------------------*/

#ifndef BasicNeChemistryModel_H
#define BasicNeChemistryModel_H

#include "basicNeChemistryModel.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class fvMesh;

/*---------------------------------------------------------------------------*\
                     class BasicNeChemistryModel Declaration
\*---------------------------------------------------------------------------*/

template<class ReactionThermo>
class BasicNeChemistryModel
:
    public basicNeChemistryModel
{
    // Private Member Functions

        //- Construct as copy (not implemented)
        BasicNeChemistryModel
        (
            const BasicNeChemistryModel<ReactionThermo>&
        ) = delete;

        //- No copy assignment
        void operator=(const BasicNeChemistryModel<ReactionThermo>&) = delete;


protected:

    // Protected data

        //- Thermo
        ReactionThermo& thermo_;


public:

    //- Runtime type information
    TypeName("BasicNeChemistryModel");


    //- Thermo type
    typedef ReactionThermo ne2TReactionThermo;


    //- Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        BasicNeChemistryModel,
        thermo,
        (ReactionThermo& thermo),
        (thermo)
    );


    // Constructors

        //- Construct from thermo
        BasicNeChemistryModel(ReactionThermo& thermo);


    //- Selector
    static autoPtr<BasicNeChemistryModel<ReactionThermo>> New
    (
        ReactionThermo& thermo
    );


    //- Destructor
    virtual ~BasicNeChemistryModel();

    // Bring inherited RR(label) overloads into scope
    using basicNeChemistryModel::RR;
    
    // Member Functions

        //- Return access to the thermo package
        inline ReactionThermo& thermo();

        //- Return const access to the thermo package
        inline const ReactionThermo& thermo() const;

        virtual const PtrList<volScalarField::Internal>& RR() const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "BasicNeChemistryModel.C"
#endif

#include "BasicNeChemistryModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
