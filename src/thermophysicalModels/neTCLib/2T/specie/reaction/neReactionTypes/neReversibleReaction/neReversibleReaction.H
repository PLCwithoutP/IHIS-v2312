/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019-2021 OpenCFD Ltd.
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
    Foam::neReversibleReaction

Group
    grpSpecieReactions

Description
    Extension of neReaction handling reversible reactions, with the backward
    rate obtained at runtime from the equilibrium constant rather than from a
    separately fitted Arrhenius expression:

        k_b(T) = k_f(T) / Keq(T)

    evaluated in log space, since at low temperature k_f and Keq individually
    underflow while their ratio remains well scaled.

    Two-temperature note
    --------------------
    The forward rate is evaluated at whatever temperature the chemistry model
    supplies -- for a Park two-temperature dissociation that is the
    rate-controlling temperature Tf = T^alpha * Tv^(1-alpha). Detailed balance,
    however, holds at the translational-rotational temperature only, so the
    backward rate must be evaluated at TTR. The chemistry model is responsible
    for passing the correct temperature to each of kf() and kr().

    Equilibrium constant
    --------------------
    Keq is read as the five-coefficient curve fit

        ln Keq = A1 + A2*Z + A3*Z^2 + A4*Z^3 + A5*Z^4,      Z = 10000/T

    with Keq in mol/cm^3, converted internally to (kmol/m^3)^dn.

    NOTE the form. This is a plain quartic in Z. It is NOT
        A1/Z + A2 + A3*ln(Z) + A4*Z + A5*Z^2
    which is a different (and for these coefficients, wrong) correlation --
    using it overestimates Keq by e^42 at 2000 K and e^366 at 300 K, which
    silently reduces the mechanism to an irreversible one.

    Verified against statistical mechanics for all five Park reactions over
    2000-20000 K, worst deviation 0.58 in ln Keq. See validate_Keq.py.

    Low-temperature branch
    ----------------------
    The quartic is a high-temperature correlation and diverges violently below
    roughly 1500 K, where the Z^4 term dominates and changes sign. Below
    KeqTmin the fit is replaced by the physically correct asymptote

        ln Keq = alpha - theta/T

    where theta is the reaction energy divided by R (dictionary entry
    KeqTheta) and alpha is set at construction so the two branches agree in
    value at KeqTmin. Accuracy of the continuation at 300 K is better than a
    factor of 5 for all five reactions.

    Required dictionary entries
    ---------------------------
        A, beta, Ta         forward Arrhenius, A in m^3/(kmol s)
        KeqCoeffs           (A1 A2 A3 A4 A5)   [mandatory]
        KeqTheta            reaction energy / R, in K   [mandatory]
        KeqTmin             optional, default 2000 K
        alphaPark           optional, forward only
        twoTemperature      optional, default false. MUST be true for
                            alphaPark to have any effect.

    KeqCoeffs and KeqTheta are read with readEntry (fatal if absent) rather
    than getOrDefault, deliberately: a silently defaulted equilibrium constant
    is indistinguishable from a working one until the heat flux is wrong.

SourceFiles
    neReversibleReaction.C

\*---------------------------------------------------------------------------*/

#ifndef neReversibleReaction_H
#define neReversibleReaction_H

#include "neReaction.H"
#include "FixedList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class neReversibleReaction Declaration
\*---------------------------------------------------------------------------*/

template
<
    template<class> class ReactionType,
    class ReactionThermo,
    class ReactionRate
>
class neReversibleReaction
:
    public ReactionType<ReactionThermo>
{
    // Private Data
    //  NOTE: C++ initialises members in declaration order regardless of the
    //  order written in the initialiser list. Keep the two in step.

        //- Forward reaction rate
        ReactionRate    k_;

        //- Park rate-controlling temperature exponent. FORWARD ONLY.
        scalar          alphaPark_;

        //- Two-temperature switch for the forward rate
        bool            twoTemperature_;

        //- A1..A5 of the Keq curve fit (Park 1985, Table 1)
        FixedList<scalar, 5> KeqCoeffs_;

        //- Reaction energy / R [K], used for the low-T continuation
        scalar          KeqTheta_;

        //- Temperature below which the fit is replaced by alpha - theta/T
        scalar          KeqTmin_;

        //- Continuation constant, set at construction so branches match
        scalar          KeqAlpha_;

        //- Net change in moles, dn = sum(rhs) - sum(lhs)
        scalar          dn_;

        //- ln of the (mol/cm^3)^dn -> (kmol/m^3)^dn factor, = dn*ln(1000)
        scalar          lnKeqUnitFactor_;


    // Private Member Functions

        //- Set dn_, lnKeqUnitFactor_ and KeqAlpha_. Called by all constructors
        //- that do not simply copy an already-initialised object.
        void initKeq();

        //- Report the constructed Keq for this reaction to the log
        void reportKeq() const;

        //- ln of the equilibrium constant on a concentration basis,
        //- units (kmol/m^3)^dn.
        //  T MUST be the translational-rotational temperature.
        inline scalar lnKeq(const scalar T) const;

        //- No copy assignment
        void operator=
        (
            const neReversibleReaction
            <ReactionType, ReactionThermo, ReactionRate>&
        ) = delete;


public:

    //- Runtime type information
    TypeName("reversible");


    // Constructors

        //- Construct from components
        neReversibleReaction
        (
            const ReactionType<ReactionThermo>& reaction,
            const ReactionRate& reactionRate
        );

        //- Construct as copy given new speciesTable
        neReversibleReaction
        (
            const neReversibleReaction
            <ReactionType, ReactionThermo, ReactionRate>&,
            const speciesTable& species
        );

        //- Construct from dictionary
        neReversibleReaction
        (
            const speciesTable& species,
            const ReactionTable<ReactionThermo>& thermoDatabase,
            const dictionary& dict
        );

        //- Construct and return a clone
        virtual autoPtr<neReaction<ReactionThermo>> clone() const
        {
            return autoPtr<neReaction<ReactionThermo>>
            (
                new neReversibleReaction
                <ReactionType, ReactionThermo, ReactionRate>(*this)
            );
        }

        //- Construct and return a clone with new speciesTable
        virtual autoPtr<neReaction<ReactionThermo>> clone
        (
            const speciesTable& species
        ) const
        {
            return autoPtr<neReaction<ReactionThermo>>
            (
                new neReversibleReaction
                <ReactionType, ReactionThermo, ReactionRate>(*this, species)
            );
        }


    //- Destructor
    virtual ~neReversibleReaction() = default;


    // Member Functions

        //- Forward rate constant.
        //  T is whatever the chemistry model supplies: for a 2T dissociation
        //  reaction that is the Park rate-controlling temperature Tf.
        virtual scalar kf
        (
            const scalar p,
            const scalar T,
            const scalarField& c
        ) const;

        //- Backward rate constant, kb = kf(TTR)/Keq(TTR).
        //  T MUST be TTR. Never Tf.
        virtual scalar kr
        (
            const scalar p,
            const scalar T,
            const scalarField& c
        ) const;

        //- Backward rate constant given a forward rate.
        //  kfwd is DELIBERATELY IGNORED: in a two-temperature mechanism the
        //  caller's kfwd was evaluated at Tf, but detailed balance holds at
        //  TTR only, so kf is recomputed here at the supplied T.
        //  Do not "optimise" this to kfwd/Keq.
        virtual scalar kr
        (
            const scalar kfwd,
            const scalar p,
            const scalar T,
            const scalarField& c
        ) const;

        //- Write
        virtual void write(Ostream&) const;

        virtual scalar alphaPark() const { return alphaPark_; }

        virtual bool twoTemperature() const { return twoTemperature_; }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "neReversibleReaction.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //