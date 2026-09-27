/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2021 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "neReversibleReaction.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template
<
    template<class> class ReactionType,
    class ReactionThermo,
    class ReactionRate
>
void Foam::neReversibleReaction<ReactionType, ReactionThermo, ReactionRate>::
initKeq()
{
    // Change in moles across the reaction
    dn_ = 0.0;

    for (const auto& sc : this->rhs())
    {
        dn_ += sc.stoichCoeff;
    }

    for (const auto& sc : this->lhs())
    {
        dn_ -= sc.stoichCoeff;
    }

    // Park's fits give Keq on a mol/cm^3 basis; OpenFOAM concentrations are
    // kmol/m^3. 1 mol/cm^3 = 1e3 kmol/m^3, hence a factor 1000^dn.
    lnKeqUnitFactor_ = dn_*Foam::log(1000.0);

    // Continuation constant: match the two branches in value at KeqTmin,
    //   fit(Tmin) = alpha - theta/Tmin
    const scalar Tm = max(KeqTmin_, SMALL);
    const scalar Zm = 10000.0/Tm;

    const scalar fitAtTmin =
        KeqCoeffs_[0]
      + KeqCoeffs_[1]*Zm
      + KeqCoeffs_[2]*Zm*Zm
      + KeqCoeffs_[3]*Zm*Zm*Zm
      + KeqCoeffs_[4]*Zm*Zm*Zm*Zm;

    KeqAlpha_ = fitAtTmin + KeqTheta_/Tm;
}


template
<
    template<class> class ReactionType,
    class ReactionThermo,
    class ReactionRate
>
void Foam::neReversibleReaction<ReactionType, ReactionThermo, ReactionRate>::
reportKeq() const
{
    // Printed once per reaction at construction. This exists so that a
    // wrong or defaulted equilibrium constant is visible in the log instead
    // of showing up months later as a wrong wall heat flux.
    const scalarField cDummy(0);

    const scalar lnkb300  = k_.logk(0.0, 300.0,  cDummy) - lnKeq(300.0);
    const scalar lnkb4000 = k_.logk(0.0, 4000.0, cDummy) - lnKeq(4000.0);

    Info<< "    " << this->name()
        << ": dn = " << dn_
        << ", lnKeq(300) = " << lnKeq(300.0)
        << ", lnKeq(4000) = " << lnKeq(4000.0)
        << ", ln kb(300) = " << lnkb300
        << ", ln kb(4000) = " << lnkb4000 << endl;
}


template
<
    template<class> class ReactionType,
    class ReactionThermo,
    class ReactionRate
>
inline Foam::scalar
Foam::neReversibleReaction<ReactionType, ReactionThermo, ReactionRate>::lnKeq
(
    const scalar T
) const
{
    const scalar Teff = max(T, SMALL);

    if (Teff >= KeqTmin_)
    {
        const scalar Z = 10000.0/Teff;

        // Quartic in Z. See the class header before changing this.
        return
            lnKeqUnitFactor_
          + KeqCoeffs_[0]
          + KeqCoeffs_[1]*Z
          + KeqCoeffs_[2]*Z*Z
          + KeqCoeffs_[3]*Z*Z*Z
          + KeqCoeffs_[4]*Z*Z*Z*Z;
    }

    // Low-temperature asymptote. The quartic diverges below ~1500 K (the Z^4
    // term dominates and flips sign); ln Keq -> alpha - theta/T is the correct
    // leading behaviour and is continuous with the fit at KeqTmin by
    // construction of KeqAlpha_.
    return lnKeqUnitFactor_ + KeqAlpha_ - KeqTheta_/Teff;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template
<
    template<class> class ReactionType,
    class ReactionThermo,
    class ReactionRate
>
Foam::neReversibleReaction<ReactionType, ReactionThermo, ReactionRate>::
neReversibleReaction
(
    const ReactionType<ReactionThermo>& reaction,
    const ReactionRate& k
)
:
    ReactionType<ReactionThermo>(reaction),
    k_(k),
    alphaPark_(-1.0),
    twoTemperature_(false),
    KeqCoeffs_(Zero),
    KeqTheta_(0.0),
    KeqTmin_(2000.0),
    KeqAlpha_(0.0),
    dn_(0.0),
    lnKeqUnitFactor_(0.0)
{
    initKeq();
}


template
<
    template<class> class ReactionType,
    class ReactionThermo,
    class ReactionRate
>
Foam::neReversibleReaction<ReactionType, ReactionThermo, ReactionRate>::
neReversibleReaction
(
    const speciesTable& species,
    const ReactionTable<ReactionThermo>& thermoDatabase,
    const dictionary& dict
)
:
    ReactionType<ReactionThermo>(species, thermoDatabase, dict),
    k_(species, dict),
    alphaPark_(dict.getOrDefault<scalar>("alphaPark", -1.0)),
    twoTemperature_(dict.getOrDefault<bool>("twoTemperature", false)),
    KeqCoeffs_(Zero),
    KeqTheta_(0.0),
    KeqTmin_(2000.0),
    KeqAlpha_(0.0),
    dn_(0.0),
    lnKeqUnitFactor_(0.0)
{
    // Mandatory. Never getOrDefault: a silently defaulted equilibrium
    // constant produces a plausible-looking but irreversible mechanism.
    dict.readEntry("KeqCoeffs", KeqCoeffs_);
    dict.readEntry("KeqTheta", KeqTheta_);

    KeqTmin_ = dict.getOrDefault<scalar>("KeqTmin", 2000.0);

    if (KeqTmin_ < SMALL)
    {
        FatalIOErrorInFunction(dict)
            << "KeqTmin must be positive for reaction " << this->name()
            << exit(FatalIOError);
    }

    if (KeqTheta_ < 0.0)
    {
        FatalIOErrorInFunction(dict)
            << "KeqTheta must be positive (reaction energy / R, in K)"
            << " for reaction " << this->name() << nl
            << "Got " << KeqTheta_
            << exit(FatalIOError);
    }

    initKeq();

    // Consistency check. The change in moles is available two independent
    // ways: from the stoichiometric coefficients, and from the reaction
    // thermo assembled by neReaction::setThermo() as Y()/W().
    const scalar nmThermo = this->Y()/this->W();

    if (mag(nmThermo - dn_) > 1e-6*max(mag(dn_), 1.0))
    {
        WarningInFunction
            << "Reaction " << this->name() << ": change in moles from"
            << " stoichiometry (" << dn_ << ") disagrees with the"
            << " reaction thermo (" << nmThermo << ")." << nl
            << "The equilibrium constant will be inconsistent." << endl;
    }

    reportKeq();
}


template
<
    template<class> class ReactionType,
    class ReactionThermo,
    class ReactionRate
>
Foam::neReversibleReaction<ReactionType, ReactionThermo, ReactionRate>::
neReversibleReaction
(
    const neReversibleReaction<ReactionType, ReactionThermo, ReactionRate>& rr,
    const speciesTable& species
)
:
    ReactionType<ReactionThermo>(rr, species),
    k_(rr.k_),
    alphaPark_(rr.alphaPark_),
    twoTemperature_(rr.twoTemperature_),
    KeqCoeffs_(rr.KeqCoeffs_),
    KeqTheta_(rr.KeqTheta_),
    KeqTmin_(rr.KeqTmin_),
    KeqAlpha_(rr.KeqAlpha_),
    dn_(rr.dn_),
    lnKeqUnitFactor_(rr.lnKeqUnitFactor_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template
<
    template<class> class ReactionType,
    class ReactionThermo,
    class ReactionRate
>
Foam::scalar
Foam::neReversibleReaction<ReactionType, ReactionThermo, ReactionRate>::kf
(
    const scalar p,
    const scalar T,
    const scalarField& c
) const
{
    return k_(p, T, c);
}


template
<
    template<class> class ReactionType,
    class ReactionThermo,
    class ReactionRate
>
Foam::scalar
Foam::neReversibleReaction<ReactionType, ReactionThermo, ReactionRate>::kr
(
    const scalar p,
    const scalar T,
    const scalarField& c
) const
{
    // Log space throughout. At a 300 K wall a dissociation kf is O(1e-149)
    // and Keq is O(1e-164); forming either explicitly and dividing invites
    // underflow, and any floor placed on Keq silently destroys recombination.
    const scalar lnkb = k_.logk(p, T, c) - lnKeq(T);

    return Foam::exp(min(lnkb, 300.0));
}


template
<
    template<class> class ReactionType,
    class ReactionThermo,
    class ReactionRate
>
Foam::scalar
Foam::neReversibleReaction<ReactionType, ReactionThermo, ReactionRate>::kr
(
    const scalar kfwd,
    const scalar p,
    const scalar T,
    const scalarField& c
) const
{
    return kr(p, T, c);   // kfwd intentionally unused, see header
}


template
<
    template<class> class ReactionType,
    class ReactionThermo,
    class ReactionRate
>
void
Foam::neReversibleReaction<ReactionType, ReactionThermo, ReactionRate>::write
(
    Ostream& os
) const
{
    ReactionType<ReactionThermo>::write(os);
    k_.write(os);

    if (alphaPark_ > 0)
    {
        os.writeEntry("alphaPark", alphaPark_);
    }

    os.writeEntry("twoTemperature", twoTemperature_);
    os.writeEntry("KeqCoeffs", KeqCoeffs_);
    os.writeEntry("KeqTheta", KeqTheta_);
    os.writeEntry("KeqTmin", KeqTmin_);
}


// ************************************************************************* //