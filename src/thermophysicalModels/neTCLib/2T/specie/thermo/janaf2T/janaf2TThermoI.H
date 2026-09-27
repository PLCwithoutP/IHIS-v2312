/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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

#include "janaf2TThermo.H"
#include "specie2T.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::janaf2TThermo<EquationOfState>::janaf2TThermo
(
    const EquationOfState& st,
    const scalar TAbsLow,
    const scalar TAbsHigh,
    const scalar TFirstInt,
    const scalar TSecondInt,
    const scalar ThetaVib,
    const typename janaf2TThermo<EquationOfState>::coeffArray& highCpCoeffs,
    const typename janaf2TThermo<EquationOfState>::coeffArray& midCpCoeffs,
    const typename janaf2TThermo<EquationOfState>::coeffArray& lowCpCoeffs,
    const bool convertCoeffs,
    const bool warnOutOfRange
)
:
    EquationOfState(st),
    TAbsLow_(TAbsLow),
    TAbsHigh_(TAbsHigh),
    TFirstInt_(TFirstInt),
    TSecondInt_(TSecondInt),
    ThetaVib_(ThetaVib),
    convertCoeffs_(convertCoeffs),
    warnOutOfRange_(warnOutOfRange)
{
    if (convertCoeffs_)
    {
        for (label coefLabel=0; coefLabel<nCoeffs_; coefLabel++)
        {
            highCpCoeffs_[coefLabel] = highCpCoeffs[coefLabel];
            midCpCoeffs_[coefLabel] = midCpCoeffs[coefLabel];
            lowCpCoeffs_[coefLabel] = lowCpCoeffs[coefLabel];//*this->R()*1e-3;
        }
    }
    else
    {
        for (label coefLabel=0; coefLabel<nCoeffs_; coefLabel++)
        {
            highCpCoeffs_[coefLabel] = highCpCoeffs[coefLabel];
            midCpCoeffs_[coefLabel] = midCpCoeffs[coefLabel];
            lowCpCoeffs_[coefLabel] = lowCpCoeffs[coefLabel];
        }
    }
}


template<class EquationOfState>
inline const typename Foam::janaf2TThermo<EquationOfState>::coeffArray&
Foam::janaf2TThermo<EquationOfState>::coeffs
(
    const scalar TTR
) const
{
    if (TTR < TFirstInt_)
    {
        return lowCpCoeffs_;
    }
    else if ((TTR >= TFirstInt_) && (TTR < TSecondInt_))
    {
        return midCpCoeffs_;
    }
    else
    {
        return highCpCoeffs_;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::janaf2TThermo<EquationOfState>::janaf2TThermo
(
    const word& name,
    const janaf2TThermo& jt
)
:
    EquationOfState(name, jt),
    TAbsLow_(jt.TAbsLow_),
    TAbsHigh_(jt.TAbsHigh_),
    TFirstInt_(jt.TFirstInt_),
    TSecondInt_(jt.TSecondInt_),
    ThetaVib_(jt.ThetaVib_),
    warnOutOfRange_(jt.warnOutOfRange_)
{
    for (label coefLabel=0; coefLabel<nCoeffs_; coefLabel++)
    {
        highCpCoeffs_[coefLabel] = jt.highCpCoeffs_[coefLabel];
        midCpCoeffs_[coefLabel] = jt.midCpCoeffs_[coefLabel];
        lowCpCoeffs_[coefLabel] = jt.lowCpCoeffs_[coefLabel];
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::scalar Foam::janaf2TThermo<EquationOfState>::limit
(
    const scalar TTR
) const
{
    if (TTR < TAbsLow_ || TTR > TAbsHigh_)
    {
        if (warnOutOfRange_)
        {
            WarningInFunction
                << "attempt to use janaf2TThermo<EquationOfState>"
                   " out of temperature range "
                << TAbsLow_ << " -> " << TAbsHigh_ << ";  T = " << TTR
                << endl;
        }
        return clamp(TTR, TAbsLow_, TAbsHigh_);
    }
    return TTR;
}


template<class EquationOfState>
inline Foam::scalar Foam::janaf2TThermo<EquationOfState>::TAbsLow() const
{
    return TAbsLow_;
}


template<class EquationOfState>
inline Foam::scalar Foam::janaf2TThermo<EquationOfState>::TAbsHigh() const
{
    return TAbsHigh_;
}

template<class EquationOfState>
inline Foam::scalar Foam::janaf2TThermo<EquationOfState>::TFirstInt() const
{
    return TFirstInt_;
}

template<class EquationOfState>
inline Foam::scalar Foam::janaf2TThermo<EquationOfState>::TSecondInt() const
{
    return TSecondInt_;
}

template<class EquationOfState>
inline Foam::scalar Foam::janaf2TThermo<EquationOfState>::ThetaVib() const
{
    return ThetaVib_;
}

template<class EquationOfState>
inline const typename Foam::janaf2TThermo<EquationOfState>::coeffArray&
Foam::janaf2TThermo<EquationOfState>::highCpCoeffs() const
{
    return highCpCoeffs_;
}

template<class EquationOfState>
inline const typename Foam::janaf2TThermo<EquationOfState>::coeffArray&
Foam::janaf2TThermo<EquationOfState>::midCpCoeffs() const
{
    return midCpCoeffs_;
}

template<class EquationOfState>
inline const typename Foam::janaf2TThermo<EquationOfState>::coeffArray&
Foam::janaf2TThermo<EquationOfState>::lowCpCoeffs() const
{
    return lowCpCoeffs_;
}


template<class EquationOfState>
inline Foam::scalar Foam::janaf2TThermo<EquationOfState>::CpTR
(
    const scalar p,
    const scalar TTR
) const
{
    const coeffArray& a = coeffs(TTR);
    return
        a[0]/(TTR*TTR) + a[1]/TTR + a[2] + (a[3] + (a[4] + (a[5] + a[6]*TTR)*TTR)*TTR)*TTR
      + EquationOfState::CpT(p, TTR);
}


template<class EquationOfState>
inline Foam::scalar Foam::janaf2TThermo<EquationOfState>::Ha
(
    const scalar p,
    const scalar TTR
) const
{
    const coeffArray& a = coeffs(TTR);
    return
    (
        ((((a[6]/5.0*TTR + a[5]/4.0)*TTR + a[4]/3.0)*TTR + a[3]/2.0)*TTR + a[2])*TTR
        + a[1]*log(TTR) - a[0]/TTR
      + a[7]
    ) + EquationOfState::HT(p, TTR);
}


template<class EquationOfState>
inline Foam::scalar Foam::janaf2TThermo<EquationOfState>::Hs
(
    const scalar p,
    const scalar TTR
) const
{
    return Ha(p, TTR) - Hc();
}


template<class EquationOfState>
inline Foam::scalar Foam::janaf2TThermo<EquationOfState>::Hc() const
{
    const coeffArray& a = lowCpCoeffs_;
    return
    (
        ((((a[6]/5.0*Tstd + a[5]/4.0)*Tstd + a[4]/3.0)*Tstd + a[3]/2.0)*Tstd + a[2])*Tstd
        + a[1]*log(Tstd) - a[0]/Tstd
        + a[7]
    );
}


template<class EquationOfState>
inline Foam::scalar Foam::janaf2TThermo<EquationOfState>::S
(
    const scalar p,
    const scalar TTR
) const
{
    const coeffArray& a = coeffs(TTR);
    return
    (
      -a[0]/(TTR*TTR*2) - a[1]/TTR + a[2]*log(TTR) + (a[3] + (a[4]/2 + (a[5]/3 + a[6]/4*TTR)*TTR)*TTR)*TTR
      + a[8]
    ) + EquationOfState::S(p, TTR);
}


template<class EquationOfState>
inline Foam::scalar Foam::janaf2TThermo<EquationOfState>::Gstd
(
    const scalar p,
    const scalar TTR
) const
{
    return
    (
        Ha(p,TTR) - S(p,TTR)*TTR
    );
}


template<class EquationOfState>
inline Foam::scalar Foam::janaf2TThermo<EquationOfState>::dCpdT
(
    const scalar p,
    const scalar TTR
) const
{
    const coeffArray& a = coeffs(TTR);
    return
        (-2)*a[0]/(TTR*TTR*TTR) - a[1]/(TTR*TTR) + 1 + a[3] + 2*a[4]*TTR + 3*a[5]*TTR*TTR + 4*a[6]*TTR*TTR*TTR;
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class EquationOfState>
inline void Foam::janaf2TThermo<EquationOfState>::operator+=
(
    const janaf2TThermo<EquationOfState>& jt
)
{
    scalar Y1 = this->Y();

    EquationOfState::operator+=(jt);

    if (mag(this->Y()) > SMALL)
    {
        Y1 /= this->Y();
        const scalar Y2 = jt.Y()/this->Y();
        {
            scalar wTheta = 0.0, sumTheta = 0.0;
            if (ThetaVib_ > SMALL)    { wTheta += Y1; sumTheta += Y1*ThetaVib_; }
            if (jt.ThetaVib_ > SMALL) { wTheta += Y2; sumTheta += Y2*jt.ThetaVib_; }
            ThetaVib_ = (wTheta > SMALL) ? sumTheta/wTheta : 0.0;
        }
        TAbsLow_ = max(TAbsLow_, jt.TAbsLow_);
        TAbsHigh_ = min(TAbsHigh_, jt.TAbsHigh_);

        if
        (
            janaf2TThermo<EquationOfState>::debug
         && !equal(TFirstInt_, jt.TFirstInt_)
        )
        {
            FatalErrorInFunction
                << "TFirstInt " << TFirstInt_ << " for "
                << (this->name().size() ? this->name() : "others")
                << " != " << jt.TFirstInt_ << " for "
                << (jt.name().size() ? jt.name() : "others")
                << exit(FatalError);
        }

        if
        (
            janaf2TThermo<EquationOfState>::debug
         && !equal(TSecondInt_, jt.TSecondInt_)
        )
        {
            FatalErrorInFunction
                << "TSecondInt " << TSecondInt_ << " for "
                << (this->name().size() ? this->name() : "others")
                << " != " << jt.TSecondInt_ << " for "
                << (jt.name().size() ? jt.name() : "others")
                << exit(FatalError);
        }

        for
        (
            label coefLabel=0;
            coefLabel<janaf2TThermo<EquationOfState>::nCoeffs_;
            coefLabel++
        )
        {
            highCpCoeffs_[coefLabel] =
                Y1*highCpCoeffs_[coefLabel]
              + Y2*jt.highCpCoeffs_[coefLabel];
            midCpCoeffs_[coefLabel] =
                Y1*midCpCoeffs_[coefLabel]
              + Y2*jt.midCpCoeffs_[coefLabel];

            lowCpCoeffs_[coefLabel] =
                Y1*lowCpCoeffs_[coefLabel]
              + Y2*jt.lowCpCoeffs_[coefLabel];
        }
    }
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::janaf2TThermo<EquationOfState> Foam::operator+
(
    const janaf2TThermo<EquationOfState>& jt1,
    const janaf2TThermo<EquationOfState>& jt2
)
{
    EquationOfState eofs = jt1;
    eofs += jt2;

    if (mag(eofs.Y()) < SMALL)
    {
        return janaf2TThermo<EquationOfState>
        (
            eofs,
            jt1.TAbsLow_,
            jt1.TAbsHigh_,
            jt1.TFirstInt_,
            jt1.TSecondInt_,
            jt1.ThetaVib_,
            jt1.highCpCoeffs_,
            jt1.midCpCoeffs_,
            jt1.lowCpCoeffs_,
            jt1.convertCoeffs_,
            jt1.warnOutOfRange_
        );
    }
    else
    {
        const scalar Y1 = jt1.Y()/eofs.Y();
        const scalar Y2 = jt2.Y()/eofs.Y();

        //const scalar thetaVibMix = jt1.ThetaVib_ * Y1 + jt2.ThetaVib_ * Y2;
        scalar wTheta = 0.0, sumTheta = 0.0;
        if (jt1.ThetaVib_ > SMALL) { wTheta += Y1; sumTheta += Y1*jt1.ThetaVib_; }
        if (jt2.ThetaVib_ > SMALL) { wTheta += Y2; sumTheta += Y2*jt2.ThetaVib_; }
        const scalar thetaVibMix = (wTheta > SMALL) ? sumTheta/wTheta : 0.0;
        
        typename janaf2TThermo<EquationOfState>::coeffArray highCpCoeffs;
        typename janaf2TThermo<EquationOfState>::coeffArray midCpCoeffs;
        typename janaf2TThermo<EquationOfState>::coeffArray lowCpCoeffs;

        for
        (
            label coefLabel=0;
            coefLabel<janaf2TThermo<EquationOfState>::nCoeffs_;
            coefLabel++
        )
        {
            highCpCoeffs[coefLabel] =
                Y1*jt1.highCpCoeffs_[coefLabel]
              + Y2*jt2.highCpCoeffs_[coefLabel];

            midCpCoeffs[coefLabel] =
                Y1*jt1.midCpCoeffs_[coefLabel]
              + Y2*jt2.midCpCoeffs_[coefLabel];

            lowCpCoeffs[coefLabel] =
                Y1*jt1.lowCpCoeffs_[coefLabel]
              + Y2*jt2.lowCpCoeffs_[coefLabel];
        }

        if
        (
            janaf2TThermo<EquationOfState>::debug
         && !equal(jt1.TFirstInt_, jt2.TFirstInt_)
        )
        {
            FatalErrorInFunction
                << "TFirstInt " << jt1.TFirstInt_ << " for "
                << (jt1.name().size() ? jt1.name() : "others")
                << " != " << jt2.TFirstInt_ << " for "
                << (jt2.name().size() ? jt2.name() : "others")
                << exit(FatalError);
        }

        if
        (
            janaf2TThermo<EquationOfState>::debug
         && !equal(jt1.TSecondInt_, jt2.TSecondInt_)
        )
        {
            FatalErrorInFunction
                << "TSecondInt " << jt1.TSecondInt_ << " for "
                << (jt1.name().size() ? jt1.name() : "others")
                << " != " << jt2.TSecondInt_ << " for "
                << (jt2.name().size() ? jt2.name() : "others")
                << exit(FatalError);
        }
        return janaf2TThermo<EquationOfState>
        (
            eofs,
            max(jt1.TAbsLow_, jt2.TAbsLow_),
            min(jt1.TAbsHigh_, jt2.TAbsHigh_),
            jt1.TFirstInt_,
            jt1.TSecondInt_,
            thetaVibMix,
            highCpCoeffs,
            midCpCoeffs,
            lowCpCoeffs,
            jt1.convertCoeffs_,
            jt1.warnOutOfRange_
        );
    }
}


template<class EquationOfState>
inline Foam::janaf2TThermo<EquationOfState> Foam::operator*
(
    const scalar s,
    const janaf2TThermo<EquationOfState>& jt
)
{
    return janaf2TThermo<EquationOfState>
    (
        s*static_cast<const EquationOfState&>(jt),
        jt.TAbsLow_,
        jt.TAbsHigh_,
        jt.TFirstInt_,
        jt.TSecondInt_,
        jt.ThetaVib_,
        jt.highCpCoeffs_,
        jt.midCpCoeffs_,
        jt.lowCpCoeffs_,
        jt.convertCoeffs_,
        jt.warnOutOfRange_
    );
}


template<class EquationOfState>
inline Foam::janaf2TThermo<EquationOfState> Foam::operator==
(
    const janaf2TThermo<EquationOfState>& jt1,
    const janaf2TThermo<EquationOfState>& jt2
)
{
    EquationOfState eofs
    (
        static_cast<const EquationOfState&>(jt1)
     == static_cast<const EquationOfState&>(jt2)
    );

    const scalar Y1 = jt1.Y()/eofs.Y();
    const scalar Y2 = jt2.Y()/eofs.Y();
    //const scalar thetaVibMix = jt1.ThetaVib_ * Y1 + jt2.ThetaVib_ * Y2;
    scalar wTheta = 0.0, sumTheta = 0.0;
    if (jt1.ThetaVib_ > SMALL) { wTheta += Y1; sumTheta += Y1*jt1.ThetaVib_; }
    if (jt2.ThetaVib_ > SMALL) { wTheta += Y2; sumTheta += Y2*jt2.ThetaVib_; }
    const scalar thetaVibMix = (wTheta > SMALL) ? sumTheta/wTheta : 0.0;
    
    typename janaf2TThermo<EquationOfState>::coeffArray highCpCoeffs;
    typename janaf2TThermo<EquationOfState>::coeffArray midCpCoeffs;
    typename janaf2TThermo<EquationOfState>::coeffArray lowCpCoeffs;

    for
    (
        label coefLabel=0;
        coefLabel<janaf2TThermo<EquationOfState>::nCoeffs_;
        coefLabel++
    )
    {
        highCpCoeffs[coefLabel] =
            Y1*jt1.highCpCoeffs_[coefLabel]
          - Y2*jt2.highCpCoeffs_[coefLabel];

        midCpCoeffs[coefLabel] =
            Y1*jt1.midCpCoeffs_[coefLabel]
          - Y2*jt2.midCpCoeffs_[coefLabel];

        lowCpCoeffs[coefLabel] =
            Y1*jt1.lowCpCoeffs_[coefLabel]
          - Y2*jt2.lowCpCoeffs_[coefLabel];
    }

    if
    (
        janaf2TThermo<EquationOfState>::debug
     && !equal(jt2.TFirstInt_, jt1.TFirstInt_)
    )
    {
        FatalErrorInFunction
            << "TFirstInt " << jt2.TFirstInt_ << " for "
            << (jt2.name().size() ? jt2.name() : "others")
            << " != " << jt1.TFirstInt_ << " for "
            << (jt1.name().size() ? jt1.name() : "others")
            << exit(FatalError);
    }

    return janaf2TThermo<EquationOfState>
    (
        eofs,
        max(jt2.TAbsLow_, jt1.TAbsLow_),
        min(jt2.TAbsHigh_, jt1.TAbsHigh_),
        jt2.TFirstInt_,
        jt2.TSecondInt_,
        thetaVibMix,
        highCpCoeffs,
        midCpCoeffs,
        lowCpCoeffs,
        jt2.convertCoeffs_,
        jt2.warnOutOfRange_
    );
}


// ************************************************************************* //
