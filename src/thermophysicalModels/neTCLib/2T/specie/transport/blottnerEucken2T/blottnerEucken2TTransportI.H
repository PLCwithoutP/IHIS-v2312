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

\*---------------------------------------------------------------------------*/

#include "specie2T.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo2T>
inline Foam::blottnerEucken2TTransport<Thermo2T>::blottnerEucken2TTransport
(
    const Thermo2T& t,
    const scalar AB,
    const scalar BB,
    const scalar CB,
    const scalar Sc
)
:
    Thermo2T(t),
    AB_(AB),
    BB_(BB),
    CB_(CB),
    Sc_(Sc)
{}


template<class Thermo2T>
inline Foam::blottnerEucken2TTransport<Thermo2T>::blottnerEucken2TTransport
(
    const word& name,
    const blottnerEucken2TTransport& st
)
:
    Thermo2T(name, st),
    AB_(st.AB_),
    BB_(st.BB_),
    CB_(st.CB_),
    Sc_(st.Sc_)
{}


template<class Thermo2T>
inline Foam::autoPtr<Foam::blottnerEucken2TTransport<Thermo2T>>
Foam::blottnerEucken2TTransport<Thermo2T>::clone() const
{
    return autoPtr<blottnerEucken2TTransport<Thermo2T>>::New(*this);
}


template<class Thermo2T>
inline Foam::autoPtr<Foam::blottnerEucken2TTransport<Thermo2T>>
Foam::blottnerEucken2TTransport<Thermo2T>::New
(
    const dictionary& dict
)
{
    return autoPtr<blottnerEucken2TTransport<Thermo2T>>::New(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo2T>
inline Foam::scalar Foam::blottnerEucken2TTransport<Thermo2T>::mu
(
    const scalar TTR
) const
{
    //Info << "Mu : " << 0.1*exp((AB_*log(TTR) + BB_)*log(TTR) + CB_) << endl;
    return 0.1*exp((AB_*log(TTR) + BB_)*log(TTR) + CB_);
}


template<class Thermo2T>
inline Foam::scalar Foam::blottnerEucken2TTransport<Thermo2T>::kappaTR
(
    const scalar p, const scalar TTR
) const
{
    // Eucken: kappa_tr = mu*(5/2*Cv_trans + Cv_rot).
    //
    // The 5/2 factor applies to TRANSLATION ONLY. Note that CvT() in
    // HtoEthermo2T.H returns CpTR - R, which for a diatomic is 2.5R --
    // translation AND rotation, not translation alone. Subtracting CvR()
    // recovers the translational part (1.5R for both atoms and molecules).
    //
    // Without the subtraction a diatomic gets mu*7.25R instead of mu*4.75R,
    // a factor of 1.526 too high. Atoms are unaffected either way, since
    // CvR() returns 0 for them.
    //
    // DO NOT "fix" this by changing CvT() itself: gamma() in thermo2TI.H is
    // CpTR/CvT and requires the full translational-rotational value.

    const scalar Cv_R_ = this->CvR(p, TTR);
    const scalar Cv_t_ = this->CvT(p, TTR) - Cv_R_;

    return mu(TTR)*(2.5*Cv_t_ + Cv_R_);
}

template<class Thermo2T>
inline Foam::scalar Foam::blottnerEucken2TTransport<Thermo2T>::kappaVib
(
    const scalar p,
    const scalar TTR,
    const scalar TVib,
    const scalar ThetaVib
) const
{
    if (ThetaVib < SMALL)
    {
        return 0.0;
    }
    
    return mu(TTR)*this->CvVib(p, TTR, TVib, ThetaVib)/Sc_;
}


template<class Thermo2T>
inline Foam::scalar Foam::blottnerEucken2TTransport<Thermo2T>::alphah
(
    const scalar p,
    const scalar TTR
) const
{

    return kappaTR(p, TTR)/this->CpTR(p, TTR);
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Thermo2T>
inline void Foam::blottnerEucken2TTransport<Thermo2T>::operator+=
(
    const blottnerEucken2TTransport<Thermo2T>& st
)
{
    scalar Y1 = this->Y();

    Thermo2T::operator+=(st);

    if (mag(this->Y()) > SMALL)
    {
        Y1 /= this->Y();
        scalar Y2 = st.Y()/this->Y();

        AB_ = Y1*AB_ + Y2*st.AB_;
        BB_ = Y1*BB_ + Y2*st.BB_;
        CB_ = Y1*CB_ + Y2*st.CB_;
        Sc_ = Y1*Sc_ + Y2*st.Sc_;
    }
}


template<class Thermo2T>
inline void Foam::blottnerEucken2TTransport<Thermo2T>::operator*=
(
    const scalar s
)
{
    Thermo2T::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Thermo2T>
inline Foam::blottnerEucken2TTransport<Thermo2T> Foam::operator+
(
    const blottnerEucken2TTransport<Thermo2T>& st1,
    const blottnerEucken2TTransport<Thermo2T>& st2
)
{
    Thermo2T t
    (
        static_cast<const Thermo2T&>(st1) + static_cast<const Thermo2T&>(st2)
    );

    if (mag(t.Y()) < SMALL)
    {
        return blottnerEucken2TTransport<Thermo2T>
        (
            t, st1.AB_, st1.BB_, st1.CB_, st1.Sc_
        );
    }
    else
    {
        scalar Y1 = st1.Y()/t.Y();
        scalar Y2 = st2.Y()/t.Y();

        return blottnerEucken2TTransport<Thermo2T>
        (
            t,
            Y1*st1.AB_ + Y2*st2.AB_,
            Y1*st1.BB_ + Y2*st2.BB_,
            Y1*st1.CB_ + Y2*st2.CB_,
            Y1*st1.Sc_ + Y2*st2.Sc_
        );
    }
}


template<class Thermo2T>
inline Foam::blottnerEucken2TTransport<Thermo2T> Foam::operator*
(
    const scalar s,
    const blottnerEucken2TTransport<Thermo2T>& st
)
{
    return blottnerEucken2TTransport<Thermo2T>
    (
        s*static_cast<const Thermo2T&>(st),
        st.AB_,
        st.BB_,
        st.CB_,
        st.Sc_
    );
}


// ************************************************************************* //
