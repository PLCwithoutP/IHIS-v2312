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

#include "thermo2T.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Thermo2T, template<class> class Type>
inline Foam::species::thermo2T<Thermo2T, Type>::thermo2T
(
    const Thermo2T& sp
)
:
    Thermo2T(sp)
{}


template<class Thermo2T, template<class> class Type>
inline Foam::scalar Foam::species::thermo2T<Thermo2T, Type>::TTR
(
    scalar f,
    scalar p,
    scalar T0,
    scalar (thermo2T<Thermo2T, Type>::*F)(const scalar, const scalar) const,
    scalar (thermo2T<Thermo2T, Type>::*dFdT)(const scalar, const scalar)
        const,
    scalar (thermo2T<Thermo2T, Type>::*limit)(const scalar) const
) const
{
    if (T0 < 0)
    {
        FatalErrorInFunction
            << "Negative initial temperature T0: " << T0
            << abort(FatalError);
    }

    scalar Test = T0;
    scalar Tnew = T0;
    scalar Ttol = T0*tol_;
    int    iter = 0;

    do
    {
        Test = Tnew;
        Tnew =
            (this->*limit)
            (Test - ((this->*F)(p, Test) - f)/(this->*dFdT)(p, Test));

        if (iter++ > maxIter_)
        {
            FatalErrorInFunction
                << "Maximum number of iterations exceeded: " << maxIter_
                << " when starting from T0:" << T0
                << " old T:" << Test << " new T:" << Tnew
                << " f:" << f
                << " p:" << p
                << " tol:" << Ttol
                << abort(FatalError);
        }

    } while (mag(Tnew - Test) > Ttol);

    return Tnew;
}

template<class Thermo2T, template<class> class Type>
inline Foam::scalar Foam::species::thermo2T<Thermo2T, Type>::TVib
(
    scalar f,
    scalar p,
    scalar TTR,
    scalar TVib0,
    scalar (thermo2T<Thermo2T, Type>::*F)(const scalar, const scalar, const scalar, const scalar) const,
    scalar (thermo2T<Thermo2T, Type>::*dFdT)(const scalar, const scalar, const scalar, const scalar)
        const,
    scalar (thermo2T<Thermo2T, Type>::*limit)(const scalar) const
) const
{
    const Foam::specie2T& sp = static_cast<const Foam::specie2T&>(*this);

    if (sp.isAtomic())
    {
        return 0;
    }
    if (TVib0 < 0)
    {
        FatalErrorInFunction
            << "Negative initial temperature TVib0: " << TVib0
            << abort(FatalError);
    }

    scalar Test      = TVib0;
    scalar Tnew      = TVib0;
    scalar Ttol      = TVib0*tol_;
    int    iter      = 0;

    do
    {
        Test = Tnew;

        const scalar fVal  = (this->*F)   (p, TTR, Test, this->ThetaVib()) - f;
        const scalar dfdT  = (this->*dFdT)(p, TTR, Test, this->ThetaVib());

        Tnew = (this->*limit)(Test - fVal/dfdT);

        // Secondary exit: function residual is smaller than what a step
        // of size Ttol would produce given the current slope.
        // Catches convergence when the function is flat (low TVib / high θ/T).
        if (mag(fVal) <= Ttol*mag(dfdT))
        {
            return Tnew;
        }

        if (iter++ > maxIter_)
        {
            FatalErrorInFunction
                << "Maximum number of iterations exceeded: " << maxIter_
                << " when starting from TVib0:" << TVib0
                << " old T:" << Test << " new T:" << Tnew
                << " f:" << f
                << " fVal:" << fVal
                << " dfdT:" << dfdT
                << " p:" << p
                << " tol:" << Ttol
                << abort(FatalError);
        }

    } while (mag(Tnew - Test) > Ttol);

    return Tnew;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo2T, template<class> class Type>
inline Foam::species::thermo2T<Thermo2T, Type>::thermo2T
(
    const word& name,
    const thermo2T& st
)
:
    Thermo2T(name, st)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo2T, template<class> class Type>
inline Foam::word
Foam::species::thermo2T<Thermo2T, Type>::hName()
{
    return Type<thermo2T<Thermo2T, Type>>::energyName();
}


template<class Thermo2T, template<class> class Type>
inline Foam::scalar
Foam::species::thermo2T<Thermo2T, Type>::CpvTR(const scalar p, const scalar TTR) const
{
    return Type<thermo2T<Thermo2T, Type>>::CpvTR(*this, p, TTR);
}


template<class Thermo2T, template<class> class Type>
inline Foam::scalar
Foam::species::thermo2T<Thermo2T, Type>::gamma(const scalar p, const scalar TTR) const
{
    #ifdef __clang__
    volatile const scalar CpTR = this->CpTR(p, TTR);
    volatile const scalar R = this->CpMCv(p, TTR);
    #else
    const scalar CpTR = this->CpTR(p, TTR);
    const scalar R = this->CpMCv(p, TTR);
    #endif
    //Info << "Gamma : " << CpTR/(CpTR - R) << nl << endl;
    return CpTR/(CpTR - R);
}


template<class Thermo2T, template<class> class Type>
inline Foam::scalar
Foam::species::thermo2T<Thermo2T, Type>::CpByCpvTR
(
    const scalar p,
    const scalar TTR
) const
{
    return Type<thermo2T<Thermo2T, Type>>::CpByCpvTR(*this, p, TTR);
}


template<class Thermo2T, template<class> class Type>
inline Foam::scalar
Foam::species::thermo2T<Thermo2T, Type>::H(const scalar p, const scalar TTR) const
{
    return Type<thermo2T<Thermo2T, Type>>::H(*this, p, TTR);
}

template<class Thermo2T, template<class> class Type>
inline Foam::scalar
Foam::species::thermo2T<Thermo2T, Type>::ET(const scalar p, const scalar TTR) const
{
    return Type<thermo2T<Thermo2T, Type>>::ET(*this, p, TTR);
}

template<class Thermo2T, template<class> class Type>
inline Foam::scalar
Foam::species::thermo2T<Thermo2T, Type>::ER(const scalar p, const scalar TTR) const
{
    return Type<thermo2T<Thermo2T, Type>>::ER(*this, p, TTR);
}

template<class Thermo2T, template<class> class Type>
inline Foam::scalar
Foam::species::thermo2T<Thermo2T, Type>::EV(const scalar p, const scalar TTR, const scalar TVib, const scalar ThetaVib) const
{
    return Type<thermo2T<Thermo2T, Type>>::EV(*this, p, TTR, TVib, ThetaVib);
}

template<class Thermo2T, template<class> class Type>
inline Foam::scalar
Foam::species::thermo2T<Thermo2T, Type>::G(const scalar p, const scalar TTR) const
{
    return this->Ha(p, TTR) - TTR*this->S(p, TTR);
}


template<class Thermo2T, template<class> class Type>
inline Foam::scalar
Foam::species::thermo2T<Thermo2T, Type>::cpTR(const scalar p, const scalar TTR) const
{
    return this->CpTR(p, TTR)*this->W();
    //return this->CpTR(p, TTR);
}


template<class Thermo2T, template<class> class Type>
inline Foam::scalar
Foam::species::thermo2T<Thermo2T, Type>::ha(const scalar p, const scalar TTR) const
{
    return this->Ha(p, TTR)*this->W();
}


template<class Thermo2T, template<class> class Type>
inline Foam::scalar
Foam::species::thermo2T<Thermo2T, Type>::hs(const scalar p, const scalar TTR) const
{
    return this->Hs(p, TTR)*this->W();
}


template<class Thermo2T, template<class> class Type>
inline Foam::scalar
Foam::species::thermo2T<Thermo2T, Type>::hc() const
{
    return this->Hc()*this->W();
}


template<class Thermo2T, template<class> class Type>
inline Foam::scalar
Foam::species::thermo2T<Thermo2T, Type>::s(const scalar p, const scalar TTR) const
{
    return this->S(p, TTR)*this->W();
}


template<class Thermo2T, template<class> class Type>
inline Foam::scalar
Foam::species::thermo2T<Thermo2T, Type>::h(const scalar p, const scalar TTR) const
{
    return this->H(p, TTR)*this->W();
}


template<class Thermo2T, template<class> class Type>
inline Foam::scalar
Foam::species::thermo2T<Thermo2T, Type>::cvT(const scalar p, const scalar TTR) const
{
    return this->CvT(p, TTR)*this->W();
}

template<class Thermo2T, template<class> class Type>
inline Foam::scalar
Foam::species::thermo2T<Thermo2T, Type>::cvR(const scalar p, const scalar TTR) const
{
    return this->CvR(p, TTR)*this->W();
}

template<class Thermo2T, template<class> class Type>
inline Foam::scalar
Foam::species::thermo2T<Thermo2T, Type>::cvVib(const scalar p, const scalar TTR, const scalar TVib, const scalar ThetaVib) const
{
    return this->CvVib(p, TTR, TVib, ThetaVib)*this->W();
}

template<class Thermo2T, template<class> class Type>
inline Foam::scalar
Foam::species::thermo2T<Thermo2T, Type>::esT(const scalar p, const scalar TTR) const
{
    return this->EsT(p, TTR)*this->W();
}

template<class Thermo2T, template<class> class Type>
inline Foam::scalar
Foam::species::thermo2T<Thermo2T, Type>::esR(const scalar p, const scalar TTR) const
{
    return this->EsR(p, TTR)*this->W();
}

template<class Thermo2T, template<class> class Type>
inline Foam::scalar
Foam::species::thermo2T<Thermo2T, Type>::esVib(const scalar p, const scalar TTR, const scalar TVib, const scalar ThetaVib) const
{
    return this->EsVib(p, TTR, TVib, ThetaVib)*this->W();
}

template<class Thermo2T, template<class> class Type>
inline Foam::scalar
Foam::species::thermo2T<Thermo2T, Type>::eaT(const scalar p, const scalar TTR) const
{
    return this->EaT(p, TTR)*this->W();
}


template<class Thermo2T, template<class> class Type>
inline Foam::scalar
Foam::species::thermo2T<Thermo2T, Type>::g(const scalar p, const scalar TTR) const
{
    return this->G(p, TTR)*this->W();
}


template<class Thermo2T, template<class> class Type>
inline Foam::scalar
Foam::species::thermo2T<Thermo2T, Type>::a(const scalar p, const scalar TTR) const
{
    return this->A(p, TTR)*this->W();
}


template<class Thermo2T, template<class> class Type>
inline Foam::scalar
Foam::species::thermo2T<Thermo2T, Type>::K(const scalar p, const scalar TTR) const
{
    scalar arg = -this->Y()*this->G(Pstd, TTR)/(RR*TTR);

    if (arg < 600)
    {
        return exp(arg);
    }
    else
    {
        return VGREAT;
    }
}


template<class Thermo2T, template<class> class Type>
inline Foam::scalar
Foam::species::thermo2T<Thermo2T, Type>::Kp(const scalar p, const scalar TTR) const
{
    return K(p, TTR);
}


template<class Thermo2T, template<class> class Type>
inline Foam::scalar
Foam::species::thermo2T<Thermo2T, Type>::Kc(const scalar p, const scalar TTR) const
{
    const scalar nm = this->Y()/this->W();

    if (equal(nm, SMALL))
    {
        return Kp(p, TTR);
    }
    else
    {
        return Kp(p, TTR)*pow(Pstd/(RR*TTR), nm);
    }
}


template<class Thermo2T, template<class> class Type>
inline Foam::scalar Foam::species::thermo2T<Thermo2T, Type>::Kx
(
    const scalar p,
    const scalar TTR
) const
{
    const scalar nm = this->Y()/this->W();

    if (equal(nm, SMALL))
    {
        return Kp(p, TTR);
    }
    else
    {
        return Kp(p, TTR)*pow(Pstd/p, nm);
    }
}


template<class Thermo2T, template<class> class Type>
inline Foam::scalar Foam::species::thermo2T<Thermo2T, Type>::Kn
(
    const scalar p,
    const scalar TTR,
    const scalar n
) const
{
    const scalar nm = this->Y()/this->W();

    if (equal(nm, SMALL))
    {
        return Kp(p, TTR);
    }
    else
    {
        return Kp(p, TTR)*pow(n*Pstd/p, nm);
    }
}

template<class Thermo2T, template<class> class Type>
inline Foam::scalar Foam::species::thermo2T<Thermo2T, Type>::TE_TR
(
    const scalar eTR,
    const scalar p,
    const scalar TTR0
) const
{
    return Type<thermo2T<Thermo2T, Type>>::TE_TR(*this, eTR, p, TTR0);
}

template<class Thermo2T, template<class> class Type>
inline Foam::scalar Foam::species::thermo2T<Thermo2T, Type>::TH_TR
(
    const scalar h,
    const scalar p,
    const scalar TTR0
) const
{
    return Type<thermo2T<Thermo2T, Type>>::TH_TR(*this, h, p, TTR0);
}

template<class Thermo2T, template<class> class Type>
inline Foam::scalar Foam::species::thermo2T<Thermo2T, Type>::TE_Vib
(
    const scalar esVib,
    const scalar p,
    const scalar TTR0,
    const scalar TVib0
) const
{
    return Type<thermo2T<Thermo2T, Type>>::TE_Vib(*this, esVib, p, TTR0, TVib0);
}


template<class Thermo2T, template<class> class Type>
inline Foam::scalar Foam::species::thermo2T<Thermo2T, Type>::THs_TR
(
    const scalar hs,
    const scalar p,
    const scalar TTR0
) const
{
    return TTR
    (
        hs,
        p,
        TTR0,
        &thermo2T<Thermo2T, Type>::Hs,
        &thermo2T<Thermo2T, Type>::CpTR,
        &thermo2T<Thermo2T, Type>::limit
    );
}


template<class Thermo2T, template<class> class Type>
inline Foam::scalar Foam::species::thermo2T<Thermo2T, Type>::THa_TR
(
    const scalar ha,
    const scalar p,
    const scalar TTR0
) const
{
    return TTR
    (
        ha,
        p,
        TTR0,
        &thermo2T<Thermo2T, Type>::Ha,
        &thermo2T<Thermo2T, Type>::CpTR,
        &thermo2T<Thermo2T, Type>::limit
    );
}


template<class Thermo2T, template<class> class Type>
inline Foam::scalar Foam::species::thermo2T<Thermo2T, Type>::TEs_TR
(
    const scalar es,
    const scalar p,
    const scalar TTR0
) const
{
    return TTR
    (
        es,
        p,
        TTR0,
        &thermo2T<Thermo2T, Type>::EsT,
        &thermo2T<Thermo2T, Type>::CvT,
        &thermo2T<Thermo2T, Type>::limit
    );
}

template<class Thermo2T, template<class> class Type>
inline Foam::scalar Foam::species::thermo2T<Thermo2T, Type>::TEs_Vib
(
    const scalar esVib,
    const scalar p,
    const scalar TVib0,
    const scalar TTR0
) const
{
    return TVib
    (
        esVib,
        p,
        TTR0,
        TVib0,
        &thermo2T<Thermo2T, Type>::EsVib,
        &thermo2T<Thermo2T, Type>::CvVib,
        &thermo2T<Thermo2T, Type>::limit
    );
}


template<class Thermo2T, template<class> class Type>
inline Foam::scalar Foam::species::thermo2T<Thermo2T, Type>::TEa_TR
(thermo2T_H
    const scalar ea,
    const scalar p,
    const scalar TTR0
) const
{
    return TTR
    (
        ea,
        p,
        TTR0,
        &thermo2T<Thermo2T, Type>::EaT,
        &thermo2T<Thermo2T, Type>::CvT,
        &thermo2T<Thermo2T, Type>::limit
    );
}


template<class Thermo2T, template<class> class Type>
inline Foam::scalar
Foam::species::thermo2T<Thermo2T, Type>::dKcdTbyKc
(
    const scalar p,
    const scalar TTR
) const
{
    const scalar dKcdTbyKc =
        (this->S(Pstd, TTR) + this->Gstd(TTR)/TTR)*this->Y()/(RR*TTR);

    const scalar nm = this->Y()/this->W();
    if (equal(nm, SMALL))
    {
        return dKcdTbyKc;
    }
    else
    {
        return dKcdTbyKc - nm/TTR;
    }
}


template<class Thermo2T, template<class> class Type>
inline Foam::scalar
Foam::species::thermo2T<Thermo2T, Type>::dcpdT(const scalar p, const scalar TTR) const
{
    return this->dCpdT(p, TTR)*this->W();
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Thermo2T, template<class> class Type>
inline void Foam::species::thermo2T<Thermo2T, Type>::operator+=
(
    const thermo2T<Thermo2T, Type>& st
)
{
    Thermo2T::operator+=(st);
}


template<class Thermo2T, template<class> class Type>
inline void Foam::species::thermo2T<Thermo2T, Type>::operator*=(const scalar s)
{
    Thermo2T::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Thermo2T, template<class> class Type>
inline Foam::species::thermo2T<Thermo2T, Type> Foam::species::operator+
(
    const thermo2T<Thermo2T, Type>& st1,
    const thermo2T<Thermo2T, Type>& st2
)
{
    return thermo2T<Thermo2T, Type>
    (
        static_cast<const Thermo2T&>(st1) + static_cast<const Thermo2T&>(st2)
    );
}


template<class Thermo2T, template<class> class Type>
inline Foam::species::thermo2T<Thermo2T, Type> Foam::species::operator*
(
    const scalar s,
    const thermo2T<Thermo2T, Type>& st
)
{
    return thermo2T<Thermo2T, Type>
    (
        s*static_cast<const Thermo2T&>(st)
    );
}


template<class Thermo2T, template<class> class Type>
inline Foam::species::thermo2T<Thermo2T, Type> Foam::species::operator==
(
    const thermo2T<Thermo2T, Type>& st1,
    const thermo2T<Thermo2T, Type>& st2
)
{
    return thermo2T<Thermo2T, Type>
    (
        static_cast<const Thermo2T&>(st1) == static_cast<const Thermo2T&>(st2)
    );
}


// ************************************************************************* //
