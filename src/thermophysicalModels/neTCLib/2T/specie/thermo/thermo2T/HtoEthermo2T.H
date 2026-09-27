inline scalar CvT
(
    const scalar p,
    const scalar TTR
) const
{
    #ifdef __clang__
    volatile const scalar cpTR = CpTR(p, TTR);
    volatile const scalar cpmcvTR = EquationOfState::CpMCv(p, TTR);
    return cpTR - cpmcvTR;
    #else
    return CpTR(p, TTR) - EquationOfState::CpMCv(p, TTR);
    #endif
}

inline scalar CvR
(
    const scalar p,
    const scalar TTR
) const
{
    #ifdef __clang__
    if (static_cast<const specie2T&>(*this).isAtomic())
    {
        return 0;
    }
    else
    {
        volatile const scalar cpmcvTR = EquationOfState::CpMCv(p, TTR);
        return cpmcvTR;
    }
    #else
    if (static_cast<const specie2T&>(*this).isAtomic())
    {
        return 0;
    }
    else
    {
        return EquationOfState::CpMCv(p, TTR);
    }
    #endif
}

inline scalar CvVib
(
    const scalar p,
    const scalar TTR,
    const scalar TVib,
    const scalar ThetaVib
) const
{
    #ifdef __clang__
    if (static_cast<const specie2T&>(*this).isAtomic())
    {
        return 0;
    }
    else
    {
        volatile const scalar cpmcvTR = EquationOfState::CpMCv(p, TTR);
        volatile const scalar x = 0.5*ThetaVib/TVib;
        volatile const scalar xOverSinh = x/sinh(x);
        return cpmcvTR*xOverSinh*xOverSinh;
    }
    #else
    if (static_cast<const specie2T&>(*this).isAtomic())
    {
        return 0;
    }
    if (TVib <= SMALL)
    {
        return 0;   // derivative → 0 as T → 0 for SHO
    }
    else
    {
        const scalar x = 0.5*ThetaVib/TVib;
        const scalar xOverSinh = x/sinh(x);
        return EquationOfState::CpMCv(p, TTR)*xOverSinh*xOverSinh;
    }
    #endif
}

inline scalar EsT
(
    const scalar p,
    const scalar TTR
) const
{
    #ifdef __clang__
    volatile const scalar hsT = Hs(p, TTR);
    volatile const scalar rho = EquationOfState::rho(p, TTR);
    return hsT - p/rho;
    #else
    return Hs(p, TTR) - p/EquationOfState::rho(p, TTR);
    #endif
}

inline scalar EsR
(
    const scalar p,
    const scalar TTR
) const
{
    #ifdef __clang__
    if (static_cast<const specie2T&>(*this).isAtomic())
    {
        return 0;
    }
    else
    {
        volatile const scalar cpmcvTR = EquationOfState::CpMCv(p, TTR);
        return cpmcvTR*TTR;
    }
    #else
    if (static_cast<const specie2T&>(*this).isAtomic())
    {
        return 0;
    }
    else
    {
        return EquationOfState::CpMCv(p, TTR)*TTR;
    }
    #endif
}

inline scalar EsVib
(
    const scalar p,
    const scalar TTR,
    const scalar TVib,
    const scalar ThetaVib
) const
{
    #ifdef __clang__
    if (static_cast<const specie2T&>(*this).isAtomic())
    {
        return 0;
    }
    else
    {
        volatile const scalar div = exp(ThetaVib/TVib) - 1;
        volatile const scalar cpmcvTR = EquationOfState::CpMCv(p, TTR);
        return cpmcvTR*ThetaVib/div;
    }
    #else
    if (static_cast<const specie2T&>(*this).isAtomic())
    {
        return 0;
    }
    if (TVib <= SMALL)
    {
        return 0;   // SHO ground state: no vibrational energy at T→0
    }
    else
    {
        return EquationOfState::CpMCv(p, TTR)*ThetaVib/(exp(ThetaVib/TVib) - 1);
    }
    #endif
}

inline scalar EaT
(
    const scalar p,
    const scalar TTR
) const
{
    #ifdef __clang__
    volatile const scalar haT = Ha(p, TTR);
    volatile const scalar rho = EquationOfState::rho(p, TTR);
    return haT - p/rho;
    #else
    return Ha(p, TTR) - p/EquationOfState::rho(p, TTR);
    #endif
}
