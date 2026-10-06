/*---------------------------------------------------------------------------*\

    IHIS: In-house Hypersonics Implicit solver

    Copyright (C) 2014-2018 Johan Heyns - CSIR, South Africa
    Copyright (C) 2014-2018 Oliver Oxtoby - CSIR, South Africa
    Copyright (C) 2020-2021 Oliver Oxtoby
    Copyright (C) 2011-2013 OpenFOAM Foundation

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

\*---------------------------------------------------------------------------*/

#include "ihisModule.H"
#include "thermodynamicConstants.H"   // constant::thermodynamic::RR [J/kmol/K]


namespace Foam
{

namespace
{

//- Molecules below this mass fraction carry no meaningful e_v,m = eps_m/Y_m (D2)
constexpr scalar YMolMin = 1e-10;


//  Closed-form inverse of the harmonic oscillator (D2)
//  e_v = R*theta/(exp(theta/Tv) - 1)  =>  Tv = theta/log1p(R*theta/e_v)
//  Exact for neTCLib's EsVib (CpMCv = R); no start value, cannot fail.
inline scalar TvHO
(
    const scalar ev,
    const scalar R,
    const scalar theta,
    const scalar TvMin
)
{
    if (ev <= VSMALL)
    {
        return TvMin;                                // no vibrational energy -> floor
    }
    return max(theta/std::log1p(R*theta/ev), TvMin);
}


//  Energy-equivalent mixture vibrational temperature at one point
//  Solves  sum_m Y_m e_v,m(T) = sum_m eps_m  over molecules with Y_m >= YMolMin.
//  f(T) is increasing and convex (c_vv grows with T), and the root lies in
//  [min_m Tv_m, max_m Tv_m]; Newton started at the maximum therefore converges
//  monotonically from the right (Ortega & Rheinboldt 1970, 13.3).
inline scalar TvMixPoint
(
    const compositionInterface& comp,
    const labelList& molI,               // species index of molecule m
    const UList<scalar>& Ym,             // Y of each molecule at this point
    const UList<scalar>& Tvm,            // Tv of each molecule at this point
    const UList<scalar>& epsm,           // eps of each molecule at this point
    const scalar p,
    const scalar TTR,
    const scalar TvMin
)
{
    scalar target = 0;
    scalar Tlo = GREAT;
    scalar Thi = -GREAT;
    label nPresent = 0;

    forAll(molI, m)
    {
        if (Ym[m] >= YMolMin)
        {
            target += epsm[m];
            Tlo = min(Tlo, Tvm[m]);
            Thi = max(Thi, Tvm[m]);
            ++nPresent;
        }
    }

    if (nPresent == 0)
    {
        return TTR;                                  // fully dissociated: no molecule
    }
    if (Thi - Tlo <= 1e-12*Thi)
    {
        return Thi;                                  // uniform Tv (incl. Nm = 1): exact
    }

    scalar T = Thi;
    for (label iter = 0; iter < 50; ++iter)
    {
        scalar f = -target;
        scalar df = 0;

        forAll(molI, m)
        {
            if (Ym[m] < YMolMin) continue;

            const label s = molI[m];
            const scalar th = comp.thetaVib(s);
            f  += Ym[m]*comp.EsVib(s, p, TTR, T, th);
            df += Ym[m]*comp.CvVib(s, p, TTR, T, th);
        }

        if (df <= VSMALL)
        {
            return max(T, TvMin);                    // vibrationally frozen point
        }

        const scalar dT = f/df;
        T -= dT;

        if (mag(dT) <= 1e-12*T)
        {
            return max(T, TvMin);
        }
    }

    FatalErrorInFunction
        << "Mixture Tv Newton did not converge in 50 iterations"
        << " (T = " << T << ", target = " << target
        << ", Tv range [" << Tlo << ", " << Thi << "])"
        << exit(FatalError);

    return T;
}

//- Forward at one point: Tv_m -> e_v,m, eps_m; returns the mixture Tv
inline scalar encodePoint
(
    const compositionInterface& comp,
    const labelList& molI,
    const UList<scalar>& Ym,
    const UList<scalar>& Tvm,
    UList<scalar>& evm,
    UList<scalar>& epsm,
    const scalar p,
    const scalar TTR,
    const scalar TvMin
)
{
    forAll(molI, m)
    {
        const label s = molI[m];
        evm[m]  = comp.EsVib(s, p, TTR, Tvm[m], comp.thetaVib(s));
        epsm[m] = Ym[m]*evm[m];
    }
    return TvMixPoint(comp, molI, Ym, Tvm, epsm, p, TTR, TvMin);
}


//- Inverse at one point: eps_m -> e_v,m, Tv_m; returns the mixture Tv.
//  Trace molecules (Y_m < YMolMin) take Tv_m = mixture Tv of the others and
//  e_v,m = e_v,m(mixture Tv); their eps_m is never altered (energy kept, D2).
inline scalar decodePoint
(
    const compositionInterface& comp,
    const labelList& molI,
    const UList<scalar>& Ym,
    const UList<scalar>& epsm,
    UList<scalar>& evm,
    UList<scalar>& Tvm,
    const scalar p,
    const scalar TTR,
    const scalar TvMin,
    label& nFloor
)
{
    forAll(molI, m)
    {
        if (Ym[m] >= YMolMin)
        {
            const label s = molI[m];
            evm[m] = max(epsm[m], scalar(0))/Ym[m];
            Tvm[m] = TvHO
            (
                evm[m],
                constant::thermodynamic::RR/comp.W(s),
                comp.thetaVib(s),
                TvMin
            );
            if (Tvm[m] <= TvMin && evm[m] > VSMALL)
            {
                ++nFloor;
            }
        }
    }

    const scalar TvMix = TvMixPoint(comp, molI, Ym, Tvm, epsm, p, TTR, TvMin);

    forAll(molI, m)
    {
        if (Ym[m] < YMolMin)
        {
            const label s = molI[m];
            Tvm[m] = TvMix;
            evm[m] = comp.EsVib(s, p, TTR, TvMix, comp.thetaVib(s));
        }
    }
    return TvMix;
}

} // End anonymous namespace


scalar ihisModule::vibTvMin() const
{
    return pThermo_->TTR().mesh().solutionDict().subDict("pseudoTime")
        .lookupOrDefault<scalar>("TMin", ROOTSMALL);
}


void ihisModule::vibEncode()
{
    ne2TReactionThermo& neThermo = pThermo_();
    const compositionInterface& comp = neThermo.composition();
    const PtrList<volScalarField>& Y = comp.Y();

    const scalarField& pI = neThermo.p().primitiveField();
    const scalarField& TI = neThermo.TTR().primitiveField();
    scalarField& TvMixI   = neThermo.TVib().primitiveFieldRef();

    const label Nm = molI_.size();
    const scalar TvMin = vibTvMin();

    // Bind each field's cell values once, outside the cell loop
    UPtrList<const scalarField> YI(Nm), TvI(Nm);
    UPtrList<scalarField> evI(Nm), epsI(Nm);
    forAll(molI_, m)
    {
        YI.set(m, &Y[molI_[m]].primitiveField());
        TvI.set(m, &TVibSpecies_[molI_[m]].primitiveField());
        evI.set(m, &evVib_[m].primitiveFieldRef());
        epsI.set(m, &epsVib_[m].primitiveFieldRef());
    }

    scalarField Ym(Nm), Tvm(Nm), evm(Nm), epsm(Nm);

    forAll(TI, celli)
    {
        forAll(molI_, m)
        {
            Ym[m]  = YI[m][celli];
            Tvm[m] = TvI[m][celli];
        }

        TvMixI[celli] = encodePoint
        (
            comp, molI_, Ym, Tvm, evm, epsm, pI[celli], TI[celli], TvMin
        );

        forAll(molI_, m)
        {
            evI[m][celli]  = evm[m];
            epsI[m][celli] = epsm[m];
        }
    }

    vibBoundaries();
}


label ihisModule::vibDecode()
{
    ne2TReactionThermo& neThermo = pThermo_();
    const compositionInterface& comp = neThermo.composition();
    const PtrList<volScalarField>& Y = comp.Y();

    const scalarField& pI = neThermo.p().primitiveField();
    const scalarField& TI = neThermo.TTR().primitiveField();
    scalarField& TvMixI   = neThermo.TVib().primitiveFieldRef();

    const label Nm = molI_.size();
    const scalar TvMin = vibTvMin();

    UPtrList<const scalarField> YI(Nm), epsI(Nm);
    UPtrList<scalarField> evI(Nm), TvI(Nm);
    forAll(molI_, m)
    {
        YI.set(m, &Y[molI_[m]].primitiveField());
        epsI.set(m, &epsVib_[m].primitiveField());
        evI.set(m, &evVib_[m].primitiveFieldRef());
        TvI.set(m, &TVibSpecies_[molI_[m]].primitiveFieldRef());
    }

    scalarField Ym(Nm), Tvm(Nm), evm(Nm), epsm(Nm);
    label nFloor = 0;

    forAll(TI, celli)
    {
        forAll(molI_, m)
        {
            Ym[m]   = YI[m][celli];
            epsm[m] = epsI[m][celli];
        }

        TvMixI[celli] = decodePoint
        (
            comp, molI_, Ym, epsm, evm, Tvm, pI[celli], TI[celli], TvMin, nFloor
        );

        forAll(molI_, m)
        {
            evI[m][celli] = evm[m];
            TvI[m][celli] = Tvm[m];
        }
    }

    vibBoundaries();

    return returnReduce(nFloor, sumOp<label>());
}


void ihisModule::vibBoundaries()
{
    ne2TReactionThermo& neThermo = pThermo_();
    const compositionInterface& comp = neThermo.composition();
    const PtrList<volScalarField>& Y = comp.Y();

    const label Nm = molI_.size();
    const scalar TvMin = vibTvMin();

    // Atoms: hold T_TR, internal and boundary (Ns-long layout for V-T/V-V, F20)
    boolList isMol(Y.size(), false);
    forAll(molI_, m)
    {
        isMol[molI_[m]] = true;
    }
    forAll(Y, s)
    {
        if (!isMol[s])
        {
            TVibSpecies_[s] == neThermo.TTR();
        }
    }

    // Boundary rule, step 2: Tv_m from its own BCs (fixed kept, zeroGradient
    // copied, processor neighbours swapped); e_v, eps swapped on coupled patches
    forAll(molI_, m)
    {
        TVibSpecies_[molI_[m]].correctBoundaryConditions();
        evVib_[m].correctBoundaryConditions();
        epsVib_[m].correctBoundaryConditions();
    }

    volScalarField& TvMix = neThermo.TVib();
    TvMix.correctBoundaryConditions();

    // Boundary rule, step 3: non-coupled faces, forward from face Tv_m and Y
    const volScalarField::Boundary& pBf = neThermo.p().boundaryField();
    const volScalarField::Boundary& TBf = neThermo.TTR().boundaryField();
    volScalarField::Boundary& TvMixBf   = TvMix.boundaryFieldRef();

    scalarField Ym(Nm), Tvm(Nm), evm(Nm), epsm(Nm);

    forAll(TvMixBf, patchi)
    {
        if (TvMixBf[patchi].coupled()) continue;

        const label nF = TvMixBf[patchi].size();
        scalarField TvMixP(nF);
        List<scalarField> evP(Nm, scalarField(nF));
        List<scalarField> epsP(Nm, scalarField(nF));

        forAll(TvMixP, facei)
        {
            forAll(molI_, m)
            {
                const label s = molI_[m];
                Ym[m]  = Y[s].boundaryField()[patchi][facei];
                Tvm[m] = TVibSpecies_[s].boundaryField()[patchi][facei];
            }

            TvMixP[facei] = encodePoint
            (
                comp, molI_, Ym, Tvm, evm, epsm,
                pBf[patchi][facei], TBf[patchi][facei], TvMin
            );

            forAll(molI_, m)
            {
                evP[m][facei]  = evm[m];
                epsP[m][facei] = epsm[m];
            }
        }

        TvMixBf[patchi] == TvMixP;
        forAll(molI_, m)
        {
            evVib_[m].boundaryFieldRef()[patchi]  == evP[m];
            epsVib_[m].boundaryFieldRef()[patchi] == epsP[m];
        }
    }

    // Mixture vibrational energy for the decode/encode of rhoE (D1, F33)
    volScalarField& eVm = eVibMix_();
    eVm == dimensionedScalar(eVm.dimensions(), Zero);
    forAll(epsVib_, m)
    {
        eVm == eVm + epsVib_[m];
    }
}

void ihisModule::vibCheckModel() const
{
    const compositionInterface& comp = pThermo_->composition();
    const scalar p = 1e5;                                  // perfect gas: p-independent
    const scalarList Tcheck({50, 300, 3000, 20000});

    Info<< "IHIS: vibrational model checks" << nl;

    forAll(molI_, m)
    {
        const label s = molI_[m];
        const word& sName = comp.species()[s];
        const scalar th = comp.thetaVib(s);
        const scalar R  = constant::thermodynamic::RR/comp.W(s);

        // (1) Round trip Tv -> EsVib -> TvHO -> Tv (guards the closed-form assumption)
        scalar maxErr = 0;
        for (const scalar T : Tcheck)
        {
            const scalar ev = comp.EsVib(s, p, T, T, th);
            const scalar Tback = TvHO(ev, R, th, 0);
            maxErr = max(maxErr, mag(Tback - T)/T);
        }
        if (maxErr > 1e-10)
        {
            FatalErrorInFunction
                << "Closed-form Tv inverse does not match EsVib for " << sName
                << " (max relative error " << maxErr << ")." << nl
                << "    The library's vibrational energy is no longer a harmonic"
                << " oscillator with R = RR/W; TvHO (ihisModuleVib.C) must be"
                << " revised." << exit(FatalError);
        }

        // (2) Eucken test: kappaVib/CvVib independent of Tv?
        const scalar TTR = 1000;
        const scalar r1 =
            comp.kappaVib(s, p, TTR, 300, th)/comp.CvVib(s, p, TTR, 300, th);
        const scalar r2 =
            comp.kappaVib(s, p, TTR, 3000, th)/comp.CvVib(s, p, TTR, 3000, th);
        const bool eucken = mag(r1 - r2) <= 1e-12*mag(r1);

        Info<< "    " << sName << ": theta_v = " << th << " K"
            << ", inverse round trip " << maxErr
            << ", kappaVib/CvVib " << (eucken ? "Eucken-type" : "NOT Eucken-type")
            << nl;
    }
    Info<< endl;
}

} // End namespace Foam