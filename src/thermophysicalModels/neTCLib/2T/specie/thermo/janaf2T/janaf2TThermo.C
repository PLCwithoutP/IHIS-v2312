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

#include "janaf2TThermo.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class EquationOfState>
void Foam::janaf2TThermo<EquationOfState>::checkInputData() const
{
    if (TAbsLow_ >= TAbsHigh_)
    {
        FatalErrorInFunction
            << "TAbsLow(" << TAbsLow_ << ") >= TAbsHigh(" << TAbsHigh_ << ')'
            << exit(FatalError);
    }

    if ((TFirstInt_ <= TAbsLow_))
    {
        FatalErrorInFunction
            << "TFirstInt(" << TFirstInt_ << ") <= TAbsLow(" << TAbsLow_ << ')'
            << exit(FatalError);
    }

    if ((TSecondInt_ <= TAbsLow_))
    {
        FatalErrorInFunction
            << "TSecondInt(" << TSecondInt_ << ") <= TAbsLow(" << TAbsLow_ << ')'
            << exit(FatalError);
    }

    if (TFirstInt_ > TAbsHigh_)
    {
        FatalErrorInFunction
            << "TFirstInt(" << TFirstInt_ << ") > TAbsHigh(" << TAbsHigh_ << ')'
            << exit(FatalError);
    }

    if (TSecondInt_ > TAbsHigh_)
    {
        FatalErrorInFunction
            << "TSecondInt(" << TSecondInt_ << ") > TAbsHigh(" << TAbsHigh_ << ')'
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class EquationOfState>
Foam::janaf2TThermo<EquationOfState>::janaf2TThermo(const dictionary& dict)
:
    EquationOfState(dict),
    TAbsLow_(dict.subDict("thermodynamics").get<scalar>("TAbsLow")),
    TAbsHigh_(dict.subDict("thermodynamics").get<scalar>("TAbsHigh")),
    TFirstInt_(dict.subDict("thermodynamics").get<scalar>("TFirstInt")),
    TSecondInt_(dict.subDict("thermodynamics").get<scalar>("TSecondInt")),
    ThetaVib_(dict.subDict("thermodynamics").get<scalar>("ThetaVibrational")),
    warnOutOfRange_(dict.subDict("thermodynamics").getOrDefault<bool>("warnOutOfRange", true)),
    highCpCoeffs_(dict.subDict("thermodynamics").lookup("highCpCoeffs")),
    midCpCoeffs_(dict.subDict("thermodynamics").lookup("midCpCoeffs")),
    lowCpCoeffs_(dict.subDict("thermodynamics").lookup("lowCpCoeffs"))
{
    // Convert coefficients to mass-basis
    for (label coefLabel=0; coefLabel<nCoeffs_; coefLabel++)
    {
        highCpCoeffs_[coefLabel] *= this->R();
        midCpCoeffs_[coefLabel] *= this->R();
        lowCpCoeffs_[coefLabel] *= this->R();
    }

    checkInputData();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class EquationOfState>
void Foam::janaf2TThermo<EquationOfState>::write(Ostream& os) const
{
    EquationOfState::write(os);

    // Convert coefficients back to dimensionless form
    coeffArray highCpCoeffs;
    coeffArray midCpCoeffs;
    coeffArray lowCpCoeffs;
    for (label coefLabel=0; coefLabel<nCoeffs_; coefLabel++)
    {
        highCpCoeffs[coefLabel] = highCpCoeffs_[coefLabel]/this->R();
        midCpCoeffs[coefLabel] = midCpCoeffs_[coefLabel]/this->R();
        lowCpCoeffs[coefLabel] = lowCpCoeffs_[coefLabel]/this->R();
    }

    // Entries in dictionary format
    {
        os.beginBlock("thermodynamics");
        os.writeEntry("TAbsLow", TAbsLow_);
        os.writeEntry("TAbsHigh", TAbsHigh_);
        os.writeEntry("TFirstInt", TFirstInt_);
        os.writeEntry("TSecondInt", TSecondInt_);
        os.writeEntry("ThetaVib", ThetaVib_);
        os.writeEntry("highCpCoeffs", highCpCoeffs);
        os.writeEntry("midCpCoeffs", midCpCoeffs);
        os.writeEntry("lowCpCoeffs", lowCpCoeffs);
        os.endBlock();
    }
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class EquationOfState>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const janaf2TThermo<EquationOfState>& jt
)
{
    jt.write(os);
    return os;
}


// ************************************************************************* //
