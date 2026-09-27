/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2015-2021 OpenCFD Ltd.
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

#include "he2TThermo.H"
#include "gradientEnergyFvPatchScalarField.H"
#include "mixedEnergyFvPatchScalarField.H"
#include "gradient2TEnergyFvPatchScalarField.H"
#include "mixed2TEnergyFvPatchScalarField.H"
#include "gradientVibEnergyFvPatchScalarField.H"
#include "mixedVibEnergyFvPatchScalarField.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Basic2TThermo, class MixtureType>
void Foam::he2TThermo<Basic2TThermo, MixtureType>::
hBoundaryCorrection(volScalarField& e)
{
    volScalarField::Boundary& eBf = e.boundaryFieldRef();

    forAll(eBf, patchi)
    {
        if (isA<gradientEnergyFvPatchScalarField>(eBf[patchi]))
        {
            refCast<gradientEnergyFvPatchScalarField>(eBf[patchi]).gradient()
                = eBf[patchi].fvPatchField::snGrad();
        }
        else if (isA<mixedEnergyFvPatchScalarField>(eBf[patchi]))
        {
            refCast<mixedEnergyFvPatchScalarField>(eBf[patchi]).refGrad()
                = eBf[patchi].fvPatchField::snGrad();
        }
        else if (isA<gradient2TEnergyFvPatchScalarField>(eBf[patchi]))
        {
            refCast<gradient2TEnergyFvPatchScalarField>(eBf[patchi]).gradient()
                = eBf[patchi].fvPatchField::snGrad();
        }
        else if (isA<mixed2TEnergyFvPatchScalarField>(eBf[patchi]))
        {
            refCast<mixed2TEnergyFvPatchScalarField>(eBf[patchi]).refGrad()
                = eBf[patchi].fvPatchField::snGrad();
        }
    }
}

template<class Basic2TThermo, class MixtureType>
void Foam::he2TThermo<Basic2TThermo, MixtureType>::
eVibBoundaryCorrection(volScalarField& e)
{
    volScalarField::Boundary& eBf = e.boundaryFieldRef();

    forAll(eBf, patchi)
    {
        if (isA<gradientEnergyFvPatchScalarField>(eBf[patchi]))
        {
            refCast<gradientEnergyFvPatchScalarField>(eBf[patchi]).gradient()
                = eBf[patchi].fvPatchField::snGrad();
        }
        else if (isA<mixedEnergyFvPatchScalarField>(eBf[patchi]))
        {
            refCast<mixedEnergyFvPatchScalarField>(eBf[patchi]).refGrad()
                = eBf[patchi].fvPatchField::snGrad();
        }
        else if (isA<gradientVibEnergyFvPatchScalarField>(eBf[patchi]))
        {
            refCast<gradientVibEnergyFvPatchScalarField>(eBf[patchi]).gradient()
                = eBf[patchi].fvPatchField::snGrad();
        }
        else if (isA<mixedVibEnergyFvPatchScalarField>(eBf[patchi]))
        {
            refCast<mixedVibEnergyFvPatchScalarField>(eBf[patchi]).refGrad()
                = eBf[patchi].fvPatchField::snGrad();
        }
    }
}

template<class Basic2TThermo, class MixtureType>
void Foam::he2TThermo<Basic2TThermo, MixtureType>::init
(
    const volScalarField& p,
    const volScalarField& TTR,
    const volScalarField& TVib,
    volScalarField& h,
    volScalarField& eTR,
    volScalarField& eT,
    volScalarField& eR,
    volScalarField& eVib
)
{
    const bool use2T = this->twoTemperature();
    const scalar theta = this->cellMixture(0).ThetaVib();
    const label nSpecie = this->Y().size();

    scalarField& hCells   = h.primitiveFieldRef();
    scalarField& eTRCells = eTR.primitiveFieldRef();    
    scalarField& eVibCells = eVib.primitiveFieldRef();
    const scalarField& pCells   = p.primitiveField();
    const scalarField& TTRCells = TTR.primitiveField();
    const scalarField& TVibCells = TVib.primitiveField();

    forAll(hCells, celli)
    {
        const scalar TVibUse = use2T ? TVibCells[celli] : TTRCells[celli];

        hCells[celli] =
            this->cellMixture(celli).H(pCells[celli], TTRCells[celli]);

        eTRCells[celli] =                               
            this->cellMixture(celli).EsT               
            (                                          
                pCells[celli],                         
                TTRCells[celli]                        
            );                                         

        scalar evMix = 0;
        for (label speciei = 0; speciei < nSpecie; ++speciei)
        {
            evMix +=
                this->Y()[speciei][celli]
              * this->EsVib
                (
                    speciei, pCells[celli], TTRCells[celli],
                    TVibUse, this->composition().thetaVib(speciei)
                );
        }
        eVibCells[celli] = evMix;
    }

    volScalarField::Boundary& hBf    = h.boundaryFieldRef();
    volScalarField::Boundary& eTRBf  = eTR.boundaryFieldRef();  
    volScalarField::Boundary& eVibBf = eVib.boundaryFieldRef();

    forAll(hBf, patchi)
    {
        const fvPatchField<scalar> TVibBfUse =
            use2T ? TVib.boundaryField()[patchi] : TTR.boundaryField()[patchi];

        hBf[patchi] == this->h
        (
            p.boundaryField()[patchi],
            TTR.boundaryField()[patchi],
            patchi
        );

        hBf[patchi].useImplicit(TTR.boundaryField()[patchi].useImplicit());

        scalarField& peTR    = eTRBf[patchi];
        const scalarField& pp    = p.boundaryField()[patchi];
        const scalarField& pTTR  = TTR.boundaryField()[patchi];

        forAll(peTR, facei)
        {
            peTR[facei] =
                this->patchFaceMixture(patchi, facei).EsT
                (
                    pp[facei],
                    pTTR[facei]
                );
        }

        eTRBf[patchi].useImplicit(TTR.boundaryField()[patchi].useImplicit());  

        eVibBf[patchi] == this->eVib
        (
            p.boundaryField()[patchi],
            TTR.boundaryField()[patchi],
            TVibBfUse,
            theta,
            patchi
        );

        eVibBf[patchi].useImplicit(TTR.boundaryField()[patchi].useImplicit());
    }

    this->hBoundaryCorrection(h);
    this->hBoundaryCorrection(eTR);     
    this->eVibBoundaryCorrection(eVib);

    if (p.nOldTimes())
    {
        init
        (
            p.oldTime(),
            TTR.oldTime(),
            TVib.oldTime(),
            h.oldTime(),
            eTR.oldTime(),              
            eT.oldTime(),
            eR.oldTime(),
            eVib.oldTime()
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Basic2TThermo, class MixtureType>
Foam::he2TThermo<Basic2TThermo, MixtureType>::he2TThermo
(
    const fvMesh& mesh,
    const word& phaseName
)
:
    Basic2TThermo(mesh, phaseName),
    MixtureType(*this, mesh, phaseName),

    h_
    (
        IOobject
        (
            Basic2TThermo::phasePropertyName
            (
                MixtureType::thermoType::hName()
            ),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimEnergy/dimMass,
        this->hBoundaryTypes(),
        this->hBoundaryBaseTypes()
    ),
    eTR_
    (
        IOobject
        (
            Basic2TThermo::phasePropertyName
            (
                word("eTR")
            ),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimEnergy/dimMass,
        this->hBoundaryTypes(),
        this->hBoundaryBaseTypes()
    ),
    eT_
    (
        IOobject
        (
            Basic2TThermo::phasePropertyName
            (
                word("eT")
            ),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimEnergy/dimMass,
        this->hBoundaryTypes(),
        this->hBoundaryBaseTypes()
    ),
    eR_
    (
        IOobject
        (
            Basic2TThermo::phasePropertyName
            (
                word("eR")
            ),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimEnergy/dimMass,
        this->hBoundaryTypes(),
        this->hBoundaryBaseTypes()
    ),
    eVib_
    (
        IOobject
        (
            Basic2TThermo::phasePropertyName
            (
                word("eVib")
            ),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimEnergy/dimMass,
        this->eVibBoundaryTypes(),
        this->eVibBoundaryBaseTypes()
    )
{
    init
    (
        this->p_, 
        this->TTR_, 
        this->TVib_, 
        h_, 
        eTR_, 
        eT_, 
        eR_, 
        eVib_
    );
}


template<class Basic2TThermo, class MixtureType>
Foam::he2TThermo<Basic2TThermo, MixtureType>::he2TThermo
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& phaseName
)
:
    Basic2TThermo(mesh, dict, phaseName),
    MixtureType(*this, mesh, phaseName),

    h_
    (
        IOobject
        (
            Basic2TThermo::phasePropertyName
            (
                MixtureType::thermoType::hName()
            ),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimEnergy/dimMass,
        this->hBoundaryTypes(),
        this->hBoundaryBaseTypes()
    ),
    eTR_
    (
        IOobject
        (
            Basic2TThermo::phasePropertyName
            (
                word("eTR")
            ),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimEnergy/dimMass,
        this->hBoundaryTypes(),
        this->hBoundaryBaseTypes()
    ),
    eT_
    (
        IOobject
        (
            Basic2TThermo::phasePropertyName
            (
                word("eT")
            ),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimEnergy/dimMass,
        this->hBoundaryTypes(),
        this->hBoundaryBaseTypes()
    ),
    eR_
    (
        IOobject
        (
            Basic2TThermo::phasePropertyName
            (
                word("eR")
            ),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimEnergy/dimMass,
        this->hBoundaryTypes(),
        this->hBoundaryBaseTypes()
    ),
    eVib_
    (
        IOobject
        (
            Basic2TThermo::phasePropertyName
            (
                word("eVib")
            ),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimEnergy/dimMass,
        this->eVibBoundaryTypes(),
        this->eVibBoundaryBaseTypes()
    )
{
    init
    (
        this->p_, 
        this->TTR_, 
        this->TVib_, 
        h_, 
        eTR_, 
        eT_, 
        eR_, 
        eVib_
    );
}


template<class Basic2TThermo, class MixtureType>
Foam::he2TThermo<Basic2TThermo, MixtureType>::he2TThermo
(
    const fvMesh& mesh,
    const word& phaseName,
    const word& dictionaryName
)
:
    Basic2TThermo(mesh, phaseName, dictionaryName),
    MixtureType(*this, mesh, phaseName),

    h_
    (
        IOobject
        (
            Basic2TThermo::phasePropertyName
            (
                MixtureType::thermoType::hName()
            ),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimEnergy/dimMass,
        this->hBoundaryTypes(),
        this->hBoundaryBaseTypes()
    ),
    eTR_
    (
        IOobject
        (
            Basic2TThermo::phasePropertyName
            (
                word("eTR")
            ),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimEnergy/dimMass,
        this->hBoundaryTypes(),
        this->hBoundaryBaseTypes()
    ),
    eT_
    (
        IOobject
        (
            Basic2TThermo::phasePropertyName
            (
                word("eT")
            ),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimEnergy/dimMass,
        this->hBoundaryTypes(),
        this->hBoundaryBaseTypes()
    ),
    eR_
    (
        IOobject
        (
            Basic2TThermo::phasePropertyName
            (
                word("eR")
            ),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimEnergy/dimMass,
        this->hBoundaryTypes(),
        this->hBoundaryBaseTypes()
    ),
    eVib_
    (
        IOobject
        (
            Basic2TThermo::phasePropertyName
            (
                word("eVib")
            ),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimEnergy/dimMass,
        this->eVibBoundaryTypes(),
        this->eVibBoundaryBaseTypes()
    )
{
    init
    (
        this->p_, 
        this->TTR_, 
        this->TVib_, 
        h_, 
        eTR_, 
        eT_, 
        eR_, 
        eVib_
    );
}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Basic2TThermo, class MixtureType>
Foam::he2TThermo<Basic2TThermo, MixtureType>::~he2TThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Basic2TThermo, class MixtureType>
Foam::tmp<Foam::volScalarField> Foam::he2TThermo<Basic2TThermo, MixtureType>::h
(
    const volScalarField& p,
    const volScalarField& TTR
) const
{
    const fvMesh& mesh = this->TTR_.mesh();

    auto the = volScalarField::New
    (
        "h",
        IOobject::NO_REGISTER,
        mesh,
        h_.dimensions()
    );
    auto& h = the.ref();

    scalarField& hCells = h.primitiveFieldRef();
    const scalarField& pCells = p;
    const scalarField& TTRCells = TTR;

    forAll(hCells, celli)
    {
        hCells[celli] =
            this->cellMixture(celli).H(pCells[celli], TTRCells[celli]);
    }

    volScalarField::Boundary& hBf = h.boundaryFieldRef();

    forAll(hBf, patchi)
    {
        scalarField& hp = hBf[patchi];
        const scalarField& pp = p.boundaryField()[patchi];
        const scalarField& TTRp = TTR.boundaryField()[patchi];

        forAll(hp, facei)
        {
            hp[facei] =
                this->patchFaceMixture(patchi, facei).H(pp[facei], TTRp[facei]);
        }
    }

    return the;
}


template<class Basic2TThermo, class MixtureType>
Foam::tmp<Foam::scalarField> Foam::he2TThermo<Basic2TThermo, MixtureType>::h
(
    const scalarField& p,
    const scalarField& TTR,
    const labelList& cells
) const
{
    auto the = tmp<scalarField>::New(TTR.size());
    auto& h = the.ref();

    forAll(TTR, celli)
    {
        h[celli] = this->cellMixture(cells[celli]).H(p[celli], TTR[celli]);
    }

    return the;
}


template<class Basic2TThermo, class MixtureType>
Foam::tmp<Foam::scalarField> Foam::he2TThermo<Basic2TThermo, MixtureType>::h
(
    const scalarField& p,
    const scalarField& TTR,
    const label patchi
) const
{
    auto the = tmp<scalarField>::New(TTR.size());
    auto& h = the.ref();

    forAll(TTR, facei)
    {
        h[facei] =
            this->patchFaceMixture(patchi, facei).H(p[facei], TTR[facei]);
    }

    return the;
}

template<class Basic2TThermo, class MixtureType>
Foam::tmp<Foam::volScalarField> Foam::he2TThermo<Basic2TThermo, MixtureType>::eTR
(
    const volScalarField& p,
    const volScalarField& TTR
) const
{
    const fvMesh& mesh = this->TTR_.mesh();

    auto test = volScalarField::New
    (
        "eT",
        IOobject::NO_REGISTER,
        mesh,
        eT_.dimensions()
    );
    auto& eT = test.ref();

    scalarField& eTCells = eT.primitiveFieldRef();
    const scalarField& pCells = p;
    const scalarField& TTRCells = TTR;

    forAll(eTCells, celli)
    {
        eTCells[celli] =
            this->cellMixture(celli).ET(pCells[celli], TTRCells[celli]);
    }

    volScalarField::Boundary& eTBf = eT.boundaryFieldRef();

    forAll(eTBf, patchi)
    {
        scalarField& esTp = eTBf[patchi];
        const scalarField& pp = p.boundaryField()[patchi];
        const scalarField& TTRp = TTR.boundaryField()[patchi];

        forAll(esTp, facei)
        {
            esTp[facei] =
                this->patchFaceMixture(patchi, facei).ET(pp[facei], TTRp[facei]);
        }
    }

    return test;
}


template<class Basic2TThermo, class MixtureType>
Foam::tmp<Foam::scalarField> Foam::he2TThermo<Basic2TThermo, MixtureType>::eTR
(
    const scalarField& p,
    const scalarField& TTR,
    const labelList& cells
) const
{
    auto test = tmp<scalarField>::New(TTR.size());
    auto& eT = test.ref();

    forAll(TTR, celli)
    {
        eT[celli] = this->cellMixture(cells[celli]).ET(p[celli], TTR[celli]);
    }

    return test;
}


template<class Basic2TThermo, class MixtureType>
Foam::tmp<Foam::scalarField> Foam::he2TThermo<Basic2TThermo, MixtureType>::eTR
(
    const scalarField& p,
    const scalarField& TTR,
    const label patchi
) const
{
    auto test = tmp<scalarField>::New(TTR.size());
    auto& eT = test.ref();

    forAll(TTR, facei)
    {
        eT[facei] =
            this->patchFaceMixture(patchi, facei).ET(p[facei], TTR[facei]);
    }

    return test;
}

template<class Basic2TThermo, class MixtureType>
Foam::tmp<Foam::volScalarField> Foam::he2TThermo<Basic2TThermo, MixtureType>::eT
(
    const volScalarField& p,
    const volScalarField& TTR
) const
{
    const fvMesh& mesh = this->TTR_.mesh();

    auto test = volScalarField::New
    (
        "eT",
        IOobject::NO_REGISTER,
        mesh,
        eT_.dimensions()
    );
    auto& eT = test.ref();

    scalarField& eTCells = eT.primitiveFieldRef();
    const scalarField& pCells = p;
    const scalarField& TTRCells = TTR;

    forAll(eTCells, celli)
    {
        eTCells[celli] =
            this->cellMixture(celli).ET(pCells[celli], TTRCells[celli]);
    }

    volScalarField::Boundary& eTBf = eT.boundaryFieldRef();

    forAll(eTBf, patchi)
    {
        scalarField& esTp = eTBf[patchi];
        const scalarField& pp = p.boundaryField()[patchi];
        const scalarField& TTRp = TTR.boundaryField()[patchi];

        forAll(esTp, facei)
        {
            esTp[facei] =
                this->patchFaceMixture(patchi, facei).ET(pp[facei], TTRp[facei]);
        }
    }

    return test;
}


template<class Basic2TThermo, class MixtureType>
Foam::tmp<Foam::scalarField> Foam::he2TThermo<Basic2TThermo, MixtureType>::eT
(
    const scalarField& p,
    const scalarField& TTR,
    const labelList& cells
) const
{
    auto test = tmp<scalarField>::New(TTR.size());
    auto& eT = test.ref();

    forAll(TTR, celli)
    {
        eT[celli] = this->cellMixture(cells[celli]).ET(p[celli], TTR[celli]);
    }

    return test;
}


template<class Basic2TThermo, class MixtureType>
Foam::tmp<Foam::scalarField> Foam::he2TThermo<Basic2TThermo, MixtureType>::eT
(
    const scalarField& p,
    const scalarField& TTR,
    const label patchi
) const
{
    auto test = tmp<scalarField>::New(TTR.size());
    auto& eT = test.ref();

    forAll(TTR, facei)
    {
        eT[facei] =
            this->patchFaceMixture(patchi, facei).ET(p[facei], TTR[facei]);
    }

    return test;
}

template<class Basic2TThermo, class MixtureType>
Foam::tmp<Foam::volScalarField> Foam::he2TThermo<Basic2TThermo, MixtureType>::eR
(
    const volScalarField& p,
    const volScalarField& TTR
) const
{
    const fvMesh& mesh = this->TTR_.mesh();

    auto tesr = volScalarField::New
    (
        "eR",
        IOobject::NO_REGISTER,
        mesh,
        eR_.dimensions()
    );
    auto& eR = tesr.ref();

    scalarField& eRCells = eR.primitiveFieldRef();
    const scalarField& pCells = p;
    const scalarField& TTRCells = TTR;

    forAll(eRCells, celli)
    {
        eRCells[celli] =
            this->cellMixture(celli).ER(pCells[celli], TTRCells[celli]);
    }

    volScalarField::Boundary& eRBf = eR.boundaryFieldRef();

    forAll(eRBf, patchi)
    {
        scalarField& esRp = eRBf[patchi];
        const scalarField& pp = p.boundaryField()[patchi];
        const scalarField& TTRp = TTR.boundaryField()[patchi];

        forAll(esRp, facei)
        {
            esRp[facei] =
                this->patchFaceMixture(patchi, facei).ER(pp[facei], TTRp[facei]);
        }
    }

    return tesr;
}


template<class Basic2TThermo, class MixtureType>
Foam::tmp<Foam::scalarField> Foam::he2TThermo<Basic2TThermo, MixtureType>::eR
(
    const scalarField& p,
    const scalarField& TTR,
    const labelList& cells
) const
{
    auto tesr = tmp<scalarField>::New(TTR.size());
    auto& eR = tesr.ref();

    forAll(TTR, celli)
    {
        eR[celli] = this->cellMixture(cells[celli]).ER(p[celli], TTR[celli]);
    }

    return tesr;
}


template<class Basic2TThermo, class MixtureType>
Foam::tmp<Foam::scalarField> Foam::he2TThermo<Basic2TThermo, MixtureType>::eR
(
    const scalarField& p,
    const scalarField& TTR,
    const label patchi
) const
{
    auto tesr = tmp<scalarField>::New(TTR.size());
    auto& eR = tesr.ref();

    forAll(TTR, facei)
    {
        eR[facei] =
            this->patchFaceMixture(patchi, facei).ER(p[facei], TTR[facei]);
    }

    return tesr;
}

template<class Basic2TThermo, class MixtureType>
Foam::tmp<Foam::volScalarField> Foam::he2TThermo<Basic2TThermo, MixtureType>::eVib
(
    const volScalarField& p,
    const volScalarField& TTR,
    const volScalarField& TVib,
    const scalar ThetaVib
) const
{
    const bool use2T = this->twoTemperature();
    const fvMesh& mesh = this->TTR_.mesh();

    auto tesvib = volScalarField::New
    (
        "eVib",
        IOobject::NO_REGISTER,
        mesh,
        eVib_.dimensions()
    );
    auto& eVib = tesvib.ref();

    scalarField& eVibCells = eVib.primitiveFieldRef();
    const scalarField& pCells = p;
    const scalarField& TTRCells = TTR;
    const scalarField& TVibCells = TVib;

    const PtrList<volScalarField>& Yvf = this->Y();
    const label nSpecie = Yvf.size();

    forAll(eVibCells, celli)
    {
        const scalar TVibUse = use2T ? TVibCells[celli] : TTRCells[celli];

        scalar evMix = 0;
        for (label speciei = 0; speciei < nSpecie; ++speciei)
        {
            evMix +=
                Yvf[speciei][celli]
              * this->EsVib
                (
                    speciei, pCells[celli], TTRCells[celli],
                    TVibUse, this->composition().thetaVib(speciei)
                );
        }
        eVibCells[celli] = evMix;
    }

    volScalarField::Boundary& eVibBf = eVib.boundaryFieldRef();

    forAll(eVibBf, patchi)
    {
        scalarField& esVibp = eVibBf[patchi];
        const scalarField& pp = p.boundaryField()[patchi];
        const scalarField& TTRp = TTR.boundaryField()[patchi];
        const scalarField& TVibp = TVib.boundaryField()[patchi];

        forAll(esVibp, facei)
        {
            const scalar TVibUse = use2T ? TVibp[facei] : TTRp[facei];

            scalar evMix = 0;
            for (label speciei = 0; speciei < nSpecie; ++speciei)
            {
                evMix +=
                    Yvf[speciei].boundaryField()[patchi][facei]
                  * this->EsVib
                    (
                        speciei, pp[facei], TTRp[facei],
                        TVibUse, this->composition().thetaVib(speciei)
                    );
            }
            esVibp[facei] = evMix;
        }
    }

    return tesvib;
}


template<class Basic2TThermo, class MixtureType>
Foam::tmp<Foam::scalarField> Foam::he2TThermo<Basic2TThermo, MixtureType>::eVib
(
    const scalarField& p,
    const scalarField& TTR,                
    const scalarField& TVib,
    const scalar ThetaVib,
    const labelList& cells
) const
{
    const bool use2T = this->twoTemperature();
    auto tesvib = tmp<scalarField>::New(TTR.size());
    auto& eVib = tesvib.ref();

    const PtrList<volScalarField>& Yvf = this->Y();
    const label nSpecie = Yvf.size();

    forAll(TTR, celli)
    {
        const scalar TVibUse = use2T ? TVib[celli] : TTR[celli];
        const label gCelli = cells[celli];

        scalar evMix = 0;
        for (label speciei = 0; speciei < nSpecie; ++speciei)
        {
            evMix +=
                Yvf[speciei][gCelli]
              * this->EsVib
                (
                    speciei, p[celli], TTR[celli],
                    TVibUse, this->composition().thetaVib(speciei)
                );
        }
        eVib[celli] = evMix;
    }

    return tesvib;
}


template<class Basic2TThermo, class MixtureType>
Foam::tmp<Foam::scalarField> Foam::he2TThermo<Basic2TThermo, MixtureType>::eVib
(
    const scalarField& p,
    const scalarField& TTR,          
    const scalarField& TVib,
    const scalar ThetaVib,
    const label patchi
) const
{
    const bool use2T = this->twoTemperature();
    auto tesvib = tmp<scalarField>::New(TTR.size());
    auto& eVib = tesvib.ref();

    const PtrList<volScalarField>& Yvf = this->Y();
    const label nSpecie = Yvf.size();

    forAll(TTR, facei)
    {
        const scalar TVibUse = use2T ? TVib[facei] : TTR[facei];

        scalar evMix = 0;
        for (label speciei = 0; speciei < nSpecie; ++speciei)
        {
            evMix +=
                Yvf[speciei].boundaryField()[patchi][facei]
              * this->EsVib
                (
                    speciei, p[facei], TTR[facei],
                    TVibUse, this->composition().thetaVib(speciei)
                );
        }
        eVib[facei] = evMix;
    }

    return tesvib;
}

template<class Basic2TThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::he2TThermo<Basic2TThermo, MixtureType>::hc() const
{
    const fvMesh& mesh = this->TTR_.mesh();

    auto thc = volScalarField::New
    (
        "hc",
        IOobject::NO_REGISTER,
        mesh,
        h_.dimensions()
    );
    auto& hcf = thc.ref();

    scalarField& hcCells = hcf.primitiveFieldRef();

    forAll(hcCells, celli)
    {
        hcCells[celli] = this->cellMixture(celli).Hc();
    }

    volScalarField::Boundary& hcfBf = hcf.boundaryFieldRef();

    forAll(hcfBf, patchi)
    {
        scalarField& hcp = hcfBf[patchi];

        forAll(hcp, facei)
        {
            hcp[facei] = this->patchFaceMixture(patchi, facei).Hc();
        }
    }

    return thc;
}


template<class Basic2TThermo, class MixtureType>
Foam::tmp<Foam::scalarField> Foam::he2TThermo<Basic2TThermo, MixtureType>::CpTR
(
    const scalarField& p,
    const scalarField& TTR,
    const label patchi
) const
{
    auto tCpTR = tmp<scalarField>::New(TTR.size());
    auto& cpTR = tCpTR.ref();

    forAll(TTR, facei)
    {
        cpTR[facei] =
            this->patchFaceMixture(patchi, facei).CpTR(p[facei], TTR[facei]);
    }

    return tCpTR;
}


template<class Basic2TThermo, class MixtureType>
Foam::tmp<Foam::scalarField>
Foam::he2TThermo<Basic2TThermo, MixtureType>::CpTR
(
    const scalarField& p,
    const scalarField& TTR,
    const labelList& cells
) const
{
    auto tCpTR = tmp<scalarField>::New(TTR.size());
    auto& CpTR = tCpTR.ref();

    forAll(cells, i)
    {
        const label celli = cells[i];
        CpTR[i] = this->cellMixture(celli).CpTR(p[i], TTR[i]);
    }

    return tCpTR;
}


template<class Basic2TThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::he2TThermo<Basic2TThermo, MixtureType>::CpTR() const
{
    const fvMesh& mesh = this->TTR_.mesh();

    auto tCpTR = volScalarField::New
    (
        "CpTR",
        IOobject::NO_REGISTER,
        mesh,
        dimEnergy/dimMass/dimTemperature
    );
    auto& cpTR = tCpTR.ref();

    forAll(this->TTR_, celli)
    {
        cpTR[celli] =
            this->cellMixture(celli).CpTR(this->p_[celli], this->TTR_[celli]);
    }

    volScalarField::Boundary& cpTRBf = cpTR.boundaryFieldRef();

    forAll(cpTRBf, patchi)
    {
        const fvPatchScalarField& pp = this->p_.boundaryField()[patchi];
        const fvPatchScalarField& pTTR = this->TTR_.boundaryField()[patchi];
        fvPatchScalarField& pCpTR = cpTRBf[patchi];

        forAll(pTTR, facei)
        {
            pCpTR[facei] =
                this->patchFaceMixture(patchi, facei).CpTR(pp[facei], pTTR[facei]);
        }
    }

    return tCpTR;
}

template<class Basic2TThermo, class MixtureType>
Foam::tmp<Foam::scalarField> Foam::he2TThermo<Basic2TThermo, MixtureType>::CvTR
(
    const scalarField& p,
    const scalarField& TTR,
    const label patchi
) const
{
    auto tCvTR = tmp<scalarField>::New(TTR.size());
    auto& cvTR = tCvTR.ref();

    forAll(TTR, facei)
    {
        cvTR[facei] =
            this->patchFaceMixture(patchi, facei).CvT(p[facei], TTR[facei]);
    }

    return tCvTR;
}


template<class Basic2TThermo, class MixtureType>
Foam::tmp<Foam::scalarField>
Foam::he2TThermo<Basic2TThermo, MixtureType>::CvTR
(
    const scalarField& p,
    const scalarField& TTR,
    const labelList& cells
) const
{
    auto tCvTR = tmp<scalarField>::New(TTR.size());
    auto& CvTR = tCvTR.ref();

    forAll(cells, i)
    {
        const label celli = cells[i];
        CvTR[i] = this->cellMixture(celli).CvT(p[i], TTR[i]);
    }

    return tCvTR;
}


template<class Basic2TThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::he2TThermo<Basic2TThermo, MixtureType>::CvTR() const
{
    const fvMesh& mesh = this->TTR_.mesh();

    auto tCvTR = volScalarField::New
    (
        "CvTR",
        IOobject::NO_REGISTER,
        mesh,
        dimEnergy/dimMass/dimTemperature
    );
    auto& cvTR = tCvTR.ref();

    forAll(this->TTR_, celli)
    {
        cvTR[celli] =
            this->cellMixture(celli).CvT(this->p_[celli], this->TTR_[celli]);
    }

    volScalarField::Boundary& cvTRBf = cvTR.boundaryFieldRef();

    forAll(cvTRBf, patchi)
    {
        const fvPatchScalarField& pp = this->p_.boundaryField()[patchi];
        const fvPatchScalarField& pTTR = this->TTR_.boundaryField()[patchi];
        fvPatchScalarField& pCvTR = cvTRBf[patchi];

        forAll(pTTR, facei)
        {
            pCvTR[facei] =
                this->patchFaceMixture(patchi, facei).CvT(pp[facei], pTTR[facei]);
        }
    }

    return tCvTR;
}


template<class Basic2TThermo, class MixtureType>
Foam::tmp<Foam::scalarField>
Foam::he2TThermo<Basic2TThermo, MixtureType>::CvT
(
    const scalarField& p,
    const scalarField& TTR,
    const label patchi
) const
{
    auto tCvT = tmp<scalarField>::New(TTR.size());
    auto& cvT = tCvT.ref();

    forAll(TTR, facei)
    {
        cvT[facei] =
            this->patchFaceMixture(patchi, facei).CvT(p[facei], TTR[facei]);
    }

    return tCvT;
}

template<class Basic2TThermo, class MixtureType>
Foam::tmp<Foam::scalarField>
Foam::he2TThermo<Basic2TThermo, MixtureType>::CvR
(
    const scalarField& p,
    const scalarField& TTR,
    const label patchi
) const
{
    auto tCvR = tmp<scalarField>::New(TTR.size());
    auto& cvR = tCvR.ref();

    forAll(TTR, facei)
    {
        cvR[facei] =
            this->patchFaceMixture(patchi, facei).CvR(p[facei], TTR[facei]);
    }

    return tCvR;
}

template<class Basic2TThermo, class MixtureType>
Foam::tmp<Foam::scalarField>
Foam::he2TThermo<Basic2TThermo, MixtureType>::CvVib
(
    const scalarField& p,
    const scalarField& TTR,
    const scalarField& TVib,
    const scalar ThetaVib,
    const label patchi
) const
{
    auto tCvVib = tmp<scalarField>::New(TVib.size());
    auto& cvVib = tCvVib.ref();

    forAll(TVib, facei)
    {
        cvVib[facei] =
            this->patchFaceMixture(patchi, facei).CvVib(p[facei], TTR[facei], TVib[facei], ThetaVib);
    }

    return tCvVib;
}


template<class Basic2TThermo, class MixtureType>
Foam::tmp<Foam::scalarField>
Foam::he2TThermo<Basic2TThermo, MixtureType>::rhoEoS
(
    const scalarField& p,
    const scalarField& TTR,
    const labelList& cells
) const
{
    auto tRho = tmp<scalarField>::New(TTR.size());
    auto& rho = tRho.ref();

    forAll(cells, i)
    {
        const label celli = cells[i];
        rho[i] = this->cellMixture(celli).rho(p[i], TTR[i]);
    }

    return tRho;
}


template<class Basic2TThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::he2TThermo<Basic2TThermo, MixtureType>::CvT() const
{
    const fvMesh& mesh = this->TTR_.mesh();

    auto tCvT = volScalarField::New
    (
        "CvT",
        IOobject::NO_REGISTER,
        mesh,
        dimEnergy/dimMass/dimTemperature
    );
    auto& cvT = tCvT.ref();

    forAll(this->TTR_, celli)
    {
        cvT[celli] =
            this->cellMixture(celli).CvT(this->p_[celli], this->TTR_[celli]);
    }

    volScalarField::Boundary& cvTBf = cvT.boundaryFieldRef();

    forAll(cvTBf, patchi)
    {
        cvTBf[patchi] = CvT
        (
            this->p_.boundaryField()[patchi],
            this->TTR_.boundaryField()[patchi],
            patchi
        );
    }

    return tCvT;
}

template<class Basic2TThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::he2TThermo<Basic2TThermo, MixtureType>::CvR() const
{
    const fvMesh& mesh = this->TTR_.mesh();

    auto tCvR = volScalarField::New
    (
        "CvR",
        IOobject::NO_REGISTER,
        mesh,
        dimEnergy/dimMass/dimTemperature
    );
    auto& cvR = tCvR.ref();

    forAll(this->TTR_, celli)
    {
        cvR[celli] =
            this->cellMixture(celli).CvR(this->p_[celli], this->TTR_[celli]);
    }

    volScalarField::Boundary& cvRBf = cvR.boundaryFieldRef();

    forAll(cvRBf, patchi)
    {
        cvRBf[patchi] = CvR
        (
            this->p_.boundaryField()[patchi],
            this->TTR_.boundaryField()[patchi],
            patchi
        );
    }

    return tCvR;
}

template<class Basic2TThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::he2TThermo<Basic2TThermo, MixtureType>::CvVib() const
{
    const fvMesh& mesh = this->TVib_.mesh();
    const scalar theta = this->cellMixture(0).ThetaVib();  // getter on Thermo2T

    auto tCvVib = volScalarField::New
    (
        "CvVib",
        IOobject::NO_REGISTER,
        mesh,
        dimEnergy/dimMass/dimTemperature
    );
    auto& cvVib = tCvVib.ref();

    forAll(this->TVib_, celli)
    {
        cvVib[celli] =
            this->cellMixture(celli).CvVib(this->p_[celli], this->TTR_[celli] ,this->TVib_[celli], theta);
    }

    volScalarField::Boundary& cvVibBf = cvVib.boundaryFieldRef();

    forAll(cvVibBf, patchi)
    {
        cvVibBf[patchi] = CvVib
        (
            this->p_.boundaryField()[patchi],
            this->TTR_.boundaryField()[patchi],
            this->TVib_.boundaryField()[patchi],
            theta,
            patchi
        );
    }

    return tCvVib;
}


template<class Basic2TThermo, class MixtureType>
Foam::tmp<Foam::scalarField> Foam::he2TThermo<Basic2TThermo, MixtureType>::gamma
(
    const scalarField& p,
    const scalarField& TTR,
    const label patchi
) const
{
    auto tgamma = tmp<scalarField>::New(TTR.size());
    auto& gamma = tgamma.ref();

    forAll(TTR, facei)
    {
        gamma[facei] =
            this->patchFaceMixture(patchi, facei).gamma(p[facei], TTR[facei]);
    }

    return tgamma;
}


template<class Basic2TThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::he2TThermo<Basic2TThermo, MixtureType>::gamma() const
{
    const fvMesh& mesh = this->TTR_.mesh();

    auto tgamma = volScalarField::New
    (
        "gamma",
        IOobject::NO_REGISTER,
        mesh,
        dimless
    );
    auto& gamma = tgamma.ref();

    forAll(this->TTR_, celli)
    {
        gamma[celli] =
            this->cellMixture(celli).gamma(this->p_[celli], this->TTR_[celli]);
    }

    volScalarField::Boundary& gammaBf = gamma.boundaryFieldRef();

    forAll(gammaBf, patchi)
    {
        const fvPatchScalarField& pp = this->p_.boundaryField()[patchi];
        const fvPatchScalarField& pTTR = this->TTR_.boundaryField()[patchi];
        fvPatchScalarField& pgamma = gammaBf[patchi];

        forAll(pTTR, facei)
        {
            pgamma[facei] = this->patchFaceMixture(patchi, facei).gamma
            (
                pp[facei],
                pTTR[facei]
            );
        }
    }

    return tgamma;
}


template<class Basic2TThermo, class MixtureType>
Foam::tmp<Foam::scalarField> Foam::he2TThermo<Basic2TThermo, MixtureType>::CpvTR
(
    const scalarField& p,
    const scalarField& TTR,
    const label patchi
) const
{
    auto tCpvTR = tmp<scalarField>::New(TTR.size());
    auto& CpvTR = tCpvTR.ref();

    forAll(TTR, facei)
    {
        CpvTR[facei] =
            this->patchFaceMixture(patchi, facei).CpvTR(p[facei], TTR[facei]);
    }

    return tCpvTR;
}


template<class Basic2TThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::he2TThermo<Basic2TThermo, MixtureType>::CpvTR() const
{
    const fvMesh& mesh = this->TTR_.mesh();

    auto tCpvTR = volScalarField::New
    (
        "CpvTR",
        IOobject::NO_REGISTER,
        mesh,
        dimEnergy/dimMass/dimTemperature
    );
    auto& CpvTR = tCpvTR.ref();

    forAll(this->TTR_, celli)
    {
        CpvTR[celli] =
            this->cellMixture(celli).CpvTR(this->p_[celli], this->TTR_[celli]);
    }

    volScalarField::Boundary& CpvTRBf = CpvTR.boundaryFieldRef();

    forAll(CpvTRBf, patchi)
    {
        const fvPatchScalarField& pp = this->p_.boundaryField()[patchi];
        const fvPatchScalarField& pTTR = this->TTR_.boundaryField()[patchi];
        fvPatchScalarField& pCpvTR = CpvTRBf[patchi];

        forAll(pTTR, facei)
        {
            pCpvTR[facei] =
                this->patchFaceMixture(patchi, facei).CpvTR(pp[facei], pTTR[facei]);
        }
    }

    return tCpvTR;
}


template<class Basic2TThermo, class MixtureType>
Foam::tmp<Foam::scalarField> Foam::he2TThermo<Basic2TThermo, MixtureType>::CpByCpvTR
(
    const scalarField& p,
    const scalarField& TTR,
    const label patchi
) const
{
    auto tCpByCpvTR = tmp<scalarField>::New(TTR.size());
    auto& CpByCpvTR = tCpByCpvTR.ref();

    forAll(TTR, facei)
    {
        CpByCpvTR[facei] =
            this->patchFaceMixture(patchi, facei).CpByCpvTR(p[facei], TTR[facei]);
    }

    return tCpByCpvTR;
}


template<class Basic2TThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::he2TThermo<Basic2TThermo, MixtureType>::CpByCpvTR() const
{
    const fvMesh& mesh = this->TTR_.mesh();

    auto tCpByCpvTR = volScalarField::New
    (
        "CpByCpvTR",
        IOobject::NO_REGISTER,
        mesh,
        dimless
    );
    auto& CpByCpvTR = tCpByCpvTR.ref();

    forAll(this->TTR_, celli)
    {
        CpByCpvTR[celli] = this->cellMixture(celli).CpByCpvTR
        (
            this->p_[celli],
            this->TTR_[celli]
        );
    }

    volScalarField::Boundary& CpByCpvTRBf =
        CpByCpvTR.boundaryFieldRef();

    forAll(CpByCpvTRBf, patchi)
    {
        const fvPatchScalarField& pp = this->p_.boundaryField()[patchi];
        const fvPatchScalarField& pTTR = this->TTR_.boundaryField()[patchi];
        fvPatchScalarField& pCpByCpvTR = CpByCpvTRBf[patchi];

        forAll(pTTR, facei)
        {
            pCpByCpvTR[facei] = this->patchFaceMixture(patchi, facei).CpByCpvTR
            (
                pp[facei],
                pTTR[facei]
            );
        }
    }

    return tCpByCpvTR;
}


template<class Basic2TThermo, class MixtureType>
Foam::tmp<Foam::scalarField> Foam::he2TThermo<Basic2TThermo, MixtureType>::TE_TR
(
    const scalarField& eTR,
    const scalarField& p,
    const scalarField& TTR0,
    const labelList& cells
) const
{
    auto tTTR = tmp<scalarField>::New(eTR.size());
    auto& TTR = tTTR.ref();

    forAll(eTR, celli)
    {
        TTR[celli] =
            this->cellMixture(cells[celli]).TE_TR(eTR[celli], p[celli], TTR0[celli]);
    }

    return tTTR;
}


template<class Basic2TThermo, class MixtureType>
Foam::tmp<Foam::scalarField> Foam::he2TThermo<Basic2TThermo, MixtureType>::TE_TR
(
    const scalarField& eTR,
    const scalarField& p,
    const scalarField& TTR0,
    const label patchi
) const
{

    auto tTTR = tmp<scalarField>::New(eTR.size());
    auto& TTR = tTTR.ref();

    forAll(eTR, facei)
    {
        TTR[facei] = this->patchFaceMixture
        (
            patchi,
            facei
        ).TE_TR(eTR[facei], p[facei], TTR0[facei]);
    }

    return tTTR;
}

template<class Basic2TThermo, class MixtureType>
Foam::tmp<Foam::scalarField> Foam::he2TThermo<Basic2TThermo, MixtureType>::TH_TR
(
    const scalarField& h,
    const scalarField& p,
    const scalarField& TTR0,
    const labelList& cells
) const
{
    auto tTTR = tmp<scalarField>::New(h.size());
    auto& TTR = tTTR.ref();

    forAll(h, celli)
    {
        TTR[celli] =
            this->cellMixture(cells[celli]).TH_TR(h[celli], p[celli], TTR0[celli]);
    }

    return tTTR;
}


template<class Basic2TThermo, class MixtureType>
Foam::tmp<Foam::scalarField> Foam::he2TThermo<Basic2TThermo, MixtureType>::TH_TR
(
    const scalarField& h,
    const scalarField& p,
    const scalarField& TTR0,
    const label patchi
) const
{

    auto tTTR = tmp<scalarField>::New(h.size());
    auto& TTR = tTTR.ref();

    forAll(h, facei)
    {
        TTR[facei] = this->patchFaceMixture
        (
            patchi,
            facei
        ).TH_TR(h[facei], p[facei], TTR0[facei]);
    }

    return tTTR;
}

template<class Basic2TThermo, class MixtureType>
Foam::tmp<Foam::scalarField> Foam::he2TThermo<Basic2TThermo, MixtureType>::TE_Vib
(
    const scalarField& h,
    const scalarField& p,
    const scalarField& TTR0,
    const scalarField& TVib0,
    const labelList& cells
) const
{

    auto tTVib = tmp<scalarField>::New(h.size());
    auto& TVib = tTVib.ref();

    forAll(h, celli)
    {
        TVib[celli] =
            this->cellMixture(cells[celli]).TE_Vib(h[celli], p[celli], TTR0[celli], TVib0[celli]);
    }

    return tTVib;
}


template<class Basic2TThermo, class MixtureType>
Foam::tmp<Foam::scalarField> Foam::he2TThermo<Basic2TThermo, MixtureType>::TE_Vib
(
    const scalarField& h,
    const scalarField& p,
    const scalarField& TTR0,
    const scalarField& TVib0,
    const label patchi
) const
{

    auto tTVib = tmp<scalarField>::New(h.size());
    auto& TVib = tTVib.ref();

    forAll(h, facei)
    {
        TVib[facei] = this->patchFaceMixture
        (
            patchi,
            facei
        ).TE_Vib(h[facei], p[facei], TTR0[facei], TVib0[facei]);
    }

    return tTVib;
}


template<class Basic2TThermo, class MixtureType>
Foam::tmp<Foam::volScalarField> Foam::he2TThermo<Basic2TThermo, MixtureType>::W
(
) const
{
    const fvMesh& mesh = this->TTR_.mesh();

    auto tW = volScalarField::New
    (
        "W",
        IOobject::NO_REGISTER,
        mesh,
        dimMass/dimMoles
    );
    auto& W = tW.ref();

    scalarField& WCells = W.primitiveFieldRef();

    forAll(WCells, celli)
    {
        WCells[celli] = this->cellMixture(celli).W();
    }

    auto& WBf = W.boundaryFieldRef();

    forAll(WBf, patchi)
    {
        scalarField& Wp = WBf[patchi];
        forAll(Wp, facei)
        {
            Wp[facei] = this->patchFaceMixture(patchi, facei).W();
        }
    }

    return tW;
}


template<class Basic2TThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::he2TThermo<Basic2TThermo, MixtureType>::kappaTR() const
{
    auto tK = Foam::tmp<Foam::volScalarField>
    (
        new Foam::volScalarField
        (
            Foam::IOobject
            (
                "kappaTR",
                this->p_.time().timeName(),
                this->p_.mesh(),
                Foam::IOobject::NO_READ,
                Foam::IOobject::NO_WRITE
            ),
            this->p_.mesh(),
            Foam::dimensionedScalar
            (
                "zero",
                dimEnergy/(dimTime*dimLength*dimTemperature),
                0.0
            )
        )
    );

    auto& K = tK.ref();

    // Internal field
    {
        const Foam::scalarField& pi = this->p_.primitiveField();
        const Foam::scalarField& Ti = this->TTR_.primitiveField();
        Foam::scalarField& Ki = K.primitiveFieldRef();

        forAll(Ki, i)
        {
            const typename MixtureType::thermoType& mixture_ =
                this->cellMixture(i);

            Ki[i] = mixture_.kappaTR(pi[i], Ti[i]);
        }
    }

    // Patch fields
    forAll(K.boundaryField(), patchi)
    {
        K.boundaryFieldRef()[patchi] = this->kappaTR(patchi);
    }



    K.rename("kappaTR");
    return tK;
}


template<class Basic2TThermo, class MixtureType>
Foam::tmp<Foam::scalarField>
Foam::he2TThermo<Basic2TThermo, MixtureType>::kappaTR
(
    const label patchi
) const
{
    const Foam::scalarField& pp = this->p_.boundaryField()[patchi];
    const Foam::scalarField& Tp = this->TTR_.boundaryField()[patchi];

    Foam::tmp<Foam::scalarField> tKp(new Foam::scalarField(Tp.size()));
    Foam::scalarField& Kp = tKp.ref();

    forAll(Kp, facei)
    {
        const typename MixtureType::thermoType& mixture_ =
            this->patchFaceMixture(patchi, facei);

        Kp[facei] = mixture_.kappaTR(pp[facei], Tp[facei]);
    }

    return tKp;
}


template<class Basic2TThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::he2TThermo<Basic2TThermo, MixtureType>::kappaVib() const
{
    auto tK = Foam::tmp<Foam::volScalarField>
    (
        new Foam::volScalarField
        (
            Foam::IOobject
            (
                "kappaVib",
                this->p_.time().timeName(),
                this->p_.mesh(),
                Foam::IOobject::NO_READ,
                Foam::IOobject::NO_WRITE
            ),
            this->p_.mesh(),
            // W/(m·K)
            Foam::dimensionedScalar
            (
                "zero",
                dimEnergy/(dimTime*dimLength*dimTemperature),
                0.0
            )
        )
    );

    auto& K = tK.ref();

    // internal field
    {
        const Foam::scalarField& pi = this->p_.primitiveField();
        const Foam::scalarField& Ti = this->TTR_.primitiveField();
        const Foam::scalarField& TVi = this->TVib_.primitiveField();
        Foam::scalarField&       Ki  = K.primitiveFieldRef();
        const PtrList<volScalarField>& Yvf = this->Y();
        const label nSpecie = Yvf.size();

        forAll(Ki, i)
        {
            
            scalar kvMix = 0;
            for (label speciei = 0; speciei < nSpecie; ++speciei)
            {
                kvMix +=
                    max(Yvf[speciei][i], scalar(0))
                  * this->composition().kappaVib
                    (
                        speciei, pi[i], Ti[i], TVi[i],
                        this->composition().thetaVib(speciei)
                    );
            }
            Ki[i] = kvMix;
        }
    }

    
    // patch fields
    forAll(K.boundaryField(), patchi)
    {
        K.boundaryFieldRef()[patchi] = this->kappaVib(patchi);
    }

    K.rename("kappaVib");
    return tK;
}


template<class Basic2TThermo, class MixtureType>
Foam::tmp<Foam::scalarField>
Foam::he2TThermo<Basic2TThermo, MixtureType>::kappaVib
(
    const label patchi
) const
{
    const Foam::scalarField& pp  = this->p_.boundaryField()[patchi];
    const Foam::scalarField& Tp  = this->TTR_.boundaryField()[patchi];
    const Foam::scalarField& TVp = this->TVib_.boundaryField()[patchi];

    Foam::tmp<Foam::scalarField> tKp(new Foam::scalarField(Tp.size()));
    Foam::scalarField& Kp = tKp.ref();

    const PtrList<volScalarField>& Yvf = this->Y();
    const label nSpecie = Yvf.size();

    forAll(Kp, facei)
    {
        scalar kvMix = 0;
        for (label speciei = 0; speciei < nSpecie; ++speciei)
        {
            kvMix +=
                max(Yvf[speciei].boundaryField()[patchi][facei], scalar(0))
              * this->composition().kappaVib
                (
                    speciei, pp[facei], Tp[facei], TVp[facei],
                    this->composition().thetaVib(speciei)
                );
        }
        Kp[facei] = kvMix;
    }

    return tKp;
}


template<class Basic2TThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::he2TThermo<Basic2TThermo, MixtureType>::alphaheTR() const
{
    auto tAlpha = Foam::tmp<Foam::volScalarField>
    (
        new Foam::volScalarField
        (
            Foam::IOobject
            (
                "alphaheTR",
                this->p_.time().timeName(),
                this->p_.mesh(),
                Foam::IOobject::NO_READ,
                Foam::IOobject::NO_WRITE
            ),
            this->p_.mesh(),
            Foam::dimensionedScalar
            (
                "zero",
                this->alpha_.dimensions(),
                0.0
            )
        )
    );

    auto& A = tAlpha.ref();

    // internal field
    {
        const Foam::scalarField& pi = this->p_.primitiveField();
        const Foam::scalarField& Ti = this->TTR_.primitiveField();
        Foam::scalarField& Ai = A.primitiveFieldRef();

        forAll(Ai, i)
        {
            const typename MixtureType::thermoType& mixture_ =
                this->cellMixture(i);

            Ai[i] =
                mixture_.kappaTR(pi[i], Ti[i])
               /mixture_.CpvTR(pi[i], Ti[i]);
        }
    }

    // patch fields
    forAll(A.boundaryField(), patchi)
    {
        A.boundaryFieldRef()[patchi] = this->alphaheTR(patchi);
    }

    A.rename("alphaheTR");
    return tAlpha;
}


template<class Basic2TThermo, class MixtureType>
Foam::tmp<Foam::scalarField>
Foam::he2TThermo<Basic2TThermo, MixtureType>::alphaheTR
(
    const label patchi
) const
{
    const Foam::scalarField& pp = this->p_.boundaryField()[patchi];
    const Foam::scalarField& Tp = this->TTR_.boundaryField()[patchi];

    Foam::tmp<Foam::scalarField> tAp(new Foam::scalarField(Tp.size()));
    Foam::scalarField& Ap = tAp.ref();

    forAll(Ap, facei)
    {
        const typename MixtureType::thermoType& mixture_ =
            this->patchFaceMixture(patchi, facei);

        Ap[facei] =
            mixture_.kappaTR(pp[facei], Tp[facei])
           /mixture_.CpvTR(pp[facei], Tp[facei]);
    }

    return tAp;
}

template<class Basic2TThermo, class MixtureType>
Foam::scalar Foam::he2TThermo<Basic2TThermo, MixtureType>::thetaVib() const
{
    return this->cellMixture(0).ThetaVib();
}


template<class Basic2TThermo, class MixtureType>
bool Foam::he2TThermo<Basic2TThermo, MixtureType>::read()
{
    if (Basic2TThermo::read())
    {
        MixtureType::read(*this);
        return true;
    }

    return false;
}


// ************************************************************************* //
