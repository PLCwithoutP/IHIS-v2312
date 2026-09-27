/*---------------------------------------------------------------------------*\
    neTCChemistryAdapter.H

    Concrete implementation of chemistryInterface backed by neTCLib chemistry.
    Wraps autoPtr<BasicNeChemistryModel<ne2TReactionThermo>>.

    This header is included ONLY from chemistryInterface.C — it is an
    implementation detail, not part of the public interface.
\*---------------------------------------------------------------------------*/

#ifndef neTCChemistryAdapter_H
#define neTCChemistryAdapter_H

#include "chemistryInterface.H"
#include "BasicNeChemistryModel.H"
#include "ne2TReactionThermo.H"

namespace Foam
{

class neTCChemistryAdapter
:
    public chemistryInterface
{
    autoPtr<BasicNeChemistryModel<ne2TReactionThermo>> model_;

public:

    neTCChemistryAdapter(ne2TReactionThermo& thermo)
    :
        model_(BasicNeChemistryModel<ne2TReactionThermo>::New(thermo))
    {}

    scalar solve(scalar deltaT) override
    {
        return model_->solve(deltaT);
    }

    const volScalarField::Internal& RR(label i) const override
    {
        return model_->RR(i);
    }

    const PtrList<volScalarField::Internal>& RR() const override
    {
        return model_->RR();
    }

    tmp<volScalarField> Qdot() const override
    {
        return model_->Qdot();
    }
};

} // End namespace Foam

#endif
