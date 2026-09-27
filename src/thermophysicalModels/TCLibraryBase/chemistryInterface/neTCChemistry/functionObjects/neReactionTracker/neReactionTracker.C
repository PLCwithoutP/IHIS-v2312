#include "neReactionTracker.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "BasicNeChemistryModel.H"
#include "ne2TReactionThermo.H"

template<class ChemistryModelType>
void Foam::functionObjects::neReactionTracker<ChemistryModelType>::writeFileHeader
(
    Ostream& os
) const
{
    writeHeader(os, "Reaction tracker");
    volRegion::writeFileHeader(*this, os);
    writeHeaderValue(os, "activityTolerance", activityTolerance_);

    writeCommented(os, "Time");
    writeTabbed(os, "ReactionName");
    writeTabbed(os, "ReactionEquation");
    writeTabbed(os, "ActiveCellCount");
    os << endl;
}


template<class ChemistryModelType>
Foam::functionObjects::neReactionTracker<ChemistryModelType>::neReactionTracker
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    volRegion(fvMeshFunctionObject::mesh_, dict),
    writeFile(obr_, name, typeName, dict),
    chemistryModel_
    (
        fvMeshFunctionObject::mesh_.lookupObject<ChemistryModelType>
        (
            "chemistryProperties"
        )
    ),
    activityTolerance_(dict.getOrDefault<scalar>("activityTolerance", 1e-20)),
    writeToLog_(dict.getOrDefault<bool>("writeToLog", true))
{
    writeFileHeader(file());
}


template<class ChemistryModelType>
bool Foam::functionObjects::neReactionTracker<ChemistryModelType>::read
(
    const dictionary& dict
)
{
    regionFunctionObject::read(dict);

    activityTolerance_ = dict.getOrDefault<scalar>("activityTolerance", 1e-20);
    writeToLog_ = dict.getOrDefault<bool>("writeToLog", true);

    return true;
}


template<class ChemistryModelType>
bool Foam::functionObjects::neReactionTracker<ChemistryModelType>::execute()
{
    return true;
}


template<class ChemistryModelType>
bool Foam::functionObjects::neReactionTracker<ChemistryModelType>::write()
{
    const label nReaction = chemistryModel_.nReaction();

    volRegion::update();

    const bool useAll = this->volRegion::useAllCells();

    for (label ri = 0; ri < nReaction; ++ri)
    {
        const word rxnName = chemistryModel_.reactionName(ri);
        const string rxnEq = chemistryModel_.reactionEquation(ri);

        tmp<volScalarField::Internal> tOmegaR
        (
            chemistryModel_.calculateReactionOmega(ri)
        );
        const volScalarField::Internal& omegaR = tOmegaR();

        label activeCount = 0;

        if (useAll)
        {
            forAll(omegaR, celli)
            {
                if (mag(omegaR[celli]) > activityTolerance_)
                {
                    ++activeCount;
                }
            }
        }
        else
        {
            const labelList& ids = cellIDs();

            forAll(ids, i)
            {
                const label celli = ids[i];

                if (mag(omegaR[celli]) > activityTolerance_)
                {
                    ++activeCount;
                }
            }
        }

        writeCurrentTime(file());
        file() << token::TAB << rxnName
               << token::TAB << rxnEq
               << token::TAB << activeCount
               << nl;

        if (writeToLog_)
        {
            Info<< rxnName << ": "
                << rxnEq
                << " ----> "
                << activeCount
                << " active cells"
                << nl;
        }
    }

    file() << nl << endl;

    return true;
}


// * * * * * * * * * * * * * * Runtime registration * * * * * * * * * * * * //

namespace Foam
{
    typedef
        functionObjects::neReactionTracker
        <
            BasicNeChemistryModel<ne2TReactionThermo>
        >
        ne2TReactionTracker;

    defineTemplateTypeNameAndDebugWithName
    (
        ne2TReactionTracker,
        "ne2TReactionTracker",
        0
    );

    addToRunTimeSelectionTable
    (
        functionObject,
        ne2TReactionTracker,
        dictionary
    );
}
