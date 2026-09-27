#ifndef functionObjects_neReactionTracker_H
#define functionObjects_neReactionTracker_H

#include "fvMeshFunctionObject.H"
#include "volRegion.H"
#include "writeFile.H"

namespace Foam
{
namespace functionObjects
{

template<class ChemistryModelType>
class neReactionTracker
:
    public fvMeshFunctionObject,
    public volRegion,
    public writeFile
{
    // Private Data

        //- Chemistry model
        const ChemistryModelType& chemistryModel_;

        //- Threshold for considering a reaction active in a cell
        scalar activityTolerance_;

        //- Also print to standard log
        bool writeToLog_;


    // Private Member Functions

        //- Write file header
        virtual void writeFileHeader(Ostream& os) const;

        //- Disallow default bitwise copy construction
        neReactionTracker(const neReactionTracker&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const neReactionTracker&) = delete;


public:

    //- Runtime type information
    TypeName("neReactionTracker");


    // Constructors

        neReactionTracker
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~neReactionTracker() = default;


    // Member Functions

        //- Read dictionary
        virtual bool read(const dictionary& dict);

        //- Execute
        virtual bool execute();

        //- Write
        virtual bool write();
};

} // End namespace functionObjects
} // End namespace Foam

#endif
