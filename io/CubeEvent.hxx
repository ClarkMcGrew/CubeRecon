#ifndef CubeEvent_hxx_seen
#define CubeEvent_hxx_seen
#include "CubeAlgorithmResult.hxx"
#include "CubeG4Hit.hxx"

#include <TRef.h>
#include <TROOT.h>

namespace Cube {
    class Event;
}

/// An event processed through the Cube reconstruction.
class Cube::Event : public Cube::AlgorithmResult {
public:
    Event(void);
    virtual ~Event();

    static Cube::Event* CurrentEvent();
    void MakeCurrentEvent() const;

    virtual void Initialize(int run = -1, int event = -1,
                            TObject* evt = NULL) {
        SetName("CubeEvent");
        SetTitle("Hits generated by ERepSim");
        fRunId = run;
        fEventId = event;
        fEDepSimEvent = TRef(evt);
        G4Hits.clear();
        Cube::AlgorithmResult::Initialize();
    }

    /// The run number for this event.
    int GetRunId() const {return fRunId;}

    /// The event number.
    int GetEventId() const {return fEventId;}

    void ls(Option_t *opt = "") const {
        TROOT::IndentLevel();
        std::cout << "%%% Event : " << fRunId << " / " << fEventId << std::endl;
        TROOT::IncreaseDirLevel();
        fEDepSimEvent.ls();
        Cube::AlgorithmResult::ls(opt);
        TROOT::DecreaseDirLevel();
    }

private:

    /// The run number
    int fRunId;

    /// The event number.
    int fEventId;

    /// The parent EDepSim event (if available).
    TRef fEDepSimEvent;

public:
    /////////////////////////////////////////////////////////////////////////
    /// The stuff below this comment "doesn't exist".  This is a cheaters way
    /// to get some of the important truth information needed to debug the
    /// reconstruction.
    ////////////////////////////////////////////////////////////////////////
    typedef std::vector<Cube::Handle<Cube::G4Hit>> G4HitVector;
    G4HitVector G4Hits;

    ClassDef(Event,1);
};
#endif

// Local Variables:
// mode:c++
// c-basic-offset:4
// compile-command:"$(git rev-parse --show-toplevel)/build/cube-build.sh force"
// End:
