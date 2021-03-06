#ifndef CubeTrackState_hxx_seen
#define CubeTrackState_hxx_seen

#include "CubeReconState.hxx"

namespace Cube {
    class TrackState;
}

/// A state holding parameters associated with a TReconTrack, and the
/// intermediate states.
class Cube::TrackState:
    public Cube::ReconState
{
public:
    TrackState();
    TrackState(const TrackState& init);
    virtual ~TrackState();
    virtual TrackState& operator=(const TrackState& rhs);

    ENERGY_DEPOSIT_STATE_DECLARATION;
    POSITION_STATE_DECLARATION;
    DIRECTION_STATE_DECLARATION;
    CURVATURE_STATE_DECLARATION;
    WIDTH_STATE_DECLARATION;

    ENERGY_DEPOSIT_STATE_PRIVATE;
    POSITION_STATE_PRIVATE;
    DIRECTION_STATE_PRIVATE;
    CURVATURE_STATE_PRIVATE;
    WIDTH_STATE_PRIVATE;

    ClassDef(TrackState,1);
};
#endif

// Local Variables:
// mode:c++
// c-basic-offset:4
// compile-command:"$(git rev-parse --show-toplevel)/build/cube-build.sh force"
// End:
