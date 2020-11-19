#ifndef ToolG4Hits_hxx_seen
#define ToolG4Hits_hxx_seen

#include <CubeHandle.hxx>
#include <CubeG4Hit.hxx>
#include <CubeHitSelection.hxx>
#include <CubeEvent.hxx>
#include <CubeReconObject.hxx>

#include <vector>

namespace Cube {
    namespace Tool {
        /// Get all the hit segments associated with this hit.  The resulting
        /// vector is not sorted, and may contain duplicates.
        std::vector<Cube::Handle<Cube::G4Hit>>
        HitG4Hits(Cube::Event& event, Cube::Handle<Cube::Hit> hit);

        /// Get all the hit segments associated with this hit selection.  The
        /// resulting vector is sorted by the hit segment time (earliest to
        /// latest).
        std::vector<Cube::Handle<Cube::G4Hit>>
        SelectionG4Hits(Cube::Event& event, Cube::HitSelection& hits);

        /// Get all the hit segments associated with this reconstruction
        /// object.  The resulting vector is sorted by the hit segment time
        /// (earliest to latest).
        std::vector<Cube::Handle<Cube::G4Hit>>
        ObjectG4Hits(Cube::Event& event, Cube::ReconObject& object);

        /// Get all the hit segments associated with this track id.  The
        /// resulting vector is sorted by the hit segment time (earliest to
        /// latest).
        std::vector<Cube::Handle<Cube::G4Hit>>
        TrajectoryG4Hits(Cube::Event& event, int trackId);
    }
}
#endif

// Local Variables:
// mode:c++
// c-basic-offset:4
// compile-command:"$(git rev-parse --show-toplevel)/build/cube-build.sh force"
// End:
