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
        std::vector<Cube::Handle<Cube::G4Hit>>
        HitG4Hits(Cube::Event& event, Cube::Handle<Cube::Hit> hit);

        std::vector<Cube::Handle<Cube::G4Hit>>
        SelectionG4Hits(Cube::Event& event, Cube::HitSelection& hits);

        std::vector<Cube::Handle<Cube::G4Hit>>
        ObjectG4Hits(Cube::Event& event, Cube::ReconObject& object);

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
