#ifndef ToolCubeTruth_hxx_seen
#define ToolCubeTruth_hxx_seen

#include <CubeEvent.hxx>
#include <CubeHit.hxx>
#include <CubeHandle.hxx>

namespace Cube {
    namespace Tool {
        /// Sum the true deposit for a cube.  This eliminates crosstalk hits.
        double CubeDeposit(Cube::Event& event, Cube::Handle<Cube::Hit> hit);

        /// Sum the true deposit for a cube coming from cross talk.
        double CubeCrossTalk(Cube::Event& event, Cube::Handle<Cube::Hit> hit);

        /// Determine the true time that the cube was hit.  The time is based
        /// on when it crosses a certain energy deposition threshold.
        double CubeTime(Cube::Event& event, Cube::Handle<Cube::Hit> hit,
                        double threshold = 0.1 /*MeV*/);
    }
}
#endif

// Local Variables:
// mode:c++
// c-basic-offset:4
// compile-command:"$(git rev-parse --show-toplevel)/build/cube-build.sh force"
// End:
