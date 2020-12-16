#ifndef ToolMainTrajectory_hxx_seen
#define ToolMainTrajectory_hxx_seen

#include <CubeEvent.hxx>
#include <CubeReconObject.hxx>

namespace Cube {
    namespace Tool {
        /// Find the trajectory that contributed most to the track.  The
        /// longest trajectory wins.
        int MainTrajectory(Cube::Event& event, Cube::ReconObject& object);

        /// Find the fraction of the track length that comes from the longest
        /// trajectory.
        double MainPurity(Cube::Event& event, Cube::ReconObject& object);

        /// Find how complete the object is.  This is the fraction of the main
        /// trajectory track length that is contained in the reconstruction
        /// objet.
        double MainCompleteness(Cube::Event& event, Cube::ReconObject& object);
    }
}
#endif

// Local Variables:
// mode:c++
// c-basic-offset:4
// compile-command:"$(git rev-parse --show-toplevel)/build/cube-build.sh force"
// End:
