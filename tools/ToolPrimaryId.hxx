#ifndef ToolPrimaryId_hxx_seen
#define ToolPrimaryId_hxx_seen

#include <CubeEvent.hxx>
namespace Cube {
    namespace Tool {
        /// Find the primary trajectory id (this is the trajectory that was
        /// started by a GEANT4 primary particle).  The trajectories are a
        /// tree of parents with their children.  This looks at the ancestors
        /// to find the originating parent trajectory.
        int PrimaryId(Cube::Event& event, int trajId);

        /// Return a vector of the trajectory id values for all of the primary
        /// particles in the event.  This is trying to find the "useful"
        /// primary trajectories, which means it should include the michel
        /// electrons, and the photons from pizero decay.
        std::vector<int> AllPrimaries(Cube::Event& event);
    }
}
#endif

// Local Variables:
// mode:c++
// c-basic-offset:4
// compile-command:"$(git rev-parse --show-toplevel)/build/cube-build.sh force"
// End:
