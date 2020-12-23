#ifndef ToolRecon_hxx_seen
#define ToolRecon_hxx_seen

#include <CubeUnits.hxx>
#include <CubeEvent.hxx>
#include <CubeReconObject.hxx>

namespace Cube {
    namespace Tool {
        // Return a positive integer if the hit is contained (i.e. not on the
        // edge), otherwise zero.  This is the number of cubes the current hit
        // is from the edge of the detector.
        double DistanceBetweenObjects(Cube::ReconObject& obj1,
                                      Cube::ReconObject& obj2);

        // Return a positive integer if the object is contained, otherwise,
        // return zero.
        bool AreNeighboringObjects(Cube::ReconObject& obj1,
                                   Cube::ReconObject& obj2,
                                   double dist = 4.0*unit::cm);
    }
}
#endif

// Local Variables:
// mode:c++
// c-basic-offset:4
// compile-command:"$(git rev-parse --show-toplevel)/build/cube-build.sh force"
// End:
