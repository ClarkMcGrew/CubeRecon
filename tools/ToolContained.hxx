#ifndef ToolContainedObject_hxx_seen
#define ToolContainedObject_hxx_seen

#include <CubeEvent.hxx>
#include <CubeReconObject.hxx>

#include <TVector3.h>

namespace Cube {
    namespace Tool {
        // Return a positive integer if the hit is contained (i.e. not on the
        // edge), otherwise zero.  This is the number of cubes the current hit
        // is from the edge of the detector.
        int ContainedHit(Cube::Hit& hit);

        // Return a positive integer if the object is contained, otherwise,
        // return zero.
        int ContainedObject(Cube::ReconObject& object);

        // Return a positive value (in mm) if the object is contained,
        // otherwise, return zero or negative.
        double ContainedPoint(TVector3 point);
    }
}
#endif

// Local Variables:
// mode:c++
// c-basic-offset:4
// compile-command:"$(git rev-parse --show-toplevel)/build/cube-build.sh force"
// End:
