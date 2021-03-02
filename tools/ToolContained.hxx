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
        int Contained3DSTHit(const Cube::Hit& hit);

        // Return a positive integer if the object is contained, otherwise,
        // return zero.
        int Contained3DSTObject(const Cube::ReconObject& object);

        // Return a positive value (in mm) if the object is contained,
        // otherwise, return zero or negative.
        double Contained3DSTPoint(const TVector3& point);
        double Contained3DSTRight(const TVector3& point);
        double Contained3DSTLeft(const TVector3& point);
        double Contained3DSTBottom(const TVector3& point);
        double Contained3DSTTop(const TVector3& point);
        double Contained3DSTUpstream(const TVector3& point);
        double Contained3DSTDownstream(const TVector3& point);
    }
}
#endif

// Local Variables:
// mode:c++
// c-basic-offset:4
// compile-command:"$(git rev-parse --show-toplevel)/build/cube-build.sh force"
// End:
