#ifndef ToolRecon_hxx_seen
#define ToolRecon_hxx_seen

#include <CubeUnits.hxx>
#include <CubeEvent.hxx>
#include <CubeReconObject.hxx>

namespace Cube {
    namespace Tool {
        // Find the distance between two reconstruction objects.  This will be
        // the distance between the two closest hits in the objects.
        double DistanceBetweenObjects(Cube::ReconObject& obj1,
                                      Cube::ReconObject& obj2);

        // Return true if the objects are neighboring.  They are neighboring
        // if the distance between is less than the third parameter (default
        // is 4cm).
        bool AreNeighboringObjects(Cube::ReconObject& obj1,
                                   Cube::ReconObject& obj2,
                                   double dist = 4.0*unit::cm);

        // Return the names of the detectors that contain the object.  This is
        // based on the hits are in the object.
        std::string ObjectDetectors(Cube::ReconObject& obj);

    }
}
#endif

// Local Variables:
// mode:c++
// c-basic-offset:4
// compile-command:"$(git rev-parse --show-toplevel)/build/cube-build.sh force"
// End:
