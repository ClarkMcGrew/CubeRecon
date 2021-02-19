#ifndef CubeRecurseGeometry_hxx_seen
#define CubeRecurseGeometry_hxx_seen
#include <iostream>

#include <TGeoManager.h>
#include <TGeoNode.h>

namespace Cube {
    class RecurseGeometry;
}

/// RecurseGeometry is a base class that provides tools to recurse through a
/// ROOT Geometry and accumulate information.  The daughter class is expected
/// to provide an "bool Action(int depth)" method.  The argument is the
/// current depth in the geometry (relative to the starting node).  This
/// method will be called for each node during a depth-first traversal.  It
/// can prevent going deeper in the tree by returning false.
///
/// Trivial Example:
///
/// \code
///    class ExampleRecurse: public Cube::RecurseGeometry {
///    public:
///        int Action(const TGeoNode* node, int depth) {
///            std::cout << "At node " << node->GetName() << std::endl;
///            return true;
///        }
///    }
///
///    gGeoManager->CdTop();
///    ExampleRecurse recurseThis;
///    recurseThis.Apply();
/// \endcode
class Cube::RecurseGeometry {
public:
    RecurseGeometry() {}
    virtual ~RecurseGeometry() {}

    /// An example method to override!  This should return false if the
    /// recursion should be stopped.  The action should act on input node, and
    /// NOT change the geometry.  The input arguments are the current node,
    /// and the current depth in the geometry (zero is the node where the
    /// recursion started).
    virtual bool Action(const TGeoNode* node, int depth) {
        // Example of how to get the user extenstion object out of the node.
        TGeoExtension* extension;
        extension = node->GetUserExtension();
        if (extension) {
            std::cout << "A user extension for node " << node->GetName()
                      << " exists at depth " << depth
                      << std::endl;
        }
        // Example of how to get a framework extension.  These are part of
        // TGeoManager, and you probably don't want to use them!.
        extension = node->GetFWExtension();
        if (extension) {
            std::cout << "A framework extension for node " << node->GetName()
                      << " exists at depth " << depth
                      << std::endl;
        }
        // An example of how to get the volume.
        TGeoVolume* volume = node->GetVolume();
        // Example of how to get the user extenstion object out of the volume.
        extension = volume->GetUserExtension();
        if (extension) {
            std::cout << "A user extension for volume " << volume->GetName()
                      << " exists at depth " << depth
                      << std::endl;
        }
        // Example of how to get the framework extenstion out of the volume.
        extension = volume->GetFWExtension();
        if (extension) {
            std::cout << "A framework extension for volume "<<volume->GetName()
                      << " exists at depth " << depth
                      << std::endl;
        }
        return true;
    }

    /// Apply the action recursively to every node starting with the current
    /// node (depth first).
    bool Apply(int depth=0) {
        if (!gGeoManager) return false; // Should NEVER happen!
        TGeoNode * node = gGeoManager->GetCurrentNode();
        if (!node) return false;        // Should NEVER happen!

        if (Action(node,depth)) {
            for (int i=0; i< node->GetNdaughters(); ++i) {
                gGeoManager->CdDown(i);
                Apply(depth+1);
                gGeoManager->CdUp();
            }
        }

        return true;
    }
};
#endif

// Local Variables:
// mode:c++
// c-basic-offset:4
// compile-command:"$(git rev-parse --show-toplevel)/build/cube-build.sh force"
// End:
