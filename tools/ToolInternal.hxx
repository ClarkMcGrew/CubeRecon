#ifndef ToolInternal_hxx_Seen
#define ToolInternal_hxx_Seen
// Obligatory Wizard of Oz Reference: Ignore the code behind the curtain

// This is a quick and dirty solution for saving information needed by all of
// the tools.  This should NOT be accessed directly outside of the tools, and
// you should pretend it doesn't exist.  I'm going to pretend I didn't
// actually do this.

// This defines a bunch of global variables that are filled in
// Cube::Tool::Intiailize, and used in the other tool routines.  Don't do this
// at home, and I shouldn't do it here.  They are defined in
// ToolInitialize.cxx
namespace Cube {
    namespace Tool {
        namespace Internal {
            extern int g3DSTPlanes;
            extern double g3DSTPlaneMin;
            extern double g3DSTPlaneMax;
            extern int g3DSTBars;
            extern double g3DSTBarMin;
            extern double g3DSTBarMax;
            extern int g3DSTCubes;
            extern double g3DSTCubeMin;
            extern double g3DSTCubeMax;
        }
    }
}
#endif

// Local Variables:
// mode:c++
// c-basic-offset:4
// compile-command:"$(git rev-parse --show-toplevel)/build/cube-build.sh force"
// End:
