#include "ToolInitialize.hxx"
#include "ToolInternal.hxx"

#include <CubeRecurseGeometry.hxx>

#include <TGeoManager.h>
#include <TGeoNode.h>
#include <TGeoVolume.h>

int Cube::Tool::Internal::g3DSTPlanes = 0;
double Cube::Tool::Internal::g3DSTPlaneMin = 1E+20;
double Cube::Tool::Internal::g3DSTPlaneMax = -1E+20;
int Cube::Tool::Internal::g3DSTBars = 0;
double Cube::Tool::Internal::g3DSTBarMin = 1E+20;
double Cube::Tool::Internal::g3DSTBarMax = -1E+20;
int Cube::Tool::Internal::g3DSTCubes = 0;
double Cube::Tool::Internal::g3DSTCubeMin = 1E+20;
double Cube::Tool::Internal::g3DSTCubeMax = -1E+20;

namespace {
    bool gInitialized = false;
    class ToolRecursion: public Cube::RecurseGeometry {
        bool Action(const TGeoNode* node, int depth) {
            std::string path(gGeoManager->GetPath());
            std::string name(node->GetName());

            if (path.find("vol3DST_PV") != std::string::npos) {
                return Handle3DST(node,name);
            }

            return true;
        }

        // Get information about the 3DST
        bool Handle3DST(const TGeoNode* node, const std::string& name)  {
            double local[3] = {0,0,0};
            double master[3];
            gGeoManager->LocalToMaster(local,master);

            if (name.find("3DSTPlane_PV") != std::string::npos) {
                ++Cube::Tool::Internal::g3DSTPlanes;
                Cube::Tool::Internal::g3DSTPlaneMin
                    = std::min(Cube::Tool::Internal::g3DSTPlaneMin,master[2]);
                Cube::Tool::Internal::g3DSTPlaneMax
                    = std::max(Cube::Tool::Internal::g3DSTPlaneMax,master[2]);
            }
            else if (name.find("3DSTBar_PV") != std::string::npos) {
                ++Cube::Tool::Internal::g3DSTBars;
                Cube::Tool::Internal::g3DSTBarMin
                    = std::min(Cube::Tool::Internal::g3DSTBarMin,master[1]);
                Cube::Tool::Internal::g3DSTBarMax
                    = std::max(Cube::Tool::Internal::g3DSTBarMax,master[1]);
            }
            else if (name.find("volcube_PV") != std::string::npos) {
                ++Cube::Tool::Internal::g3DSTCubes;
                Cube::Tool::Internal::g3DSTCubeMin
                    = std::min(Cube::Tool::Internal::g3DSTCubeMin,master[0]);
                Cube::Tool::Internal::g3DSTCubeMax
                    = std::max(Cube::Tool::Internal::g3DSTCubeMax,master[0]);
            }

            // Short circuit.  Only look at the daughters for the "_0"
            // plane and bar.
            if (name.find("3DSTPlane_PV") != std::string::npos
                && name.find("_0") == std::string::npos) {
                return false;
            }

            if (name.find("3DSTBar_PV") != std::string::npos
                && name.find("_0") == std::string::npos) {
                return false;
            }

            return true;
        }

    };
}

void Cube::Tool::Initialize() {
    if (gInitialized) return;
    gGeoManager->CdTop();
    ToolRecursion recurse;
    recurse.Apply();
    gInitialized = true;
    std::cout << "Tools initialized" << std::endl;
}

// Local Variables:
// mode:c++
// c-basic-offset:4
// compile-command:"$(git rev-parse --show-toplevel)/build/cube-build.sh force"
// End:
