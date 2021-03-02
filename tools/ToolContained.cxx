#include "ToolContained.hxx"
#include "ToolInitialize.hxx"
#include "ToolInternal.hxx"

#include <CubeEvent.hxx>
#include <CubeReconObject.hxx>
#include <CubeHit.hxx>
#include <CubeHandle.hxx>
#include <CubeInfo.hxx>
#include <CubeUnits.hxx>

#include <TVector3.h>

#include <vector>

int Cube::Tool::Contained3DSTHit(const Cube::Hit& hit) {
    Cube::Tool::Initialize();
    int id = hit.GetIdentifier();
    int num = Cube::Info::CubeNumber(id);
    int bar = Cube::Info::CubeBar(id);
    int pln = Cube::Info::CubePlane(id);
    int low = std::min(std::min(num,bar),pln);
    // This should be carried in the data, but for now will need to be hard
    // coded.
    int hiNum = Cube::Tool::Internal::g3DSTCubes - num - 1;
    int hiBar = Cube::Tool::Internal::g3DSTBars - bar - 1;
    int hiPln = Cube::Tool::Internal::g3DSTPlanes - pln - 1;
    int hi = std::min(std::min(hiNum,hiBar),hiPln);
    if (low < 0) return 0;
    if (hi < 0) return 0;
    return std::min(low,hi);
}

int Cube::Tool::Contained3DSTObject(const Cube::ReconObject& object) {
    Cube::Handle<Cube::HitSelection> hits = object.GetHitSelection();
    if (!hits) return 0;
    if (hits->empty()) return 0;
    int minDist = 1E+6;
    for (Cube::HitSelection::iterator h = hits->begin();
         h != hits->end(); ++h) {
        int dist = Cube::Tool::Contained3DSTHit(*(*h));
        minDist = std::min(dist,minDist);
    }
    return minDist;
}

// Return the distance inside of the 3DST.
double Cube::Tool::Contained3DSTPoint(const TVector3& point) {
    double minDist = 1E+20;
    minDist = std::min(Cube::Tool::Contained3DSTRight(point),minDist);
    minDist = std::min(Cube::Tool::Contained3DSTLeft(point),minDist);
    minDist = std::min(Cube::Tool::Contained3DSTBottom(point),minDist);
    minDist = std::min(Cube::Tool::Contained3DSTTop(point),minDist);
    minDist = std::min(Cube::Tool::Contained3DSTUpstream(point),minDist);
    minDist = std::min(Cube::Tool::Contained3DSTDownstream(point),minDist);
    return minDist;
}

double Cube::Tool::Contained3DSTUpstream(const TVector3& point) {
    Cube::Tool::Initialize();
    return point.Z()-Cube::Tool::Internal::g3DSTPlaneMin;
}

double Cube::Tool::Contained3DSTDownstream(const TVector3& point) {
    Cube::Tool::Initialize();
    return Cube::Tool::Internal::g3DSTPlaneMax-point.Z();
}

double Cube::Tool::Contained3DSTBottom(const TVector3& point) {
    Cube::Tool::Initialize();
    return point.Y()-Cube::Tool::Internal::g3DSTBarMin;
}

double Cube::Tool::Contained3DSTTop(const TVector3& point) {
    Cube::Tool::Initialize();
    return Cube::Tool::Internal::g3DSTBarMax-point.Y();
}

double Cube::Tool::Contained3DSTRight(const TVector3& point) {
    Cube::Tool::Initialize();
    return point.X()-Cube::Tool::Internal::g3DSTCubeMin;
}

double Cube::Tool::Contained3DSTLeft(const TVector3& point) {
    Cube::Tool::Initialize();
    return Cube::Tool::Internal::g3DSTCubeMax-point.X();
}

// Local Variables:
// mode:c++
// c-basic-offset:4
// compile-command:"$(git rev-parse --show-toplevel)/build/cube-build.sh force"
// End:
