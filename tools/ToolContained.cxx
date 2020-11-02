#include "ToolContained.hxx"

#include <CubeEvent.hxx>
#include <CubeReconObject.hxx>
#include <CubeHit.hxx>
#include <CubeHandle.hxx>
#include <CubeInfo.hxx>

#include <TVector3.h>

#include <vector>

int Cube::Tool::ContainedHit(Cube::Hit& hit) {
    int id = hit.GetIdentifier();
    int num = Cube::Info::CubeNumber(id);
    int bar = Cube::Info::CubeBar(id);
    int pln = Cube::Info::CubePlane(id);
    int low = std::min(std::min(num,bar),pln);
    // This should be carried in the data, but for now will need to be hard
    // coded.
    int maxNum = 252;
    int maxBar = 236;
    int maxPln = 200;
    int hiNum = maxNum - num - 1;
    int hiBar = maxBar - bar - 1;
    int hiPln = maxPln - pln - 1;
    int hi = std::min(std::min(hiNum,hiBar),hiPln);
    if (low < 0) return 0;
    if (hi < 0) return 0;
    return std::min(low,hi);
}

int Cube::Tool::ContainedObject(Cube::ReconObject& object) {
    Cube::Handle<Cube::HitSelection> hits = object.GetHitSelection();
    if (!hits) return 0;
    if (hits->empty()) return 0;
    int minDist = 1E+6;
    for (Cube::HitSelection::iterator h = hits->begin();
         h != hits->end(); ++h) {
        int dist = Cube::Tool::ContainedHit(*(*h));
        minDist = std::min(dist,minDist);
    }
    return minDist;
}

double Cube::Tool::ContainedPoint(TVector3 point) {
    double minDist = 1E+10;
    double xBoundary = 1260.0;
    minDist = std::min(xBoundary-point.X(),minDist);
    minDist = std::min(point.X()+xBoundary,minDist);
    double yBoundary = 1180;
    minDist = std::min(yBoundary-point.Y(),minDist);
    minDist = std::min(point.Y()+yBoundary,minDist);
    double zMin = 3510.0;
    double zMax = 5510.0;
    minDist = std::min(zMax-point.Z(),minDist);
    minDist = std::min(point.Z()-zMin,minDist);
    return minDist;
}

// Local Variables:
// mode:c++
// c-basic-offset:4
// compile-command:"$(git rev-parse --show-toplevel)/build/cube-build.sh force"
// End:
