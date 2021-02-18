#include "ToolRecon.hxx"

#include <CubeHandle.hxx>
#include <CubeHitSelection.hxx>
#include <CubeHit.hxx>
#include <CubeInfo.hxx>

double Cube::Tool::DistanceBetweenObjects(Cube::ReconObject& obj1,
                                          Cube::ReconObject& obj2) {
    Cube::Handle<Cube::HitSelection> hit1 = obj1.GetHitSelection();
    Cube::Handle<Cube::HitSelection> hit2 = obj2.GetHitSelection();
    if (!hit1 || !hit2) return 1E+20;

    double minDist = 1E+20;
    for (Cube::HitSelection::iterator h1 = hit1->begin();
         h1 != hit1->end(); ++h1) {
        for (Cube::HitSelection::iterator h2 = hit2->begin();
             h2 != hit2->end(); ++h2) {
            double d = ((*h1)->GetPosition()-(*h2)->GetPosition()).Mag();
            if (d < minDist) minDist = d;
        }
    }

    return minDist;
}

bool Cube::Tool::AreNeighboringObjects(Cube::ReconObject& obj1,
                                       Cube::ReconObject& obj2,
                                       double dist) {
    double d = DistanceBetweenObjects(obj1,obj2);
    return (d < dist);
}

std::string Cube::Tool::ObjectDetectors(Cube::ReconObject& obj) {
    std::string result;
    Cube::Handle<Cube::HitSelection> hits = obj.GetHitSelection();
    bool has3DST = false;
    bool hasDSTPC = false;
    bool hasTopTPC = false;
    bool hasBotTPC = false;
    bool hasECal = false;
    for (Cube::HitSelection::iterator h = hits->begin();
         h != hits->end(); ++h) {
        if (Cube::Info::Is3DST((*h)->GetIdentifier())) {
            has3DST = true;
        }
        else if (Cube::Info::IsDownstreamTPC((*h)->GetIdentifier())) {
            hasDSTPC = true;
        }
        else if (Cube::Info::IsTopTPC((*h)->GetIdentifier())) {
            hasTopTPC = true;
        }
        else if (Cube::Info::IsBottomTPC((*h)->GetIdentifier())) {
            hasBotTPC = true;
        }
        else if (Cube::Info::IsECal((*h)->GetIdentifier())) {
            hasECal = true;
        }
    }
    result = ":";
    if (has3DST) {
        result += "3DST:";
    }
    if (hasDSTPC) {
        result += "TPC-Downstream:";
    }
    if (hasTopTPC) {
        result += "TPC-Top:";
    }
    if (hasBotTPC) {
        result += "TPC-Bottom:";
    }
    if (hasECal) {
        result += "ECal:";
    }

    return result;
}

// Local Variables:
// mode:c++
// c-basic-offset:4
// compile-command:"$(git rev-parse --show-toplevel)/build/cube-build.sh force"
// End:
