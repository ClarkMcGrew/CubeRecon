#include "ToolRecon.hxx"

#include <CubeHandle.hxx>
#include <CubeHitSelection.hxx>
#include <CubeHit.hxx>

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

// Local Variables:
// mode:c++
// c-basic-offset:4
// compile-command:"$(git rev-parse --show-toplevel)/build/cube-build.sh force"
// End:
