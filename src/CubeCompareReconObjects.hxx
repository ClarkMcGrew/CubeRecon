#ifndef CubeCompareReconObjects_hxx_seen
#define CubeCompareReconObjects_hxx_seen

#include <CubeReconNode.hxx>
#include <CubeReconCluster.hxx>
#include <CubeReconTrack.hxx>
#include <CubeReconVertex.hxx>
#include <CubeHandle.hxx>

namespace Cube {

    /// This is a predicate to order recon objects.  This puts the vertices,
    /// followed by tracks and clusters.  Within any type the objects are
    /// ordered by decreasing size (i.e. largest first).  If everything about
    /// two objects is equal, they are ordered by the value of the pointers.
    struct CompareReconObjects {
        bool operator () (Cube::Handle<Cube::ReconObject> lhs,
                          Cube::Handle<Cube::ReconObject> rhs) {
            Cube::Handle<Cube::ReconVertex> lv = lhs;
            Cube::Handle<Cube::ReconVertex> rv = rhs;
            Cube::Handle<Cube::ReconTrack> lt = lhs;
            Cube::Handle<Cube::ReconTrack> rt = rhs;
            Cube::Handle<Cube::ReconCluster> lc = lhs;
            Cube::Handle<Cube::ReconCluster> rc = rhs;

            if (lv && rv) {
                if (lv->GetConstituents()&&!rv->GetConstituents()) return true;
                if (!lv->GetConstituents()&&rv->GetConstituents()) return false;
                if (lv->GetConstituents() && rv->GetConstituents()) {
                    if (lv->GetConstituents()->size()
                        != rv->GetConstituents()->size()) {
                        return (lv->GetConstituents()->size()
                                > rv->GetConstituents()->size());
                    }
                    int s1 = -1;
                    if (lv->GetHitSelection()) s1=lv->GetHitSelection()->size();
                    int s2 = -1;
                    if (rv->GetHitSelection()) s2=rv->GetHitSelection()->size();
                    return (s1 > s2);
                }
            }
            if (lv && !rv) return true;
            if (!lv && rv) return false;

            if (lt && rt) return lt->GetNodes().size() > rt->GetNodes().size();
            if (lt && !rt) return true;
            if (!lt && rt) return false;

            if (lc && rc) {
                return (lc->GetHitSelection()->size()
                        > rc->GetHitSelection()->size());
            }
            if (lc && !rc) return true;
            if (!lc && rc) return false;

            return Cube::GetPointer(lhs) < Cube::GetPointer(rhs);
        }
    };
};
#endif

// Local Variables:
// mode:c++
// c-basic-offset:4
// compile-command:"$(git rev-parse --show-toplevel)/build/cube-build.sh force"
// End:
