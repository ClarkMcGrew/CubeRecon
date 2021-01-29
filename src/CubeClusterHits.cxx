#include "CubeClusterHits.hxx"
#include "CubeMakeUsed.hxx"
#include "CubeClusterManagement.hxx"
#include "CubeInfo.hxx"

#include <CubeLog.hxx>
#include <CubeHandle.hxx>
#include <CubeReconCluster.hxx>
#include <CubeAlgorithmResult.hxx>
#include <TTmplDensityCluster.hxx>

namespace {
    // Determine the proximity of hits.  The proximity is defined as the
    // number of "cube to cube" jumps needed to go between the cubes.
    // Diagonal steps are allowed.
    struct CubeProximity {
        double operator()(Cube::Handle<Cube::Hit> lhs,
                          Cube::Handle<Cube::Hit> rhs) {
            if (Cube::Info::Is3DST(lhs->GetIdentifier())
                && Cube::Info::Is3DST(rhs->GetIdentifier())) {
                double dx = Cube::Info::CubeNumber(lhs->GetIdentifier())
                    - Cube::Info::CubeNumber(rhs->GetIdentifier());

                double dy = Cube::Info::CubeBar(lhs->GetIdentifier())
                    - Cube::Info::CubeBar(rhs->GetIdentifier());

                double dz = Cube::Info::CubePlane(lhs->GetIdentifier())
                    - Cube::Info::CubePlane(rhs->GetIdentifier());

                double d = std::max(std::abs(dx), std::abs(dy));
                return std::max(d,std::abs(dz));
            }
            if (Cube::Info::IsTPC(lhs->GetIdentifier())
                && Cube::Info::IsTPC(rhs->GetIdentifier())) {
                TVector3 dist
                    = lhs->GetPosition() - rhs->GetPosition();
                for (int i = 0; i < 3; ++i) {
                    dist[i] = std::abs(dist[i])
                        - lhs->GetSize()[i]
                        - rhs->GetSize()[i];
                    if (dist[i] < 0.0) dist[i] = 0.0;
                }
                return dist.Mag()/10.0;
            }
            if (Cube::Info::IsECal(lhs->GetIdentifier())
                && Cube::Info::IsECal(rhs->GetIdentifier())) {
                TVector3 dist
                    = lhs->GetPosition() - rhs->GetPosition();
                for (int i = 0; i < 3; ++i) {
                    dist[i] = std::abs(dist[i])
                        - lhs->GetSize()[i]
                        - rhs->GetSize()[i];
                    if (dist[i] < 0.0) dist[i] = 0.0;
                }
                return dist.Mag()/40.0;
            }
            return 1E+20;
        }
    };
}

Cube::ClusterHits::ClusterHits()
    : Cube::Algorithm("ClusterHits") {
    // We need simply connected hits, so this isn't actually tunable, but can
    // be overridden with the setters
    fCubeNeighborhood = 1;
    fMinimumPoints = 1;
}

Cube::Handle<Cube::AlgorithmResult>
Cube::ClusterHits::Process(const Cube::AlgorithmResult& input,
                           const Cube::AlgorithmResult&,
                           const Cube::AlgorithmResult&) {
    Cube::Handle<Cube::HitSelection> inputHits = input.GetHitSelection();
    CUBE_LOG(2) << "Process Cube::ClusterHits" << std::endl;

    // Create the result for this algorithm.
    Cube::Handle<Cube::AlgorithmResult> result = CreateResult();

    // Check that we have hits!
    if (!inputHits) {
        return result;
    }
    if (inputHits->empty()) {
        return result;
    }

    // Create the container for the final objects.
    Cube::Handle<Cube::ReconObjectContainer>
        finalObjects(new Cube::ReconObjectContainer("final"));
    result->AddObjectContainer(finalObjects);

    // Apply the DBScan algorithm.
    // Work around template parsing bug in some GCC versions...
    typedef Cube::Handle<Cube::Hit> Arg;

    // Make a typedef for the ClusterAlgorithm.
    typedef TTmplDensityCluster<Arg, CubeProximity> ClusterCubes;

    // Make sure that we don't ever run into roundoff issues.
    const double maxDist = fCubeNeighborhood + 0.5;
    ClusterCubes clusterCubes(fMinimumPoints, maxDist);

    clusterCubes.Cluster(inputHits->begin(), inputHits->end());

    // Transfer all of the clusters to the final objects.
    for (int i=0; i<clusterCubes.GetClusterCount(); ++i) {
        const ClusterCubes::Points& points = clusterCubes.GetCluster(i);
        Cube::Handle<Cube::ReconCluster> cluster
            = Cube::CreateCluster("clusterHits",points.begin(),points.end());
        finalObjects->push_back(cluster);
    }

    // Build the hit selections.
    Cube::MakeUsed makeUsed(*inputHits);
    result = makeUsed(result);

    return result;
}

// Local Variables:
// mode:c++
// c-basic-offset:4
// compile-command:"$(git rev-parse --show-toplevel)/build/cube-build.sh force"
// End:
