#include "CubeFindKinks.hxx"
#include "CubeMakeUsed.hxx"
#include "CubeClusterManagement.hxx"

#include <CubeLog.hxx>
#include <CubeHandle.hxx>
#include <CubeReconCluster.hxx>
#include <CubeAlgorithmResult.hxx>
#include <CubeUnits.hxx>

#include <list>

Cube::FindKinks::FindKinks()
    : Cube::Algorithm("FindKinks") {

    fScanLength = 7;

    fMinimumScanLength = 4;

    fKinkThreshold = 1.7*unit::cm;

    fLengthFraction = 0.6;

}

Cube::Handle<Cube::AlgorithmResult>
Cube::FindKinks::Process(const Cube::AlgorithmResult& input,
                         const Cube::AlgorithmResult&,
                         const Cube::AlgorithmResult&) {
    Cube::Handle<Cube::HitSelection> inputHits = input.GetHitSelection();
    CUBE_LOG(2) << "Cube::FindKinks Process" << std::endl;

    // Create the result for this algorithm.
    Cube::Handle<Cube::AlgorithmResult> result = CreateResult();

    // Check that we have hits!
    if (!inputHits) {
        return result;
    }
    if (inputHits->empty()) {
        return result;
    }

    // Sanity check
    if (fScanLength < fMinimumScanLength) {
        fScanLength = fMinimumScanLength;
    }

    // Create the container for the final objects.
    Cube::Handle<Cube::ReconObjectContainer>
        finalObjects(new Cube::ReconObjectContainer("final"));
    result->AddObjectContainer(finalObjects);

    std::list<Cube::Handle<Cube::ReconCluster>> clusterList;

    // Copy the input clusters into a list.
    Cube::Handle<Cube::ReconObjectContainer> inputObjects
        = input.GetObjectContainer("final");
    for (Cube::ReconObjectContainer::iterator o = inputObjects->begin();
         o != inputObjects->end(); ++ o) {
        Cube::Handle<Cube::HitSelection> objectHits = (*o)->GetHitSelection();
        if (!objectHits) continue;

        if (objectHits->empty()) {
            CUBE_ERROR << "empty object" << std::endl;
            continue;
        }

        clusterList.push_back((*o));
    }

    int splits = 0;
    // Check if the clusters should be split.
    for (std::list<Cube::Handle<Cube::ReconCluster>>::iterator
             o = clusterList.begin(); o != clusterList.end();) {
        int split = -1;
        while (true) {
            Cube::Handle<Cube::HitSelection> hits = (*o)->GetHitSelection();
            int scanLength = fScanLength;
            if (hits->size() < scanLength) scanLength = hits->size();
            if (scanLength < fMinimumScanLength) break;
            double length = (hits->front()->GetPosition()
                             - hits->back()->GetPosition()).Mag();
            int trialSplit = -1;
            double biggestKink = 0.0;
            for (int i = scanLength; i < hits->size(); ++i) {
                // The kink distance calculation is really clumsy!
                TVector3 front = hits->at(i-scanLength)->GetPosition();
                TVector3 back = hits->at(i)->GetPosition();
                TVector3 dir = (back - front).Unit();
                for (int j = i-scanLength+1; j<i; ++j) {
                    TVector3 middle = hits->at(j)->GetPosition();
                    double dist = dir*(middle-front);
                    double kink = ((middle-front) - dist*dir).Mag();
                    if (biggestKink < kink) {
                        biggestKink = kink;
                        trialSplit = j;
                    }
                }
            }
            if (biggestKink > length*fLengthFraction) {
                split = trialSplit;
                break;
            }
            if (biggestKink > fKinkThreshold) {
                split = trialSplit;
                break;
            }
            break;
        }
        if (split<0) {
            ++o;
            continue;
        }

        ++splits;
        // Save the cluster to split.
        Cube::Handle<Cube::ReconCluster> splitCluster = (*o);
        Cube::Handle<Cube::HitSelection> splitHits
            = splitCluster->GetHitSelection();

        // Break the cluster in two and push both onto the list.  The order
        // doesn't matter.  REMEMBER, the last hit of the first and the first
        // hit of the last are shared.
        Cube::Handle<Cube::ReconCluster> cluster1
            = Cube::CreateCluster("findKinks",
                                  splitHits->begin(),
                                  splitHits->begin()+split+1);
        clusterList.push_back(cluster1);
        Cube::Handle<Cube::ReconCluster> cluster2
            = Cube::CreateCluster("findKinks",
                                  splitHits->begin()+split,
                                  splitHits->end());
        clusterList.push_back(cluster2);

        // Remove the cluster that was split from the list.
        o = clusterList.erase(o);
    }

    // Copy the clusters into the finalObjects.
    for (std::list<Cube::Handle<Cube::ReconCluster>>::iterator
             o = clusterList.begin();
         o != clusterList.end(); ++o) {
        finalObjects->push_back(*o);
    }

    CUBE_LOG(1) << "Cube::FindKinks splits " << splits << std::endl;

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
