#include "CubeGrowClusters.hxx"
#include "CubeMakeUsed.hxx"
#include "CubeClusterManagement.hxx"
#include "CubeTrackFit.hxx"
#include "CubeCreateTrack.hxx"
#include "CubeSafeLine.hxx"

#include <CubeLog.hxx>
#include <CubeHandle.hxx>
#include <CubeReconCluster.hxx>
#include <CubeHitSelection.hxx>
#include <CubeAlgorithmResult.hxx>
#include <CubeUnits.hxx>

#include <list>
#include <set>

Cube::GrowClusters::GrowClusters()
    : Cube::Algorithm("GrowClusters") {
    fMaxLineHits = 6;
    // = Cube::TOARuntimeParams::Get().GetParameterI(
    // "sfgRecon.GrowClusters.MaxLineHits");

    fChi2Threshold = 8.0;
    // = Cube::TOARuntimeParams::Get().GetParameterD(
    // "sfgRecon.GrowClusters.Chi2Threshold");

    fChargePerHitThreshold = 15.0;
    // = Cube::TOARuntimeParams::Get().GetParameterD(
    // "sfgRecon.GrowClusters.MinChargePerHit");

}

namespace {
    bool clusterSize(const Cube::Handle<Cube::ReconCluster>& lhs,
                     const Cube::Handle<Cube::ReconCluster>& rhs) {
        return lhs->GetHitSelection()->size() > rhs->GetHitSelection()->size();
    }
}

Cube::Handle<Cube::AlgorithmResult>
Cube::GrowClusters::Process(const Cube::AlgorithmResult& input,
                     const Cube::AlgorithmResult&,
                     const Cube::AlgorithmResult&) {
    CUBE_LOG(0) << "Process Cube::GrowClusters" << std::endl;

    Cube::Handle<Cube::HitSelection> inputHits = input.GetHitSelection();

    // Create the result for this algorithm.
    Cube::Handle<Cube::AlgorithmResult> result = CreateResult();

    // Check that we have hits!
    if (!inputHits) {
        CUBE_ERROR << "No hits" << std::endl;
        return result;
    }
    if (inputHits->empty()) {
        CUBE_ERROR << "Hits empty" << std::endl;
        return result;
    }

    // Create the container for the final objects.
    Cube::Handle<Cube::ReconObjectContainer>
        finalObjects(new Cube::ReconObjectContainer("final"));
    result->AddObjectContainer(finalObjects);

    typedef std::list<Cube::Handle<Cube::ReconCluster>> ClusterList;
    ClusterList clusterList;

    // Copy the input clusters into a list, with error checking...
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

    clusterList.sort(clusterSize);
    double clusterSizeCut = std::log(clusterList.size()+1.0)/2.3;

    CUBE_LOG(0) << "Grow tracks from " << clusterList.size() << " clusters"
                << " w/ minimum seed size of " << clusterSizeCut
                << std::endl;

    // A set of hits that are already interior to a track.
    std::set<Cube::Handle<Cube::Hit>> interiorHits;

    // Iterate over the clusters until we stop growing any clusters.
    Cube::HitSelection combinedHits;
    typedef std::pair<ClusterList::iterator,ClusterList::iterator> ClusterPair;
    int growingClusters = 0;
    int iterations = 0;
    do {
        growingClusters = 0;
        // Save the best matching pair of clusters based on a heuristic for
        // what is meant by the best pair.
        ClusterPair bestPair
            = std::make_pair(clusterList.end(), clusterList.end());
        double bestMatch = -1.0;
        for (ClusterList::iterator c1 = clusterList.begin();
             c1 != clusterList.end(); ++c1) {
            if ((*c1)->GetHitSelection()->size() < clusterSizeCut) break;
            double c1PerHit
                = (*c1)->GetEDeposit() / (*c1)->GetHitSelection()->size();
            if (c1PerHit < fChargePerHitThreshold) continue;
            for (ClusterList::iterator c2 = std::next(c1);
                 c2 != clusterList.end(); ++c2) {
                double c2PerHit
                    = (*c2)->GetEDeposit()/(*c2)->GetHitSelection()->size();
                if (c2PerHit < fChargePerHitThreshold) continue;

                // Are the two iterators for clusters that are neighbors?
                if (!Cube::AreNeighbors(*(*c1), *(*c2))) continue;

                // If we get here, the clusters could plausibly be joined.
                // Check if it's a good idea.
                Cube::HitSelection::iterator sharedHit
                    = Cube::CombineNeighbors((*c1)->GetHitSelection()->begin(),
                                             (*c1)->GetHitSelection()->end(),
                                             (*c2)->GetHitSelection()->begin(),
                                             (*c2)->GetHitSelection()->end(),
                                             combinedHits);

                // Check to see if the shared hit is already an interior hit.
                if (interiorHits.find(*sharedHit) != interiorHits.end()) {
                    continue;
                }

                // Find the fit range for the fits.
                Cube::HitSelection::iterator s1 = combinedHits.begin();
                Cube::HitSelection::iterator e1 = sharedHit+1;
                Cube::HitSelection::iterator s2 = sharedHit;
                Cube::HitSelection::iterator e2 = combinedHits.end();

                // Compress the range to something reasonable...
                if (sharedHit - combinedHits.begin() > fMaxLineHits) {
                    s1 = sharedHit - fMaxLineHits;
                }
                if (combinedHits.end() - sharedHit > fMaxLineHits) {
                    e2 = s2 + fMaxLineHits;
                }

                // Protect against acute angles.
                double ang
                    = ((*sharedHit)->GetPosition() - (*s1)->GetPosition())
                    * ((*(e2-1))->GetPosition() - (*sharedHit)->GetPosition());
                if (ang <= 0.0) continue;

                // Temporaries used to find the line fit for the segment
                // of line before the shared hit.
                TVector3 pos1, dir1;
                double chi1, ndof1;
                Cube::SafeLine(s1,e1,pos1,dir1,chi1,ndof1);

                // Temporaries used to find the line fit for the segment
                // of line after the shared hit.
                TVector3 pos2, dir2;
                double chi2, ndof2;
                Cube::SafeLine(s2,e2,pos2,dir2,chi2,ndof2);

                // Temporaries used to find the line fit for the combination
                // of segements (including the shared hit).
                TVector3 pos3, dir3;
                double chi3, ndof3;
                Cube::SafeLine(s1,e2,pos3,dir3,chi3,ndof3);

                // Find the chi-squared for the new fit relative to the two
                // segments.  The dChi is adding one degree of freedom, so
                // follows the chi-squared for a single degree of freedom.
                double dChi = chi3 - chi2 - chi1;

                // Check to see if the fit is good enough to combine.  If not,
                // just skip to the next pair of clusters.
                if (fChi2Threshold < dChi) continue;

                // Calculate the heuristic that defines the order that
                // clusters are combined into larger clusters.  This favors
                // large input clusters with a "good" line fit.  Notice that
                // the values of n1 and n2 are capped by the value of
                // fMaxLineHits so that once clusters get over that size, the
                // decision is based on the quality of the line fit.

                // Favor combining larger clusters.
                double n1 = 1.0*(e1 - s1);
                double n2 = 1.0*(e2 - s2);
                double sizeHeuristic = n1*n2;
                double sizeWeight = 1.0;

                // Favor combining high charge clusters.
                double q1 = (*s1)->GetCharge();
                for (Cube::HitSelection::iterator h = s1; h != e1; ++h) {
                    q1 = std::min(q1,(*h)->GetCharge());
                }
                double q2 = (*s2)->GetCharge();
                for (Cube::HitSelection::iterator h = s2; h != e2; ++h) {
                    q2 = std::min(q2,(*h)->GetCharge());
                }
                double chargeHeuristic = q1*q2;
                // A mip is about 40.  This is 20^2.
                double chargeWeight = 1.0/400.0;

                // Favor good fits
                double goodnessHeuristic = - dChi;
                double goodnessWeight = 1.0;

                double heuristic = sizeWeight*sizeHeuristic;
                heuristic += chargeWeight*chargeHeuristic;
                heuristic += goodnessWeight*goodnessHeuristic;

                if (heuristic > bestMatch) {
                    // If we get here, there is a pair of clusters to combine,
                    // and this might be the one.
                    ++growingClusters;
                    bestMatch = heuristic;
                    bestPair = std::make_pair(c1,c2);
                }
            }
        }

        CUBE_LOG(3) << "Iteration: " << ++iterations
                    << " " << clusterList.size()
                    << " " << growingClusters
                    << std::endl;

        // Check if we found any clusters to grow and stop if we don't.
        if (growingClusters<1) break;

        // Combine the hits from the clusters.
        Cube::HitSelection newClusterHits;
        Cube::HitSelection::iterator sharedHit =
            Cube::CombineNeighbors(
                (*bestPair.first)->GetHitSelection()->begin(),
                (*bestPair.first)->GetHitSelection()->end(),
                (*bestPair.second)->GetHitSelection()->begin(),
                (*bestPair.second)->GetHitSelection()->end(),
                newClusterHits);

        // Add the shared hit to the set of hits that can't be added back to
        // a track.  (This prevents crossing tracks).
        interiorHits.insert(*sharedHit);

        // Make the new cluster
        Cube::Handle<Cube::ReconCluster> newCluster
            = Cube::CreateCluster("clusterGrow",
                                     newClusterHits.begin(),
                                     newClusterHits.end());

        // Remove the two clusters being combined from the list (after
        // combining them so the iterators are valid for that step).
        clusterList.erase(bestPair.first);
        clusterList.erase(bestPair.second);

        // Add the new cluster to the list.
        clusterList.push_front(newCluster);

    } while (growingClusters > 0);

    // Copy the remaining clusters into the finalObjects.  Everything in the
    // clusterList has already been combined as much as possible, and is a
    // candidate "track".
    Cube::TrackFit fitter;
    for (ClusterList::iterator c1 = clusterList.begin();
         c1 != clusterList.end(); ++c1) {
        if ((*c1)->GetHitSelection()->size() < 4) {
            finalObjects->push_back((*c1));
            continue;
        }
        if (((*c1)->GetEDeposit()/(*c1)->GetHitSelection()->size())
            < fChargePerHitThreshold) {
            finalObjects->push_back((*c1));
            continue;
        }
        // Make sure that the tracks "go" in the positive Z direction.  This
        // is a convienence so that typical tracks are all headed in the same
        // directions.  If the track is perpendicular to Z, the direction is
        // left unchanged.
        if ((*c1)->GetHitSelection()->front()->GetPosition().Z()
            > (*c1)->GetHitSelection()->back()->GetPosition().Z()) {
            std::reverse((*c1)->GetHitSelection()->begin(),
                         (*c1)->GetHitSelection()->end());
        }
        Cube::Handle<Cube::ReconTrack> track
            = Cube::CreateTrackFromHits("growClusters",
                                        (*c1)->GetHitSelection()->begin(),
                                        (*c1)->GetHitSelection()->end());
        track = fitter(track);
        finalObjects->push_back(track);
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
