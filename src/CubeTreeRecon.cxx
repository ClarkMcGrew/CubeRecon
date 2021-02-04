#include "CubeTreeRecon.hxx"
#include "CubeHitUtilities.hxx"
#include "CubeMakeUsed.hxx"
#include "CubeClusterHits.hxx"
#include "CubeSpanningTree.hxx"
#include "CubeFindKinks.hxx"
#include "CubeGrowClusters.hxx"
#include "CubeGrowTracks.hxx"
#include "CubeMergeXTalk.hxx"
#include "CubeBuildPairwiseVertices.hxx"

#include <CubeLog.hxx>
#include <CubeInfo.hxx>
#include <CubeHandle.hxx>
#include <CubeAlgorithmResult.hxx>

#include <algorithm>

#ifdef EXPLORE_RUN_TIME
#include <time.h>
#include <TH2F.h>
TH2F* bogo_treeReconTime = NULL;
#endif

Cube::TreeRecon::TreeRecon(): Cube::Algorithm("TreeRecon") { }

Cube::Handle<Cube::AlgorithmResult>
Cube::TreeRecon::Process(const Cube::AlgorithmResult& input,
                         const Cube::AlgorithmResult&,
                         const Cube::AlgorithmResult&) {
    CUBE_LOG(2) << "Process TreeRecon" << std::endl;

    // The input hits are assumed to come from a single time-slice, and have
    // had any ghost hit removal already done.  Every hit in the input
    // *should* be used in the output.  Every hit in the input *should* be a
    // cube hit.
    Cube::Handle<Cube::HitSelection> inputHits = input.GetHitSelection();

#ifdef EXPLORE_RUN_TIME
    double ticks;
    ticks = clock();
#endif

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

    // Check if this is for the tpc.
    bool isTPC = false;
    for (Cube::HitSelection::iterator h = inputHits->begin();
         h != inputHits->end(); ++h) {
        if (!Cube::Info::IsTPC((*h)->GetIdentifier())) continue;
        isTPC = true;
        break;
    }

    // Create the container for the final objects.
    Cube::Handle<Cube::ReconObjectContainer>
        finalObjects(new Cube::ReconObjectContainer("final"));

    // Apply all of the algorithms to the input.  The "do/while" construct
    // implements finalization.  The final objects from currentResult will be
    // copied into the final recon container.
    Cube::Handle<Cube::AlgorithmResult> currentResult;;
    do {
        // Cluster the 3D hits so they are all simply connected.
        Cube::Handle<Cube::AlgorithmResult> clusterHits
            = Run<Cube::ClusterHits>(*inputHits);
        if (!clusterHits) break;
        currentResult = clusterHits;
        result->AddAlgorithmResult(currentResult);

        // Further break the clusters based on the spanning tree.
        Cube::Handle<Cube::AlgorithmResult> spanningTree
            = Run<Cube::SpanningTree>(*currentResult);
        if (!spanningTree) break;
        currentResult = spanningTree;
        result->AddAlgorithmResult(currentResult);

        Cube::Handle<Cube::AlgorithmResult> findKinks;
        if (!isTPC)
            findKinks = Run<Cube::FindKinks>(*currentResult);
        else {
            // Further break the clusters based on kinks.
            std::unique_ptr<Cube::FindKinks> ptr(new Cube::FindKinks);
            ptr->SetScanLength(50);
            ptr->SetMinimumScanLength(30);
            findKinks = ptr->Process(*currentResult);
        }
        if (!findKinks) break;
        currentResult = findKinks;
        result->AddAlgorithmResult(currentResult);

        // Grow the clusters into track like objects.
        Cube::Handle<Cube::AlgorithmResult> growClusters
            = Run<Cube::GrowClusters>(*currentResult);
        if (!growClusters) break;
        currentResult = growClusters;
        result->AddAlgorithmResult(currentResult);

        // Grow the tracks to prevent gaps.
        Cube::Handle<Cube::AlgorithmResult> growTracks
            = Run<Cube::GrowTracks>(*currentResult);
        if (!growTracks) break;
        currentResult = growTracks;
        result->AddAlgorithmResult(currentResult);

        // Grow the tracks to prevent gaps.
        Cube::Handle<Cube::AlgorithmResult> mergeXTalk
            = Run<Cube::MergeXTalk>(*currentResult);
        if (!mergeXTalk) break;
        currentResult = mergeXTalk;
        result->AddAlgorithmResult(currentResult);

    } while (false);            // Always stops...

    // Copy the best current result to the output
    if (currentResult) {
        Cube::Handle<Cube::ReconObjectContainer> objects
            = currentResult->GetObjectContainer("final");
        if (objects) {
            for (Cube::ReconObjectContainer::iterator o
                     = objects->begin();
                 o != objects->end(); ++o) {
                Cube::Handle<Cube::ReconTrack> t = *o;
                if (!t) {
                    CUBE_ERROR << "Only tracks please!!"
                               << std::endl;
                    continue;
                }
                finalObjects->push_back(t);
            }
        }
    }

    // Now build vertices.  This should be after all of the tracks and showers
    // are built.   The tracks in the TPC are not included.
    if (!isTPC) {
        Cube::Handle<Cube::AlgorithmResult> buildVertices
            = Run<Cube::BuildPairwiseVertices>(*currentResult);
        currentResult = buildVertices;
        if (currentResult) {
            result->AddAlgorithmResult(currentResult);
            Cube::Handle<Cube::ReconObjectContainer> objects
                = currentResult->GetObjectContainer("final");
            std::copy(objects->begin(), objects->end(),
                      std::back_inserter(*finalObjects));
        }
    }

    /// Apply a backstop algorithm for any hits that didn't make into the
    /// final tracking result.  This includes very small clusters of hits, and
    /// also very big clusters.  The small clusters are below the tracking
    /// threshold.  The big clusters take to long to process, and so are not
    /// tracked.
    do {
        // Take all of the hits in the final objects, subtract them from the
        // input hits, and then recluster the result.
        Cube::Handle<Cube::HitSelection> finalHits
            = Cube::AllHitSelection(*finalObjects);

        // Make a local copy of the input hits so they can be sorted.
        Cube::HitSelection allHits;
        std::copy(inputHits->begin(), inputHits->end(),
                  std::back_inserter(allHits));

        // Sort the hits so set_difference can be used to remove hits that are
        // in final.
        std::sort(allHits.begin(), allHits.end());
        std::sort(finalHits->begin(), finalHits->end());

        // Take the set difference between all hits and the final hits.
        Cube::HitSelection needsClustering;
        std::set_difference(allHits.begin(),allHits.end(),
                            finalHits->begin(),finalHits->end(),
                            std::back_inserter(needsClustering));

        ///////////////////////////////////////////////////////////////
        // Handle the left over hits.
        ///////////////////////////////////////////////////////////////
        if (needsClustering.size() < 1) break;

        // Cluster the leftover hits.
        std::unique_ptr<Cube::ClusterHits>
            reclusterAlgo(new Cube::ClusterHits);
        reclusterAlgo->SetName("ReclusterHits");
        reclusterAlgo->SetCubeNeighborhood(2);
        reclusterAlgo->SetCubeCount(1);
        Cube::Handle<Cube::AlgorithmResult> recluster
            = reclusterAlgo->Process(needsClustering);
        if (!recluster) break;
        currentResult = recluster;
        Cube::Handle<Cube::ReconObjectContainer> newObjects
            = currentResult->GetObjectContainer("final");
        if (!newObjects) break;

        result->AddAlgorithmResult(currentResult);
        std::copy(newObjects->begin(), newObjects->end(),
                  std::back_inserter(*finalObjects));

    } while (false);

    // Save the final objects last.
    result->AddObjectContainer(finalObjects);

    // Build the selections of used and unused hits.
    Cube::MakeUsed makeUsed(*inputHits);
    result = makeUsed(result);

#ifdef EXPLORE_RUN_TIME
    double newTicks = clock();
    newTicks -= ticks;
    newTicks /= 1000000.0;
    if (newTicks < 0) newTicks += 72.0*60.0;
    std::cout << "**************** TICKS " << newTicks << std::endl;
    if (!bogo_treeReconTime) {
        std::cout << "MAKE HISTOGRAM BOGO_TREERECONTIME" << std::endl;
        bogo_treeReconTime = new TH2F("treeReconTime", "Time vs Hits",
                                      100, 0.0, std::log10(20000.0),
                                      100, 0.0, std::log10(10000));
    }
    double hitCount = inputHits->size();
    std::cout << hitCount << " " << newTicks << std::endl;
    if (hitCount < 1) hitCount = 1;
    if (newTicks < 1) newTicks = 1;
    hitCount = std::log10(hitCount);
    newTicks = std::log10(newTicks);
    if (newTicks < 0.01) newTicks = 0.01;
    bogo_treeReconTime->Fill(hitCount,newTicks);
#endif

    return result;
}

// Local Variables:
// mode:c++
// c-basic-offset:4
// compile-command:"$(git rev-parse --show-toplevel)/build/cube-build.sh force"
// End:
