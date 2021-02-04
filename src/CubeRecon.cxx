#include "CubeRecon.hxx"
#include "CubeCleanHits.hxx"
#include "CubeTimeSlice.hxx"
#include "CubeTreeRecon.hxx"
#include "CubeMakeUsed.hxx"
#include "CubeClusterManagement.hxx"
#include "CubeCompareReconObjects.hxx"

#include <CubeLog.hxx>
#include <CubeInfo.hxx>
#include <CubeHandle.hxx>
#include <CubeReconCluster.hxx>
#include <CubeAlgorithmResult.hxx>

#include <sstream>
#include <iomanip>

Cube::Recon::Recon(): Cube::Algorithm("CubeRecon") {
    fOversizeCut = 20000;
}

Cube::Handle<Cube::AlgorithmResult>
Cube::Recon::Process(const Cube::AlgorithmResult& input,
                     const Cube::AlgorithmResult&,
                     const Cube::AlgorithmResult&) {
    CUBE_LOG(0) << "Process Recon" << std::endl;
    Cube::Handle<Cube::HitSelection> inputHits = input.GetHitSelection();

    // Create the result for this algorithm.
    Cube::Handle<Cube::AlgorithmResult> result = CreateResult();

    // Check that we have hits!
    if (!inputHits) {
        return result;
    }
    if (inputHits->empty()) {
        return result;
    }

    // Clean out junk hits.
    Cube::Handle<Cube::AlgorithmResult> cleanHits
        = Run<Cube::CleanHits>(*inputHits);
    if (!cleanHits) return Cube::Handle<Cube::AlgorithmResult>();
    result->AddAlgorithmResult(cleanHits);

    // Slice the event up by time.
    Cube::Handle<Cube::AlgorithmResult> timeSlice
        = Run<Cube::TimeSlice>(*cleanHits);
    if (!timeSlice) return Cube::Handle<Cube::AlgorithmResult>();
    result->AddAlgorithmResult(timeSlice);

    // Create the container for the final objects.
    Cube::Handle<Cube::ReconObjectContainer>
        finalObjects(new Cube::ReconObjectContainer("final"));

    // Process each time slice as a separate cube reconstruction.
    Cube::Handle<Cube::ReconObjectContainer> slices
        = timeSlice->GetObjectContainer("final");
    int count = 0;
    for (Cube::ReconObjectContainer::iterator object = slices->begin();
         object != slices->end(); ++object) {
        Cube::Handle<Cube::HitSelection> hits = (*object)->GetHitSelection();
        if (!hits || hits->empty()) {
            CUBE_ERROR << "No hits is slice" << std::endl;
            continue;
        }

        // Protect against really large time slices.
        if (hits->size() > fOversizeCut) {
            CUBE_ERROR << "Skip oversize slice with "
                       << hits->size() << " hits"
                       << std::endl;
            Cube::Handle<Cube::ReconCluster> cluster
                = Cube::CreateCluster("cubeReconLarge",
                                      hits->begin(), hits->end());
            finalObjects->push_back(cluster);
            continue;
        }

        // We need to make a local copy to allow the HitSelection to be
        // translated into a AlgorithmResult.
        Cube::HitSelection local("cubes");
        std::copy(hits->begin(), hits->end(), std::back_inserter(local));
        Cube::Handle<Cube::AlgorithmResult>
            treeRecon = Run<Cube::TreeRecon>(local);
        if (!treeRecon) {
            CUBE_ERROR << "Unsuccessful cube reconstruction" << std::endl;
            continue;
        }

        // Build the name and add it as a sub algorithm.
        std::ostringstream nm;
        nm << "Slice" << std::setfill('0') << std::setw(2) << count++;
        treeRecon->SetName(nm.str().c_str());
        result->AddAlgorithmResult(treeRecon);
        // Copy the results to final.
        Cube::Handle<Cube::ReconObjectContainer> objects
            = treeRecon->GetObjectContainer("final");
        if (objects) {
            std::copy(objects->begin(), objects->end(),
                      std::back_inserter(*finalObjects));
        }
    }

    std::sort(finalObjects->begin(), finalObjects->end(),
              Cube::CompareReconObjects());

    // Save the final objects last.
    result->AddObjectContainer(finalObjects);

    // Collect all of the hits into the used and unused hit selections.
    Cube::MakeUsed makeUsed(*inputHits);
    result = makeUsed(result);

    return result;
}

// Local Variables:
// mode:c++
// c-basic-offset:4
// compile-command:"$(git rev-parse --show-toplevel)/build/cube-build.sh force"
// End:
