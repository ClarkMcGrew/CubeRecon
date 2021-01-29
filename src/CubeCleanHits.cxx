#include "CubeCleanHits.hxx"
#include "CubeMakeUsed.hxx"
#include "CubeClusterManagement.hxx"

#include <CubeAlgorithmResult.hxx>
#include <CubeLog.hxx>
#include <CubeInfo.hxx>
#include <CubeHandle.hxx>
#include <CubeHit.hxx>
#include <CubeHitSelection.hxx>

Cube::CleanHits::CleanHits()
    : Cube::Algorithm("CleanHits") {
    fGhostHitThreshold = 9.5;
}

Cube::Handle<Cube::AlgorithmResult>
Cube::CleanHits::Process(const Cube::AlgorithmResult& input,
                         const Cube::AlgorithmResult&,
                         const Cube::AlgorithmResult&) {
    Cube::Handle<Cube::HitSelection> inputHits = input.GetHitSelection();
    CUBE_LOG(0) << "Process CleanHits" << std::endl;

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

    // Create the containers for the final objects and hits.
    Cube::Handle<Cube::HitSelection> unusedHits(
        new Cube::HitSelection("unused"));
    result->AddHitSelection(unusedHits);

    Cube::Handle<Cube::HitSelection> usedHits(
        new Cube::HitSelection("used"));
    result->AddHitSelection(usedHits);

    // The input hits possibly a include mix of simple and composite hits.
    for (Cube::HitSelection::iterator hit = inputHits->begin();
         hit != inputHits->end(); ++hit) {
        if (Cube::Info::Is3DST((*hit)->GetIdentifier())
            && (*hit)->GetCharge() < fGhostHitThreshold) {
            unusedHits->push_back(*hit);
            continue;
        }
        usedHits->push_back(*hit);
    }

    CUBE_LOG(0) << "CleanHits: Hits rejected " << unusedHits->size()
                << " leaving " << usedHits->size() << " hits" << std::endl;

    return result;
}

// Local Variables:
// mode:c++
// c-basic-offset:4
// compile-command:"$(git rev-parse --show-toplevel)/build/cube-build.sh force"
// End:
