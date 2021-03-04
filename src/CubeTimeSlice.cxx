#include "CubeTimeSlice.hxx"
#include "CubeClusterManagement.hxx"

#include "CubeUnits.hxx"
#include "CubeInfo.hxx"

#include <CubeReconCluster.hxx>

#include <cmath>
#include <iostream>
#include <algorithm>

namespace {
    struct timeSort {
        bool operator() (const Cube::Handle<Cube::Hit>& lhs,
                         const Cube::Handle<Cube::Hit>& rhs) {
            // Order the times by detector and then time.
            int lhsDet = 0;
            if (Cube::Info::IsTPC(lhs->GetIdentifier())) lhsDet = 1;
            else if (Cube::Info::IsECal(lhs->GetIdentifier())) lhsDet = 2;
            else if (Cube::Info::Is3DST(lhs->GetIdentifier())) lhsDet = 3;
            else lhsDet = 4;
            int rhsDet = 0;
            if (Cube::Info::IsTPC(rhs->GetIdentifier())) rhsDet = 1;
            else if (Cube::Info::IsECal(rhs->GetIdentifier())) rhsDet = 2;
            else if (Cube::Info::Is3DST(rhs->GetIdentifier())) rhsDet = 3;
            else rhsDet = 4;
            if (lhsDet != rhsDet) return lhsDet < rhsDet;
            return lhs->GetTime() < rhs->GetTime();
        }
    };
};

Cube::TimeSlice::TimeSlice()
    : Cube::Algorithm("TimeSlice") {
    fGapCut = 40.0 * unit::ns;
}

Cube::TimeSlice::~TimeSlice() { }

Cube::Handle<Cube::AlgorithmResult>
Cube::TimeSlice::Process(const Cube::AlgorithmResult& in,
                         const Cube::AlgorithmResult&,
                         const Cube::AlgorithmResult&) {
    CUBE_LOG(1) << "TimeSlice::Process" << std::endl;

    // Get the hits and check that they exist.
    Cube::Handle<Cube::HitSelection> hits = in.GetHitSelection();
    if (!hits) return Cube::Handle<Cube::AlgorithmResult>();
    if (hits->empty()) return Cube::Handle<Cube::AlgorithmResult>();

    // Create the new result
    Cube::Handle<Cube::AlgorithmResult> result(CreateResult());

    // Create the container for all of the hits used in this result.
    Cube::Handle<Cube::HitSelection>  usedHits(new Cube::HitSelection("used"));
    result->AddHitSelection(usedHits);

    // Create the container for the final objects.
    Cube::Handle<Cube::ReconObjectContainer> finalObjects(
        new Cube::ReconObjectContainer("final"));
    result->AddObjectContainer(finalObjects);

    // Copy all of the input hits to the output
    std::copy(hits->begin(), hits->end(), std::back_inserter(*usedHits));
    std::sort(usedHits->begin(), usedHits->end(), timeSort());

    // Check for any gaps in the hits and build into separate clusters.
    Cube::HitSelection::iterator first = usedHits->begin();
    Cube::HitSelection::iterator last = usedHits->begin();
    int lastDet = 0;
    while (last != usedHits->end()) {
        Cube::HitSelection::iterator prev = last++;
        int prevDet = lastDet;
        if (last != usedHits->end()) {
            if (Cube::Info::IsTPC((*last)->GetIdentifier())) lastDet = 1;
            else if (Cube::Info::IsECal((*last)->GetIdentifier())) lastDet = 2;
            else if (Cube::Info::Is3DST((*last)->GetIdentifier())) lastDet = 3;
            else lastDet = 4;
        }
        if (last==usedHits->end()
            || (lastDet != 0 && lastDet != prevDet)
            || ((*last)->GetTime()-(*prev)->GetTime())>fGapCut) {
            CUBE_LOG(2) << "Time slice with " << last-first << " `hits "
                        << " from " << (*first)->GetTime()
                        << " to " << (*prev)->GetTime()
                        << std::endl;

            // Build a new cluster.
            Cube::Handle<Cube::ReconCluster> timeCluster
                = Cube::CreateCluster("timeSlice",
                                      first,last);
            finalObjects->push_back(timeCluster);

            // Start a new cluster
            first = last;
        }
    }

    return result;
}

// Local Variables:
// mode:c++
// c-basic-offset:4
// compile-command:"$(git rev-parse --show-toplevel)/build/cube-build.sh force"
// End:
